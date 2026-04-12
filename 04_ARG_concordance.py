#!/usr/bin/env python3
import os
import re
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Attach ARG contigs to MAG bins and backbone plasmid clusters")
parser.add_argument("--sample", required=True, help="Sample ID")
parser.add_argument("--basedir", default=".", help="Base directory containing per-sample subdirectories")
args = parser.parse_args()

SAMPLE = args.sample.strip()

# =========================
# Paths
# =========================

sample_dir = os.path.join(args.basedir, SAMPLE)

cleaned_edges_fp   = os.path.join(sample_dir, f"{SAMPLE}_cleaned.txt")
contig_reads_fp    = os.path.join(sample_dir, f"{SAMPLE}_contig_reads.csv")
binmap_fp          = os.path.join(sample_dir, f"contig_bin_mapping_{SAMPLE}.txt")
cluster_fp         = os.path.join(sample_dir, f"{SAMPLE}_cluster_membership_filtered.csv")
norm_fp            = os.path.join(sample_dir, f"{SAMPLE}_normalized_filtered.csv")
amrfinder_fp       = os.path.join(sample_dir, f"{SAMPLE}_amrfinder_AMR.tsv")
plasmid_summary_fp = os.path.join(sample_dir, f"{SAMPLE}_final.contigs_plasmid_summary.tsv")

out_dir = os.path.join(sample_dir, "ARG_mapping_out")
os.makedirs(out_dir, exist_ok=True)

# =========================
# Helpers
# =========================
def pick_first_existing_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise ValueError(f"None of these columns exist: {candidates}. Found: {list(df.columns)[:50]}")

def safe_mean(series):
    s = pd.to_numeric(series, errors="coerce")
    s = s.replace([np.inf, -np.inf], np.nan).dropna()
    return float(s.mean()) if len(s) else np.nan

def add_fractions(df, key_cols, value_col):
    # per ARG, compute fraction across targets
    totals = df.groupby(key_cols[0])[value_col].sum().rename("total").reset_index()
    out = df.merge(totals, on=key_cols[0], how="left")
    out["fraction"] = np.where(out["total"] > 0, out[value_col] / out["total"], np.nan)
    return out.drop(columns=["total"])

def topcall_table(df, arg_col, target_col, score_col, prefix):
    # returns one row per ARG: top target + top score + top fraction + top1/top2 ratio
    df = df.copy()
    df = df.sort_values([arg_col, score_col], ascending=[True, False])
    # top1
    top1 = df.groupby(arg_col).head(1).copy()
    top1 = top1[[arg_col, target_col, score_col, "fraction"]].rename(columns={
        target_col: f"{prefix}_top",
        score_col: f"{prefix}_top_score",
        "fraction": f"{prefix}_top_fraction"
    })
    # top2 score for ratio
    top2 = df.groupby(arg_col).nth(1).reset_index()
    if len(top2):
        top2 = top2[[arg_col, score_col]].rename(columns={score_col: f"{prefix}_second_score"})
        top1 = top1.merge(top2, on=arg_col, how="left")
        top1[f"{prefix}_top1_top2_ratio"] = np.where(
            top1[f"{prefix}_second_score"] > 0,
            top1[f"{prefix}_top_score"] / top1[f"{prefix}_second_score"],
            np.nan
        )
        top1 = top1.drop(columns=[f"{prefix}_second_score"])
    else:
        top1[f"{prefix}_top1_top2_ratio"] = np.nan
    return top1

# =========================
# 1) Load contig reads (len + mapped)
# =========================
reads = pd.read_csv(contig_reads_fp, header=None)
# Handle possible separators if needed
if reads.shape[1] == 1:
    reads = pd.read_csv(contig_reads_fp, sep="\t", header=None)
if reads.shape[1] == 1:
    reads = pd.read_csv(contig_reads_fp, sep=r"\s+", header=None)

if reads.shape[1] < 3:
    raise ValueError(f"Unexpected contig_reads format: {contig_reads_fp} has {reads.shape[1]} columns.")

reads = reads.iloc[:, :4]  # keep first 4 if more
reads.columns = ["contig_id", "len", "mapped", "unmapped"]
reads["contig_id"] = reads["contig_id"].astype(str).str.strip()
reads["len"] = pd.to_numeric(reads["len"], errors="coerce")
reads["mapped"] = pd.to_numeric(reads["mapped"], errors="coerce")

mean_len = safe_mean(reads["len"])
mean_mapped = safe_mean(reads.loc[reads["mapped"] > 0, "mapped"])  # avoid tons of zeros
if np.isnan(mean_mapped):
    mean_mapped = safe_mean(reads["mapped"])

print(f"[{SAMPLE}] contigs in reads: {len(reads):,} | mean_len={mean_len:.2f} | mean_mapped={mean_mapped:.2f}")

# =========================
# 2) Load edges (unfiltered)
# =========================
edges = pd.read_csv(cleaned_edges_fp, sep=r"\s+|\t+", engine="python", header=None)
if edges.shape[1] < 3:
    raise ValueError(f"Unexpected cleaned edge format: {cleaned_edges_fp} has {edges.shape[1]} columns.")
edges = edges.iloc[:, :3]
edges.columns = ["A", "B", "V3"]
edges["A"] = edges["A"].astype(str).str.strip()
edges["B"] = edges["B"].astype(str).str.strip()
edges["V3"] = pd.to_numeric(edges["V3"], errors="coerce").fillna(0.0)

print(f"[{SAMPLE}] edges loaded: {len(edges):,}")

# Attach contig metadata to both ends
reads_small = reads[["contig_id", "len", "mapped"]].copy()
edges = edges.merge(reads_small, left_on="A", right_on="contig_id", how="left").drop(columns=["contig_id"])
edges = edges.rename(columns={"len": "A_len", "mapped": "A_mapped"})
edges = edges.merge(reads_small, left_on="B", right_on="contig_id", how="left").drop(columns=["contig_id"])
edges = edges.rename(columns={"len": "B_len", "mapped": "B_mapped"})

# Comprehensive contig-level normalized weight (both ends)
# W = V3 * (mean_mapped/A_mapped)*(mean_mapped/B_mapped) * (mean_len/A_len)*(mean_len/B_len)
def compute_w(row):
    a_m, b_m = row["A_mapped"], row["B_mapped"]
    a_l, b_l = row["A_len"], row["B_len"]
    if pd.isna(a_m) or pd.isna(b_m) or pd.isna(a_l) or pd.isna(b_l):
        return np.nan
    if a_m <= 0 or b_m <= 0 or a_l <= 0 or b_l <= 0:
        return np.nan
    return row["V3"] * (mean_mapped / a_m) * (mean_mapped / b_m) * (mean_len / a_l) * (mean_len / b_l)

edges["W_contig_norm"] = edges.apply(compute_w, axis=1)

# =========================
# Load geNomad plasmid summary (plasmid contig set)
# =========================
plasmid_contigs = set()
if os.path.exists(plasmid_summary_fp):
    ps = pd.read_csv(plasmid_summary_fp, sep="\t")
    if "seq_name" not in ps.columns:
        raise ValueError(f"'seq_name' column not found in {plasmid_summary_fp}. Columns: {list(ps.columns)}")
    plasmid_contigs = set(ps["seq_name"].astype(str).str.strip().unique())
    print(f"[{SAMPLE}] plasmid contigs (geNomad): {len(plasmid_contigs):,}")
else:
    print(f"[{SAMPLE}] WARNING: plasmid summary not found: {plasmid_summary_fp} (ARG_location will be unknown)")


# =========================
# 3) Load AMR contigs (ARG contigs)
# =========================
amr = pd.read_csv(amrfinder_fp, sep="\t", engine="python")
# Try to find a contig column robustly
contig_col = pick_first_existing_col(amr, ["Contig id", "Contig_ID", "Contig", "Contig.id", "contig", "contig_id"])
amr_contigs = set(amr[contig_col].astype(str).str.strip().unique())
print(f"[{SAMPLE}] AMR contigs found: {len(amr_contigs):,}")

# optional gene/subclass columns if present
gene_col = None
subc_col = None
for c in ["Gene symbol", "Gene_symbol", "gene_symbol", "Gene"]:
    if c in amr.columns:
        gene_col = c
        break
for c in ["Subclass", "subclass", "Class", "class"]:
    if c in amr.columns:
        subc_col = c
        break

amr_anno = amr[[contig_col]].copy()
amr_anno.columns = ["ARG_contig"]
amr_anno["ARG_contig"] = amr_anno["ARG_contig"].astype(str).str.strip()
if gene_col:
    amr_anno["ARG_gene"] = amr[gene_col].astype(str)
if subc_col:
    amr_anno["ARG_subclass"] = amr[subc_col].astype(str)
amr_anno = amr_anno.drop_duplicates()

# Label ARG contigs by geNomad plasmid call (plasmid vs non_plasmid vs unknown)
if len(plasmid_contigs) > 0:
    amr_anno["ARG_location"] = np.where(amr_anno["ARG_contig"].isin(plasmid_contigs), "plasmid", "non_plasmid")
else:
    amr_anno["ARG_location"] = "unknown"

# Only plasmid-labeled (and optionally unknown) are eligible for plasmid-cluster attachment claims
# Choose policy: plasmid_only or plasmid_or_unknown
POLICY = "plasmid_only"   # change to "plasmid_or_unknown" if you want to include unknowns

if POLICY == "plasmid_only":
    amr_anno["eligible_for_plasmid_attachment"] = (amr_anno["ARG_location"] == "plasmid")
else:
    amr_anno["eligible_for_plasmid_attachment"] = (amr_anno["ARG_location"] != "non_plasmid")


# Sanity: how many ARG contigs appear in edges?
arg_in_edges = len(set(edges["A"]).intersection(amr_contigs) | set(edges["B"]).intersection(amr_contigs))
print(f"[{SAMPLE}] AMR contigs that appear in edge list: {arg_in_edges:,}")

# =========================
# 4) Load contig -> Bin mapping
# =========================
binmap = pd.read_csv(binmap_fp, sep=r"\s+|\t+", engine="python")
contig_id_col = pick_first_existing_col(binmap, ["Contig_ID", "ContigID", "contig_id", "contig"])
bin_col       = pick_first_existing_col(binmap, ["Bin", "bin", "MAG", "mag"])
binmap = binmap[[contig_id_col, bin_col]].copy()
binmap.columns = ["contig", "Bin"]
binmap["contig"] = binmap["contig"].astype(str).str.strip()
binmap["Bin"] = binmap["Bin"].astype(str).str.strip()

# =========================
# 5) Load backbone cluster membership (contig -> cluster)
# =========================
clu = pd.read_csv(cluster_fp)
node_col = pick_first_existing_col(clu, ["Node", "node", "contig", "contig_id"])
clu_col  = pick_first_existing_col(clu, ["Cluster", "cluster", "membership"])
clu = clu[[node_col, clu_col]].copy()
clu.columns = ["contig", "Cluster"]
clu["contig"] = clu["contig"].astype(str).str.strip()
clu["Cluster"] = clu["Cluster"].astype(str).str.strip()

# =========================
# 6) Use existing normalized_filtered.csv for TotalMappedMAG(Bin)
# =========================
norm = pd.read_csv(norm_fp)
bin_norm_col = pick_first_existing_col(norm, ["Bin", "bin"])
tm_col = pick_first_existing_col(norm, ["Total_Mapped_Reads_MAG", "Total_Mapped_Reads", "TotalMappedMAG"])
cluster_norm_col = pick_first_existing_col(norm, ["Cluster", "cluster"])
nhic_col = pick_first_existing_col(norm, ["Normalized_HiC_Contacts", "Log_Normalized_HiC_Contacts", "V3"])

norm_bin = norm[[bin_norm_col, tm_col]].copy()
norm_bin.columns = ["Bin", "TotalMappedMAG"]
norm_bin["Bin"] = norm_bin["Bin"].astype(str).str.strip()
norm_bin["TotalMappedMAG"] = pd.to_numeric(norm_bin["TotalMappedMAG"], errors="coerce")
# One value per bin (take median to be robust)
norm_bin = norm_bin.groupby("Bin", as_index=False)["TotalMappedMAG"].median()
mean_TotalMappedMAG = safe_mean(norm_bin["TotalMappedMAG"])
print(f"[{SAMPLE}] bins with TotalMappedMAG: {len(norm_bin):,} | mean_TotalMappedMAG={mean_TotalMappedMAG:.2f}")

# Also compute Bin->Cluster strength for triangulation
norm_bc = norm[[bin_norm_col, cluster_norm_col, nhic_col]].copy()
norm_bc.columns = ["Bin", "Cluster", "BinClusterStrength"]
norm_bc["Bin"] = norm_bc["Bin"].astype(str).str.strip()
norm_bc["Cluster"] = norm_bc["Cluster"].astype(str).str.strip()
norm_bc["BinClusterStrength"] = pd.to_numeric(norm_bc["BinClusterStrength"], errors="coerce").fillna(0.0)
norm_bc = norm_bc.groupby(["Bin", "Cluster"], as_index=False)["BinClusterStrength"].sum()
# top cluster per bin
norm_bc = norm_bc.sort_values(["Bin", "BinClusterStrength"], ascending=[True, False])
bin_top_cluster = norm_bc.groupby("Bin").head(1).rename(columns={
    "Cluster": "Bin_top_Cluster",
    "BinClusterStrength": "Bin_top_Cluster_strength"
})[["Bin", "Bin_top_Cluster", "Bin_top_Cluster_strength"]]

# =========================
# 7) Build oriented ARG-edge tables
# =========================
# A is ARG, partner is B
argA = edges[edges["A"].isin(amr_contigs)].copy()
argA = argA.rename(columns={"A": "ARG_contig", "B": "partner"})
argA["partner_side"] = "B"

# B is ARG, partner is A
argB = edges[edges["B"].isin(amr_contigs)].copy()
argB = argB.rename(columns={"B": "ARG_contig", "A": "partner"})
argB["partner_side"] = "A"

arg_edges = pd.concat([argA, argB], ignore_index=True)
arg_edges["ARG_contig"] = arg_edges["ARG_contig"].astype(str)
arg_edges["partner"] = arg_edges["partner"].astype(str)

# Keep only needed columns + attach partner contig reads for sanity (optional)
# (We already have A_len/A_mapped/B_len/B_mapped, but after orientation it’s messy;
#  we keep the precomputed W_contig_norm which already used both ends.)
arg_edges = arg_edges[["ARG_contig", "partner", "V3", "W_contig_norm"]].copy()

# =========================
# 8) ARG -> Bin (host)
# =========================
arg_bin = arg_edges.merge(binmap, left_on="partner", right_on="contig", how="left").drop(columns=["contig"])
arg_bin = arg_bin.dropna(subset=["Bin"])

# Raw support
arg_bin_raw = arg_bin.groupby(["ARG_contig", "Bin"], as_index=False)["V3"].sum().rename(columns={"V3": "support_raw"})
arg_bin_raw = add_fractions(arg_bin_raw, ["ARG_contig", "Bin"], "support_raw")
arg_bin_raw = arg_bin_raw.merge(amr_anno, on="ARG_contig", how="left")

# Contig-normalized support
arg_bin_norm = arg_bin.groupby(["ARG_contig", "Bin"], as_index=False)["W_contig_norm"].sum().rename(columns={"W_contig_norm": "support_contig_norm"})
arg_bin_norm = add_fractions(arg_bin_norm, ["ARG_contig", "Bin"], "support_contig_norm")
arg_bin_norm = arg_bin_norm.merge(amr_anno, on="ARG_contig", how="left")

# Contig + Bin abundance normalization (use TotalMappedMAG(Bin))
arg_bin_norm2 = arg_bin_norm.merge(norm_bin, on="Bin", how="left")
arg_bin_norm2["support_contig_bin_norm"] = arg_bin_norm2["support_contig_norm"] * (mean_TotalMappedMAG / arg_bin_norm2["TotalMappedMAG"])
arg_bin_norm2 = arg_bin_norm2.drop(columns=["fraction"])  # recompute fraction on new support
arg_bin_norm2 = add_fractions(arg_bin_norm2, ["ARG_contig", "Bin"], "support_contig_bin_norm")

# Save outputs
arg_bin_raw.to_csv(os.path.join(out_dir, f"{SAMPLE}_ARG_to_Bin_raw.csv"), index=False)
arg_bin_norm.to_csv(os.path.join(out_dir, f"{SAMPLE}_ARG_to_Bin_norm_contig.csv"), index=False)
arg_bin_norm2.to_csv(os.path.join(out_dir, f"{SAMPLE}_ARG_to_Bin_norm_contig_bin.csv"), index=False)

print(f"[{SAMPLE}] ARG->Bin rows (raw): {len(arg_bin_raw):,} | (contig-norm): {len(arg_bin_norm):,} | (contig+bin): {len(arg_bin_norm2):,}")

# =========================
# 9) ARG -> Cluster (plasmid lineage via backbone contigs)
# =========================
arg_clu = arg_edges.merge(clu, left_on="partner", right_on="contig", how="left").drop(columns=["contig"])
arg_clu = arg_clu.dropna(subset=["Cluster"])

# Raw
arg_clu_raw = arg_clu.groupby(["ARG_contig", "Cluster"], as_index=False)["V3"].sum().rename(columns={"V3": "support_raw"})
arg_clu_raw = add_fractions(arg_clu_raw, ["ARG_contig", "Cluster"], "support_raw")
arg_clu_raw = arg_clu_raw.merge(amr_anno, on="ARG_contig", how="left")

# Contig-normalized
arg_clu_norm = arg_clu.groupby(["ARG_contig", "Cluster"], as_index=False)["W_contig_norm"].sum().rename(columns={"W_contig_norm": "support_contig_norm"})
arg_clu_norm = add_fractions(arg_clu_norm, ["ARG_contig", "Cluster"], "support_contig_norm")
arg_clu_norm = arg_clu_norm.merge(amr_anno, on="ARG_contig", how="left")

arg_clu_raw.to_csv(os.path.join(out_dir, f"{SAMPLE}_ARG_to_Cluster_raw.csv"), index=False)
arg_clu_norm.to_csv(os.path.join(out_dir, f"{SAMPLE}_ARG_to_Cluster_norm_contig.csv"), index=False)

print(f"[{SAMPLE}] ARG->Cluster rows (raw): {len(arg_clu_raw):,} | (contig-norm): {len(arg_clu_norm):,}")

# =========================
# 10) Triangulation (bins: keep fraction>=0.2; clusters: keep all)
# =========================

BIN_FRAC_THRESH = 0.0  # your request

# ---- Build Bin -> Cluster ranking table from your existing normalized table ----
rank_tbl = norm_bc.copy()
rank_tbl = rank_tbl.sort_values(["Bin", "BinClusterStrength"], ascending=[True, False])
rank_tbl["bin_cluster_rank"] = rank_tbl.groupby("Bin").cumcount() + 1
rank_tbl["bin_max_strength"] = rank_tbl.groupby("Bin")["BinClusterStrength"].transform("max")
rank_tbl["bin_rel_strength"] = np.where(
    rank_tbl["bin_max_strength"] > 0,
    rank_tbl["BinClusterStrength"] / rank_tbl["bin_max_strength"],
    np.nan
)

# -------------------------
# A) ARG -> Bin candidates (keep ALL bins with fraction >= 0.2)
#    Use contig+bin normalized mapping: arg_bin_norm2
# -------------------------
arg_bin_candidates = arg_bin_norm2.copy()
arg_bin_candidates["fraction"] = pd.to_numeric(arg_bin_candidates["fraction"], errors="coerce")
arg_bin_candidates = arg_bin_candidates[arg_bin_candidates["fraction"] >= BIN_FRAC_THRESH].copy()

# If an ARG has no bins passing threshold, keep placeholder to mark "missing"
all_args = pd.DataFrame({"ARG_contig": sorted(amr_anno["ARG_contig"].astype(str).unique())})
has_any_bin = arg_bin_candidates[["ARG_contig"]].drop_duplicates()
missing_bin_args = all_args.merge(has_any_bin, on="ARG_contig", how="left", indicator=True)
missing_bin_args = missing_bin_args[missing_bin_args["_merge"] == "left_only"][["ARG_contig"]]

if len(missing_bin_args):
    tmp = missing_bin_args.copy()
    tmp["Bin"] = np.nan
    tmp["support_contig_bin_norm"] = np.nan
    tmp["fraction"] = np.nan
    # carry annotation columns if present
    for c in ["ARG_gene", "ARG_subclass", "ARG_location", "eligible_for_plasmid_attachment", "TotalMappedMAG", "support_contig_norm"]:
        if c in arg_bin_candidates.columns and c not in tmp.columns:
            tmp[c] = np.nan
    arg_bin_candidates = pd.concat([arg_bin_candidates, tmp], ignore_index=True)

arg_bin_candidates.to_csv(
    os.path.join(out_dir, f"{SAMPLE}_ARG_to_Bin_candidates_frac_ge_{BIN_FRAC_THRESH}.csv"),
    index=False
)

print(f"[{SAMPLE}] ARG->Bin candidate rows (fraction>={BIN_FRAC_THRESH}): {len(arg_bin_candidates):,}")

# -------------------------
# B) ARG -> Cluster list (keep ALL clusters)
#    Use contig-normalized mapping: arg_clu_norm
# -------------------------
arg_cluster_all = arg_clu_norm.copy()
if "fraction" in arg_cluster_all.columns:
    arg_cluster_all["fraction"] = pd.to_numeric(arg_cluster_all["fraction"], errors="coerce")

# Add placeholder rows for ARGs with no cluster links
has_any_cluster = arg_cluster_all[["ARG_contig"]].drop_duplicates()
missing_cluster_args = all_args.merge(has_any_cluster, on="ARG_contig", how="left", indicator=True)
missing_cluster_args = missing_cluster_args[missing_cluster_args["_merge"] == "left_only"][["ARG_contig"]]

if len(missing_cluster_args):
    tmp = missing_cluster_args.copy()
    tmp["Cluster"] = np.nan
    tmp["support_contig_norm"] = np.nan
    tmp["fraction"] = np.nan
    for c in ["ARG_gene", "ARG_subclass", "ARG_location", "eligible_for_plasmid_attachment"]:
        if c in arg_cluster_all.columns and c not in tmp.columns:
            tmp[c] = np.nan
    arg_cluster_all = pd.concat([arg_cluster_all, tmp], ignore_index=True)

arg_cluster_all.to_csv(os.path.join(out_dir, f"{SAMPLE}_ARG_to_Cluster_all.csv"), index=False)

# -------------------------
# C) Per (ARG, Cluster): test if ANY candidate Bin supports that cluster
#     (Bins kept earlier with fraction>=BIN_FRAC_THRESH; clusters = all)
# -------------------------

# 1) Rename columns to avoid merge suffix conflicts
#    - cluster side:
arg_cluster_all = arg_cluster_all.rename(columns={
    "fraction": "cluster_fraction",
    "support_contig_norm": "cluster_support_contig_norm"
})
#    - bin side:
arg_bin_candidates = arg_bin_candidates.rename(columns={
    "fraction": "bin_fraction",
    "support_contig_bin_norm": "bin_support_contig_bin_norm"
})

# 2) Build long table: (ARG, Cluster) × (candidate Bins for that ARG)
tri_long = arg_cluster_all.merge(
    arg_bin_candidates[["ARG_contig", "Bin", "bin_fraction", "bin_support_contig_bin_norm"]],
    on="ARG_contig",
    how="left"
)

# 3) Attach Bin→Cluster rank/rel_strength evidence (from normalized MAG↔cluster table)
tri_long = tri_long.merge(
    rank_tbl[["Bin", "Cluster", "bin_cluster_rank", "bin_rel_strength", "BinClusterStrength"]],
    on=["Bin", "Cluster"],
    how="left"
)

# 4) Bin has this cluster at all?
tri_long["bin_has_ARG_cluster"] = ~tri_long["bin_cluster_rank"].isna()

# 5) Classify at the bin-level for each (ARG, Cluster, Bin) evidence row
def classify_binlevel(row):
    # missing if no cluster or no bin
    if pd.isna(row.get("Cluster")) or pd.isna(row.get("Bin")):
        return "missing"
    if not bool(row.get("bin_has_ARG_cluster")):
        return "no_support"
    r = row.get("bin_cluster_rank")
    rs = row.get("bin_rel_strength")
    if (pd.notna(r) and r <= 3) or (pd.notna(rs) and rs >= 0.2):
        return "supported"
    return "weak_support"

tri_long["triangulation_class_binlevel"] = tri_long.apply(classify_binlevel, axis=1)

# (optional) save the full long evidence table
tri_long.to_csv(
    os.path.join(out_dir, f"{SAMPLE}_ARG_triangulation_long_binsFracGE{BIN_FRAC_THRESH}.csv"),
    index=False
)

# -------------------------
# D) Collapse to per-(ARG, Cluster): if ANY candidate bin supported -> supported
# -------------------------
def collapse_classes(series):
    classes = [c for c in series.dropna().tolist() if c != ""]
    if "supported" in classes:
        return "supported"
    if "weak_support" in classes:
        return "weak_support"
    if "no_support" in classes:
        return "no_support"
    return "missing"

tri_collapsed = tri_long.groupby(["ARG_contig", "Cluster"], as_index=False).agg(
    arg_cluster_support_contignorm=("cluster_support_contig_norm", "first"),
    arg_cluster_fraction=("cluster_fraction", "first"),
    n_candidate_bins=("Bin", lambda x: x.notna().sum()),
    triangulation_class=("triangulation_class_binlevel", collapse_classes),
)

# best evidence among candidate bins (helpful diagnostics)
best_rank = tri_long.groupby(["ARG_contig", "Cluster"])["bin_cluster_rank"].min().reset_index().rename(
    columns={"bin_cluster_rank": "best_bin_rank_for_cluster"}
)
best_rel = tri_long.groupby(["ARG_contig", "Cluster"])["bin_rel_strength"].max().reset_index().rename(
    columns={"bin_rel_strength": "best_bin_rel_strength_for_cluster"}
)
tri_collapsed = tri_collapsed.merge(best_rank, on=["ARG_contig", "Cluster"], how="left")
tri_collapsed = tri_collapsed.merge(best_rel, on=["ARG_contig", "Cluster"], how="left")

# add ARG annotations (gene/subclass/location/eligibility) for interpretation
tri_collapsed = tri_collapsed.merge(amr_anno, on="ARG_contig", how="left")

# binary concordance: "supported" or "weak_support" -> concordant, others -> non-concordant
tri_collapsed["concordance"] = tri_collapsed["triangulation_class"].apply(
    lambda x: "concordant" if x in ("supported", "weak_support") else "non-concordant"
)

tri_out = os.path.join(out_dir, f"{SAMPLE}_ARG_triangulation_allclusters_binsFracGE{BIN_FRAC_THRESH}.csv")
tri_collapsed.to_csv(tri_out, index=False)

print(f"[{SAMPLE}] Wrote: {tri_out}")
print("[ALL CLUSTERS] triangulation_class counts:")
print(tri_collapsed["triangulation_class"].value_counts(dropna=False))

# -------------------------
# E) Optional: ARG-level summary (how many clusters supported)
# -------------------------
arg_summary = tri_collapsed.groupby("ARG_contig", as_index=False).agg(
    n_clusters_total=("Cluster", lambda x: x.notna().sum()),
    n_clusters_supported=("triangulation_class", lambda x: (x == "supported").sum()),
    n_clusters_weak=("triangulation_class", lambda x: (x == "weak_support").sum()),
    n_clusters_no=("triangulation_class", lambda x: (x == "no_support").sum()),
    n_clusters_missing=("triangulation_class", lambda x: (x == "missing").sum()),
)
arg_summary = arg_summary.merge(amr_anno, on="ARG_contig", how="left")

arg_out = os.path.join(out_dir, f"{SAMPLE}_ARG_triangulation_ARGlevel_summary_binsFracGE{BIN_FRAC_THRESH}.csv")
arg_summary.to_csv(arg_out, index=False)
print(f"[{SAMPLE}] Wrote: {arg_out}")

# -------------------------
# F) Plasmid-eligible-only views (does NOT delete rows; just saves filtered copies)
# -------------------------
eligible_mask = (tri_collapsed.get("eligible_for_plasmid_attachment") == True)
tri_collapsed[eligible_mask].to_csv(
    os.path.join(out_dir, f"{SAMPLE}_ARG_triangulation_allclusters_ELIGIBLEONLY_binsFracGE{BIN_FRAC_THRESH}.csv"),
    index=False
)
arg_summary[arg_summary.get("eligible_for_plasmid_attachment") == True].to_csv(
    os.path.join(out_dir, f"{SAMPLE}_ARG_triangulation_ARGlevel_summary_ELIGIBLEONLY_binsFracGE{BIN_FRAC_THRESH}.csv"),
    index=False
)

# =========================
# Extra summaries: plasmid-eligible only
# =========================
summary_fp = os.path.join(out_dir, f"{SAMPLE}_ARG_summary_topcalls.csv")
sumdf = pd.read_csv(summary_fp)

# count how many ARG contigs are eligible
eligible = sumdf[sumdf["eligible_for_plasmid_attachment"] == True].copy()
print(f"[{SAMPLE}] eligible ARG contigs for plasmid-attachment (POLICY={POLICY}): {eligible['ARG_contig'].nunique():,}")

# among eligible, how many have any ARG->Cluster mapping?
has_cluster = eligible[eligible["ARG_to_Cluster_top"].notna()]
print(f"[{SAMPLE}] eligible ARG contigs with a top cluster: {has_cluster['ARG_contig'].nunique():,}")

# distribution of triangulation_class among eligible only
if "triangulation_class" in eligible.columns:
    print("[eligible only] triangulation_class:")
    print(eligible["triangulation_class"].value_counts(dropna=False))

eligible.to_csv(os.path.join(out_dir, f"{SAMPLE}_ARG_summary_topcalls_eligible_only.csv"), index=False)

print(tri_collapsed["triangulation_class"].value_counts(dropna=False))
print(f"[{SAMPLE}] Concordance: {(tri_collapsed['concordance']=='concordant').sum()} concordant, "
      f"{(tri_collapsed['concordance']=='non-concordant').sum()} non-concordant")
print(f"[{SAMPLE}] Wrote updated triangulation outputs to: {out_dir}")


