#!/usr/bin/env python3
import os, re, argparse
import pandas as pd
from collections import defaultdict

# Mobility ranking if we need a "best" mobility across nodes in a lineage
MOB_ORDER = {"no_marker": 1, "mobilizable": 2, "mobilizable/conjugative": 3}

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def parse_week_num(x):
    # "Week 1" -> 1
    if pd.isna(x):
        return None
    m = re.search(r"(\d+)", str(x))
    return int(m.group(1)) if m else None

def site_to_idx(site):
    s = str(site).lower()
    if "influent" in s:
        return "1"
    if "activated" in s or "sludge" in s:
        return "2"
    if "effluent" in s:
        return "3"
    return None

def city_to_letter(city):
    """Map facility name to single-letter code.
    Edit this mapping to match your dataset."""
    CITY_MAP = {}  # e.g., {"facility_1": "A", "facility_2": "B", "facility_3": "C"}
    c = str(city).lower().strip()
    if c in CITY_MAP:
        return CITY_MAP[c]
    return str(city).strip()[0].upper() if str(city).strip() else None

def taxonomy_level(tax, level="family"):
    # expects GTDB-like: ...;f__X;g__Y;s__Z
    if pd.isna(tax):
        return None
    parts = str(tax).split(";")
    pref = {"family": "f__", "genus": "g__", "species": "s__"}[level]
    for p in parts:
        p = p.strip()
        if p.startswith(pref):
            val = p.replace(pref, "")
            return val if val else None
    return None

def best_mobility(mobs):
    mobs = [m for m in mobs if pd.notna(m)]
    if not mobs:
        return None
    return max(mobs, key=lambda m: MOB_ORDER.get(m, 0))

def node_to_comp_cluster(node):
    # "Inf_C22" -> ("Inf", 22)
    m = re.match(r"^(Inf|AS|Eff)_C(\d+)$", str(node))
    if not m:
        return None, None
    return m.group(1), int(m.group(2))

def comp_to_sample(group, comp):
    # group "X1" => Inf=X11, AS=X12, Eff=X13
    return f"{group}{ {'Inf':'1','AS':'2','Eff':'3'}[comp] }"

def clean_numeric_cluster(df, col="Cluster", drop_no_cluster_assigned=True, label="table"):
    """
    Robustly coerce Cluster column to numeric int, dropping non-numeric entries.
    Handles values like 'No Cluster Assigned', '', '33.0', etc.
    """
    if col not in df.columns:
        raise SystemExit(f"[ERROR] {label}: missing required column '{col}'")

    df = df.copy()
    df[f"{col}_raw"] = df[col].astype(str)

    if drop_no_cluster_assigned:
        df = df[df[f"{col}_raw"].str.strip().str.lower().ne("no cluster assigned")].copy()

    df[col] = pd.to_numeric(df[f"{col}_raw"], errors="coerce")
    before = len(df)
    df = df.dropna(subset=[col]).copy()
    after = len(df)
    df[col] = df[col].astype(int)

    dropped = before - after
    if dropped > 0:
        print(f"[INFO] {label}: dropped {dropped:,} rows with non-numeric Cluster values")
    return df

def main():
    ap = argparse.ArgumentParser(description="Annotate persistence lineages with mobility + host taxa continuity.")
    ap.add_argument("--persistence_dir", required=True,
                    help="Folder with group subfolders containing step2_lineage_nodes.tsv and step2_lineage_summary.tsv")
    ap.add_argument("--mobility_categories", required=True,
                    help="CSV with Sample, Cluster, cluster_mobility (sample×cluster mobility categories)")
    ap.add_argument("--host_table", required=True,
                    help="Bin-Cluster host association table with columns: Sample (or City/Week/Site), Cluster, Bin, Taxonomy, and host strength column")
    ap.add_argument("--host_strength_col", default="Log_Normalized_HiC_Contacts",
                    help="Column used as host-link strength")
    ap.add_argument("--host_level", default="family", choices=["family", "genus", "species"],
                    help="Taxonomic level for host transfer (family/genus/species)")
    ap.add_argument("--host_mode", default="frac", choices=["frac", "topk"],
                    help="How to define host set per cluster: 'frac' keeps taxa with frac>=min_host_frac; "
                         "'topk' keeps top K taxa")
    ap.add_argument("--min_host_frac", type=float, default=0.0,
                    help="Host fraction threshold if host_mode=frac (default: 0.0, keep all nonzero hosts)")
    ap.add_argument("--topk_hosts", type=int, default=3,
                    help="Top K hosts if host_mode=topk")
    ap.add_argument("--outdir", required=True,
                    help="Output folder for annotated tables")
    args = ap.parse_args()

    ensure_dir(args.outdir)

    # -------------------------
    # Mobility categories
    # -------------------------
    mob = pd.read_csv(args.mobility_categories)
    need_m = {"Sample", "Cluster", "cluster_mobility"}
    miss = need_m - set(mob.columns)
    if miss:
        raise SystemExit(f"[ERROR] mobility_categories missing columns: {sorted(miss)}")

    mob = clean_numeric_cluster(mob, col="Cluster", drop_no_cluster_assigned=True, label="mobility_categories")

    mob_key = mob[["Sample", "Cluster", "cluster_mobility"]].drop_duplicates()

    # -------------------------
    # Host table
    # -------------------------
    host = pd.read_csv(args.host_table)

    # If Sample is missing, reconstruct from City/Week/Site
    if "Sample" not in host.columns:
        need_h = {"City", "Week", "Site"}
        missh = need_h - set(host.columns)
        if missh:
            raise SystemExit(f"[ERROR] host_table missing Sample and also missing {sorted(missh)} to reconstruct Sample.")
        wk = host["Week"].apply(parse_week_num)
        ct = host["City"].apply(city_to_letter)
        st = host["Site"].apply(site_to_idx)

        if wk.isna().any():
            print("[WARN] Some Week values could not be parsed; these rows will likely drop later.")
        if ct.isna().any():
            print("[WARN] Some City values could not be parsed; these rows will likely drop later.")
        if st.isna().any():
            print("[WARN] Some Site values could not be parsed; these rows will likely drop later.")

        host["Sample"] = ct.astype(str) + wk.astype("Int64").astype(str) + st.astype(str)

    need_h2 = {"Sample", "Cluster", "Bin", "Taxonomy", args.host_strength_col}
    miss2 = need_h2 - set(host.columns)
    if miss2:
        raise SystemExit(f"[ERROR] host_table missing columns: {sorted(miss2)}")

    host = clean_numeric_cluster(host, col="Cluster", drop_no_cluster_assigned=True, label="host_table")

    # Parse taxonomy to desired level
    host["host_taxon"] = host["Taxonomy"].apply(lambda x: taxonomy_level(x, args.host_level))
    host = host.dropna(subset=["host_taxon"]).copy()

    # Coerce host strength numeric
    host["host_strength"] = pd.to_numeric(host[args.host_strength_col], errors="coerce")
    host = host.dropna(subset=["host_strength"]).copy()

    # Aggregate host strength at (Sample, Cluster, host_taxon)
    host_agg = (host
                .groupby(["Sample", "Cluster", "host_taxon"], as_index=False)["host_strength"]
                .sum())

    # Build host set per (Sample, Cluster)
    host_sets = {}
    host_top1 = {}
    host_list_str = {}

    for (s, c), sub in host_agg.groupby(["Sample", "Cluster"]):
        sub = sub.sort_values("host_strength", ascending=False).copy()
        tot = sub["host_strength"].sum()
        if tot <= 0:
            continue

        if args.host_mode == "topk":
            keep = sub.head(args.topk_hosts)
        else:
            sub["frac"] = sub["host_strength"] / tot
            keep = sub[sub["frac"] >= args.min_host_frac]
            if keep.empty:
                keep = sub.head(1)

        hs = list(keep["host_taxon"].astype(str))
        host_sets[(s, int(c))] = set(hs)
        host_top1[(s, int(c))] = str(sub.iloc[0]["host_taxon"])
        host_list_str[(s, int(c))] = ";".join(hs)

    # -------------------------
    # Iterate group folders
    # -------------------------
    groups = sorted([d for d in os.listdir(args.persistence_dir)
                     if os.path.isdir(os.path.join(args.persistence_dir, d))])

    all_nodes = []
    all_lineage = []

    for g in groups:
        gdir = os.path.join(args.persistence_dir, g)
        fn_nodes = os.path.join(gdir, "step2_lineage_nodes.tsv")
        fn_sum = os.path.join(gdir, "step2_lineage_summary.tsv")
        if not (os.path.exists(fn_nodes) and os.path.exists(fn_sum)):
            continue

        nodes = pd.read_csv(fn_nodes, sep="\t")
        summ = pd.read_csv(fn_sum, sep="\t")

        nodes["group"] = g
        summ["group"] = g

        # parse node -> compartment + cluster_id
        comps, cls = zip(*nodes["node"].apply(node_to_comp_cluster))
        nodes["compartment"] = comps
        nodes["cluster_id"] = cls

        # construct Sample and numeric Cluster
        nodes["Sample"] = nodes["compartment"].apply(lambda comp: comp_to_sample(g, comp))
        nodes["Cluster"] = nodes["cluster_id"].astype(int)

        # join mobility
        nodes = nodes.merge(mob_key, on=["Sample", "Cluster"], how="left")

        # attach host annotations
        nodes["host_top1"] = nodes.apply(lambda r: host_top1.get((r["Sample"], int(r["Cluster"]))), axis=1)
        nodes["host_list"] = nodes.apply(lambda r: host_list_str.get((r["Sample"], int(r["Cluster"]))), axis=1)

        all_nodes.append(nodes)

        # lineage-level rollups
        for lid, sub in nodes.groupby("lineage_id"):
            hs = {}
            top = {}
            mob_comp = {}

            for comp, subc in sub.groupby("compartment"):
                s = subc.iloc[0]["Sample"]
                c = int(subc.iloc[0]["Cluster"])
                hs[comp] = host_sets.get((s, c), set())
                top[comp] = host_top1.get((s, c))
                mob_comp[comp] = subc.iloc[0]["cluster_mobility"]

            lineage_mob = best_mobility(list(mob_comp.values()))

            def shared(a, b):
                if a not in hs or b not in hs:
                    return None
                return len(hs[a].intersection(hs[b])) > 0

            row = {
                "group": g,
                "lineage_id": lid,
                "persistence_class": sub.iloc[0].get("persistence_class", None),
                "lineage_mobility_best": lineage_mob,
                "mobility_Inf": mob_comp.get("Inf"),
                "mobility_AS": mob_comp.get("AS"),
                "mobility_Eff": mob_comp.get("Eff"),
                "topHost_Inf": top.get("Inf"),
                "topHost_AS": top.get("AS"),
                "topHost_Eff": top.get("Eff"),
                "sharedHost_Inf_AS": shared("Inf", "AS"),
                "sharedHost_AS_Eff": shared("AS", "Eff"),
                "sharedHost_Inf_Eff": shared("Inf", "Eff"),
                "hostList_Inf": ";".join(sorted(hs.get("Inf", set()))) if "Inf" in hs else None,
                "hostList_AS": ";".join(sorted(hs.get("AS", set()))) if "AS" in hs else None,
                "hostList_Eff": ";".join(sorted(hs.get("Eff", set()))) if "Eff" in hs else None,
                "dominantShift_Inf_AS": (top.get("Inf") != top.get("AS")) if (top.get("Inf") and top.get("AS")) else None,
                "dominantShift_AS_Eff": (top.get("AS") != top.get("Eff")) if (top.get("AS") and top.get("Eff")) else None,
                "dominantShift_Inf_Eff": (top.get("Inf") != top.get("Eff")) if (top.get("Inf") and top.get("Eff")) else None,
            }
            all_lineage.append(row)

    nodes_all = pd.concat(all_nodes, ignore_index=True) if all_nodes else pd.DataFrame()
    lin_all = pd.DataFrame(all_lineage)

    # -------------------------
    # Write outputs
    # -------------------------
    nodes_all.to_csv(os.path.join(args.outdir, "lineage_nodes_annotated.tsv"), sep="\t", index=False)
    lin_all.to_csv(os.path.join(args.outdir, "lineage_summary_mobility_hosts.tsv"), sep="\t", index=False)

    # Stratified counts: persistence class × mobility
    if not lin_all.empty:
        c1 = (lin_all
              .groupby(["persistence_class", "lineage_mobility_best"], as_index=False)
              .size()
              .rename(columns={"size": "n_lineages"}))
        c1.to_csv(os.path.join(args.outdir, "counts_persistence_by_mobility.tsv"), sep="\t", index=False)

        # Host continuity by mobility
        def mean_bool(s):
            s = s.dropna()
            return float(s.mean()) if len(s) else float("nan")

        c2 = (lin_all
              .groupby(["lineage_mobility_best"], as_index=False)
              .agg(
                  n_lineages=("lineage_id", "count"),
                  frac_shared_Inf_AS=("sharedHost_Inf_AS", mean_bool),
                  frac_shared_AS_Eff=("sharedHost_AS_Eff", mean_bool),
                  frac_shared_Inf_Eff=("sharedHost_Inf_Eff", mean_bool),
              ))
        c2.to_csv(os.path.join(args.outdir, "host_continuity_by_mobility.tsv"), sep="\t", index=False)

    print(f"[DONE] Wrote outputs to {args.outdir}")

if __name__ == "__main__":
    main()

