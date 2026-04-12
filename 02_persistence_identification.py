#!/usr/bin/env python3
import os, re, csv, sys, argparse
from collections import defaultdict
import pandas as pd
import networkx as nx

FA_EXTS = (".fa", ".fna", ".fasta", ".fas")

def die(msg):
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(1)

def ensure_dir(p):
    os.makedirs(p, exist_ok=True)

def sniff_delim(path):
    with open(path, "r", newline="") as f:
        line = f.readline()
    if "\t" in line and line.count("\t") >= line.count(","):
        return "\t"
    if "," in line:
        return ","
    return "\t"

def fasta_lengths(path):
    lens = {}
    cur = None
    cur_len = 0
    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    lens[cur] = cur_len
                hdr = line[1:].strip()
                cur = hdr.split()[0]
                cur_len = 0
            else:
                cur_len += len(line.strip())
    if cur is not None:
        lens[cur] = cur_len
    return lens

def read_membership_csv(path):
    contig_to_cluster = {}
    with open(path, "r", newline="") as f:
        r = csv.DictReader(f)
        if "Node" not in r.fieldnames or "Cluster" not in r.fieldnames:
            die(f"Membership {path} missing Node/Cluster. Found: {r.fieldnames}")
        for row in r:
            node = (row.get("Node") or "").strip()
            cl = (row.get("Cluster") or "").strip()
            if not node or not cl:
                continue
            try:
                cl_int = int(float(cl))
            except:
                continue
            contig_to_cluster[node] = cl_int
    return contig_to_cluster

def find_fasta(indir, sample_id):
    cands = [fn for fn in os.listdir(indir) if fn.startswith(sample_id) and fn.endswith(FA_EXTS)]
    if not cands:
        return None
    plasmid = [x for x in cands if "plasmid" in x.lower()]
    pick = plasmid[0] if plasmid else sorted(cands)[0]
    return os.path.join(indir, pick)

def union_len(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    s, e = intervals[0]
    tot = 0
    for s2, e2 in intervals[1:]:
        if s2 <= e + 1:
            e = max(e, e2)
        else:
            tot += (e - s + 1)
            s, e = s2, e2
    tot += (e - s + 1)
    return tot

def orient_blast(hits_path, mapA, mapB, delim, max_check=20000):
    q_in_A = q_in_B = s_in_A = s_in_B = 0
    n = 0
    with open(hits_path, "r", newline="") as f:
        r = csv.reader(f, delimiter=delim)
        for row in r:
            if len(row) < 12:
                continue
            q, s = row[0].strip(), row[1].strip()
            if q in mapA: q_in_A += 1
            if q in mapB: q_in_B += 1
            if s in mapA: s_in_A += 1
            if s in mapB: s_in_B += 1
            n += 1
            if n >= max_check:
                break
    swap = (q_in_B > q_in_A) and (s_in_A > s_in_B)
    return swap

def compute_cluster_bp(map_contig_to_cluster, contig_len):
    bp = defaultdict(int)
    for c, cl in map_contig_to_cluster.items():
        L = contig_len.get(c)
        if L is None:
            continue
        bp[cl] += L
    return bp

def optionC_cluster_similarity(hits_path, mapA, mapB, lenA, lenB,
                              min_pid=95.0, min_aln=500, max_evalue=1e-10,
                              min_edge_hits=1):
    delim = sniff_delim(hits_path)
    swap = orient_blast(hits_path, mapA, mapB, delim)

    A_intervals = defaultdict(list)  # (Acontig, clB) -> intervals on A contig coords
    B_intervals = defaultdict(list)  # (Bcontig, clA) -> intervals on B contig coords
    edge_hits = defaultdict(int)
    edge_maxpid = defaultdict(float)

    def parse(row):
        try:
            q = row[0].strip(); s = row[1].strip()
            pid = float(row[2]); aln = int(float(row[3]))
            qs = int(row[6]); qe = int(row[7])
            ss = int(row[8]); se = int(row[9])
            ev = float(row[10])
        except:
            return None
        if pid < min_pid or aln < min_aln or ev > max_evalue:
            return None
        qs, qe = (qs, qe) if qs <= qe else (qe, qs)
        ss, se = (ss, se) if ss <= se else (se, ss)
        return q, s, pid, qs, qe, ss, se

    with open(hits_path, "r", newline="") as f:
        r = csv.reader(f, delimiter=delim)
        for row in r:
            if len(row) < 12:
                continue
            pr = parse(row)
            if pr is None:
                continue
            q, s, pid, qs, qe, ss, se = pr

            if not swap:
                clA = mapA.get(q); clB = mapB.get(s)
                if clA is None or clB is None:
                    continue
                A_intervals[(q, clB)].append((qs, qe))
                B_intervals[(s, clA)].append((ss, se))
            else:
                clA = mapA.get(s); clB = mapB.get(q)
                if clA is None or clB is None:
                    continue
                A_intervals[(s, clB)].append((ss, se))
                B_intervals[(q, clA)].append((qs, qe))

            edge_hits[(clA, clB)] += 1
            if pid > edge_maxpid[(clA, clB)]:
                edge_maxpid[(clA, clB)] = pid

    cluster_bp_A = compute_cluster_bp(mapA, lenA)
    cluster_bp_B = compute_cluster_bp(mapB, lenB)

    aligned_A_to_B = defaultdict(int)  # (clA, clB)
    aligned_B_to_A = defaultdict(int)  # (clA, clB)

    for (contigA, clB), ivals in A_intervals.items():
        clA = mapA.get(contigA)
        if clA is None:
            continue
        aligned_A_to_B[(clA, clB)] += union_len(ivals)

    for (contigB, clA), ivals in B_intervals.items():
        clB = mapB.get(contigB)
        if clB is None:
            continue
        aligned_B_to_A[(clA, clB)] += union_len(ivals)

    rows = []
    for (clA, clB), hits in edge_hits.items():
        if hits < min_edge_hits:
            continue
        totA = cluster_bp_A.get(clA, 0)
        totB = cluster_bp_B.get(clB, 0)
        alnA = aligned_A_to_B.get((clA, clB), 0)
        alnB = aligned_B_to_A.get((clA, clB), 0)
        covA = (alnA / totA) if totA > 0 else 0.0
        covB = (alnB / totB) if totB > 0 else 0.0
        rows.append({
            "cluster_A": clA,
            "cluster_B": clB,
            "hits": hits,
            "max_pid": edge_maxpid.get((clA, clB), 0.0),
            "aligned_bp_A_to_B": alnA,
            "total_bp_cluster_A": totA,
            "cov_A_to_B": covA,
            "aligned_bp_B_to_A": alnB,
            "total_bp_cluster_B": totB,
            "cov_B_to_A": covB,
            "sym_cov_min": min(covA, covB)
        })
    return pd.DataFrame(rows)

def node_label(etype, clA, clB):
    if etype == "Inf-AS":
        return f"Inf_C{clA}", f"AS_C{clB}"
    if etype == "AS-Eff":
        return f"AS_C{clA}", f"Eff_C{clB}"
    if etype == "Inf-Eff":
        return f"Inf_C{clA}", f"Eff_C{clB}"
    raise ValueError(etype)

def make_inf_eff_pairs(df_inf_eff):
    if df_inf_eff is None or df_inf_eff.empty:
        return pd.DataFrame()
    out = pd.DataFrame({
        "Inf_node": df_inf_eff["cluster_A"].astype(int).map(lambda x: f"Inf_C{x}"),
        "Eff_node": df_inf_eff["cluster_B"].astype(int).map(lambda x: f"Eff_C{x}"),
        "sym_cov_min": df_inf_eff["sym_cov_min"],
        "cov_Inf_to_Eff": df_inf_eff["cov_A_to_B"],
        "cov_Eff_to_Inf": df_inf_eff["cov_B_to_A"],
        "hits": df_inf_eff["hits"],
        "max_pid": df_inf_eff["max_pid"],
        "total_bp_Inf_cluster": df_inf_eff["total_bp_cluster_A"],
        "total_bp_Eff_cluster": df_inf_eff["total_bp_cluster_B"],
    }).sort_values("sym_cov_min", ascending=False)
    return out

def step2_persistence(edges, min_sym=0.2, min_hits=1, min_cluster_bp=0):
    filtered = {}
    for etype, df in edges.items():
        if df is None or df.empty:
            filtered[etype] = pd.DataFrame()
            continue
        keep = df[(df["sym_cov_min"] >= min_sym) & (df["hits"] >= min_hits)].copy()
        if min_cluster_bp > 0:
            keep = keep[(keep["total_bp_cluster_A"] >= min_cluster_bp) & (keep["total_bp_cluster_B"] >= min_cluster_bp)]
        filtered[etype] = keep

    # strict Inf-AS-Eff triples (bridge via AS cluster id)
    triples = pd.DataFrame()
    if not filtered["Inf-AS"].empty and not filtered["AS-Eff"].empty:
        a = filtered["Inf-AS"].rename(columns={"cluster_A":"inf_cluster","cluster_B":"as_cluster"})
        b = filtered["AS-Eff"].rename(columns={"cluster_A":"as_cluster","cluster_B":"eff_cluster"})
        m = a.merge(b, on="as_cluster", suffixes=("_inf_as","_as_eff"))
        if not m.empty:
            triples = pd.DataFrame({
                "Inf_node": m["inf_cluster"].apply(lambda x: f"Inf_C{int(x)}"),
                "AS_node":  m["as_cluster"].apply(lambda x: f"AS_C{int(x)}"),
                "Eff_node": m["eff_cluster"].apply(lambda x: f"Eff_C{int(x)}"),
                "sym_inf_as": m["sym_cov_min_inf_as"],
                "sym_as_eff": m["sym_cov_min_as_eff"],
                "sym_min_chain": m[["sym_cov_min_inf_as","sym_cov_min_as_eff"]].min(axis=1),
                "hits_inf_as": m["hits_inf_as"],
                "hits_as_eff": m["hits_as_eff"],
                "maxpid_inf_as": m["max_pid_inf_as"],
                "maxpid_as_eff": m["max_pid_as_eff"],
            }).sort_values("sym_min_chain", ascending=False)

    # connected components = lineages
    G = nx.Graph()
    edge_records = []
    for etype, df in filtered.items():
        if df.empty:
            continue
        for _, r in df.iterrows():
            u, v = node_label(etype, int(r["cluster_A"]), int(r["cluster_B"]))
            w = float(r["sym_cov_min"])
            G.add_edge(u, v, weight=w, edge_type=etype)
            edge_records.append((u, v, etype))

    if G.number_of_nodes() == 0:
        return filtered, triples, pd.DataFrame(), pd.DataFrame()

    comps = list(nx.connected_components(G))
    node_to_lineage = {}
    for i, nodeset in enumerate(comps, start=1):
        for n in nodeset:
            node_to_lineage[n] = i

    def compartment(n):
        if n.startswith("Inf_"): return "Inf"
        if n.startswith("AS_"): return "AS"
        if n.startswith("Eff_"): return "Eff"
        return "UNK"

    lineage_edge_types = defaultdict(set)
    for u, v, etype in edge_records:
        lid = node_to_lineage[u]
        lineage_edge_types[lid].add(etype)

    node_rows = []
    sum_rows = []
    for lid, nodeset in enumerate(comps, start=1):
        comps_present = set(compartment(n) for n in nodeset)
        etypes = lineage_edge_types.get(lid, set())
        has_inf_as = "Inf-AS" in etypes
        has_as_eff = "AS-Eff" in etypes
        has_inf_eff = "Inf-Eff" in etypes

        if ("Inf" in comps_present) and ("AS" in comps_present) and ("Eff" in comps_present) and has_inf_as and has_as_eff:
            cls = "Inf-AS-Eff"
        elif ("Inf" in comps_present) and ("Eff" in comps_present) and has_inf_eff:
            cls = "Inf-Eff"
        elif ("Inf" in comps_present) and ("AS" in comps_present) and has_inf_as:
            cls = "Inf-AS"
        elif ("AS" in comps_present) and ("Eff" in comps_present) and has_as_eff:
            cls = "AS-Eff"
        else:
            cls = "single_or_ambiguous"

        for n in sorted(nodeset):
            node_rows.append({"lineage_id": lid, "node": n, "compartment": compartment(n), "persistence_class": cls})

        sum_rows.append({
            "lineage_id": lid,
            "persistence_class": cls,
            "compartments": ";".join(sorted(comps_present)),
            "n_nodes": len(nodeset),
            "has_InfAS_edge": has_inf_as,
            "has_ASEff_edge": has_as_eff,
            "has_InfEff_edge": has_inf_eff
        })

    return filtered, triples, pd.DataFrame(node_rows), pd.DataFrame(sum_rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, help="Folder containing BLAST + membership + FASTA")
    ap.add_argument("--outdir", required=True, help="Output folder")

    # Step1 thresholds
    ap.add_argument("--min-pid", type=float, default=95.0)
    ap.add_argument("--min-aln", type=int, default=500)
    ap.add_argument("--max-evalue", type=float, default=1e-10)
    ap.add_argument("--min-edge-hits", type=int, default=1)

    # Step2 thresholds
    ap.add_argument("--min-sym", type=float, default=0.2, help="Threshold on sym_cov_min for cluster-edge acceptance")
    ap.add_argument("--min-hits", type=int, default=1)
    ap.add_argument("--min-cluster-bp", type=int, default=0, help="Require both clusters >= this bp in Step2 filtering")

    # Group detection
    ap.add_argument("--group-regex", default=r"^inf_vs_AS_(.+)\.csv$", help="Regex to extract group id from inf_vs_AS filenames")
    args = ap.parse_args()

    indir = args.indir
    outdir = args.outdir
    ensure_dir(outdir)

    rgx = re.compile(args.group_regex)
    groups = sorted({rgx.match(fn).group(1) for fn in os.listdir(indir) if rgx.match(fn)})
    if not groups:
        die(f"No groups found in {indir} with regex {args.group_regex}")

    all_counts = []
    all_triples = []
    all_inf_eff_pairs = []

    for g in groups:
        inf_id = f"{g}1"
        as_id  = f"{g}2"
        eff_id = f"{g}3"

        mem_inf = os.path.join(indir, f"{inf_id}_cluster_membership_filtered.csv")
        mem_as  = os.path.join(indir, f"{as_id}_cluster_membership_filtered.csv")
        mem_eff = os.path.join(indir, f"{eff_id}_cluster_membership_filtered.csv")
        if not (os.path.exists(mem_inf) and os.path.exists(mem_as) and os.path.exists(mem_eff)):
            print(f"[WARN] {g}: missing membership. Skipping.")
            continue

        fa_inf = find_fasta(indir, inf_id)
        fa_as  = find_fasta(indir, as_id)
        fa_eff = find_fasta(indir, eff_id)
        if not (fa_inf and fa_as and fa_eff):
            print(f"[WARN] {g}: missing FASTA for {inf_id}/{as_id}/{eff_id}. Skipping.")
            continue

        blast_inf_as = os.path.join(indir, f"inf_vs_AS_{g}.csv")
        blast_as_eff = os.path.join(indir, f"AS_vs_eff_{g}.csv")
        blast_inf_eff = os.path.join(indir, f"inf_vs_eff_{g}.csv")  # direct Inf-Eff

        if not (os.path.exists(blast_inf_as) and os.path.exists(blast_as_eff)):
            print(f"[WARN] {g}: missing BLAST inf-vs-AS or AS-vs-eff. Skipping.")
            continue

        gout = os.path.join(outdir, g)
        ensure_dir(gout)

        map_inf = read_membership_csv(mem_inf)
        map_as  = read_membership_csv(mem_as)
        map_eff = read_membership_csv(mem_eff)
        len_inf = fasta_lengths(fa_inf)
        len_as  = fasta_lengths(fa_as)
        len_eff = fasta_lengths(fa_eff)

        # Step1
        df_inf_as = optionC_cluster_similarity(
            blast_inf_as, map_inf, map_as, len_inf, len_as,
            min_pid=args.min_pid, min_aln=args.min_aln, max_evalue=args.max_evalue,
            min_edge_hits=args.min_edge_hits
        )
        df_as_eff = optionC_cluster_similarity(
            blast_as_eff, map_as, map_eff, len_as, len_eff,
            min_pid=args.min_pid, min_aln=args.min_aln, max_evalue=args.max_evalue,
            min_edge_hits=args.min_edge_hits
        )
        df_inf_eff = None
        if os.path.exists(blast_inf_eff):
            df_inf_eff = optionC_cluster_similarity(
                blast_inf_eff, map_inf, map_eff, len_inf, len_eff,
                min_pid=args.min_pid, min_aln=args.min_aln, max_evalue=args.max_evalue,
                min_edge_hits=args.min_edge_hits
            )
        else:
            print(f"[INFO] {g}: no direct {os.path.basename(blast_inf_eff)} found; Inf-Eff will be absent unless you create this BLAST.")

        df_inf_as.to_csv(os.path.join(gout, "step1_cluster_similarity_inf_as.tsv"), sep="\t", index=False)
        df_as_eff.to_csv(os.path.join(gout, "step1_cluster_similarity_as_eff.tsv"), sep="\t", index=False)
        if df_inf_eff is not None:
            df_inf_eff.to_csv(os.path.join(gout, "step1_cluster_similarity_inf_eff.tsv"), sep="\t", index=False)

        # Step2 filtering + lineages
        edges = {
            "Inf-AS": df_inf_as,
            "AS-Eff": df_as_eff,
            "Inf-Eff": df_inf_eff if df_inf_eff is not None else pd.DataFrame()
        }
        filtered, triples, node_tbl, sum_tbl = step2_persistence(
            edges,
            min_sym=args.min_sym, min_hits=args.min_hits, min_cluster_bp=args.min_cluster_bp
        )

        filtered["Inf-AS"].to_csv(os.path.join(gout, "step2_edges_inf_as.filtered.tsv"), sep="\t", index=False)
        filtered["AS-Eff"].to_csv(os.path.join(gout, "step2_edges_as_eff.filtered.tsv"), sep="\t", index=False)
        if df_inf_eff is not None:
            filtered["Inf-Eff"].to_csv(os.path.join(gout, "step2_edges_inf_eff.filtered.tsv"), sep="\t", index=False)

        triples.to_csv(os.path.join(gout, "step2_strict_triples_inf_as_eff.tsv"), sep="\t", index=False)
        node_tbl.to_csv(os.path.join(gout, "step2_lineage_nodes.tsv"), sep="\t", index=False)
        sum_tbl.to_csv(os.path.join(gout, "step2_lineage_summary.tsv"), sep="\t", index=False)

        # ALSO: write Inf-Eff pair list (direct Inf-Eff only)
        if df_inf_eff is not None:
            inf_eff_pairs = make_inf_eff_pairs(filtered["Inf-Eff"])
            inf_eff_pairs.to_csv(os.path.join(gout, "step2_pairs_inf_eff.tsv"), sep="\t", index=False)
            if not inf_eff_pairs.empty:
                tmp = inf_eff_pairs.copy()
                tmp.insert(0, "group", g)
                all_inf_eff_pairs.append(tmp)

        # group counts
        if not sum_tbl.empty:
            c = sum_tbl["persistence_class"].value_counts().rename_axis("persistence_class").reset_index(name="n_lineages")
        else:
            c = pd.DataFrame({"persistence_class": [], "n_lineages": []})
        c.insert(0, "group", g)
        all_counts.append(c)

        if not triples.empty:
            tt = triples.copy()
            tt.insert(0, "group", g)
            all_triples.append(tt)

        print(f"[INFO] {g}: wrote Step1+Step2 to {gout}")

    # combined summaries
    pd.concat(all_counts, ignore_index=True).to_csv(os.path.join(outdir, "ALL_groups_lineage_counts.tsv"), sep="\t", index=False) if all_counts else pd.DataFrame().to_csv(os.path.join(outdir, "ALL_groups_lineage_counts.tsv"), sep="\t", index=False)
    pd.concat(all_triples, ignore_index=True).to_csv(os.path.join(outdir, "ALL_groups_strict_triples.tsv"), sep="\t", index=False) if all_triples else pd.DataFrame().to_csv(os.path.join(outdir, "ALL_groups_strict_triples.tsv"), sep="\t", index=False)
    pd.concat(all_inf_eff_pairs, ignore_index=True).to_csv(os.path.join(outdir, "ALL_groups_inf_eff_pairs.tsv"), sep="\t", index=False) if all_inf_eff_pairs else pd.DataFrame().to_csv(os.path.join(outdir, "ALL_groups_inf_eff_pairs.tsv"), sep="\t", index=False)

    print(f"[DONE] Pipeline complete. Outputs in: {outdir}")

if __name__ == "__main__":
    main()

