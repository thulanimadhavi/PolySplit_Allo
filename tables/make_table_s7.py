#!/usr/bin/env python3
"""Generate Supplementary Table S7 (contig misplacement, per dataset and read chemistry).

For each dataset x chemistry it reads three files produced by the PolySplit-Allo run:
  <run_dir>/contigs_to_ref.paf             contigs aligned to the chromosome-anchored reference
  <run_dir>/all_contig_labels_repaired.tsv PolySplit final contig -> subgenome labels
  <chrom_map>                              reference chromosome -> subgenome (TSV: chrom<TAB>subg)

Each contig's TRUE subgenome is its majority under positional best-identity assignment (per query
base, the subgenome of the highest-identity alignment covering it). A contig is MISPLACED when the
subgenome PolySplit assigned it (cluster mapped to truth by best one-to-one) differs from that
majority. Contigs with no anchored alignment are excluded. Prints the Table S7 markdown.

Usage:  make_table_s7.py datasets.tsv
        datasets.tsv columns (tab-separated, one row per dataset x chemistry, header line ignored):
        name<TAB>reads<TAB>run_dir<TAB>chrom_map
"""
import sys, collections, numpy as np

def load_map(fn):
    d = {}
    for ln in open(fn):
        p = ln.rstrip("\n").split("\t")
        if len(p) >= 2:
            d[p[0]] = p[1]
    return d

def analyze(run_dir, chrom_map):
    m = load_map(chrom_map)
    subs = sorted(set(m.values())); si = {s: i for i, s in enumerate(subs)}; K = len(subs)
    aln = collections.defaultdict(list); qlen = {}
    for ln in open(run_dir + "/contigs_to_ref.paf"):
        p = ln.split("\t")
        if len(p) < 12 or p[5] not in m:
            continue
        qlen[p[0]] = int(p[1])
        aln[p[0]].append((int(p[2]), int(p[3]), si[m[p[5]]], int(p[9]) / max(1, int(p[10]))))
    comp = {}
    for c, lst in aln.items():
        L = qlen[c]; sc = np.zeros(L, np.float32); sv = np.full(L, -1, np.int8)
        for qs, qe, s, dens in lst:
            qs = max(0, qs); qe = min(L, qe); msk = dens > sc[qs:qe]
            sub = sv[qs:qe]; sub[msk] = s; sc[qs:qe][msk] = dens
        comp[c] = [int((sv == k).sum()) for k in range(K)]
    lab = {}
    for ln in open(run_dir + "/all_contig_labels_repaired.tsv"):
        if ln.startswith("contig\t"):
            continue
        p = ln.rstrip("\n").split("\t"); lab[p[0]] = p[1]
    agg = collections.defaultdict(lambda: [0] * K)
    for c, L in lab.items():
        if c in comp:
            for k in range(K):
                agg[L][k] += comp[c][k]
    cl2t = {L: int(np.argmax(v)) for L, v in agg.items()}
    nmap = nmis = misbp = totbp = nlarge = largebp = 0
    for c, L in lab.items():
        if c not in comp:
            continue
        bp = comp[c]; tot = sum(bp)
        if tot == 0:
            continue
        nmap += 1; totbp += tot; maj = int(np.argmax(bp)); a = cl2t.get(L, maj)
        if a != maj:
            nmis += 1; misbp += tot
            if qlen[c] > 100000:
                nlarge += 1; largebp += tot
    return nmap, nmis, misbp, totbp, nlarge, largebp

def main(cfg):
    rows = []
    for ln in open(cfg):
        p = ln.rstrip("\n").split("\t")
        if len(p) < 4 or p[0].lower() in ("name", "#name"):
            continue
        rows.append(p[:4])
    print("| Genome | reads | contigs (anchored) | misplaced | misplaced (Mb) | placement bp acc. % | large (>100 kb) misplaced |")
    print("|---|---|---|---|---|---|---|")
    for name, reads, run_dir, cmap in rows:
        nmap, nmis, misbp, totbp, nlarge, largebp = analyze(run_dir, cmap)
        pur = 100 * (1 - misbp / max(1, totbp))
        print("| %s | %s | %d | %d | %.2f | %.2f | %d (%.2f Mb) |" % (
            name, reads, nmap, nmis, misbp / 1e6, pur, nlarge, largebp / 1e6))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: make_table_s7.py datasets.tsv")
    main(sys.argv[1])
