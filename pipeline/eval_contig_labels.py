#!/usr/bin/env python3
"""Build per-contig subgenome truth from a contigs->reference PAF + a chrom->subgenome map,
then score reference-free contig labels by the best one-to-one (label->subgenome) assignment.
Also writes the 6-column per-contig truth TSV (contig,length,Lbp,Rbp,minor_frac,class) reusable
as TRUTH_PER_CONTIG.  Usage:
  eval_contig_labels.py CHROM_SUBG.tsv contigs_to_ref.paf assembly.fasta.fai pred_labels.tsv truth_out.tsv"""
import sys
from collections import defaultdict, Counter
from itertools import permutations

CHROM_SUBG, PAF, FAI, PRED, TRUTH_OUT = sys.argv[1:6]
NONLAB = {"", "ambiguous", "unassigned", "?", "chimeric", "label", "NA", "unmapped"}

sg = {}
for ln in open(CHROM_SUBG):
    f = ln.split()
    if len(f) >= 2: sg[f[0]] = f[1]
truth_labels = sorted(set(sg.values()))

qlen = {}
for ln in open(FAI):
    p = ln.split('\t'); qlen[p[0]] = int(p[1])

bp = defaultdict(lambda: defaultdict(int))
for ln in open(PAF):
    p = ln.split('\t')
    if len(p) < 10: continue
    s = sg.get(p[5])
    if s:
        try: bp[p[0]][s] += int(p[9])
        except ValueError: pass

truth = {}
with open(TRUTH_OUT, 'w') as o:
    o.write("contig\tlength\t%s\t%s\tminor_frac\tclass\n" % (truth_labels[0], truth_labels[1]))
    for q in qlen:
        d = bp.get(q, {})
        tot = sum(d.values())
        a = d.get(truth_labels[0], 0); b = d.get(truth_labels[1], 0)
        if tot == 0:
            kls, mf = "unmapped", 0.0
        else:
            maj = max(d, key=d.get); mf = (tot - d[maj]) / tot
            kls = "pure_%s" % maj if mf < 0.05 else "mixed"
        o.write(f"{q}\t{qlen[q]}\t{a}\t{b}\t{mf:.3f}\t{kls}\n")
        if kls.startswith("pure_"): truth[q] = kls.split("_", 1)[1]

pred = {}
with open(PRED) as f:
    next(f)
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) >= 2 and p[1] not in NONLAB: pred[p[0]] = p[1]
pred_labels = sorted(set(pred.values()))
common = [c for c in pred if c in truth]

def acc_bp(mp):
    ok = tot = 0
    for c in common:
        L = qlen[c]; tot += L
        if mp.get(pred[c]) == truth[c]: ok += L
    return ok, tot

best = max((dict(zip(pred_labels, perm)) for perm in permutations(truth_labels, len(pred_labels))),
           key=lambda m: acc_bp(m)[0])
ok, tot = acc_bp(best)
npure = sum(1 for _ in truth)
print(f"[truth] {npure:,} pure contigs ; {len(qlen)-npure:,} mixed/unmapped (excluded)")
print(f"[map ] best 1:1  {best}")
print(f"[ACC ] contig-bp accuracy {100*ok/tot:.2f}%  ({ok/1e6:.1f}/{tot/1e6:.1f} Mb over {len(common)} contigs)")
conf = Counter()
for c in common:
    conf[(truth[c], best.get(pred[c]))] += qlen[c]
print("confusion (true -> predicted, % of true-subgenome bp):")
for t in truth_labels:
    row = {p: conf.get((t, p), 0) for p in truth_labels}; tt = sum(row.values()) or 1
    print(f"  true {t}: " + "  ".join(f"->{p} {100*row[p]/tt:5.1f}%" for p in truth_labels) + f"   ({tt/1e6:.1f} Mb)")
