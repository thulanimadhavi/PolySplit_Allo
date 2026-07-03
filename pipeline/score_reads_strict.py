import sys, argparse
from itertools import permutations

ap = argparse.ArgumentParser()
ap.add_argument("--labels", required=True, help="read<TAB>label (label may be a subgenome/ambiguous/unassigned)")
ap.add_argument("--paf", required=True, help="reads_to_ref.paf (truth)")
ap.add_argument("--chrom-subg", required=True)
a = ap.parse_args()

cmap = {}
for ln in open(a.chrom_subg):
    p = ln.split()
    if len(p) >= 2:
        cmap[p[0]] = p[1]
best = {}
for ln in open(a.paf):
    p = ln.split("\t")
    if len(p) < 11:
        continue
    sub = cmap.get(p[5])
    if sub is None:
        continue
    sc = int(p[9])
    if p[0] not in best or sc > best[p[0]][0]:
        best[p[0]] = (sc, sub)
truth = {r: s for r, (sc, s) in best.items()}

lab = {}
f = open(a.labels)
first = f.readline().split("\t")
if len(first) >= 2 and first[1].strip() not in ("", "ambiguous", "unassigned", "?", "label"):
    lab[first[0]] = first[1].strip()
for ln in f:
    p = ln.rstrip("\n").split("\t")
    if len(p) >= 2 and p[1] not in ("", "ambiguous", "unassigned", "?"):
        lab[p[0]] = p[1]

tl = sorted(set(truth.values())); pl = sorted(set(lab.values()))
cb = [(lab[r], truth[r]) for r in lab if r in truth]
m = max((dict(zip(pl, perm)) for perm in permutations(tl, min(len(pl), len(tl)))),
        key=lambda mp: sum(mp.get(x) == y for x, y in cb), default={})
correct = sum(1 for r, tv in truth.items() if r in lab and m.get(lab[r]) == tv)
print(f"{100 * correct / len(truth):.1f}")
