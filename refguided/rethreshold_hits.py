import sys, argparse
from itertools import permutations

ap = argparse.ArgumentParser()
ap.add_argument("--hits", required=True, help="read_labels_refguided.tsv (read_id, *_hits cols, label)")
ap.add_argument("--paf", required=True, help="reads_to_ref.paf (truth)")
ap.add_argument("--chrom-subg", required=True)
a = ap.parse_args()

f = open(a.hits)
hdr = f.readline().rstrip("\n").split("\t")
idx = [i for i, h in enumerate(hdr) if h.endswith("_hits")]
sgs = [hdr[i][:-5] for i in idx]
hits = {}
for ln in f:
    p = ln.rstrip("\n").split("\t")
    if len(p) <= max(idx):
        continue
    try:
        hits[p[0]] = tuple(int(p[i]) for i in idx)
    except ValueError:
        pass

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


def classify(h, t, rho):
    order = sorted(range(len(h)), key=lambda k: h[k], reverse=True)
    top = order[0]; rest = max((h[k] for k in order[1:]), default=0)
    if h[top] >= t and h[top] >= rho * (rest + 1):
        return sgs[top]
    return None


def acc(t, rho):
    lab = {r: c for r, h in hits.items() if (c := classify(h, t, rho)) is not None}
    tl = sorted(set(truth.values())); pl = sorted(set(lab.values()))
    cb = [(lab[r], truth[r]) for r in lab if r in truth]
    m = max((dict(zip(pl, perm)) for perm in permutations(tl, min(len(pl), len(tl)))),
            key=lambda mp: sum(mp.get(x) == y for x, y in cb), default={})
    correct = sum(1 for r, tv in truth.items() if r in lab and m.get(lab[r]) == tv)
    return 100 * correct / len(truth)


print(f"# K={len(sgs)} {sgs}; truth reads={len(truth):,}; hit reads={len(hits):,}")
print("## t sweep (rho=3.0):")
for t in [1, 2, 3, 5, 10, 20, 50, 100]:
    print(f"  t={t:<4} acc={acc(t, 3.0):.1f}")
print("## rho sweep (t=3):")
for rho in [1.5, 2, 3, 5, 10, 20, 50, 100]:
    print(f"  rho={rho:<5} acc={acc(3, rho):.1f}")
