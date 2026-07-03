import argparse, gzip, random
import numpy as np

BASE = [-1] * 256
for ch, v in zip(b"ACGTacgt", [0, 1, 2, 3, 0, 1, 2, 3]):
    BASE[ch] = v


def canon_codes(seq, k):
    fwd = 0; rev = 0; mask = (1 << (2 * k)) - 1; sh = 2 * (k - 1); valid = 0
    for ch in seq.encode():
        c = BASE[ch]
        if c < 0:
            valid = 0; fwd = 0; rev = 0; continue
        fwd = ((fwd << 2) | c) & mask
        rev = (rev >> 2) | ((3 - c) << sh)
        valid += 1
        if valid >= k:
            yield fwd if fwd < rev else rev


def iter_reads(path):
    f = gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")
    first = f.readline()
    if not first:
        return
    if first[0] == ">":
        rid = first[1:].split()[0]; chunks = []
        for line in f:
            if line[0] == ">":
                yield rid, "".join(chunks); rid = line[1:].split()[0]; chunks = []
            else:
                chunks.append(line.strip())
        yield rid, "".join(chunks)
    else:
        h = first
        while True:
            seq = f.readline(); f.readline(); f.readline()
            if not seq:
                return
            yield h[1:].split()[0], seq.strip()
            h = f.readline()
            if not h:
                return


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads", required=True)
    ap.add_argument("--band", required=True, help="kmc_dump of the single-copy count band")
    ap.add_argument("--paf", required=True, help="reads_to_ref.paf (truth)")
    ap.add_argument("--chrom-subg", required=True)
    ap.add_argument("--k", type=int, default=21)
    ap.add_argument("--nsg", type=int, default=2)
    ap.add_argument("--max-feats", type=int, default=3_000_000)
    ap.add_argument("--sample-reads", type=int, default=80_000)
    ap.add_argument("--min-len", type=int, default=5000)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--svd", type=int, default=40)
    ap.add_argument("--min-feat-reads", type=int, default=4)
    ap.add_argument("--max-feat-frac", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)

    feats = []; n_seen = 0
    with open(a.band) as f:
        for line in f:
            km = line.split("\t", 1)[0] if "\t" in line else line.split()[0]
            n_seen += 1
            if len(feats) < a.max_feats:
                feats.append(km)
            else:
                j = rng.randint(0, n_seen - 1)
                if j < a.max_feats:
                    feats[j] = km
    code_of = {}
    for j, km in enumerate(feats):
        cs = list(canon_codes(km, a.k))
        if cs:
            code_of[cs[0]] = j
    print(f"[features] sampled {len(feats):,} of {n_seen:,} band k-mers", flush=True)

    rows = []; cols = []; read_ids = []
    for rid, seq in iter_reads(a.reads):
        if len(seq) < a.min_len:
            continue
        present = set()
        for i, code in enumerate(canon_codes(seq, a.k)):
            if a.stride > 1 and (i % a.stride):
                continue
            j = code_of.get(code)
            if j is not None:
                present.add(j)
        r = len(read_ids)
        for j in present:
            rows.append(r); cols.append(j)
        read_ids.append(rid)
        if len(read_ids) >= a.sample_reads:
            break
    nread = len(read_ids); nfeat = len(feats); nnz = len(rows)
    print(f"[matrix] {nread:,} reads x {nfeat:,} features, {nnz:,} nonzeros", flush=True)

    from scipy.sparse import csr_matrix
    from sklearn.decomposition import TruncatedSVD
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import normalize
    M = csr_matrix((np.ones(nnz, dtype=np.float32), (rows, cols)), shape=(nread, nfeat))
    colsum = np.asarray(M.sum(0)).ravel()
    keep = np.where((colsum >= a.min_feat_reads) & (colsum <= a.max_feat_frac * nread))[0]
    M = M[:, keep]
    print(f"[matrix] {M.shape[1]:,} features survive the >= {a.min_feat_reads}-reads filter")
    if M.shape[1] < 5:
        print("[RESULT] no feature recurs across reads; clustering cannot separate subgenomes")
        return

    ncomp = min(a.svd, M.shape[1] - 1, nread - 1)
    X = normalize(TruncatedSVD(n_components=ncomp, random_state=a.seed).fit_transform(M))
    lab = KMeans(n_clusters=a.nsg, n_init=10, random_state=a.seed).fit_predict(X)
    sizes = [int((lab == c).sum()) for c in range(a.nsg)]
    print(f"[cluster] sizes {sizes}")

    cmap = {}
    with open(a.chrom_subg) as f:
        for ln in f:
            p = ln.split()
            if len(p) >= 2:
                cmap[p[0]] = p[1]
    sampled = set(read_ids)
    best = {}
    with open(a.paf) as f:
        for ln in f:
            p = ln.split("\t")
            if len(p) < 11 or p[0] not in sampled:
                continue
            sub = cmap.get(p[5])
            if sub is None:
                continue
            sc = int(p[9])
            if p[0] not in best or sc > best[p[0]][0]:
                best[p[0]] = (sc, sub)
    truth = {rid: s for rid, (sc, s) in best.items()}
    idx_truth = [i for i, rid in enumerate(read_ids) if rid in truth]
    print(f"[truth] {len(idx_truth):,} of {nread:,} sampled reads have a chromosome truth")
    if not idx_truth:
        print("[RESULT] no truth overlap; cannot score")
        return

    from itertools import permutations
    subs = sorted({truth[read_ids[i]] for i in idx_truth})
    clusters = sorted(set(int(lab[i]) for i in idx_truth))
    best_acc = 0.0
    for perm in permutations(subs, min(len(clusters), len(subs))):
        m = dict(zip(clusters, perm))
        correct = sum(1 for i in idx_truth if m.get(int(lab[i])) == truth[read_ids[i]])
        best_acc = max(best_acc, correct / len(idx_truth))
    print(f"[RESULT] read-clustering accuracy (best 1:1) = {100*best_acc:.1f}%   "
          f"(chance = {100.0/a.nsg:.1f}%)")


if __name__ == "__main__":
    main()
