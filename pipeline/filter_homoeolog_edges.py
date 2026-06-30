#!/usr/bin/env python3
import sys
from collections import defaultdict

inp, out = sys.argv[1], sys.argv[2]
FRAC = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5

edges = []
mx = defaultdict(int)
with open(inp) as f:
    hdr = f.readline()
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) < 3:
            continue
        try:
            w = int(p[2])
        except ValueError:
            continue
        edges.append((p[0], p[1], w))
        if w > mx[p[0]]:
            mx[p[0]] = w
        if w > mx[p[1]]:
            mx[p[1]] = w

with open(out, "w") as o:
    o.write(hdr)
    kept = 0
    for a, b, w in edges:
        if w >= FRAC * mx[a] and w >= FRAC * mx[b]:
            o.write(f"{a}\t{b}\t{w}\n")
            kept += 1
sys.stderr.write(f"[filter] kept {kept}/{len(edges)} edges (FRAC={FRAC})\n")
