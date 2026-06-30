#!/usr/bin/env python3
"""Filter homoeolog edges to each contig's strongest partner(s), removing spurious
within-subgenome paralog/repeat edges that contaminate the subgenome-contrast labelling on
highly contiguous assemblies (e.g. HiFi). A contig's true homoeolog is its dominant partner,
far stronger than within-subgenome paralog edges, so we keep an edge (a,b,w) only when w is at
least FRAC of the maximum partner weight at BOTH endpoints. This needs no truth and leaves
fragmented (e.g. ONT) assemblies unchanged, where the spurious edges are already sub-threshold.

Usage: filter_homoeolog_edges.py <in_edges.tsv> <out_edges.tsv> [frac=0.5]
"""
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
