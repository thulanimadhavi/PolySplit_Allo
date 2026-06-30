#!/usr/bin/env bash
set -uo pipefail
BASE=$DATA/camelina/Hexaploid_data/polycracker_baseline
PIPE=$POLYSPLIT/pipeline
PY=python3
WH=$DATA/camelina/Hexaploid_data/polysplit_run_hifi
CHROM=$POLYSPLIT/chrom_subg.hexaploid.tsv
TOTAL=925873

[ -s "$BASE/polycracker_contig_labels.tsv" ] || { echo "!! run make_polycracker_contig_labels.py first (run_eval.sh)"; exit 2; }

echo "### propagate polyCRACKER contig labels -> reads (ident, min-conf 0.6) ###"
"$PY" "$PIPE/propagate_to_reads.py" --paf "$WH/reads_to_contigs.paf" \
      --contig-labels "$BASE/polycracker_contig_labels.tsv" --min-conf 0.6 --weight ident \
      --out "$BASE/read_polycracker.tsv"

echo "### evaluate vs chromosome-anchored truth (same as PolySplit / ref-guided) ###"
"$PY" "$PIPE/allread_eval.py" --labels "$BASE/read_polycracker.tsv" \
      --ref-paf "$WH/reads_to_ref.paf" --total-reads "$TOTAL" --chrom-subg "$CHROM"
echo
echo "TABLE: report T2 (chromosome-truth; ambiguous+unassigned as errors), consistent with"
echo "PolySplit HiFi 96.2 and reference-guided 90.8 for the hexaploid."
