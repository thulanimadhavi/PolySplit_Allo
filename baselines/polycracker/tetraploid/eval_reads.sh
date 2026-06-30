#!/usr/bin/env bash
# Read-level eval: polyCRACKER per-contig labels -> PolySplit's IDENTICAL Step-7 vote
# (propagate_to_reads.py --weight ident --min-conf 0.6) -> allread_eval.py vs the SAME truth.
set -uo pipefail
BASE=$DATA/camelina/Tetraploid_data/polycracker_baseline
PIPE=$POLYSPLIT/pipeline
PY=python3
W=$DATA/camelina/Tetraploid_data/polysplit_run
CHROM=$POLYSPLIT/chrom_subg.tetraploid.tsv
TOTAL=1647060

[ -s "$BASE/polycracker_contig_labels.tsv" ] || { echo "!! run make_polycracker_contig_labels.py first"; exit 2; }

echo "### propagate polyCRACKER contig labels -> reads (ident, min-conf 0.6) ###"
"$PY" "$PIPE/propagate_to_reads.py" --paf "$W/reads_to_contigs.paf" \
      --contig-labels "$BASE/polycracker_contig_labels.tsv" --min-conf 0.6 --weight ident \
      --out "$BASE/read_polycracker.tsv"

echo "### evaluate vs chromosome-anchored truth (same as PolySplit / ref-guided) ###"
"$PY" "$PIPE/allread_eval.py" --labels "$BASE/read_polycracker.tsv" \
      --ref-paf "$W/reads_to_ref.paf" --total-reads "$TOTAL" --chrom-subg "$CHROM"
echo
echo "TABLE: use T2 (chromosome-truth, ambiguous+unassigned as errors), consistent with"
echo "PolySplit 97.2 and ref-guided 92.3. (napus polyCRACKER was 54.4% read / 56% contig.)"
