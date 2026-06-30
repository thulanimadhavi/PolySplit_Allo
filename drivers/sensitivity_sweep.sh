#!/usr/bin/env bash
set -uo pipefail
WORK=$1; NSG=$2; CHROM=$3; TOTAL=$4
PIPE=$POLYSPLIT/pipeline; PY=python3
for BETA in 0.50 0.55 0.60 0.65 0.70; do
  "$PY" "$PIPE/propagate_to_reads.py" --paf "$WORK/reads_to_contigs.paf" \
        --contig-labels "$WORK/all_contig_labels_repaired.tsv" --min-conf "$BETA" --weight ident \
        --out "$WORK/sweep_read_b$BETA.tsv" >/dev/null 2>&1
  echo "### beta=$BETA ###"
  "$PY" "$PIPE/allread_eval.py" --labels "$WORK/sweep_read_b$BETA.tsv" \
        --ref-paf "$WORK/reads_to_ref.paf" --total-reads "$TOTAL" --chrom-subg "$CHROM" \
    | grep -E "T2 \(over|reads: chromosome-truth"
done
