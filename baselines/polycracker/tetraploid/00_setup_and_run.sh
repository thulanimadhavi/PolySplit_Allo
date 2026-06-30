#!/usr/bin/env bash
# ============================================================================
# Step 0: stage the IDENTICAL tetraploid Flye assembly into pc_work and run the
# REAL polyCRACKER tool via the existing, tested napus driver (same image, same
# config: k=26, SpectralClustering, tsne, n_subgenomes=2, chunk=100kb). This
# keeps the tetraploid polyCRACKER number directly comparable to the napus one.
#
# RUN on a DOCKER/udocker-capable server (the napus baseline used udocker).
#   nohup bash 00_setup_and_run.sh > 00_polycracker.log 2>&1 &
# polyCRACKER on a ~365 Mb assembly is heavy (blast/kmer steps); napus 921 Mb took ~41 min.
# ============================================================================
set -euo pipefail
BASE=$DATA/camelina/Tetraploid_data/polycracker_baseline
DRIVER=$DATA/polycracker_compare/run_polycracker_napus.sh
FLYE=$DATA/camelina/Tetraploid_data/polysplit_run/flye_out/assembly.fasta
GENOME=Cmi4x_flye.fa
WORK="$BASE/pc_work"
mkdir -p "$WORK"

# stage the SAME assembly PolySplit used (real file, not symlink -- udocker mounts $WORK)
[ -s "$WORK/$GENOME" ] || cp "$FLYE" "$WORK/$GENOME"
echo "staged $GENOME ($(du -h "$WORK/$GENOME" | cut -f1)) in $WORK"
md5sum "$FLYE" "$WORK/$GENOME"   # confirm identical to PolySplit's input

# run polyCRACKER: genome, n_subgenomes=2, chunk=100kb  (identical config to napus)
WORK="$WORK" bash "$DRIVER" "$GENOME" 2 100000

echo "DONE -> $WORK/analysisOutputs_${GENOME%.fa}_n2/"
echo "next: bash run_eval.sh"
