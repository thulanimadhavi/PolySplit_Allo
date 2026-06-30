#!/usr/bin/env bash
# ============================================================================
# polyCRACKER baseline, Camelina HEXAPLOID CmiT1 (K=3). Stages the IDENTICAL HiFi
# assembly PolySplit used (the non-collapsed 591 Mb one; the ONT assembly is
# collapsed, so HiFi is the fair substrate for both methods) and runs the REAL
# polyCRACKER tool via the same tested napus driver, with n_subgenomes=3.
#
# RUN on a DOCKER/udocker-capable server (the napus/tetraploid baselines used udocker).
#   nohup bash 00_setup_and_run.sh > 00_polycracker.log 2>&1 &
# polyCRACKER on ~591 Mb is heavy (kmer/blast/spectral); expect ~1 h+ (napus 921 Mb ~41 min).
# ============================================================================
set -euo pipefail
BASE=$DATA/camelina/Hexaploid_data/polycracker_baseline
DRIVER=$DATA/polycracker_compare/run_polycracker_napus.sh
ASM=$DATA/camelina/Hexaploid_data/hifi_asm_flye/assembly.fasta
GENOME=CmiT1_hifi.fa
WORK="$BASE/pc_work"
mkdir -p "$WORK"

# stage the SAME assembly PolySplit used (real file -- udocker mounts $WORK, no symlinks)
[ -s "$WORK/$GENOME" ] || cp "$ASM" "$WORK/$GENOME"
echo "staged $GENOME ($(du -h "$WORK/$GENOME" | cut -f1)) in $WORK"
md5sum "$ASM" "$WORK/$GENOME"     # confirm identical to PolySplit's HiFi input

# run polyCRACKER: genome, n_subgenomes=3 (hexaploid), chunk=100kb (same config as napus/tetraploid)
WORK="$WORK" bash "$DRIVER" "$GENOME" 3 100000

echo "DONE -> $WORK/analysisOutputs_${GENOME%.fa}_n3/"
echo "next: bash run_eval.sh"
