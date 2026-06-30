#!/usr/bin/env bash
BASE=$DATA/camelina/Hexaploid_data/polycracker_baseline
PY=python3
PIPE=$POLYSPLIT/pipeline
HX=$DATA/camelina/Hexaploid_data
WH="$HX/polysplit_run_hifi"
FAI="$HX/hifi_asm_flye/assembly.fasta.fai"
CHROM=$POLYSPLIT/chrom_subg.hexaploid.tsv

CLUSTER="$(find "$BASE/pc_work" -type d -name clusterResults 2>/dev/null | head -1)"
if [ -z "$CLUSTER" ] || [ -z "$(ls "$CLUSTER"/subgenome_*.txt 2>/dev/null)" ]; then
  echo "!! no clusterResults/subgenome_*.txt under $BASE/pc_work -- polyCRACKER not finished yet."
  ls -d "$BASE"/pc_work/analysisOutputs_* 2>/dev/null || echo "   (no outputs yet)"
  exit 2
fi
echo "using clusterResults: $CLUSTER  ($(ls "$CLUSTER"/subgenome_*.txt | wc -l) clusters)"

echo "############ per-contig S1/S2/S3 labels (K=3) ############"
"$PY" "$BASE/make_polycracker_contig_labels.py" "$CLUSTER" "$BASE/polycracker_contig_labels.tsv" \
      "$WH/wg_purity.per_contig.tsv"

echo
echo "############ CONTIG-level accuracy (best-1:1, same eval as PolySplit) ############"
"$PY" "$PIPE/eval_contig_labels.py" "$CHROM" "$WH/contigs_to_ref.paf" "$FAI" \
      "$BASE/polycracker_contig_labels.tsv" "$BASE/wg_purity.polycracker.tsv"

echo
echo "############ READ-level accuracy ############"
bash "$BASE/eval_reads.sh"
