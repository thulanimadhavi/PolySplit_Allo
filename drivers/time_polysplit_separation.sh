#!/usr/bin/env bash
# Time the PolySplit SEPARATION stages on existing inputs -- excludes Flye (shared assembly) and the
# Hi-C / read alignment (shared prerequisites, reused by the baselines). Times the PolySplit-specific
# compute: (1) homoeolog edges (gsufsort GSA+LCP + LCP walk), (2) blocks/de-chimerize/label/repair,
# (3) propagate to reads. Re-runs to a temp dir so existing results are untouched.
#   bash time_polysplit_separation.sh <WORK_dir> <assembly.fasta> <NSG> [reads_to_contigs.paf]
set -uo pipefail
WORK=$1; ASM=$2; NSG=${3:-2}; R2C=${4:-$WORK/reads_to_contigs.paf}
PIPE=$POLYSPLIT/pipeline
PY=python3
GSUFSORT=gsufsort
T=$WORK/timing_sep; mkdir -p "$T"
secs(){ date +%s; }
echo "== PolySplit separation timing (excl. Flye + alignments) | WORK=$WORK NSG=$NSG =="

s=$(secs)
"$GSUFSORT" "$ASM" --fasta --gsa 4 8 --output "$T/k33" >/dev/null 2>&1
"$GSUFSORT" "$ASM" --fasta --lcp 1  --output "$T/k33" >/dev/null 2>&1
"$PY" "$PIPE/homoeolog_graph_from_lcp.py" --gsa "$T/k33.4.8.gsa" --lcp "$T/k33.1.lcp" \
      --fa "$ASM" --k 33 --min-copy 2 --max-copy 6 --out "$T/edges.tsv" >/dev/null 2>&1
e=$(secs); T_EDGES=$((e-s)); echo "  (1) homoeolog edges (gsufsort+LCP): ${T_EDGES}s"

s=$(secs)
export POLYSPLIT_NSG="$NSG" POLYSPLIT_FASTA="$ASM" POLYSPLIT_CONTACTS="$WORK/contacts.pkl" \
       POLYSPLIT_EDGES="$T/edges.tsv" POLYSPLIT_TRUTH=/none
( cd "$T"
  "$PY" "$PIPE/dechimerize_structural.py" >/dev/null 2>&1
  POLYSPLIT_OUT="$T/v2.tsv" "$PY" "$PIPE/label_small_contigs_v2.py" "$T/decloud_structural_contig_labels.tsv" >/dev/null 2>&1
  POLYSPLIT_LABELS_IN="$T/v2.tsv" POLYSPLIT_OUT="$T/repaired.tsv" "$PY" "$PIPE/homoeolog_repair.py" >/dev/null 2>&1 )
e=$(secs); T_LABEL=$((e-s)); echo "  (2) blocks/de-chimerize/label/repair: ${T_LABEL}s"

s=$(secs)
"$PY" "$PIPE/propagate_to_reads.py" --paf "$R2C" --contig-labels "$T/repaired.tsv" \
      --min-conf 0.6 --weight ident --out "$T/read_subg.tsv" >/dev/null 2>&1
e=$(secs); T_PROP=$((e-s)); echo "  (3) propagate to reads: ${T_PROP}s"

TOT=$((T_EDGES+T_LABEL+T_PROP))
printf "TOTAL PolySplit separation: %ds = %dm %ds\n" "$TOT" "$((TOT/60))" "$((TOT%60))"
