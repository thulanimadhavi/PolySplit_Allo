# PolySplit-Allo

Reference-free subgenome separation of allopolyploid long reads by assembly-first
homoeolog phasing. This repository accompanies the submitted manuscript and is provided
for review.

<p align="center"><img src="figures/Denovo_assembly_polysplit_V2.png" alt="Full PolySplit pipeline schematic" width="100%"></p>

## Layout
- `pipeline/`  core PolySplit stages
  - `homoeolog_graph_from_lcp.py`  shared-33mer homoeolog edges from a gsufsort GSA+LCP
  - `dechimerize_structural.py`    Hi-C blocks + structural de-chimerization (imports `decloud_blocks_v2.py`)
  - `label_small_contigs_v2.py`    repeat-composition contrast labelling + small-contig recovery
  - `homoeolog_repair.py`          homoeolog-consistency repair
  - `propagate_to_reads.py`        identity-weighted contig-label -> read propagation
  - `allread_eval.py`, `eval_contig_labels.py`  read- and contig-level evaluation
- `refguided/`  reference-guided signature baseline (K=2 and K=3)
- `drivers/`    `run_polysplit.sh` (generalized: `--type ont|hifi`, `--nsg K`) and `run_separated_assembly.sh`; plus per-dataset scripts that reproduce the paper runs (PolySplit, reference-guided, sensitivity, timing)
- `baselines/`  wrappers + cluster-to-subgenome mapping + read-eval for the comparison methods
  - `polycracker/{napus,tetraploid,hexaploid}/`  run polyCRACKER, map clusters to subgenomes (best 1:1), score
  - `subphaser/{napus,tetraploid,hexaploid}/`     YaHS scaffold + SubPhaser (scaffold-first), score
  These call the external tools polyCRACKER and SubPhaser, which must be installed separately.
- `readdirect/`  assembly-free control: cluster raw reads by single-copy k-mer incidence and score vs truth (shows the assembly step is necessary)

## Dependencies
Python 3 (numpy, scikit-learn, networkx) and the external tools Flye, bwa, samtools,
minimap2, gsufsort (32- and 64-bit builds), KMC, YaHS, and SubPhaser (baseline only).

## Running
Set the paths in `config.sh` for your environment and `source` it. To separate mixed long reads
into `K` subgenomes (`--type ont` or `--type hifi`) and then assemble each subgenome:
```
source config.sh
bash drivers/run_polysplit.sh --reads reads.fq.gz --hic1 hic_R1.fq.gz --hic2 hic_R2.fq.gz \
     --type ont --nsg 2 --work out/
bash drivers/run_separated_assembly.sh --labels out/read_subg.tsv --reads reads.fq.gz \
     --type ont --work out/ --hic1 hic_R1.fq.gz --hic2 hic_R2.fq.gz
```
`run_polysplit.sh` writes per-read subgenome labels to `out/read_subg.tsv`; add
`--ref REF --chrom-subg CS` for the optional accuracy evaluation against a reference. The
per-dataset `run_polysplit_<name>.sh` scripts reproduce the exact runs reported in the paper.
Every stage is idempotent: it skips if its output already exists.

## Data
Sequencing inputs and the evaluation reference are not redistributed here; see the
manuscript's data-availability statement for accessions.
