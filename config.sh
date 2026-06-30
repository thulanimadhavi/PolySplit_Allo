#!/usr/bin/env bash
# Edit these to your environment, then `source config.sh` before running a driver.
# All tool binaries (flye, bwa, samtools, minimap2, gsufsort/gsufsort-64, kmc, yahs)
# are assumed to be on $PATH; set the paths below only if they are not.
export DATA=/path/to/data            # root holding reads, Hi-C, and the eval reference
export POLYSPLIT=/path/to/PolySplit  # this repository
export SUBPHASER=/path/to/SubPhaser  # only for the scaffold-first baseline
export FLYE_BIN=/path/to/flye/bin    # dir containing flye + flye-minimap2 + flye-samtools
