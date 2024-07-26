#!/bin/bash

#conda activate simpleaf

ulimit -n 2048

# simpleaf configuration
export ALEVIN_FRY_HOME="/home/mkoder/scprocess_data/alevin_fry_home"
PISCEM="/home/mkoder/software/piscem-x86_64-unknown-linux-gnu/piscem"
simpleaf set-paths

# change working directory to somewhere not crazy
cd $ALEVIN_FRY_HOME

# set up this build
DATA_DIR="/home/mkoder/scprocess_data"
REF_DIR="/home/mkoder/reference_genomes/refdata-gex-GRCh38-2024-A"
IDX_DIR="$DATA_DIR/alevin_fry_home/human_2024-A_splici"

# simpleaf index
simpleaf index \
  --output $IDX_DIR \
  --fasta $REF_DIR/genome.fa \
  --gtf $REF_DIR/genes.gtf \
 # --decoy-paths $REF_DIR/genome.fa \
  --use-piscem \
  --rlen 91 \
  --threads 6

#conda deactivate
