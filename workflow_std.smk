# load modules
import yaml
import pandas as pd
import os
import re
import warnings
import glob
from snakemake.utils import validate, min_version

from scprocess_func import *

# get project parameters
PROJ_DIR, FASTQ_DIR, SHORT_TAG, FULL_TAG, YOUR_NAME, AFFILIATION, METADATA_F, METADATA_VARS, EXC_SAMPLES, SAMPLES, DATE_STAMP = \
  get_project_parameters(config)

# get setup parameters
SCPROCESS_DATA_DIR, INDEX_DIR_NAME, FASTA, GTF, USE_DECOYS = \
  get_setup_parameters(config)

# exclude any samples without fastq files
SAMPLES  = exclude_samples_without_fastq_files(PROJ_DIR, SAMPLES)

# one rule to rule them all
# rule all:
#   input:
#     expand(
#       [
#       # alevin_fry
#       af_dir    + '/af_{sample}/af_quant/',
#       af_dir    + '/af_{sample}/af_quant/alevin/quants_mat.mtx',
#       af_dir    + '/af_{sample}/af_quant/alevin/quants_mat_cols.txt',
#       af_dir    + '/af_{sample}/af_quant/alevin/quants_mat_rows.txt',
#       af_dir    + '/af_{sample}/af_counts_mat.h5',
#       af_dir    + '/af_{sample}/knee_plot_data_{sample}_' + DATE_STAMP + '.txt.gz',
#       # cellbender
#       cb_dir    + '/bender_{sample}/bender_{sample}_' + DATE_STAMP + '.h5',
#       cb_dir    + '/bender_{sample}/bender_{sample}_' + DATE_STAMP + '_filtered.h5',
#       cb_dir    + '/bender_{sample}/bender_qc_metrics_{sample}_' + DATE_STAMP + '.txt.gz',
#       # doublet_id
#       dbl_dir   + '/dbl_{sample}/scDblFinder_{sample}_outputs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
#       dbl_dir   + '/dbl_{sample}/scDblFinder_{sample}_dimreds_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
#       ],
#       sample = SAMPLES)


# define rules that are needed
include: "rules/setup.smk"
