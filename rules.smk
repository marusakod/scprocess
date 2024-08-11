# load modules
import yaml
import pandas as pd
import os
import re
import warnings
import glob
from snakemake.utils import validate, min_version

from scprocess_utils import *

SCPROCESS_DATA_DIR = os.getenv('SCPROCESS_DATA_DIR')
# get parameters
PROJ_DIR, FASTQ_DIR, SHORT_TAG, FULL_TAG, YOUR_NAME, AFFILIATION, METADATA_F, METADATA_VARS, EXC_SAMPLES, SAMPLES, DATE_STAMP = \
  get_project_parameters(config)
SPECIES, AF_MITO_STR, AF_HOME_DIR, AF_INDEX_DIR, AF_GTF_DT_F, CHEMISTRY_F = \
  get_alevin_parameters(config, SCPROCESS_DATA_DIR)
CELLBENDER_IMAGE, CELLBENDER_PROP_MAX_KEPT, DO_CELLBENDER, CUSTOM_CELLBENDER_PARAMS_F, \
FORCE_EXPECTED_CELLS, FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_LOW_COUNT_THRESHOLD, \
CELLBENDER_LEARNING_RATE = get_cellbender_parameters(config)

# specify locations
fastqs_dir    = f"{PROJ_DIR}/data/fastqs"
code_dir      = f"{PROJ_DIR}/code"
af_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_alevin_fry"
cb_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_cellbender"
sce_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_make_sce"
dbl_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_doublet_id"
qc_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_qc"
int_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_integration"
mkr_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_marker_genes"
lbl_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_label_celltypes"
meta_dir      = f"{PROJ_DIR}/output/{SHORT_TAG}_metacells"
pb_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_pseudobulk"
empty_dir     = f"{PROJ_DIR}/output/{SHORT_TAG}_empties"
zoom_dir      = f"{PROJ_DIR}/output/{SHORT_TAG}_zoom"
rmd_dir       = f"{PROJ_DIR}/analysis"
docs_dir      = f"{PROJ_DIR}/public"

# exclude any samples without fastq files
SAMPLES       = exclude_samples_without_fastq_files(FASTQ_DIR, SAMPLES)
# join all samples into a single string
SAMPLE_STR = ','.join(SAMPLES)

# one rule to rule them all
rule all:
  input:
    #chemistry
    af_dir + '/chemistry_stats.csv', 
    expand(
      [
      # alevin_fry
      af_dir    + '/af_{sample}/af_quant/',
      af_dir    + '/af_{sample}/af_quant/alevin/quants_mat.mtx',
      af_dir    + '/af_{sample}/af_quant/alevin/quants_mat_cols.txt',
      af_dir    + '/af_{sample}/af_quant/alevin/quants_mat_rows.txt',
      af_dir    + '/af_{sample}/af_counts_mat.h5',
      af_dir    + '/af_{sample}/knee_plot_data_{sample}_' + DATE_STAMP + '.txt.gz',
      af_dir    + '/af_{sample}/bender_params_{sample}_' + DATE_STAMP + '.yaml',
      # cellbender
      cb_dir    + '/bender_{sample}/bender_{sample}_' + DATE_STAMP + '.h5',
      cb_dir    + '/bender_{sample}/bender_{sample}_' + DATE_STAMP + '_filtered.h5',
      cb_dir    + '/bender_{sample}/bender_qc_metrics_{sample}_' + DATE_STAMP + '.txt.gz',
      ], sample = SAMPLES), 
      # find bender bad samples
      cb_dir    + '/bender_bad_samples_' + DATE_STAMP + '.txt'


# define rules that are needed
include: "rules/alevin_fry.smk"
include: "rules/cellbender.smk"

