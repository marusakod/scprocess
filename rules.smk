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
CELLBENDER_IMAGE, CELLBENDER_PROP_MAX_KEPT, AMBIENT_METHOD, CUSTOM_PARAMS_F, CELL_CALLS_METHOD, \
FORCE_EXPECTED_CELLS, FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE = \
  get_ambient_parameters(config)
SCE_BENDER_PROB = \
  get_make_sce_parameters(config)
DBL_MIN_FEATS = \
  get_doublet_id_parameters(config)
QC_HARD_MIN_COUNTS, QC_HARD_MIN_FEATS, QC_HARD_MAX_MITO, QC_MIN_COUNTS, QC_MIN_FEATS, \
  QC_MIN_MITO, QC_MAX_MITO, QC_MIN_SPLICE, QC_MAX_SPLICE, QC_MIN_CELLS, QC_FILTER_BENDER = \
  get_qc_parameters(config)
INT_EXC_REGEX, INT_N_HVGS, INT_N_DIMS, INT_DBL_RES, INT_DBL_CL_PROP, INT_THETA, INT_RES_LS, INT_SEL_RES = \
  get_integration_parameters(config, AF_MITO_STR)
MKR_GSEA_DIR, MKR_MIN_CL_SIZE, MKR_MIN_CELLS, MKR_NOT_OK_RE, MKR_MIN_CPM_MKR, MKR_MIN_CPM_GO, MKR_MAX_ZERO_P, MKR_GSEA_CUT, MKR_CANON_F = \
  get_marker_genes_parameters(config, SPECIES, SCPROCESS_DATA_DIR)
LBL_XGB_F, LBL_XGB_CLS_F, LBL_GENE_VAR, LBL_SEL_RES_CL, LBL_MIN_PRED, LBL_MIN_CL_PROP, LBL_MIN_CL_SIZE, LBL_SCE_SUBSETS = \
  get_label_celltypes_parameters(config, SPECIES, SCPROCESS_DATA_DIR)

# specify locations
fastqs_dir    = f"{PROJ_DIR}/data/fastqs"
code_dir      = f"{PROJ_DIR}/code"
af_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_alevin_fry"
amb_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_ambient"
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
      af_dir    + '/af_{sample}/ambient_params_{sample}_' + DATE_STAMP + '.yaml',
      # ambient (cellbender, decontx or nothing)
      amb_dir + '/ambient_{sample}/ambient_{sample}_' + DATE_STAMP + '_output_paths.yaml',
      # barcode qc metrics
      amb_dir + '/ambient_{sample}/barcodes_qc_metrics_{sample}_' + DATE_STAMP + '.txt.gz',
      # doublet id
      dbl_dir   + '/dbl_{sample}/scDblFinder_{sample}_outputs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
      dbl_dir   + '/dbl_{sample}/scDblFinder_{sample}_dimreds_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
      ], sample = SAMPLES), 
      # ambient sample statistics
      amb_dir + '/ambient_sample_statistics_' + DATE_STAMP + '.txt',  
      # make sce input df
      sce_dir + '/sce_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
      # make sce
      sce_dir + '/sce_cells_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
      # doublet_id
      dbl_dir   + '/doublet_id_files_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
      dbl_dir   + '/scDblFinder_combined_outputs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
      dbl_dir   + '/scDblFinder_combined_dimreds_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
      # qc
      qc_dir    + '/qc_dt_'   + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
      qc_dir    + '/keep_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
      # integration
      int_dir   + '/sce_clean_'           + FULL_TAG + '_' + DATE_STAMP + '.rds',
      int_dir   + '/integrated_dt_'       + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
      int_dir   + '/harmony_hvgs_'        + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
      # marker genes
      mkr_dir   + '/pb_'              + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.rds',
      mkr_dir   + '/pb_marker_genes_' + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.txt.gz',
      mkr_dir   + '/pb_hvgs_'         + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.txt.gz',
      mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'go_bp_' + DATE_STAMP + '.txt.gz',
      mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'go_cc_' + DATE_STAMP + '.txt.gz',
      mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'go_mf_' + DATE_STAMP + '.txt.gz',
      mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'paths_' + DATE_STAMP + '.txt.gz',
      mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'hlmk_' + DATE_STAMP + '.txt.gz',
      # code
      code_dir  + '/utils.R',
      code_dir  + '/ambient.R',
      code_dir  + '/qc.R', 
      code_dir  + '/make_sce.R',
      code_dir  + '/doublet_id.R',
      code_dir  + '/integration.R', 
      # code_dir  + '/label_celltypes.R',
      code_dir  + '/marker_genes.R',
      # markdowns
      rmd_dir   + '/' + SHORT_TAG + '_alevin_fry.Rmd',
      rmd_dir   + '/' + SHORT_TAG + '_ambient.Rmd', 
      rmd_dir   + '/' + SHORT_TAG + '_qc.Rmd', 
      rmd_dir   + '/' + SHORT_TAG + '_integration.Rmd', 
      #rmd_dir   + '/' + SHORT_TAG + '_label_celltypes.Rmd', 
      rmd_dir   + '/' + SHORT_TAG + f'_marker_genes_{INT_SEL_RES}.Rmd'
      # reports
      #docs_dir  + '/' + SHORT_TAG + '_alevin_fry.html',
      #docs_dir  + '/' + SHORT_TAG + '_ambient.html'


rule label_and_subset:
  input:
    lbl_dir + '/hvg_mat_for_labelling_' + LBL_GENE_VAR + '_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    lbl_dir + '/xgboost_guesses_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    lbl_dir   + '/sce_subset_specifications_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    expand([
      lbl_dir   +'/sce_subset_' + FULL_TAG + '_{s}_' + DATE_STAMP + '.rds'
      ], s = [] if LBL_SCE_SUBSETS is None else [*LBL_SCE_SUBSETS] )

# define rules that are needed
include: "rules/alevin_fry.smk"
include: "rules/ambient.smk"
include: "rules/make_sce.smk"
include: "rules/doublet_id.smk"
include: "rules/qc.smk"
include: "rules/integration.smk"
include: "rules/marker_genes.smk"
include: "rules/render_htmls.smk"
include: "rules/label_and_subset.smk"

