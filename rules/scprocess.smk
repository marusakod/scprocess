# load modules
import os
import sys
import re
import glob
import pandas as pd
import warnings
import yaml
from snakemake.utils import validate, min_version

# import utils
sys.path.append('scripts')
from scprocess_utils import *

SCPROCESS_DATA_DIR = os.getenv('SCPROCESS_DATA_DIR')

# get parameters
PROJ_DIR, FASTQ_DIR, SHORT_TAG, FULL_TAG, YOUR_NAME, AFFILIATION, METADATA_F, METADATA_VARS, EXC_SAMPLES, SAMPLES, DATE_STAMP, CUSTOM_SAMPLE_PARAMS_F, SPECIES, DEMUX_TYPE, HTO_FASTQ_DIR, FEATURE_REF, DEMUX_F, BATCH_VAR, EXC_POOLS, POOL_IDS, SAMPLE_VAR, SAMPLE_MAPPING = get_project_parameters(config, SCPROCESS_DATA_DIR)
AF_MITO_STR, AF_HOME_DIR, AF_INDEX_DIR, AF_GTF_DT_F, CHEMISTRY = \
  get_alevin_parameters(config, SCPROCESS_DATA_DIR, SPECIES)
CELLBENDER_IMAGE, CELLBENDER_VERSION, CELLBENDER_PROP_MAX_KEPT, AMBIENT_METHOD, CELL_CALLS_METHOD, \
FORCE_EXPECTED_CELLS, FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE, CELLBENDER_POSTERIOR_BATCH_SIZE = \
  get_ambient_parameters(config)
QC_HARD_MIN_COUNTS, QC_HARD_MIN_FEATS, QC_HARD_MAX_MITO, QC_MIN_COUNTS, QC_MIN_FEATS, QC_MIN_MITO, QC_MAX_MITO, QC_MIN_SPLICE, QC_MAX_SPLICE, QC_MIN_CELLS, DBL_MIN_FEATS, EXCLUDE_MITO = \
  get_qc_parameters(config)
HVG_METHOD, HVG_GROUP_VAR, HVG_CHUNK_SIZE, NUM_CHUNKS, GROUP_NAMES, CHUNK_NAMES, N_HVGS, EXCLUDE_AMBIENT_GENES = \
  get_hvg_parameters(config, METADATA_F, AF_GTF_DT_F)
INT_CL_METHOD, INT_REDUCTION, INT_N_DIMS, INT_THETA, INT_RES_LS, INT_DBL_RES, INT_DBL_CL_PROP = \
  get_integration_parameters(config)
MKR_SEL_RES, MKR_GSEA_DIR, MKR_MIN_CL_SIZE, MKR_MIN_CELLS, MKR_NOT_OK_RE, MKR_MIN_CPM_MKR, MKR_MIN_CPM_GO, MKR_MAX_ZERO_P, MKR_GSEA_CUT, CUSTOM_MKR_NAMES, CUSTOM_MKR_PATHS = \
  get_marker_genes_parameters(config, PROJ_DIR, SCPROCESS_DATA_DIR)
LBL_XGB_F, LBL_XGB_CLS_F, LBL_GENE_VAR, LBL_SEL_RES_CL, LBL_MIN_PRED, LBL_MIN_CL_PROP, LBL_MIN_CL_SIZE, LBL_TISSUE = \
  get_label_celltypes_parameters(config, SPECIES, SCPROCESS_DATA_DIR)
AMBIENT_GENES_LOGFC_THR, AMBIENT_GENES_FDR_THR = \
  get_pb_empties_parameters(config)
RETRIES, MB_RUN_MAPPING, MB_SAVE_ALEVIN_TO_H5, \
  MB_RUN_AMBIENT, MB_GET_BARCODE_QC_METRICS, \
  MB_RUN_QC, MB_RUN_HVGS, \
  MB_RUN_INTEGRATION, MB_MAKE_CLEAN_SCES, \
  MB_RUN_MARKER_GENES, MB_RENDER_HTMLS, \
  MB_LABEL_CELLTYPES, \
  MB_PB_MAKE_PBS, MB_PB_CALC_EMPTY_GENES, MB_MAKE_HTO_SCE_OBJECTS, MB_MAKE_SUBSET_SCES = \
  get_resource_parameters(config)

# specify locations
benchmark_dir = f"{PROJ_DIR}/.resources"
fastqs_dir    = f"{PROJ_DIR}/data/fastqs"
code_dir      = f"{PROJ_DIR}/code"
af_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_mapping"
af_rna_dir    = 'rna/' if DEMUX_TYPE == 'af' else ''
amb_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_ambient"
demux_dir     = f"{PROJ_DIR}/output/{SHORT_TAG}_demultiplexing"
dbl_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_doublet_id"
qc_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_qc"
hvg_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_hvg"
int_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_integration"
mkr_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_marker_genes"
lbl_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_label_celltypes"
meta_dir      = f"{PROJ_DIR}/output/{SHORT_TAG}_metacells"
pb_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_pseudobulk"
empty_dir     = f"{PROJ_DIR}/output/{SHORT_TAG}_empties"
rmd_dir       = f"{PROJ_DIR}/analysis"
docs_dir      = f"{PROJ_DIR}/public"


# exclude all samples without fastq files
if DEMUX_TYPE != "none":
  POOL_IDS = exclude_samples_without_fastq_files(FASTQ_DIR, POOL_IDS, HTO=False)
else:
  SAMPLES  = exclude_samples_without_fastq_files(FASTQ_DIR, SAMPLES, HTO=False)

# exclude all samples without hto fastq files
if DEMUX_TYPE == "af":
  POOL_IDS = exclude_samples_without_fastq_files(HTO_FASTQ_DIR, POOL_IDS, HTO=True)

# exclude pools and samples from sample mapping dictionary
SAMPLE_MAPPING = filter_sample_mapping(SAMPLE_MAPPING, POOL_IDS, SAMPLES)

# join all samples into a single string
POOL_STR   = ','.join(POOL_IDS)
SAMPLE_STR = ','.join(SAMPLES)

runs = POOL_IDS if DEMUX_TYPE != "none" else SAMPLES
RUNS_STR = ','.join(runs)

# scripts
r_scripts = [
  code_dir  + '/utils.R',
  code_dir  + '/mapping.R',
  code_dir  + '/ambient.R',
  code_dir  + '/qc.R', 
  code_dir  + '/hvgs.R', 
  code_dir  + '/integration.R', 
  code_dir  + '/marker_genes.R'
  ]
if DEMUX_TYPE == 'af':
  r_scripts.append(code_dir + '/multiplexing.R')

# alevin hto index outputs (optional)
hto_index_outs = [
    af_dir + '/hto.tsv',
    af_dir + '/t2g_hto.tsv',
    af_dir + '/hto_index/ref_indexing.log'
  ] if DEMUX_TYPE == "af" else []

# alevin hto quantification outputs (optional)
hto_af_outs = expand(
  [
  af_dir    + '/af_{run}/hto/af_quant/',
  af_dir    + '/af_{run}/hto/af_quant/alevin/quants_mat.mtx',
  af_dir    + '/af_{run}/hto/af_quant/alevin/quants_mat_cols.txt',
  af_dir    + '/af_{run}/hto/af_quant/alevin/quants_mat_rows.txt',
  af_dir    + '/af_{run}/hto/af_hto_counts_mat.h5',
  af_dir    + '/af_{run}/hto/knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz'
  ], run = runs
) if DEMUX_TYPE == 'af' else []

# seurat demultiplexing outputs (optional)
hto_sce_fs = expand(
  demux_dir + '/sce_cells_htos_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
  run = runs
  ) if DEMUX_TYPE == "af" else []

# multiplexing report (optional)
hto_rmd_f  = (rmd_dir   + '/' + SHORT_TAG + '_demultiplexing.Rmd') if DEMUX_TYPE == 'af' else []
hto_html_f = (docs_dir  + '/' + SHORT_TAG + '_demultiplexing.html') if DEMUX_TYPE == 'af' else []

# fgsea outputs (optional)
fgsea_outs = [
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{MKR_SEL_RES}_' + 'go_bp_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{MKR_SEL_RES}_' + 'go_cc_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{MKR_SEL_RES}_' + 'go_mf_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{MKR_SEL_RES}_' + 'paths_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{MKR_SEL_RES}_' + 'hlmk_' + DATE_STAMP + '.txt.gz', 
] if SPECIES in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020'] else []

# cellbender report (optional)
ambient_rmd_f  = (rmd_dir  + '/' + SHORT_TAG + '_ambient.Rmd')
ambient_html_f = (docs_dir + '/' + SHORT_TAG + '_ambient.html')

# one rule to rule them all
rule all:
  input:
    # hto outputs
    hto_index_outs,
    hto_af_outs,
    expand(
      [
      # mapping
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_quant/',
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_quant/alevin/quants_mat.mtx',
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_quant/alevin/quants_mat_cols.txt',
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_quant/alevin/quants_mat_rows.txt',
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_counts_mat.h5',
      af_dir    + '/af_{run}/' + af_rna_dir + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz',
      af_dir    + '/af_{run}/' + af_rna_dir + 'ambient_params_{run}_' + DATE_STAMP + '.yaml',
      # ambient (cellbender, decontx or nothing)
      amb_dir   + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml',
      amb_dir   + '/ambient_{run}/barcodes_qc_metrics_{run}_' + DATE_STAMP + '.txt.gz',
      # doublet id
      dbl_dir   + '/dbl_{run}/scDblFinder_{run}_outputs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
      dbl_dir   + '/dbl_{run}/scDblFinder_{run}_dimreds_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
      ], run =  runs),
    # ambient sample statistics
    amb_dir + '/ambient_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    # demultiplexing
    hto_sce_fs,
    # qc
    qc_dir  + '/qc_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    qc_dir  + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    # pseudobulks and empties
    pb_dir  + '/af_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    pb_dir  + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
    pb_dir  + '/pb_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    empty_dir + '/edger_empty_genes_' + FULL_TAG + '_all_' + DATE_STAMP + '.txt.gz', 
    # hvgs
    hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    hvg_dir + '/standardized_variance_stats_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    hvg_dir + '/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    hvg_dir + '/top_hvgs_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
    hvg_dir + '/top_hvgs_doublet_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
    # integration
    int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    expand(int_dir + '/sce_cells_clean_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.rds', sample = SAMPLES),
    int_dir  + '/sce_clean_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml', 
    # marker genes
    mkr_dir + '/pb_' + FULL_TAG + f'_{MKR_SEL_RES}_' + DATE_STAMP + '.rds',
    mkr_dir + '/pb_marker_genes_' + FULL_TAG + f'_{MKR_SEL_RES}_' + DATE_STAMP + '.txt.gz',
    mkr_dir + '/pb_hvgs_' + FULL_TAG + f'_{MKR_SEL_RES}_' + DATE_STAMP + '.txt.gz',
    # fgsea outputs
    fgsea_outs, 
    # code
    r_scripts,
    # markdowns
    rmd_dir   + '/' + SHORT_TAG + '_mapping.Rmd',
    rmd_dir   + '/' + SHORT_TAG + '_ambient.Rmd',
    rmd_dir   + '/' + SHORT_TAG + '_qc.Rmd', 
    rmd_dir   + '/' + SHORT_TAG + '_hvgs.Rmd',
    rmd_dir   + '/' + SHORT_TAG + '_integration.Rmd', 
    rmd_dir   + '/' + SHORT_TAG + f'_marker_genes_{MKR_SEL_RES}.Rmd', 
    hto_rmd_f,
    # reports
    docs_dir  + '/' + SHORT_TAG + '_mapping.html', 
    docs_dir  + '/' + SHORT_TAG + '_ambient.html', 
    docs_dir  + '/' + SHORT_TAG + '_qc.html',
    docs_dir  + '/' + SHORT_TAG + '_hvgs.html',
    docs_dir  + '/' + SHORT_TAG + '_integration.html',
    docs_dir  + '/' + SHORT_TAG + f'_marker_genes_{MKR_SEL_RES}.html',
    hto_html_f 


rule mapping:
  input:
    hto_index_outs, 
    hto_af_outs,
    expand( 
      [
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_quant/',
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_quant/alevin/quants_mat.mtx',
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_quant/alevin/quants_mat_cols.txt',
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_quant/alevin/quants_mat_rows.txt',
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_counts_mat.h5',
      af_dir    + '/af_{run}/' + af_rna_dir + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz',
      af_dir    + '/af_{run}/' + af_rna_dir + 'ambient_params_{run}_' + DATE_STAMP + '.yaml'
      ],
      run = runs), 
    rmd_dir   + '/' + SHORT_TAG + '_mapping.Rmd',
    docs_dir  + '/' + SHORT_TAG + '_mapping.html'


rule demux:
  input: 
    hto_sce_fs,
    hto_rmd_f, 
    hto_html_f


rule ambient:
  input: 
    expand(amb_dir   + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml', run = runs), 
    expand(amb_dir + '/ambient_{run}/barcodes_qc_metrics_{run}_' + DATE_STAMP + '.txt.gz', run = runs),
    amb_dir + '/ambient_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    ambient_rmd_f, 
    ambient_html_f


rule qc:
  params:
    mito_str        = AF_MITO_STR,
    exclude_mito    = EXCLUDE_MITO,
    hard_min_counts = QC_HARD_MIN_COUNTS,
    hard_min_feats  = QC_HARD_MIN_FEATS,
    hard_max_mito   = QC_HARD_MAX_MITO,
    min_counts      = QC_MIN_COUNTS,
    min_feats       = QC_MIN_FEATS,
    min_mito        = QC_MIN_MITO,
    max_mito        = QC_MAX_MITO,
    min_splice      = QC_MIN_SPLICE,
    max_splice      = QC_MAX_SPLICE,
    sample_var      = SAMPLE_VAR,
    demux_type      = DEMUX_TYPE,
    dbl_min_feats   = DBL_MIN_FEATS
  input:     
    expand(
      [
    # doublet id
      dbl_dir + '/dbl_{run}/scDblFinder_{run}_outputs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
      dbl_dir + '/dbl_{run}/scDblFinder_{run}_dimreds_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
      ], 
      run = runs),
    qc_dir  + '/qc_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    qc_dir  + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    qc_dir  + '/sce_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml', 
    rmd_dir   + '/' + SHORT_TAG + '_qc.Rmd',
    docs_dir  + '/' + SHORT_TAG + '_qc.html'


rule hvg:
  params:
    hvg_method = HVG_METHOD,
    n_hvgs = N_HVGS,
    exclude_ambient_genes = EXCLUDE_AMBIENT_GENES
  input:
    hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    hvg_dir + '/standardized_variance_stats_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    hvg_dir + '/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    hvg_dir + '/top_hvgs_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
    hvg_dir + '/top_hvgs_doublet_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5',
    rmd_dir   + '/' + SHORT_TAG + '_hvgs.Rmd',
    docs_dir  + '/' + SHORT_TAG + '_hvgs.html'


rule integration:
  input:
    int_dir   + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    rmd_dir   + '/' + SHORT_TAG + '_integration.Rmd', 
    docs_dir  + '/' + SHORT_TAG + '_integration.html'


rule marker_genes:
  input:     
    mkr_dir + '/pb_' + FULL_TAG + f'_{MKR_SEL_RES}_' + DATE_STAMP + '.rds',
    mkr_dir + '/pb_marker_genes_' + FULL_TAG + f'_{MKR_SEL_RES}_' + DATE_STAMP + '.txt.gz',
    mkr_dir + '/pb_hvgs_' + FULL_TAG + f'_{MKR_SEL_RES}_' + DATE_STAMP + '.txt.gz',
    # fgsea outputs
    fgsea_outs, 
    rmd_dir   + '/' + SHORT_TAG + f'_marker_genes_{MKR_SEL_RES}.Rmd',
    docs_dir  + '/' + SHORT_TAG + f'_marker_genes_{MKR_SEL_RES}.html'


rule label_celltypes:
  input:
    lbl_dir + '/hvg_mat_for_labelling_' + LBL_GENE_VAR + '_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    lbl_dir + '/cell_annotations_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    code_dir  + '/label_celltypes.R',
    rmd_dir   + '/' + SHORT_TAG + '_label_celltypes.Rmd', 
    docs_dir  + '/' + SHORT_TAG + '_label_celltypes.html'


# define rules that are needed
include: "mapping.smk"
include: "ambient.smk"
include: "make_hto_sce.smk"
include: "qc.smk"
include: "pb_empties.smk"
include: "hvgs.smk"
include: "integration.smk"
include: "marker_genes.smk"
include: "render_htmls.smk"
include: "label_celltypes.smk"
