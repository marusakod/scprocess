# load modules
import os
import sys
import re
import glob
import pandas as pd
import warnings
import yaml
import json
import pathlib
sys.path.append('scripts')
from scprocess_utils import *

# define some things
scprocess_dir = pathlib.Path(config.pop('scprocess_dir'))
schema_f      = scprocess_dir / "resources/schemas/config.schema.json"
scdata_dir    = pathlib.Path(os.getenv('SCPROCESS_DATA_DIR'))

# check config
config        = check_config(config, schema_f, scdata_dir, scprocess_dir)

# get lists of parameters
RUN_PARAMS, RUN_VAR = get_run_parameters(config, scdata_dir)
RUNS                = list(RUN_PARAMS.keys())
SAMPLE_PARAMS       = get_sample_parameters(config, RUNS, scdata_dir)
SAMPLES             = list(SAMPLE_PARAMS.keys())
SAMPLES_TO_RUNS     = get_samples_to_runs(config, RUNS, SAMPLES)

# unpack some variables that we use a lot
PROJ_DIR        = config['project']['proj_dir']
FULL_TAG        = config['project']['full_tag']
SHORT_TAG       = config['project']['short_tag']
DATE_STAMP      = config['project']['date_stamp']

# specify locations
benchmark_dir = f"{PROJ_DIR}/.resources"
code_dir      = f"{PROJ_DIR}/code"
af_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_mapping"
af_rna_dir    = 'rna/' if config['multiplexing']['demux_type'] == "hto" else ''
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
if config['multiplexing']['demux_type'] == "hto":
  r_scripts.append(code_dir + '/multiplexing.R')

# alevin hto index outputs (optional)
hto_index_outs = [
    af_dir + '/hto.tsv',
    af_dir + '/t2g_hto.tsv',
    af_dir + '/hto_index/ref_indexing.log'
  ] if config['multiplexing']['demux_type'] == "hto" else []

# alevin hto quantification outputs (optional)
hto_af_outs = expand(
  [
  af_dir    + '/af_{run}/hto/af_quant/',
  af_dir    + '/af_{run}/hto/af_quant/alevin/quants_mat.mtx',
  af_dir    + '/af_{run}/hto/af_quant/alevin/quants_mat_cols.txt',
  af_dir    + '/af_{run}/hto/af_quant/alevin/quants_mat_rows.txt',
  af_dir    + '/af_{run}/hto/af_hto_counts_mat.h5',
  af_dir    + '/af_{run}/hto/knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz'
  ], run = RUNS
) if config['multiplexing']['demux_type'] == "hto" else []

# seurat demultiplexing outputs (optional)
hto_sce_fs = expand(
  demux_dir + '/sce_cells_htos_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
  run = RUNS
  ) if config['multiplexing']['demux_type'] == "hto" else []

# multiplexing report (optional)
hto_rmd_f  = (rmd_dir   + '/' + SHORT_TAG + '_demultiplexing.Rmd') if config['multiplexing']['demux_type'] == "hto" else []
hto_html_f = (docs_dir  + '/' + SHORT_TAG + '_demultiplexing.html') if config['multiplexing']['demux_type'] == "hto" else []

# fgsea outputs (optional)
fgsea_outs = [
    mkr_dir   + '/fgsea_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + 'go_bp_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + 'go_cc_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + 'go_mf_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + 'paths_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + 'hlmk_' + DATE_STAMP + '.txt.gz', 
] if config['project']['species'] in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020'] else []

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
      dbl_dir   + '/dbl_{run}/scDblFinder_{run}_outputs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
      # dbl_dir   + '/dbl_{run}/scDblFinder_{run}_dimreds_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
      ], run =  RUNS),
    # ambient sample statistics
    amb_dir + '/ambient_run_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    # demultiplexing
    hto_sce_fs,
    # qc
    qc_dir  + '/qc_thresholds_by_sample_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    qc_dir  + '/qc_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', 
    qc_dir  + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', 
    qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', 
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
    mkr_dir + '/pb_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.rds',
    mkr_dir + '/pb_marker_genes_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.txt.gz',
    mkr_dir + '/pb_hvgs_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.txt.gz',
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
    rmd_dir   + '/' + SHORT_TAG + f'_marker_genes_{config['marker_genes']['mkr_sel_res']}.Rmd', 
    hto_rmd_f,
    # reports
    docs_dir  + '/' + SHORT_TAG + '_mapping.html', 
    docs_dir  + '/' + SHORT_TAG + '_ambient.html', 
    docs_dir  + '/' + SHORT_TAG + '_qc.html',
    docs_dir  + '/' + SHORT_TAG + '_hvgs.html',
    docs_dir  + '/' + SHORT_TAG + '_integration.html',
    docs_dir  + '/' + SHORT_TAG + f'_marker_genes_{config['marker_genes']['mkr_sel_res']}.html',
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
      run = RUNS), 
    rmd_dir   + '/' + SHORT_TAG + '_mapping.Rmd',
    docs_dir  + '/' + SHORT_TAG + '_mapping.html'


rule demux:
  input: 
    hto_sce_fs,
    hto_rmd_f, 
    hto_html_f


rule ambient:
  input: 
    expand(amb_dir   + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml', run = RUNS), 
    expand(amb_dir + '/ambient_{run}/barcodes_qc_metrics_{run}_' + DATE_STAMP + '.txt.gz', run = RUNS),
    amb_dir + '/ambient_run_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    ambient_rmd_f, 
    ambient_html_f


rule qc:
  params:
    mito_str        = config['mapping']['af_mito_str'],
    exclude_mito    = config['qc']['exclude_mito'],
    hard_min_counts = config['qc']['qc_hard_min_counts'],
    hard_min_feats  = config['qc']['qc_hard_min_feats'],
    hard_max_mito   = config['qc']['qc_hard_max_mito'],
    min_counts      = config['qc']['qc_min_counts'],
    min_feats       = config['qc']['qc_min_feats'],
    min_mito        = config['qc']['qc_min_mito'],
    max_mito        = config['qc']['qc_max_mito'],
    min_splice      = config['qc']['qc_min_splice'],
    max_splice      = config['qc']['qc_max_splice'],
    run_var         = RUN_VAR,
    demux_type      = config['multiplexing']['demux_type'],
    dbl_min_feats   = config['qc']['dbl_min_feats']
  input:     
    expand([
      # doublet id
      dbl_dir + '/dbl_{run}/scDblFinder_{run}_outputs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
      # dbl_dir + '/dbl_{run}/scDblFinder_{run}_dimreds_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
      ], run = RUNS),
    qc_dir    + '/qc_thresholds_by_sample_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    qc_dir    + '/qc_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', 
    qc_dir    + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', 
    qc_dir    + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', 
    qc_dir    + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    rmd_dir   + '/' + SHORT_TAG + '_qc.Rmd',
    docs_dir  + '/' + SHORT_TAG + '_qc.html'


rule hvg:
  params:
    hvg_method  = config['hvg']['hvg_method'],
    n_hvgs      = config['hvg']['hvg_n_hvgs'],
    exclude_ambient_genes = config['hvg']['hvg_exclude_ambient_genes']
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
    mkr_dir + '/pb_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.rds',
    mkr_dir + '/pb_marker_genes_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.txt.gz',
    mkr_dir + '/pb_hvgs_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.txt.gz',
    # fgsea outputs
    fgsea_outs, 
    rmd_dir   + '/' + SHORT_TAG + f'_marker_genes_{config['marker_genes']['mkr_sel_res']}.Rmd',
    docs_dir  + '/' + SHORT_TAG + f'_marker_genes_{config['marker_genes']['mkr_sel_res']}.html'


# rule label_celltypes:
#   input:
#     lbl_dir + '/hvg_mat_for_labelling_' + 'gene_id' + '_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
#     lbl_dir + '/cell_annotations_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
#     code_dir  + '/label_celltypes.R',
#     rmd_dir   + '/' + SHORT_TAG + '_label_celltypes.Rmd', 
#     docs_dir  + '/' + SHORT_TAG + '_label_celltypes.html'


# define rules that are needed
include: "mapping.smk"
if config['multiplexing']['demux_type'] == "hto":
  include: "hto.smk"
include: "ambient.smk"
include: "qc.smk"
include: "pb_empties.smk"
include: "hvgs.smk"
include: "integration.smk"
include: "marker_genes.smk"
include: "render_htmls.smk"
# include: "label_celltypes.smk"
