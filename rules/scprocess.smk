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

lm_f          = scprocess_dir / "resources/resources_lm_params_2025-11-27.csv"

# check config
config        = check_config(config, schema_f, scdata_dir, scprocess_dir)

# get lists of parameters
RUN_PARAMS, RUN_VAR     = get_run_parameters(config, scdata_dir)
RUNS                    = list(RUN_PARAMS.keys())
BATCH_PARAMS, BATCH_VAR = get_batch_parameters(config, RUNS, scdata_dir)
BATCHES                 = list(BATCH_PARAMS.keys())
RUNS_TO_BATCHES         = get_runs_to_batches(config, RUNS, BATCHES, BATCH_VAR)
LABELLER_PARAMS         = get_labeller_parameters(config, schema_f, scdata_dir)

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
  #code_dir  + '/qc.R', 
  code_dir  + '/hvgs.R', 
  code_dir  + '/integration.R', 
  code_dir  + '/marker_genes.R'
  ]
if config['multiplexing']['demux_type'] == "hto":
  r_scripts.append(code_dir + '/multiplexing.R')

# alevin hto index outputs (optional)
hto_index_outs = [
  f'{af_dir}/hto.tsv',
  f'{af_dir}/t2g_hto.tsv',
  f'{af_dir}/hto_index/ref_indexing.log'
  ] if config['multiplexing']['demux_type'] == "hto" else []

# alevin hto quantification outputs (optional)
hto_af_outs = expand(
  [
<<<<<<< HEAD
    f'{af_dir}/af_{{run}}/hto/af_quant/',
    f'{af_dir}/af_{{run}}/hto/af_quant/alevin/quants_mat.mtx',
    f'{af_dir}/af_{{run}}/hto/af_quant/alevin/quants_mat_cols.txt',
    f'{af_dir}/af_{{run}}/hto/af_quant/alevin/quants_mat_rows.txt',
    f'{af_dir}/af_{{run}}/hto/af_hto_counts_mat.h5',
    f'{af_dir}/af_{{run}}/hto/knee_plot_data_{{run}}_{DATE_STAMP}.txt.gz'
  ], run = RUNS) if config['multiplexing']['demux_type'] == "hto" else []
=======
  af_dir    + '/af_{run}/hto/af_quant/',
  af_dir    + '/af_{run}/hto/af_quant/alevin/quants_mat.mtx',
  af_dir    + '/af_{run}/hto/af_quant/alevin/quants_mat_cols.txt',
  af_dir    + '/af_{run}/hto/af_quant/alevin/quants_mat_rows.txt',
  af_dir    + '/af_{run}/hto/af_hto_counts_mat.h5',
  af_dir    + '/af_{run}/hto/knee_plot_data_{run}_' + DATE_STAMP + '.csv.gz'
  ], run = RUNS
) if config['multiplexing']['demux_type'] == "hto" else []
>>>>>>> refs/remotes/origin/dev-int

# seurat demultiplexing outputs (optional)
hto_sce_fs = expand(
  f'{demux_dir}/sce_cells_htos_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds',
  run = RUNS
  ) if config['multiplexing']['demux_type'] == "hto" else []

# multiplexing report (optional)
hto_rmd_f  = (f'{rmd_dir}/{SHORT_TAG}_demultiplexing.Rmd') if config['multiplexing']['demux_type'] == "hto" else []
hto_html_f = (f'{docs_dir}/{SHORT_TAG}_demultiplexing.html') if config['multiplexing']['demux_type'] == "hto" else []

# fgsea outputs (optional)
fgsea_outs = [
  f'{mkr_dir}/fgsea_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_go_bp_{DATE_STAMP}.txt.gz',
  f'{mkr_dir}/fgsea_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_go_cc_{DATE_STAMP}.txt.gz',
  f'{mkr_dir}/fgsea_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_go_mf_{DATE_STAMP}.txt.gz'
] if (config['project']['species'] in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020']) & config['marker_genes']['mkr_do_gsea'] else []

# one rule to rule them all
rule all:
  input:
    # hto outputs
    hto_index_outs,
    hto_af_outs,
    expand(
      [
      # mapping
\      f'{af_dir}/af_{{run}}/{af_rna_dir}af_quant/',
      f'{af_dir}/af_{{run}}/{af_rna_dir}af_quant/alevin/quants_mat.mtx',
      f'{af_dir}/af_{{run}}/{af_rna_dir}af_quant/alevin/quants_mat_cols.txt',
      f'{af_dir}/af_{{run}}/{af_rna_dir}af_quant/alevin/quants_mat_rows.txt',
      f'{af_dir}/af_{{run}}/{af_rna_dir}af_counts_mat.h5',
      f'{af_dir}/af_{{run}}/{af_rna_dir}knee_plot_data_{{run}}_{DATE_STAMP}.txt.gz',
      f'{af_dir}/af_{{run}}/{af_rna_dir}ambient_params_{{run}}_{DATE_STAMP}.yaml',
      # ambient (cellbender, decontx or nothing)
      f'{amb_dir}/ambient_{{run}}/ambient_{{run}}_{DATE_STAMP}_output_paths.yaml',
      f'{amb_dir}/ambient_{{run}}/barcodes_qc_metrics_{{run}}_{DATE_STAMP}.txt.gz',
      # doublet id
      f'{dbl_dir}/dbl_{{run}}/scDblFinder_{{run}}_outputs_{FULL_TAG}_{DATE_STAMP}.csv.gz'
      ], run =  RUNS),
    # ambient sample statistics
    f'{amb_dir}/ambient_run_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    f'{amb_dir}/paths_h5_filtered_{FULL_TAG}_{DATE_STAMP}.yaml',
    # demultiplexing
    hto_sce_fs,
    # qc
    f'{qc_dir}/qc_thresholds_by_{BATCH_VAR}_{FULL_TAG}_{DATE_STAMP}.csv',
    f'{qc_dir}/qc_all_samples_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    # pseudobulks and empties
    f'{pb_dir}/af_paths_{FULL_TAG}_{DATE_STAMP}.csv', 
    f'{pb_dir}/pb_empties_{FULL_TAG}_{DATE_STAMP}.rds', 
    f'{pb_dir}/pb_cells_all_{FULL_TAG}_{DATE_STAMP}.rds',
    f'{empty_dir}/edger_empty_genes_all_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    # hvgs
    f'{hvg_dir}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv',
    f'{hvg_dir}/standardized_variance_stats_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    f'{hvg_dir}/hvg_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    f'{hvg_dir}/top_hvgs_counts_{FULL_TAG}_{DATE_STAMP}.h5', 
    f'{hvg_dir}/top_hvgs_doublet_counts_{FULL_TAG}_{DATE_STAMP}.h5',
    #f'{hvg_dir}/chunked_counts_h5_sizes_{FULL_TAG}_{DATE_STAMP}.csv', 
    # integration
    f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    expand(f'{int_dir}/sce_cells_clean_{{batch}}_{FULL_TAG}_{DATE_STAMP}.rds', batch = BATCHES),
    f'{int_dir}/sce_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml', 
    # marker genes
    f'{mkr_dir}/pb_{FULL_TAG}_{ config['marker_genes']['mkr_sel_res'] }_{DATE_STAMP}.rds',
    f'{mkr_dir}/pb_marker_genes_{FULL_TAG}_{ config['marker_genes']['mkr_sel_res'] }_{DATE_STAMP}.txt.gz',
    f'{mkr_dir}/pb_hvgs_{FULL_TAG}_{ config['marker_genes']['mkr_sel_res'] }_{DATE_STAMP}.txt.gz',
    # fgsea outputs
    fgsea_outs, 
    # code
    r_scripts,
    # markdowns
    f'{rmd_dir}/{SHORT_TAG}_mapping.Rmd',
    f'{rmd_dir}/{SHORT_TAG}_ambient.Rmd',
    f'{rmd_dir}/{SHORT_TAG}_qc.Rmd', 
    f'{rmd_dir}/{SHORT_TAG}_hvgs.Rmd',
    f'{rmd_dir}/{SHORT_TAG}_integration.Rmd', 
    f'{rmd_dir}/{SHORT_TAG}_marker_genes_{config['marker_genes']['mkr_sel_res']}.Rmd', 
    hto_rmd_f,
    # reports
    f'{docs_dir}/{SHORT_TAG}_mapping.html', 
    f'{docs_dir}/{SHORT_TAG}_ambient.html',
    f'{docs_dir}/{SHORT_TAG}_qc.html',
    f'{docs_dir}/{SHORT_TAG}_hvgs.html',
    f'{docs_dir}/{SHORT_TAG}_integration.html',
    f'{docs_dir}/{SHORT_TAG}_marker_genes_{config['marker_genes']['mkr_sel_res']}.html',
    hto_html_f 


rule mapping:
  input:
    hto_index_outs, 
    hto_af_outs,
    expand([
      f'{af_dir}/af_{{run}}/{af_rna_dir}af_quant/',
      f'{af_dir}/af_{{run}}/{af_rna_dir}af_quant/alevin/quants_mat.mtx',
      f'{af_dir}/af_{{run}}/{af_rna_dir}af_quant/alevin/quants_mat_cols.txt',
      f'{af_dir}/af_{{run}}/{af_rna_dir}af_quant/alevin/quants_mat_rows.txt',
      f'{af_dir}/af_{{run}}/{af_rna_dir}af_counts_mat.h5',
      f'{af_dir}/af_{{run}}/{af_rna_dir}knee_plot_data_{{run}}_{DATE_STAMP}.txt.gz',
      f'{af_dir}/af_{{run}}/{af_rna_dir}ambient_params_{{run}}_{DATE_STAMP}.yaml'
      ], run = RUNS),
    f'{rmd_dir}/{SHORT_TAG}_mapping.Rmd',
    f'{docs_dir}/{SHORT_TAG}_mapping.html'


rule demux:
  input: 
    hto_sce_fs,
    hto_rmd_f, 
    hto_html_f


rule ambient:
  input: 
    expand(f'{amb_dir}/ambient_{{run}}/ambient_{{run}}_{DATE_STAMP}_output_paths.yaml', run = RUNS), 
    expand(f'{amb_dir}/ambient_{{run}}/barcodes_qc_metrics_{{run}}_{DATE_STAMP}.txt.gz', run = RUNS),
    f'{amb_dir}/ambient_run_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    f'{amb_dir}/paths_h5_filtered_{FULL_TAG}_{DATE_STAMP}.yaml',
    f'{rmd_dir}/{SHORT_TAG}_mapping.Rmd',
    f'{docs_dir}/{SHORT_TAG}_mapping.html',
    f'{rmd_dir}/{SHORT_TAG}_ambient.Rmd',
    f'{docs_dir}/{SHORT_TAG}_ambient.html'


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
      f'{dbl_dir}/dbl_{{run}}/scDblFinder_{{run}}_outputs_{FULL_TAG}_{DATE_STAMP}.csv.gz'
      ], run = RUNS),
    f'{qc_dir}/qc_thresholds_by_{BATCH_VAR}_{FULL_TAG}_{DATE_STAMP}.csv',
    f'{qc_dir}/qc_all_samples_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv', 
    f'{rmd_dir}/{SHORT_TAG}_mapping.Rmd',
    f'{docs_dir}/{SHORT_TAG}_mapping.html',
    f'{rmd_dir}/{SHORT_TAG}_ambient.Rmd',
    f'{docs_dir}/{SHORT_TAG}_ambient.html',
    f'{rmd_dir}/{SHORT_TAG}_qc.Rmd',
    f'{docs_dir}/{SHORT_TAG}_qc.html'


rule hvg:
  params:
    hvg_method  = config['hvg']['hvg_method'],
    n_hvgs      = config['hvg']['hvg_n_hvgs'],
    exclude_ambient_genes = config['hvg']['hvg_exclude_ambient_genes']
  input:
    f'{hvg_dir}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv',
    f'{hvg_dir}/standardized_variance_stats_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    f'{hvg_dir}/hvg_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    f'{hvg_dir}/top_hvgs_counts_{FULL_TAG}_{DATE_STAMP}.h5', 
    f'{hvg_dir}/top_hvgs_doublet_counts_{FULL_TAG}_{DATE_STAMP}.h5',
    f'{rmd_dir}/{SHORT_TAG}_mapping.Rmd',
    f'{docs_dir}/{SHORT_TAG}_mapping.html',
    f'{rmd_dir}/{SHORT_TAG}_ambient.Rmd',
    f'{docs_dir}/{SHORT_TAG}_ambient.html',
    f'{rmd_dir}/{SHORT_TAG}_qc.Rmd',
    f'{docs_dir}/{SHORT_TAG}_qc.html',
    f'{rmd_dir}/{SHORT_TAG}_hvgs.Rmd',
    f'{docs_dir}/{SHORT_TAG}_hvgs.html'


rule integration:
  input:
    expand(f'{int_dir}/sce_cells_clean_{{batch}}_{FULL_TAG}_{DATE_STAMP}.rds', batch = BATCHES),
    f'{int_dir}/sce_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml', 
    f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    f'{rmd_dir}/{SHORT_TAG}_mapping.Rmd',
    f'{docs_dir}/{SHORT_TAG}_mapping.html',
    f'{rmd_dir}/{SHORT_TAG}_ambient.Rmd',
    f'{docs_dir}/{SHORT_TAG}_ambient.html',
    f'{rmd_dir}/{SHORT_TAG}_qc.Rmd',
    f'{docs_dir}/{SHORT_TAG}_qc.html',
    f'{rmd_dir}/{SHORT_TAG}_hvgs.Rmd',
    f'{docs_dir}/{SHORT_TAG}_hvgs.html',
    f'{rmd_dir}/{SHORT_TAG}_integration.Rmd',
    f'{docs_dir}/{SHORT_TAG}_integration.html'


rule marker_genes:
  input:     
    f'{mkr_dir}/pb_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}.rds',
    f'{mkr_dir}/pb_marker_genes_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}.txt.gz',
    f'{mkr_dir}/pb_hvgs_{FULL_TAG}_{config['marker_genes']['mkr_sel_res']}_{DATE_STAMP}.txt.gz',
    # fgsea outputs
    fgsea_outs, 
    f'{rmd_dir}/{SHORT_TAG}_marker_genes_{config['marker_genes']['mkr_sel_res']}.Rmd',
    f'{docs_dir}/{SHORT_TAG}_marker_genes_{config['marker_genes']['mkr_sel_res']}.html'


rule label_celltypes:
  input:
    expand(f'{lbl_dir}/labels_{{labeller}}_model_{{model}}_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
      zip, 
        labeller  = [ entry['labeller'] for entry in LABELLER_PARAMS],
        model     = [ entry['model']    for entry in LABELLER_PARAMS]
      ),
    code_dir  + '/label_celltypes.R',
    f'{rmd_dir}/{SHORT_TAG}_label_celltypes.Rmd', 
    f'{docs_dir}/{SHORT_TAG}_label_celltypes.html'

rule index:
  input:
    f'{docs_dir}/index.html'


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
include: "label_celltypes.smk"

