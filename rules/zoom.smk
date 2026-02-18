# load modules
import os
import sys
import re
import glob
import pandas as pd
import polars as pl
import yaml
import warnings
from snakemake.utils import validate, min_version

# import utils
sys.path.append('scripts')
from scprocess_utils import *
from zoom import *


# define some things
scprocess_dir   = pathlib.Path(config.pop('scprocess_dir'))
proj_schema_f   = scprocess_dir / "resources/schemas/config.schema.json"
zoom_schema_f   = scprocess_dir / "resources/schemas/zoom.schema.json"
scdata_dir      = pathlib.Path(os.getenv('SCPROCESS_DATA_DIR'))
lm_f            = scprocess_dir / "resources/snakemake/resources_lm_params_2025-12-16.csv"

# check config
config          = check_config(config, proj_schema_f, scdata_dir, scprocess_dir)

# get lists of parameters
RUN_PARAMS, RUN_VAR = get_run_parameters(config, scdata_dir)
RUNS                = list(RUN_PARAMS.keys())
BATCH_PARAMS, BATCH_VAR, SAMPLES = get_batch_parameters(config, RUNS, scdata_dir)
BATCHES             = list(BATCH_PARAMS.keys())
RUNS_TO_BATCHES, RUNS_TO_SAMPLES = get_runs_to_batches(config, RUNS, BATCHES, BATCH_VAR)
RESOURCE_PARAMS     = prep_resource_params(config, proj_schema_f, lm_f, RUN_PARAMS, BATCHES)

# get zoom parameters
ZOOM_PARAMS         = get_zoom_parameters(config, zoom_schema_f, scdata_dir)
ZOOMS               = list(ZOOM_PARAMS.keys())

# unpack some variables that we use a lot
PROJ_DIR        = config['project']['proj_dir']
FULL_TAG        = config['project']['full_tag']
SHORT_TAG       = config['project']['short_tag']
DATE_STAMP      = config['project']['date_stamp']

# specify locations
benchmark_dir = f"{PROJ_DIR}/.resources"
code_dir      = f"{PROJ_DIR}/code"
amb_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_ambient"
qc_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_qc"
int_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_integration"
pb_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_pseudobulk"
hvg_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_hvg"
rmd_dir       = f"{PROJ_DIR}/analysis"
docs_dir      = f"{PROJ_DIR}/public"
zoom_dir      = f"{PROJ_DIR}/output/{SHORT_TAG}_zoom"

# define zoom marker outputs
zoom_mkr_report_outs = [
  file
  for zoom_name, mkr_sel_res, do_gsea in zip(
    ZOOMS,
    [ZOOM_PARAMS[zoom]['marker_genes']['mkr_sel_res'] for zoom in ZOOMS],
    [ZOOM_PARAMS[zoom]['marker_genes']['mkr_do_gsea'] for zoom in ZOOMS]
  )
  for file in (
    [
      '%s/%s/pb_%s_%s_%s.rds' % (zoom_dir, zoom_name, FULL_TAG, mkr_sel_res, DATE_STAMP),
      '%s/%s/pb_marker_genes_%s_%s_%s.csv.gz' % (zoom_dir, zoom_name, FULL_TAG, mkr_sel_res, DATE_STAMP),
      '%s/%s/pb_hvgs_%s_%s_%s.csv.gz' % (zoom_dir, zoom_name, FULL_TAG, mkr_sel_res, DATE_STAMP),
      '%s/%s_zoom_%s_%s.Rmd' % (rmd_dir, SHORT_TAG, zoom_name, mkr_sel_res),
      '%s/%s_zoom_%s_%s.html' % (docs_dir, SHORT_TAG, zoom_name, mkr_sel_res)
    ] + (
      [
      '%s/%s/fgsea_%s_%s_go_bp_%s.csv.gz' % (zoom_dir, zoom_name, FULL_TAG, mkr_sel_res, DATE_STAMP),
      '%s/%s/fgsea_%s_%s_go_cc_%s.csv.gz' % (zoom_dir, zoom_name, FULL_TAG, mkr_sel_res, DATE_STAMP),
      '%s/%s/fgsea_%s_%s_go_mf_%s.csv.gz' % (zoom_dir, zoom_name, FULL_TAG, mkr_sel_res, DATE_STAMP)
      ] if do_gsea and (config['project']['ref_txome'] in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020'])
        else []
      )
  )
]

zooms_to_save_sce     = [ zoom_name for zoom_name in ZOOMS if ZOOM_PARAMS[zoom_name]['zoom']['save_subset_sces']]
zooms_to_save_anndata = [ zoom_name for zoom_name in ZOOMS if ZOOM_PARAMS[zoom_name]['zoom']['save_subset_anndata']]

ZOOM_OUT_MAP = {}
for name in ZOOMS:
  ZOOM_OUT_MAP[name] = {}
  for b in BATCHES:
    ZOOM_OUT_MAP[name][b] = {}
    if name in zooms_to_save_sce:
      ZOOM_OUT_MAP[name][b]["sce"]   = f"{zoom_dir}/{name}/sce_objects/sce_cells_clean_{name}_{b}_{FULL_TAG}_{DATE_STAMP}.rds"
    if name in zooms_to_save_anndata:
      ZOOM_OUT_MAP[name][b]["adata"] = f"{zoom_dir}/{name}/anndata_objects/anndata_cells_clean_{name}_{b}_{FULL_TAG}_{DATE_STAMP}.h5ad"

zoom_all_subset_fs = [
  path 
  for name_dict in ZOOM_OUT_MAP.values() 
  for b_dict in name_dict.values() 
  for path in b_dict.values()
]

rule zoom:
  input:
    # zoom sample qc
    expand('%s/{zoom_name}/zoom_%s_statistics_%s_%s.csv' % \
      (zoom_dir, BATCH_VAR, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS),
    # zoom pseudobulks and empties
    expand('%s/{zoom_name}/pb_cells_{zoom_name}_%s_%s.rds' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS),
    expand('%s/{zoom_name}/edger_empty_genes_%s_%s.csv.gz' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS), 
    # zoom hvgs
    expand('%s/{zoom_name}/hvg_paths_%s_%s.csv' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS),  
    expand('%s/{zoom_name}/standardized_variance_stats_%s_%s.csv.gz' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS), 
    expand('%s/{zoom_name}/hvg_dt_%s_%s.csv.gz' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS), 
    expand('%s/{zoom_name}/top_hvgs_counts_%s_%s.h5' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS), 
    # zoom integration
    expand('%s/{zoom_name}/integrated_dt_%s_%s.csv.gz' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS),
    # zoom marker genes, fgsea and html report
    zoom_mkr_report_outs, 
    # zoom sce and anndata subsets (optional)
    zoom_all_subset_fs


localrules: zoom_make_tmp_pb_cells_df, zoom_make_hvg_df, zoom_merge_group_mean_var, zoom_merge_group_std_var_stats, zoom_merge_stats_for_std_variance

rule get_zoom_sample_statistics:
  input:
    qc_stats_f      = f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv'
  output:
    zoom_stats_f    = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv'
  params:
    zoom_lbls_f     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_f'],
    zoom_lbls_col   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'], 
    zoom_lbls       = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels'],
    batch_var       = BATCH_VAR,
    batches         = BATCHES,
    zoom_min_n_smpl = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc']['qc_min_cells'],
    ambient_method  = config['ambient']['ambient_method'],
  run:
    zoom_stats_df   = extract_zoom_sample_statistics(input.qc_stats_f, params.zoom_lbls_f, 
      params.zoom_lbls_col, params.zoom_lbls, params.batches, params.batch_var, 
      params.zoom_min_n_smpl, params.ambient_method)
    zoom_stats_df.write_csv(output.zoom_stats_f)


# pseudobulks and empties
rule zoom_make_one_pb_cells:
  input:
    batch_lu_f  = f'{pb_dir}/runs_to_batches_{FULL_TAG}_{DATE_STAMP}.csv',
    qc_stats_f  = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    h5_paths_f  = f'{hvg_dir}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv',
    coldata_f   = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    pb_cells_f  = temp(f'{zoom_dir}/{{zoom_name}}/tmp_pb_cells_{{zoom_name}}_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds')
  params:
    run_var       = RUN_VAR,
    batch_var     = BATCH_VAR, 
    zoom_lbls_f   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_f'],
    zoom_lbls_col = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'],
    zoom_lbls     = lambda wildcards: ','.join(ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels'])
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_make_one_pb_cells', 'memory', attempt, wildcards.run),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_make_one_pb_cells', 'time', attempt, wildcards.run)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_zoom/zoom_make_one_pb_cells_{{zoom_name}}_{{run}}_{DATE_STAMP}.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    make_pb_cells(
      sel_run     = '{wildcards.run}',
      batch_lu_f  = '{input.batch_lu_f}',
      qc_stats_f  = '{input.qc_stats_f}',
      h5_paths_f  = '{input.h5_paths_f}', 
      coldata_f   = '{input.coldata_f}',
      run_var     = '{params.run_var}',
      batch_var   = '{params.batch_var}',
      subset_f    = '{params.zoom_lbls_f}',
      subset_col  = '{params.zoom_lbls_col}',
      subset_str  = '{params.zoom_lbls}', 
      pb_cells_f  = '{output.pb_cells_f}'
    )"
    """

rule zoom_make_tmp_pb_cells_df:
  input:
    pb_cells_fs = lambda wildcards: expand(
        f'{zoom_dir}/{{zoom_name}}/tmp_pb_cells_{{zoom_name}}_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds',
        run=RUNS,
        zoom_name=[wildcards.zoom_name]  
    )
  output:
    cells_paths_f = temp(f'{zoom_dir}/{{zoom_name}}/tmp_pb_cells_paths_{FULL_TAG}_{DATE_STAMP}.csv')
  params:
    run_var = RUN_VAR,
    runs    = RUNS
  run:
    import os
    import polars as pl

    paths_df = pl.DataFrame({
        params.run_var: params.runs,
        "pb_path": input.pb_cells_fs
    })
    paths_df = paths_df.filter(
        pl.col("pb_path").map_elements(os.path.getsize, return_dtype=pl.Int64) > 0
    )
    paths_df.write_csv(output.cells_paths_f)


rule zoom_merge_pb_cells:
  input:
    cells_paths_f = f'{zoom_dir}/{{zoom_name}}/tmp_pb_cells_paths_{FULL_TAG}_{DATE_STAMP}.csv',
    pb_cells_fs   = lambda wildcards: expand(
        f'{zoom_dir}/{{zoom_name}}/tmp_pb_cells_{{zoom_name}}_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds',
        run=RUNS,
        zoom_name=[wildcards.zoom_name]  # Restrict to current zoom_name
    ),
    rowdata_f = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    pb_cells_f = f'{zoom_dir}/{{zoom_name}}/pb_cells_{{zoom_name}}_{FULL_TAG}_{DATE_STAMP}.rds'
  params:
    batch_var  = BATCH_VAR
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_merge_pb_cells', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_merge_pb_cells', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_zoom/zoom_merge_pb_cells_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    merge_pbs_cells( \
      cells_paths_f = '{input.cells_paths_f}', 
      rowdata_f     = '{input.rowdata_f}',
      batch_var     = '{params.batch_var}',
      pb_cells_f    = '{output.pb_cells_f}'
    )"
    """


rule zoom_calculate_ambient_genes:
  input:
    pb_empty_f      = f'{pb_dir}/pb_empties_{FULL_TAG}_{DATE_STAMP}.rds', 
    zoom_pb_f       = f'{zoom_dir}/{{zoom_name}}/pb_cells_{{zoom_name}}_{FULL_TAG}_{DATE_STAMP}.rds'
  output:
    zoom_empty_gs_f = f'{zoom_dir}/{{zoom_name}}/edger_empty_genes_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  params:
    zoom_fdr_thr    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['pb_empties']['ambient_genes_fdr_thr'],
    zoom_logfc_thr  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['pb_empties']['ambient_genes_logfc_thr']
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_calculate_ambient_genes', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_calculate_ambient_genes', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_zoom/zoom_calculate_ambient_genes_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    calc_empty_genes(
      pb_cells_f = '{input.zoom_pb_f}',
      pb_empty_f = '{input.pb_empty_f}',
      fdr_thr    = {params.zoom_fdr_thr}, 
      logfc_thr  = {params.zoom_logfc_thr},
      empty_gs_f = '{output.zoom_empty_gs_f}'
    )"
    """

# highly variable genes
rule zoom_make_hvg_df:
  input:
    amb_yaml_fs = expand([f'{amb_dir}/ambient_{{run}}/ambient_{{run}}_{DATE_STAMP}_output_paths.yaml'], run=RUNS)
  output:
    hvg_paths_f = f'{zoom_dir}/{{zoom_name}}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv'
  params:
    demux_type  = config['multiplexing']['demux_type'],
    run_var     = RUN_VAR,
    runs        = RUNS,
    mapping     = RUNS_TO_BATCHES,
    batch_var   = BATCH_VAR
  run:
    hvg_df = make_hvgs_input_df( 
      params.runs, input.amb_yaml_fs, params.run_var, params.batch_var, params.mapping,
      params.demux_type, FULL_TAG, DATE_STAMP, f"{zoom_dir}/{wildcards.zoom_name}"
      )
    hvg_df.write_csv(output.hvg_paths_f)


rule zoom_make_tmp_csr_matrix:
  input:
    hvg_paths_f     = f'{zoom_dir}/{{zoom_name}}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv', 
    smpl_stats_f    = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv', 
    rowdata_f       = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    clean_h5_f      = temp(expand([
      f'{zoom_dir}/{{zoom_name}}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5'
      ], batch = BATCHES, allow_missing = True))
  params: 
    zoom_lbls_f     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_f'],
    zoom_lbls_col   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'],
    zoom_lbls       = lambda wildcards: ','.join(ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels']), 
    run_var         = RUN_VAR,
    batch_var       = BATCH_VAR,
    demux_type      = config['multiplexing']['demux_type'],
    zoom_chunk_size = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_chunk_size'],
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_make_tmp_csr_matrix', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_make_tmp_csr_matrix', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_zoom/zoom_make_tmp_csr_matrix_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    python3 scripts/hvgs.py get_csr_counts \
      {input.hvg_paths_f} \
      {params.zoom_lbls_f} \
      "{params.zoom_lbls_col}" \
      "{input.smpl_stats_f}" \
      {input.rowdata_f} \
      {params.run_var} \
      {params.batch_var} \
      {params.demux_type} \
      --keep_vals_str "{params.zoom_lbls}" \
      --chunksize {params.zoom_chunk_size} \
      --ncores {threads}
    """


rule zoom_get_stats_for_std_variance_for_sample:
  input: 
    clean_h5_f    = f'{zoom_dir}/{{zoom_name}}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5', 
    smpl_stats_f  = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv', 
    rowdata_f     = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    std_var_stats_f = temp(zoom_dir + '/{zoom_name}' + '/tmp_std_var_stats_{batch}_sample_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz')
  params:
    batch_var = BATCH_VAR
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_stats_for_std_variance_for_sample', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_stats_for_std_variance_for_sample', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_zoom/zoom_get_stats_for_std_variance_for_sample_{{zoom_name}}_{{batch}}_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    python3 scripts/hvgs.py calculate_std_var_stats_for_sample \
      {wildcards.batch} \
      {params.batch_var} \
      {input.smpl_stats_f} \
      {input.clean_h5_f} \
      {input.rowdata_f} \
      {output.std_var_stats_f}
    """


rule zoom_get_mean_var_for_group:
  input:
    clean_h5_f    = expand(
      f'{zoom_dir}/{{zoom_name}}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5',
      batch = BATCHES, allow_missing = True),
    hvg_paths_f   = f'{zoom_dir}/{{zoom_name}}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv',
    rowdata_f     = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    smpl_stats_f  = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv'
  output: 
    mean_var_f    = temp(f'{zoom_dir}/{{zoom_name}}/tmp_mean_var_{{group}}_group_chunk_{{chunk}}_{FULL_TAG}_{DATE_STAMP}.csv.gz')
  params:
    metadata_f          = config['project']['sample_metadata'],
    batch_var           = BATCH_VAR,
    zoom_hvg_method     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_method'], 
    zoom_hvg_chunk_size = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_chunk_size'], 
    zoom_group_var      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_metadata_split_var']
  threads: 8
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_mean_var_for_group', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_mean_var_for_group', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_zoom/zoom_get_mean_var_for_group_{{zoom_name}}_{{group}}_chunk_{{chunk}}_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    GROUPVAR_FLAG=""
    if [ "{params.zoom_hvg_method}" = "groups" ]; then
      GROUPVAR_FLAG="--groupvar {params.zoom_group_var}"
    fi

    python3 scripts/hvgs.py calculate_mean_var_for_chunk \
      {input.hvg_paths_f} \
      {input.rowdata_f} \
      {params.metadata_f} \
      {input.smpl_stats_f} \
      {output.mean_var_f} \
      {wildcards.chunk} \
      {params.zoom_hvg_method} \
      {params.batch_var} \
      --chunksize {params.zoom_hvg_chunk_size} \
      --group {wildcards.group} \
      --ncores {threads} \
      $GROUPVAR_FLAG
    """


rule zoom_merge_group_mean_var:
  input:         
    mean_var_f    = lambda wildcards: get_zoom_raw_mean_var_files(wildcards.zoom_name, zoom_dir, ZOOM_PARAMS, FULL_TAG, DATE_STAMP)
  output:
    mean_var_merged_f = temp(f'{zoom_dir}/{{zoom_name}}/means_variances_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz')
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_merge_group_mean_var', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_merge_group_mean_var', 'time', attempt)
  run:
    merge_tmp_files(input.mean_var_f, output.mean_var_merged_f)


rule zoom_get_estimated_variances:
  input:
    mean_var_merged_f = f'{zoom_dir}/{{zoom_name}}/means_variances_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    estim_vars_f      = temp(f'{zoom_dir}/{{zoom_name}}/estimated_variances_{FULL_TAG}_{DATE_STAMP}.csv.gz')
  params: 
    zoom_hvg_method   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_method'],
    batch_var         = BATCH_VAR
  threads: 1
  retries: config['resources']['retries']
  conda:
    '../envs/hvgs.yaml'
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_estimated_variances', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_estimated_variances', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_zoom/zoom_get_estimated_variances_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  shell: """
    python3 scripts/hvgs.py calculate_estimated_vars \
      {output.estim_vars_f} \
      {params.zoom_hvg_method} \
      {params.batch_var} \
      {input.mean_var_merged_f}
    """


rule zoom_get_stats_for_std_variance_for_group:
  input: 
    clean_h5_fs   = expand(
      f'{zoom_dir}/{{zoom_name}}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5',
      batch = BATCHES, allow_missing = True
    ),
    estim_vars_f  = f'{zoom_dir}/{{zoom_name}}/estimated_variances_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    hvg_paths_f   = f'{zoom_dir}/{{zoom_name}}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv',
    rowdata_f     = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    smpl_stats_f  = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv'
  output:
    std_var_stats_f = temp(f'{zoom_dir}/{{zoom_name}}/tmp_std_var_stats_{{group}}_group_chunk_{{chunk}}_{FULL_TAG}_{DATE_STAMP}.csv.gz')
  params:
    metadata_f          = config['project']['sample_metadata'],
    batch_var       = BATCH_VAR,
    zoom_hvg_method     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_method'], 
    zoom_hvg_chunk_size = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_chunk_size'], 
    zoom_hvg_group_var  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_metadata_split_var']
  threads: 8
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_stats_for_std_variance_for_group', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_stats_for_std_variance_for_group', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_zoom/zoom_get_stats_for_std_variance_for_group_{{zoom_name}}_{{group}}_chunk_{{chunk}}_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    python3 scripts/hvgs.py calculate_std_var_stats_for_chunk \
      {input.hvg_paths_f} \
      {input.rowdata_f} \
      {params.metadata_f} \
      {input.smpl_stats_f} \
      {output.std_var_stats_f} \
      {input.estim_vars_f} \
      {wildcards.chunk} \
      {params.zoom_hvg_method} \
      {params.batch_var} \
      --chunksize {params.zoom_hvg_chunk_size} \
      --group {wildcards.group} \
      --groupvar {params.zoom_hvg_group_var} \
      --ncores {threads} 
    """


rule zoom_merge_stats_for_std_variance:
  input:
    tmp_std_var_stats_fs = lambda wildcards: get_zoom_std_var_stats_files(wildcards.zoom_name, \
      zoom_dir, ZOOM_PARAMS, FULL_TAG, DATE_STAMP, BATCHES)
  output:
    std_var_stats_merged_f= f'{zoom_dir}/{{zoom_name}}/standardized_variance_stats_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  run:
    merge_tmp_files(input.tmp_std_var_stats_fs, output.std_var_stats_merged_f)

        
rule zoom_get_highly_variable_genes:
  input:
    std_var_stats_f = f'{zoom_dir}/{{zoom_name}}/standardized_variance_stats_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    empty_gs_fs     = f'{zoom_dir}/{{zoom_name}}/edger_empty_genes_{FULL_TAG}_{DATE_STAMP}.csv.gz' 
  output:
    hvg_f =  f'{zoom_dir}/{{zoom_name}}/hvg_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  threads: 1
  retries: config['resources']['retries']
  params:
    zoom_hvg_method   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_method'],
    zoom_n_hvgs       = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_n_hvgs'],
    batch_var         = BATCH_VAR,
    zoom_exc_gs_f     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_exclude_from_file'],
    zoom_exc_ambient  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_exclude_ambient_genes']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_highly_variable_genes', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_highly_variable_genes', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_zoom/zoom_get_highly_variable_genes_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    NOAMBIENT_FLAG=""
    if [ "{params.zoom_exc_ambient}" = "True" ]; then
      NOAMBIENT_FLAG="--noambient"
    fi
    EXC_GS_F_FLAG=""
    if [ "{params.zoom_exc_gs_f}" != "None" ]; then
      EXC_GS_F_FLAG="--exc_gs_f {params.zoom_exc_gs_f}"
    fi

    python3 scripts/hvgs.py calculate_hvgs \
      {input.std_var_stats_f} \
      {output.hvg_f} \
      {input.empty_gs_fs} \
      {params.zoom_hvg_method} \
      {params.batch_var} \
      {params.zoom_n_hvgs} \
      $NOAMBIENT_FLAG \
      $EXC_GS_F_FLAG
    """


rule zoom_create_hvg_matrix:
  input: 
    clean_h5_f    = expand(
      f'{zoom_dir}/{{zoom_name}}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5',
      zoom_name = ZOOMS, batch = BATCHES
    ),
    smpl_stats_f  = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv', 
    hvg_paths_f   = f'{zoom_dir}/{{zoom_name}}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv', 
    hvg_f         = f'{zoom_dir}/{{zoom_name}}/hvg_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    hvg_mat_f     = f'{zoom_dir}/{{zoom_name}}/top_hvgs_counts_{FULL_TAG}_{DATE_STAMP}.h5'
  params:
    demux_type    = config['multiplexing']['demux_type'], 
    batch_var     = BATCH_VAR
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_create_hvg_matrix', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_create_hvg_matrix', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_zoom/zoom_create_hvg_matrix_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    python3 scripts/hvgs.py create_hvg_matrix \
      {input.smpl_stats_f} \
      {input.hvg_paths_f} \
      {input.hvg_f} \
      {output.hvg_mat_f} \
      {params.demux_type} \
      {params.batch_var} \
      --ncores {threads}
    """


rule zoom_run_integration:
  input:
    hvg_mat_f     = f'{zoom_dir}/{{zoom_name}}/top_hvgs_counts_{FULL_TAG}_{DATE_STAMP}.h5', 
    sample_qc_f   = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv', 
    coldata_f     = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    integration_f = f'{zoom_dir}/{{zoom_name}}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  params:
    demux_type            = config['multiplexing']['demux_type'],
    exclude_mito          = config['qc']['exclude_mito'],
    zoom_int_embedding    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_embedding'],
    zoom_int_n_dims       = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_n_dims'],
    zoom_int_cl_method    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_cl_method'],
    zoom_int_theta        = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_theta'],
    zoom_int_use_paga     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_use_paga'],
    zoom_int_paga_cl_res  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_paga_cl_res'],
    zoom_int_res_ls       = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_res_ls'],
    zoom_int_use_gpu      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_use_gpu'],
    batch_var             = BATCH_VAR,
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_run_integration', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_run_integration', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_zoom/zoom_run_integration_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  conda: 
    '../envs/integration.yaml'
  shell: """
    set +u
    # set use_gpu flag based on config and on whether available
    USE_GPU_FLAG=""
    if [ "{params.zoom_int_use_gpu}" == "True" && -n "$CUDA_VISIBLE_DEVICES" ]; then
       USE_GPU_FLAG="--use-gpu"
    fi
    set -u
    
    python3 scripts/integration.py run_zoom_integration \
      --hvg_mat_f     {input.hvg_mat_f} \
      --sample_qc_f   {input.sample_qc_f} \
      --coldata_f     {input.coldata_f} \
      --demux_type    {params.demux_type} \
      --exclude_mito  "{params.exclude_mito}" \
      --embedding     {params.zoom_int_embedding} \
      --n_dims        {params.zoom_int_n_dims} \
      --cl_method     {params.zoom_int_cl_method} \
      --theta         {params.zoom_int_theta} \
      --res_ls_concat "{params.zoom_int_res_ls}" \
      --integration_f {output.integration_f} \
      --batch_var     {params.batch_var} \
      $(if [ "{params.zoom_int_use_paga}" == "True" ]; then echo "--use-paga"; fi) \
      $(if [ "{params.zoom_int_use_paga}" == "True" ]; then echo "--paga-cl-res {params.zoom_int_paga_cl_res}"; fi) \
      $USE_GPU_FLAG
    """


rule zoom_run_marker_genes:
  input:
    h5ads_yaml_f   = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml',
    integration_f  = f'{zoom_dir}/{{zoom_name}}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    pb_f      = f'{zoom_dir}/{{zoom_name}}/pb_{FULL_TAG}_{{mkr_sel_res}}_{DATE_STAMP}.rds',
    mkrs_f    = f'{zoom_dir}/{{zoom_name}}/pb_marker_genes_{FULL_TAG}_{{mkr_sel_res}}_{DATE_STAMP}.csv.gz',
    pb_hvgs_f = f'{zoom_dir}/{{zoom_name}}/pb_hvgs_{FULL_TAG}_{{mkr_sel_res}}_{DATE_STAMP}.csv.gz'
  params:
    ref_txome           = config['project']['ref_txome'],
    af_gtf_dt_f         = config['mapping']['af_gtf_dt_f'],
    batch_var           = BATCH_VAR,
    zoom_mkr_sel_res     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_sel_res'],
    zoom_mkr_min_cl_size = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cl_size'], 
    zoom_mkr_min_cells   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cells']
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_run_marker_genes', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_run_marker_genes', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_zoom/zoom_run_marker_genes_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/marker_genes.R'); calculate_marker_genes(
      integration_f = '{input.integration_f}', 
      h5ads_yaml_f  = '{input.h5ads_yaml_f}',
      pb_f          = '{output.pb_f}',
      mkrs_f        = '{output.mkrs_f}',
      pb_hvgs_f     = '{output.pb_hvgs_f}',
      gtf_dt_f      = '{params.af_gtf_dt_f}',
      sel_res       = '{params.zoom_mkr_sel_res}',
      min_cl_size   =  {params.zoom_mkr_min_cl_size},
      min_cells     =  {params.zoom_mkr_min_cells},
      zoom          = 'True', 
      batch_var     = '{params.batch_var}', 
      n_cores       =  {threads})"
    """


rule zoom_run_fgsea:
  input:
    mkrs_f        = f'{zoom_dir}/{{zoom_name}}/pb_marker_genes_{FULL_TAG}_{{mkr_sel_res}}_{DATE_STAMP}.csv.gz'
  output:
    fgsea_go_bp_f = f'{zoom_dir}/{{zoom_name}}/fgsea_{FULL_TAG}_{{mkr_sel_res}}_go_bp_{DATE_STAMP}.csv.gz', 
    fgsea_go_cc_f = f'{zoom_dir}/{{zoom_name}}/fgsea_{FULL_TAG}_{{mkr_sel_res}}_go_cc_{DATE_STAMP}.csv.gz',
    fgsea_go_mf_f = f'{zoom_dir}/{{zoom_name}}/fgsea_{FULL_TAG}_{{mkr_sel_res}}_go_mf_{DATE_STAMP}.csv.gz'
  params:
    ref_txome            = config['project']['ref_txome'],
    mkr_gsea_dir         = config['marker_genes']['mkr_gsea_dir'],
    zoom_mkr_min_cpm_go  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cpm_go'],
    zoom_mkr_max_zero_p  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_max_zero_p'],
    zoom_mkr_gsea_cut    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_gsea_cut'], 
    zoom_mkr_not_ok_re   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_not_ok_re'],
    zoom_mkr_gsea_var    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_gsea_var']
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_run_fgsea', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_run_fgsea', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_zoom/zoom_run_fgsea_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.benchmark.txt'
  conda: '../envs/rlibs.yaml'
  shell:"""
    Rscript -e "source('scripts/utils.R'); source('scripts/fgsea.R'); run_fgsea(
      mkrs_f        = '{input.mkrs_f}', 
      fgsea_go_bp_f = '{output.fgsea_go_bp_f}', 
      fgsea_go_cc_f = '{output.fgsea_go_cc_f}', 
      fgsea_go_mf_f = '{output.fgsea_go_mf_f}', 
      ref_txome     = '{params.ref_txome}', 
      gsea_dir      = '{params.mkr_gsea_dir}', 
      min_cpm_go    = {params.zoom_mkr_min_cpm_go}, 
      max_zero_p    = {params.zoom_mkr_max_zero_p},
      gsea_cut      = {params.zoom_mkr_gsea_cut},
      not_ok_re     = '{params.zoom_mkr_not_ok_re}',
      gsea_var      = '{params.zoom_mkr_gsea_var}',
      n_cores       =  {threads})"
    """

rule zoom_make_subsets:
  input:
    integration_f = f'{zoom_dir}/{{zoom_name}}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    smpl_stats_f  = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    h5ads_yaml_f  = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml'
  output:
    f"{zoom_dir}/{{zoom_name}}/{{obj_type}}_objects/{{prefix}}_cells_clean_{{zoom_name}}_{{batch}}_{FULL_TAG}_{DATE_STAMP}.{{ext}}"
  params:
    batch_var     = BATCH_VAR,
    zoom_lbls_f   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_f'],
    zoom_lbls_col = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'],
    zoom_lbls     = lambda wildcards: ','.join(ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels']),
    save_sce      = lambda wildcards: "TRUE" if wildcards.zoom_name in zooms_to_save_sce else "FALSE",
    save_adata    = lambda wildcards: "TRUE" if wildcards.zoom_name in zooms_to_save_anndata else "FALSE",
    sce_path      = lambda wildcards: ZOOM_OUT_MAP[wildcards.zoom_name][wildcards.batch].get("sce", ""),
    adata_path    = lambda wildcards: ZOOM_OUT_MAP[wildcards.zoom_name][wildcards.batch].get("adata", "")
  threads: 1
  resources:
    mem_mb  = lambda w, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'zoom_make_subsets', 'memory', attempt),
    runtime = lambda w, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'zoom_make_subsets', 'time', attempt)
  conda: '../envs/rlibs.yaml'
  shell:"""
    Rscript -e "source('scripts/zoom.R');
    make_subset_objects(
      sel_b         = '{wildcards.batch}',
      batch_var     = '{params.batch_var}',
      smpl_stats_f  = '{input.smpl_stats_f}',
      h5ads_yaml_f  = '{input.h5ads_yaml_f}',
      subset_f      = '{params.zoom_lbls_f}',
      subset_col    = '{params.zoom_lbls_col}',
      subset_str    = '{params.zoom_lbls}',
      integration_f = '{input.integration_f}',
      save_sce      = {params.save_sce},
      subset_sce_f  = '{params.sce_path}',
      save_adata    = {params.save_adata},
      subset_h5ad_f = '{params.adata_path}'
    )"
    """

# render_html_zoom
rule render_html_zoom:
  input:
    r_utils_f           = f'{code_dir}/utils.R',
    r_hvgs_f            = f'{code_dir}/hvgs.R', 
    r_int_f             = f'{code_dir}/integration.R',
    r_mkr_f             = f'{code_dir}/marker_genes.R',
    qc_f                = f'{qc_dir}/qc_all_samples_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    zoom_cell_hvgs_f    = f'{zoom_dir}/{{zoom_name}}/hvg_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    zoom_int_f          = f'{zoom_dir}/{{zoom_name}}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    zoom_pb_f           = f'{zoom_dir}/{{zoom_name}}/pb_{FULL_TAG}_{{mkr_sel_res}}_{DATE_STAMP}.rds',
    zoom_pb_hvgs_f      = f'{zoom_dir}/{{zoom_name}}/pb_hvgs_{FULL_TAG}_{{mkr_sel_res}}_{DATE_STAMP}.csv.gz',
    zoom_mkrs_f         = f'{zoom_dir}/{{zoom_name}}/pb_marker_genes_{FULL_TAG}_{{mkr_sel_res}}_{DATE_STAMP}.csv.gz', 
    zoom_empty_gs_f     = f'{zoom_dir}/{{zoom_name}}/edger_empty_genes_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    pb_empty_f          = f'{pb_dir}/pb_empties_{FULL_TAG}_{DATE_STAMP}.rds', 
    fgsea_files         =lambda wildcards: list(get_zoom_conditional_fgsea_files(
        config['project']['ref_txome'],
        zoom_dir,
        FULL_TAG,
        DATE_STAMP,
        ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_do_gsea']
    ).values())
  output:
    rmd_f       = f'{rmd_dir}/{SHORT_TAG}_zoom_{{zoom_name}}_{{mkr_sel_res}}.Rmd',
    html_f      = f'{docs_dir}/{SHORT_TAG}_zoom_{{zoom_name}}_{{mkr_sel_res}}.html'
  params:
    your_name             = config['project']['your_name'],
    affiliation           = config['project']['affiliation'],
    short_tag             = config['project']['short_tag'],
    date_stamp            = config['project']['date_stamp'],
    proj_dir              = config['project']['proj_dir'],
    ref_txome             = config['project']['ref_txome'],
    metadata_f            = config['project']['sample_metadata'], 
    zoom_dir              = zoom_dir,
    batch_var             = BATCH_VAR,
    meta_vars             = ','.join(config['project']['metadata_vars']), 
    fgsea_args            = lambda wildcards, input: ", ".join([
        f"fgsea_go_bp_f = '{input.get('fgsea_go_bp_f', '')}'",
        f"fgsea_go_cc_f = '{input.get('fgsea_go_cc_f', '')}'",
        f"fgsea_go_mf_f = '{input.get('fgsea_go_mf_f', '')}',"
    ]), 
    af_gtf_dt_f           = config['mapping']['af_gtf_dt_f'],
    zoom_int_res_ls       = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_res_ls'], 
    zoom_mkr_sel_res      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_sel_res'],
    zoom_mkr_min_cpm_mkr  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cpm_mkr'], 
    zoom_mkr_min_cells    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cells'],  
    zoom_mkr_not_ok_re    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_not_ok_re'], 
    zoom_mkr_do_gsea      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_do_gsea'], 
    zoom_mkr_gsea_var     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_gsea_var'],
    zoom_mkr_gsea_cut     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_gsea_cut'], 
    zoom_custom_mkr_names = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['custom_mkr_names'], 
    zoom_custom_mkr_paths = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['custom_mkr_paths']
  threads: 1
  retries: config['resources']['retries']
  conda: 
    '../envs/rlibs.yaml'
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'render_html_zoom', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'render_html_zoom', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_zoom/render_html_zoom_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.benchmark.txt'
  shell: """
    template_f=$(realpath resources/rmd_templates/zoom.Rmd.template)
    rule="zoom"

    Rscript --vanilla -e "source('scripts/render_htmls.R'); \
    render_html(
      rule_name         = '$rule', 
      proj_dir          = '{PROJ_DIR}', 
      temp_f            =  '$template_f', 
      rmd_f             = '{output.rmd_f}', 
      your_name         = '{params.your_name}', 
      affiliation       = '{params.affiliation}', 
      short_tag         = '{params.short_tag}', 
      date_stamp        = '{params.date_stamp}', 
      threads           =  {threads},
      zoom_dir          = '{params.zoom_dir}', 
      zoom_name         = '{wildcards.zoom_name}', 
      metadata_f        = '{params.metadata_f}', 
      meta_vars_ls      = '{params.meta_vars}', 
      gtf_dt_f          = '{params.af_gtf_dt_f}', 
      qc_f              = '{input.qc_f}', 
      cell_hvgs_f       = '{input.zoom_cell_hvgs_f}',
      int_f             = '{input.zoom_int_f}', 
      pb_f              = '{input.zoom_pb_f}', 
      mkrs_f            = '{input.zoom_mkrs_f}', 
      pb_hvgs_f         = '{input.zoom_pb_hvgs_f}',
      empty_gs_f        = '{input.zoom_empty_gs_f}', 
      pb_empty_f        = '{input.pb_empty_f}', 
      {params.fgsea_args}
      int_res_ls        = '{params.zoom_int_res_ls}', 
      custom_mkr_names  = '{params.zoom_custom_mkr_names}',
      custom_mkr_paths  = '{params.zoom_custom_mkr_paths}',
      mkr_sel_res       = '{params.zoom_mkr_sel_res}', 
      mkr_not_ok_re     = '{params.zoom_mkr_not_ok_re}', 
      mkr_min_cpm_mkr   =  {params.zoom_mkr_min_cpm_mkr}, 
      mkr_min_cells     =  {params.zoom_mkr_min_cells}, 
      mkr_gsea_var      = '{params.zoom_mkr_gsea_var}',
      mkr_gsea_cut      =  {params.zoom_mkr_gsea_cut}, 
      ref_txome         = '{params.ref_txome}',
      batch_var         = '{params.batch_var}',
      do_gsea           = '{params.zoom_mkr_do_gsea}'
    )"
    """

