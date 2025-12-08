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

# check config
config          = check_config(config, proj_schema_f, scdata_dir, scprocess_dir)

# get lists of parameters
RUN_PARAMS, RUN_VAR     = get_run_parameters(config, scdata_dir)
RUNS                    = list(RUN_PARAMS.keys())
BATCH_PARAMS, BATCH_VAR = get_batch_parameters(config, RUNS, scdata_dir)
BATCHES                 = list(BATCH_PARAMS.keys())
RUNS_TO_BATCHES, BATCHES_TO_RUNS = get_batches_to_runs(config, RUNS, BATCHES)

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
rmd_dir       = f"{PROJ_DIR}/analysis"
docs_dir      = f"{PROJ_DIR}/public"
zoom_dir      = f"{PROJ_DIR}/output/{SHORT_TAG}_zoom"

# define zoom marker outputs
zoom_mkr_report_outs = expand(
  [
  '%s/{zoom_name}/pb_%s_{mkr_sel_res}_%s.rds' % (zoom_dir, FULL_TAG, DATE_STAMP), 
  '%s/{zoom_name}/pb_marker_genes_%s_{mkr_sel_res}_%s.txt.gz' % (zoom_dir, FULL_TAG, DATE_STAMP), 
  '%s/{zoom_name}/pb_hvgs_%s_{mkr_sel_res}_%s.txt.gz' % (zoom_dir, FULL_TAG, DATE_STAMP), 
  '%s/%s_zoom_{zoom_name}_{mkr_sel_res}.Rmd' % (rmd_dir, SHORT_TAG), 
  '%s/%s_zoom_{zoom_name}_{mkr_sel_res}.html' % (docs_dir, SHORT_TAG)
  ] + (
    [
    '%s/{zoom_name}/fgsea_%s_{mkr_sel_res}_go_bp_%s.txt.gz' % (zoom_dir, FULL_TAG, DATE_STAMP), 
    '%s/{zoom_name}/fgsea_%s_{mkr_sel_res}_go_cc_%s.txt.gz' % (zoom_dir, FULL_TAG, DATE_STAMP), 
    '%s/{zoom_name}/fgsea_%s_{mkr_sel_res}_go_mf_%s.txt.gz' % (zoom_dir, FULL_TAG, DATE_STAMP), 
    '%s/{zoom_name}/fgsea_%s_{mkr_sel_res}_paths_%s.txt.gz' % (zoom_dir, FULL_TAG, DATE_STAMP),  
    '%s/{zoom_name}/fgsea_%s_{mkr_sel_res}_hlmk_%s.txt.gz' % (zoom_dir, FULL_TAG, DATE_STAMP)
    ]
    if (config['project']['species'] in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020']) & config['marker_genes']['mkr_do_gsea']
    else []
   ),
  zip,
  zoom_name = ZOOMS,
  mkr_sel_res = [ZOOM_PARAMS[zoom_name]['marker_genes']["mkr_sel_res"] for zoom_name in ZOOMS]
)
zooms_to_save = [ zoom_name for zoom_name in ZOOMS if ZOOM_PARAMS[zoom_name]['zoom']['save_subset_sces'] ]

zoom_sce_outs = (
  expand(
    '%s/{zoom_name}/sce_objects/sce_cells_clean_{zoom_name}_{batch}_%s_%s.rds' % \
    (zoom_dir, FULL_TAG, DATE_STAMP),
    zoom_name = zooms_to_save,
    sample    = BATCHES
  ) if len(zooms_to_save) > 0 else []
)


rule zoom:
  input:
    # zoom sample qc
    expand('%s/{zoom_name}/zoom_sample_statistics_%s_%s.csv' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS),
    # zoom pseudobulks and empties
    expand('%s/{zoom_name}/pb_{zoom_name}_%s_%s.rds' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS),
    expand('%s/{zoom_name}/edger_empty_genes_%s_%s.csv.gz' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS), 
    # zoom hvgs
    expand('%s/{zoom_name}/hvg_paths_%s_%s.csv' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS),  
    expand('%s/{zoom_name}/standardized_variance_stats_%s_%s.txt.gz' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS), 
    expand('%s/{zoom_name}/hvg_dt_%s_%s.txt.gz' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS), 
    expand('%s/{zoom_name}/top_hvgs_counts_%s_%s.h5' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS), 
    # zoom integration
    expand('%s/{zoom_name}/integrated_dt_%s_%s.csv.gz' % \
      (zoom_dir, FULL_TAG, DATE_STAMP), zoom_name = ZOOMS),
    # zoom marker genes and html report
    zoom_mkr_report_outs, 
    # zoom sce subsets (optional)
    zoom_sce_outs


rule get_zoom_sample_statistics:
  input:
    qc_stats_f      = qc_dir + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    zoom_stats_f    = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  params:
    zoom_lbls_f     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_f'],
    zoom_lbls_col   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'], 
    zoom_lbls       = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels'],
    zoom_min_n_smpl = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc']['qc_min_cells'],
    ambient_method  = config['ambient']['ambient_method']
  run:
    zoom_stats_df   = extract_zoom_sample_statistics(input.qc_stats_f, BATCHES, 
      params.zoom_lbls_f, params.zoom_lbls_col, params.zoom_lbls, params.zoom_min_n_smpl, 
      params.ambient_method)
    zoom_stats_df.to_csv(output.zoom_stats_f, index = False)


# pseudobulks and empties
rule zoom_make_pb_subset:
  input:
    sces_yaml_f   = int_dir  + '/sce_clean_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml', 
    zoom_stats_f  = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    zoom_pb_f     = zoom_dir + '/{zoom_name}/pb_{zoom_name}_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  params:
    zoom_lbls_f   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_f'],
    zoom_lbls_col = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'],
    zoom_lbls     = lambda wildcards: ','.join(ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels'])
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('zoom_make_pb_subset', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('zoom_make_pb_subset', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_zoom/zoom_make_pb_subset_{zoom_name}_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    make_pb_cells( \
      sce_fs_yaml = '{input.sces_yaml_f}',
      qc_stats_f  = '{input.zoom_stats_f}',
      subset_f    = '{params.zoom_lbls_f}',
      subset_col  = '{params.zoom_lbls_col}', 
      subset_str  = '{params.zoom_lbls}', 
      pb_f        = '{output.zoom_pb_f}',
      n_cores     =  {threads}
    )"
    """


rule zoom_calculate_ambient_genes:
  input:
    pb_empty_f      = pb_dir + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
    zoom_pb_f       = zoom_dir + '/{zoom_name}/pb_{zoom_name}_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  output:
    zoom_empty_gs_f = zoom_dir + '/{zoom_name}/edger_empty_genes_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  params:
    zoom_fdr_thr    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['pb_empties']['ambient_genes_fdr_thr'],
    zoom_logfc_thr  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['pb_empties']['ambient_genes_logfc_thr']
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('zoom_calculate_ambient_genes', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('zoom_calculate_ambient_genes', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_zoom/zoom_calculate_ambient_genes_{zoom_name}_' + DATE_STAMP + '.benchmark.txt'
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
    amb_yaml_fs = expand([amb_dir + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml'], run=RUNS)
  output:
    hvg_paths_f = zoom_dir + '/{zoom_name}/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  params:
    demux_type  = config['multiplexing']['demux_type'],
    run_var     = RUN_VAR,
    runs        = RUNS,
    mapping     = RUNS_TO_BATCHS
  run:
    hvg_df = make_hvgs_input_df(
      params.demux_type, params.run_var, params.runs, input.amb_yaml_fs,
      params.mapping, FULL_TAG, DATE_STAMP, f"{zoom_dir}/{wildcards.zoom_name}"
      )
    hvg_df.to_csv(output.hvg_paths_f, index=False)


rule zoom_make_tmp_csr_matrix:
  input:
    hvg_paths_f     = zoom_dir + '/{zoom_name}/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    smpl_stats_f    = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    rowdata_f       = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  output:
    clean_h5_f      = temp(expand([
      zoom_dir + '/{zoom_name}' + '/chunked_counts_{batch}_' + FULL_TAG + '_' + DATE_STAMP + '.h5'
      ], zoom_name = '{zoom_name}', batch = BATCHES))
  params: 
    zoom_lbls_f     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_f'],
    zoom_lbls_col   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'],
    zoom_lbls       = lambda wildcards: ','.join(ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels']), 
    run_var         = RUN_VAR,
    demux_type      = config['multiplexing']['demux_type'],
    zoom_chunk_size = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_chunk_size'],
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('zoom_make_tmp_csr_matrix', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('zoom_make_tmp_csr_matrix', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_zoom/zoom_make_tmp_csr_matrix_{zoom_name}_' + DATE_STAMP + '.benchmark.txt'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    python3 scripts/hvgs.py get_csr_counts \
      "{input.hvg_paths_f}" \
      "{params.zoom_lbls_f}" \
      "{params.zoom_lbls_col}" \
      "{params.zoom_lbls}" \
      "{input.smpl_stats_f}" \
      "{input.rowdata_f}" \
      "{params.run_var}" \
      "{params.demux_type}" \
      --chunksize {params.zoom_chunk_size} \
      --ncores {threads}
    """


rule zoom_get_stats_for_std_variance_for_sample:
  input: 
    clean_h5_f    = zoom_dir + '/{zoom_name}' + '/chunked_counts_{batch}_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
    smpl_stats_f  = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    rowdata_f     = qc_dir   + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  output:
    std_var_stats_f = temp(zoom_dir + '/{zoom_name}' + '/tmp_std_var_stats_{batch}_sample_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('zoom_get_stats_for_std_variance_for_sample', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('zoom_get_stats_for_std_variance_for_sample', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_zoom/zoom_get_stats_for_std_variance_for_sample_{zoom_name}_{batch}_' + DATE_STAMP + '.benchmark.txt'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    python3 scripts/hvgs.py calculate_std_var_stats_for_sample \
      {wildcards.sample} \
      {input.smpl_stats_f} \
      {input.clean_h5_f} \
      {input.rowdata_f} \
      {output.std_var_stats_f}
    """


rule zoom_get_mean_var_for_group:
  input:
    clean_h5_f    = expand(
      zoom_dir + '/{zoom_name}' + '/chunked_counts_{batch}_' + FULL_TAG + '_' + DATE_STAMP + '.h5',
      zoom_name = ZOOMS, batch = BATCHES
    ),
    hvg_paths_f   = zoom_dir + '/{zoom_name}' + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    rowdata_f     = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', 
    smpl_stats_f  = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output: 
    mean_var_f    = temp(zoom_dir + '/{zoom_name}' + '/tmp_mean_var_{group}_group_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
  params:
    zoom_hvg_method     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_method'], 
    zoom_hvg_chunk_size = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_chunk_size'], 
    zoom_group_var      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_split_var']
  threads: 8
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('zoom_get_mean_var_for_group', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('zoom_get_mean_var_for_group', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_zoom/zoom_get_mean_var_for_group_{zoom_name}_{group}_chunk_{chunk}' + DATE_STAMP + '.benchmark.txt'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    python3 scripts/hvgs.py calculate_mean_var_for_chunk \
      {input.hvg_paths_f} \
      {input.rowdata_f} \
      {METADATA_F} \
      {input.smpl_stats_f} \
      {output.mean_var_f} \
      {wildcards.chunk} \
      {params.zoom_hvg_method} \
      {params.zoom_hvg_chunk_size} \
      --group {wildcards.group} \
      --groupvar {params.zoom_group_var} \
      --ncores {threads} 
    """


rule zoom_merge_group_mean_var:
  input:         
    mean_var_f    = lambda wildcards: get_zoom_raw_mean_var_files(wildcards.zoom_name, ZOOM_PARAMS, FULL_TAG, DATE_STAMP)
  output:
    mean_var_merged_f = temp(zoom_dir + '/{zoom_name}'  + '/means_variances_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('zoom_merge_group_mean_var', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('zoom_merge_group_mean_var', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  run:
    merge_tmp_files(input.mean_var_f, output.mean_var_merged_f)


rule zoom_get_estimated_variances:
  input:
    mean_var_merged_f = zoom_dir + '/{zoom_name}' + '/means_variances_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    estim_vars_f      = temp(zoom_dir + '/{zoom_name}' + '/estimated_variances_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
  params: 
    zoom_hvg_method   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_method']
  threads: 1
  retries: config['resources']['retries']
  conda:
    '../envs/hvgs.yaml'
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('zoom_get_estimated_variances', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('zoom_get_estimated_variances', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_zoom/zoom_get_estimated_variances_{zoom_name}_' + DATE_STAMP + '.benchmark.txt'
  shell: """
    python3 scripts/hvgs.py calculate_estimated_vars \
      {output.estim_vars_f} \
      {params.zoom_hvg_method} \
      {input.mean_var_merged_f}
    """


rule zoom_get_stats_for_std_variance_for_group:
  input: 
    clean_h5_fs   = expand(
      zoom_dir + '/{zoom_name}' + '/chunked_counts_{batch}_' + FULL_TAG + '_' + DATE_STAMP + '.h5',
      zoom_name = ZOOMS, batch = BATCHES
    ),
    estim_vars_f  = zoom_dir + '/{zoom_name}'  + '/estimated_variances_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    hvg_paths_f   = zoom_dir + '/{zoom_name}'  + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    rowdata_f     = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', 
    smpl_stats_f  = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    std_var_stats_f = temp(zoom_dir + '/{zoom_name}' + '/tmp_std_var_stats_{group}_group_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
  params:
    metadata_f          = config['project']['sample_metadata'],
    zoom_hvg_method     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_method'], 
    zoom_hvg_chunk_size = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_chunk_size'], 
    zoom_hvg_group_var  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_split_var']
  threads: 8
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('zoom_get_stats_for_std_variance_for_group', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('zoom_get_stats_for_std_variance_for_group', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_zoom/zoom_get_stats_for_std_variance_for_group_{zoom_name}_{group}_chunk_{chunk}' + DATE_STAMP + '.benchmark.txt'
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
    std_var_stats_merged_f= zoom_dir + '/{zoom_name}/standardized_variance_stats_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  run:
    merge_tmp_files(input.tmp_std_var_stats_fs, output.std_var_stats_merged_f)

        
rule zoom_get_highly_variable_genes:
  input:
    std_var_stats_f = zoom_dir + '/{zoom_name}/standardized_variance_stats_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    empty_gs_fs     = zoom_dir + '/{zoom_name}/edger_empty_genes_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz' 
  output:
    hvg_f =  zoom_dir + '/{zoom_name}/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  threads: 1
  retries: config['resources']['retries']
  params:
    zoom_hvg_method = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_method'],
    zoom_n_hvgs     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_n_hvgs'],
    zoom_exclude_ambient_genes = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_exclude_ambient_genes']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('zoom_get_highly_variable_genes', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('zoom_get_highly_variable_genes', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_zoom/zoom_get_highly_variable_genes_{zoom_name}_' + DATE_STAMP + '.benchmark.txt'
  conda:
    '../envs/hvgs.yaml'
  shell: """
     NOAMBIENT_FLAG=""
     if [ "{params.zoom_exclude_ambient_genes}" = "True" ]; then
       NOAMBIENT_FLAG="--noambient"
     fi

     python3 scripts/hvgs.py calculate_hvgs \
      {input.std_var_stats_f} \
      {output.hvg_f} \
      {input.empty_gs_fs} \
      {params.zoom_hvg_method} \
      {params.zoom_n_hvgs} \
      $NOAMBIENT_FLAG
     """


rule zoom_create_hvg_matrix:
  input: 
    clean_h5_f    = expand(
      zoom_dir + '/{zoom_name}' + '/chunked_counts_{batch}_' + FULL_TAG + '_' + DATE_STAMP + '.h5',
      zoom_name = ZOOMS, batch = BATCHES
    ),
    smpl_stats_f  = zoom_dir  + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    hvg_paths_f   = zoom_dir  + '/{zoom_name}/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    hvg_f         = zoom_dir  + '/{zoom_name}/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    hvg_mat_f     = zoom_dir  + '/{zoom_name}/top_hvgs_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5'
  params:
    demux_type    = config['multiplexing']['demux_type']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('zoom_create_hvg_matrix', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('zoom_create_hvg_matrix', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_zoom/zoom_create_hvg_matrix_{zoom_name}_' + DATE_STAMP + '.benchmark.txt'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    python3 scripts/hvgs.py read_top_genes \
      {input.smpl_stats_f} \
      {input.hvg_paths_f} \
      {input.hvg_f} \
      {output.hvg_mat_f} \
      {params.demux_type}
    """


rule zoom_run_integration:
  input:
    hvg_mat_f     = zoom_dir + '/{zoom_name}/top_hvgs_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
    smpl_stats_f  = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    coldata_f     = qc_dir   + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  output:
    integration_f = zoom_dir + '/{zoom_name}/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  params:
    demux_type          = config['multiplexing']['demux_type'],
    exclude_mito        = config['qc']['exclude_mito'],
    zoom_int_embedding  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_embedding'],
    zoom_int_n_dims     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_n_dims'],
    zoom_int_cl_method  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_cl_method'], 
    zoom_int_theta      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_theta'],
    zoom_int_res_ls     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_res_ls'],
    zoom_int_batch_var  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_batch_var']
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('zoom_run_integration', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('zoom_run_integration', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_zoom/zoom_run_integration_{zoom_name}_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    # run harmony
    Rscript -e "source('scripts/integration.R'); source('scripts/utils.R'); source('scripts/zoom.R'); 
      run_zoom_integration( 
        hvg_mat_f        = '{input.hvg_mat_f}', 
        smpl_stats_f     = '{input.smpl_stats_f}', 
        coldata_f        = '{input.coldata_f}', 
        demux_type       = '{params.demux_type}', 
        exclude_mito     = '{params.exclude_mito}', 
        reduction        = '{params.zoom_int_embedding}',
        n_dims           =  {params.zoom_int_n_dims}, 
        cl_method        = '{params.zoom_int_cl_method}', 
        theta            =  {params.zoom_int_theta}, 
        res_ls_concat    = '{params.zoom_int_res_ls}', 
        integration_f    = '{output.integration_f}', 
        batch_var        = '{params.zoom_int_batch_var}', 
        n_cores          =  {threads})"
    """


rule zoom_run_marker_genes:
  input:
    sces_yaml_f    = int_dir  + '/sce_clean_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml',
    integration_f  = zoom_dir + '/{zoom_name}/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  output:
    pb_f      = zoom_dir + '/{zoom_name}/pb_' + FULL_TAG + '_' + '{mkr_sel_res}_' + DATE_STAMP + '.rds',
    mkrs_f    = zoom_dir + '/{zoom_name}/pb_marker_genes_' + FULL_TAG + '_' + '{mkr_sel_res}_' + DATE_STAMP + '.txt.gz',
    pb_hvgs_f = zoom_dir + '/{zoom_name}/pb_hvgs_' + FULL_TAG + '_' + '{mkr_sel_res}_' + DATE_STAMP + '.txt.gz', 
    **get_zoom_conditional_outputs(config['project']['species'], zoom_dir, FULL_TAG, DATE_STAMP)
  params:
    species             = config['project']['species'],
    af_gtf_dt_f         = config['mapping']['af_gtf_dt_f'],
    mkr_gsea_dir        = config['marker_genes']['mkr_gsea_dir'],
    fgsea_args = lambda wildcards, output: ", ".join([
        f"fgsea_go_bp_f = '{output.get('fgsea_go_bp_f', '')}'",
        f"fgsea_go_cc_f = '{output.get('fgsea_go_cc_f', '')}'",
        f"fgsea_go_mf_f = '{output.get('fgsea_go_mf_f', '')}'",
        f"fgsea_paths_f = '{output.get('fgsea_paths_f', '')}'",
        f"fgsea_hlmk_f = '{output.get('fgsea_hlmk_f', '')}',"
    ]), 
    zoom_mkr_sel_res     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_sel_res'],
    zoom_mkr_min_cl_size = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cl_size'], 
    zoom_mkr_min_cells   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cells'], 
    zoom_mkr_not_ok_re   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_not_ok_re'],
    zoom_mkr_min_cpm_go  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cpm_go'],
    zoom_mkr_max_zero_p  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_max_zero_p'], 
    zoom_mkr_do_gsea     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_do_gsea'], 
    zoom_mkr_gsea_cut    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_gsea_cut']
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('zoom_run_marker_genes', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('zoom_run_marker_genes', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_zoom/zoom_run_marker_genes_{zoom_name}_{mkr_sel_res}_' + DATE_STAMP + '.benchmark.txt'
  conda:
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/marker_genes.R'); calculate_marker_genes(
      integration_f = '{input.integration_f}', 
      sces_yaml_f   = '{input.sces_yaml_f}',
      pb_f          = '{output.pb_f}',
      mkrs_f        = '{output.mkrs_f}',
      pb_hvgs_f     = '{output.pb_hvgs_f}',
      {params.fgsea_args}
      species       = '{params.species}',
      do_gsea       = '{params.zoom_mkr_do_gsea}', 
      gtf_dt_f      = '{params.af_gtf_dt_f}',
      gsea_dir      = '{params.mkr_gsea_dir}',
      sel_res       = '{params.zoom_mkr_sel_res}',
      min_cl_size   =  {params.zoom_mkr_min_cl_size},
      min_cells     =  {params.zoom_mkr_min_cells},
      not_ok_re     = '{params.zoom_mkr_not_ok_re}',
      min_cpm_go    =  {params.zoom_mkr_min_cpm_go},
      max_zero_p    =  {params.zoom_mkr_max_zero_p},
      gsea_cut      =  {params.zoom_mkr_gsea_cut},
      zoom          = 'True', 
      n_cores       =  {threads})"
    """


rule zoom_make_subset_sces:
  input:
    integration_f = zoom_dir + '/{zoom_name}/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', 
    smpl_stats_f  = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    sces_yaml_f   = int_dir + '/sce_clean_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml'
  output:
    clean_sce_f = zoom_dir + '/{zoom_name}/sce_objects/sce_cells_clean_{batch}_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  params:
    zoom_lbls_f     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_f'],
    zoom_lbls_col   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'],
    zoom_lbls       = lambda wildcards: ','.join(ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels'])
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('zoom_make_subset_sces', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('zoom_make_subset_sces', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_zoom/zoom_make_subset_sces_{zoom_name}_{batch}_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/zoom.R');
     make_subset_sces(
     sel_s          = '{wildcards.sample}',
     clean_sce_f    = '{output.clean_sce_f}',
     integration_f  = '{input.integration_f}',
     smpl_stats_f   = '{input.smpl_stats_f}',
     sces_yaml_f    = '{input.sces_yaml_f}',
     subset_f       = '{params.zoom_lbls_f}',
     subset_col     = '{params.zoom_lbls_col}',
     subset_str     = '{params.zoom_lbls}')"
    """


# render_html_zoom
rule render_html_zoom:
  input:
    r_utils_f           = code_dir + '/utils.R',
    r_hvgs_f            = code_dir + '/hvgs.R', 
    r_int_f             = code_dir + '/integration.R',
    r_mkr_f             = code_dir + '/marker_genes.R',
    qc_f                = qc_dir  + '/qc_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz',
    zoom_cell_hvgs_f    = zoom_dir + '/{zoom_name}/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    zoom_int_f          = zoom_dir + '/{zoom_name}/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz',
    zoom_pb_f           = zoom_dir + '/{zoom_name}/pb_' + FULL_TAG + '_' + '{mkr_sel_res}_' + DATE_STAMP + '.rds',
    zoom_pb_hvgs_f      = zoom_dir + '/{zoom_name}/pb_hvgs_' + FULL_TAG + '_' + '{mkr_sel_res}_' + DATE_STAMP + '.txt.gz',
    zoom_mkrs_f         = zoom_dir + '/{zoom_name}/pb_marker_genes_' + FULL_TAG + '_' + '{mkr_sel_res}_' + DATE_STAMP + '.txt.gz',
    zoom_empty_gs_f     = zoom_dir + '/{zoom_name}/edger_empty_genes_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    pb_empty_f          = pb_dir + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
    **get_zoom_conditional_outputs(config['project']['species'], zoom_dir, FULL_TAG, DATE_STAMP)
  output:
    rmd_f       = rmd_dir + '/' + SHORT_TAG + '_zoom' + '_{zoom_name}_{mkr_sel_res}.Rmd',
    html_f      = docs_dir + '/' + SHORT_TAG + '_zoom' + '_{zoom_name}_{mkr_sel_res}.html'
  params:
    your_name             = config['project']['your_name'],
    affiliation           = config['project']['affiliation'],
    short_tag             = config['project']['short_tag'],
    date_stamp            = config['project']['date_stamp'],
    proj_dir              = config['project']['proj_dir'],
    species               = config['project']['species'],
    metadata_f            = config['project']['sample_metadata'], 
    zoom_dir              = zoom_dir,
    meta_vars             = ','.join(config['project']['metadata_vars']), 
    fgsea_args            = lambda wildcards, input: ", ".join([
        f"fgsea_go_bp_f = '{input.get('fgsea_go_bp_f', '')}'",
        f"fgsea_go_cc_f = '{input.get('fgsea_go_cc_f', '')}'",
        f"fgsea_go_mf_f = '{input.get('fgsea_go_mf_f', '')}'",
        f"fgsea_paths_f = '{input.get('fgsea_paths_f', '')}'",
        f"fgsea_hlmk_f = '{input.get('fgsea_hlmk_f', '')}',"
    ]), 
    af_gtf_dt_f           = config['mapping']['af_gtf_dt_f'],
    zoom_int_res_ls       = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_res_ls'], 
    zoom_mkr_sel_res      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_sel_res'],
    zoom_mkr_min_cpm_mkr  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cpm_mkr'], 
    zoom_mkr_min_cells    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cells'],  
    zoom_mkr_not_ok_re    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_not_ok_re'], 
    zoom_mkr_do_gsea      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_do_gsea'], 
    zoom_mkr_gsea_cut     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_gsea_cut'], 
    zoom_custom_mkr_names = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['custom_mkr_names'], 
    zoom_custom_mkr_paths = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['custom_mkr_paths']
  threads: 1
  retries: config['resources']['retries']
  conda: 
    '../envs/rlibs.yaml'
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('render_html_zoom', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('render_html_zoom', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_zoom/render_html_zoom_{zoom_name}_{mkr_sel_res}_' + DATE_STAMP + '.benchmark.txt'
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
      mkr_gsea_cut      =  {params.zoom_mkr_gsea_cut}, 
      species           = '{params.species}',
      do_gsea           = '{params.zoom_mkr_do_gsea}'
    )"
    """

