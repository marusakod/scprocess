# Prepare selected source tables and sample-aware HVGs, then reuse join rules
# for matrix assembly, PCA, integration, pseudobulks, markers and reports.
import os
import pathlib
import sys

scprocess_dir = pathlib.Path(config.pop('scprocess_dir'))
sys.path.insert(0, str(scprocess_dir / 'scripts'))
from scprocess_utils import (get_join_zoom_spec, check_join_zoom_config,
  get_join_zoom_parameters, get_join_zoom_source_selection_parameters)

# Validate the requested zoom and inherit its parent join settings.
scdata_dir    = pathlib.Path(os.environ['SCPROCESS_DATA_DIR'])
ZOOM_NAME     = config.pop('join_zoom_name')
parent_config = config
ZOOM_SPEC     = get_join_zoom_spec(parent_config, ZOOM_NAME, scprocess_dir)
config        = check_join_zoom_config(parent_config, ZOOM_SPEC, scprocess_dir, scdata_dir)

# JOIN_ZOOM_PARAMS must exist before including join.smk: the shared rules use
# these paths/batches to distinguish the zoom from its parent analysis.
JOIN_ZOOM_PARAMS             = get_join_zoom_parameters(parent_config, config, ZOOM_SPEC)
ZOOM_SOURCE_SELECTION_PARAMS = get_join_zoom_source_selection_parameters(config, JOIN_ZOOM_PARAMS)
config['scprocess_dir']      = str(scprocess_dir)
include: 'join.smk'

# Shared join setup supplies directories, tags and resource settings above.
ZOOM_BATCHES            = list(JOIN_ZOOM_PARAMS['batches'])
JOIN_ZOOM_SCRIPT_F      = str(scprocess_dir / 'scripts/join_zoom.py')
HVGS_SCRIPT_F           = str(scprocess_dir / 'scripts/hvgs.py')
zoom_stats_pattern     = f"{join_int_dir}/tmp_raw_gene_count_statistics_{{batch}}_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
zoom_clip_pattern      = f"{join_int_dir}/tmp_clipped_gene_count_statistics_{{batch}}_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
zoom_trend_f           = f"{join_int_dir}/estimated_variances_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
zoom_batch_trend_pattern = f"{join_int_dir}/tmp_estimated_variances_{{batch}}_{JOIN_TAG}_{DATE_STAMP}.csv.gz"


rule join_zoom_prepare:
  input:
    integration_f      = JOIN_ZOOM_PARAMS['paths']['integration'],
    coldata_f          = JOIN_ZOOM_PARAMS['paths']['coldata'],
    metadata_f         = JOIN_ZOOM_PARAMS['paths']['metadata'],
    qc_f               = JOIN_ZOOM_PARAMS['paths']['qc'],
    labels_f           = JOIN_ZOOM_PARAMS['paths']['labels'],
    h5ads_yaml_f       = JOIN_ZOOM_PARAMS['paths']['manifest'],
    h5ad_fs            = [batch['path'] for batch in JOIN_ZOOM_PARAMS['batches'].values()],
    join_zoom_script_f = JOIN_ZOOM_SCRIPT_F
  output:
    cells_f        = JOIN_ZOOM_PARAMS['selected_f'],
    summary_f      = JOIN_ZOOM_PARAMS['summary_f'],
    batch_cells_fs = [batch['cells_f'] for batch in JOIN_ZOOM_PARAMS['batches'].values()],
    integration_fs = INTEGRATED_FS,
    sample_meta_fs = SAMPLE_META_FS,
    qc_fs          = QC_ALL_FS,
    h5ads_yaml_fs  = H5ADS_YAML_FS
  params:
    selection     = ZOOM_SOURCE_SELECTION_PARAMS['selection'],
    qc            = ZOOM_SOURCE_SELECTION_PARAMS['qc'],
    batches       = ZOOM_SOURCE_SELECTION_PARAMS['batches'],
    project_ids   = ZOOM_SOURCE_SELECTION_PARAMS['project_ids'],
    int_n_dims    = ZOOM_SOURCE_SELECTION_PARAMS['int_n_dims'],
    int_batch_var = ZOOM_SOURCE_SELECTION_PARAMS['int_batch_var']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'join_zoom_prepare', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'join_zoom_prepare', 'time', attempt)
  log:
    f"{logs_dir}/join_zoom_prepare_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_zoom_prepare_{DATE_STAMP}.txt"
  conda:
    '../envs/hvgs.yaml'

  # Carry large H5AD lookups through the job object rather than shell arguments.
  script: '../scripts/join_zoom.py'


rule join_zoom_calculate_raw_gene_count_statistics:
  input:
    cells_f            = lambda wildcards: JOIN_ZOOM_PARAMS['batches'][wildcards.batch]['cells_f'],
    h5ad_f             = lambda wildcards: JOIN_ZOOM_PARAMS['batches'][wildcards.batch]['path'],
    join_zoom_script_f = JOIN_ZOOM_SCRIPT_F
  output:
    statistics_f = zoom_stats_pattern
  params:
    chunk_size = JOIN_ZOOM_PARAMS['hvg']['hvg_chunk_size']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'join_zoom_calculate_raw_gene_count_statistics', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'join_zoom_calculate_raw_gene_count_statistics', 'time', attempt)
  log:
    f"{logs_dir}/join_zoom_calculate_raw_gene_count_statistics_{{batch}}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_zoom_calculate_raw_gene_count_statistics_{{batch}}_{DATE_STAMP}.txt"
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 {input.join_zoom_script_f:q} stats \
      --h5ad       {input.h5ad_f:q} \
      --cells      {input.cells_f:q} \
      --chunk-size {params.chunk_size} \
      --output     {output.statistics_f:q}
    """


rule join_zoom_fit_expected_gene_variances:
  input:
    statistics_fs      = expand(zoom_stats_pattern, batch = ZOOM_BATCHES),
    join_zoom_script_f = JOIN_ZOOM_SCRIPT_F,
    hvgs_script_f      = HVGS_SCRIPT_F
  output:
    trend_f        = zoom_trend_f,
    batch_trend_fs = expand(zoom_batch_trend_pattern, batch = ZOOM_BATCHES)
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'join_zoom_fit_expected_gene_variances', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'join_zoom_fit_expected_gene_variances', 'time', attempt)
  log:
    f"{logs_dir}/join_zoom_fit_expected_gene_variances_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_zoom_fit_expected_gene_variances_{DATE_STAMP}.txt"
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 {input.join_zoom_script_f:q} trend \
      --inputs        {input.statistics_fs:q} \
      --output        {output.trend_f:q} \
      --batch-outputs {output.batch_trend_fs:q}
    """


rule join_zoom_calculate_clipped_gene_count_statistics:
  input:
    cells_f            = lambda wildcards: JOIN_ZOOM_PARAMS['batches'][wildcards.batch]['cells_f'],
    h5ad_f             = lambda wildcards: JOIN_ZOOM_PARAMS['batches'][wildcards.batch]['path'],
    trend_f            = zoom_batch_trend_pattern,
    join_zoom_script_f = JOIN_ZOOM_SCRIPT_F
  output:
    statistics_f = zoom_clip_pattern
  params:
    chunk_size = JOIN_ZOOM_PARAMS['hvg']['hvg_chunk_size']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'join_zoom_calculate_clipped_gene_count_statistics', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'join_zoom_calculate_clipped_gene_count_statistics', 'time', attempt)
  log:
    f"{logs_dir}/join_zoom_calculate_clipped_gene_count_statistics_{{batch}}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_zoom_calculate_clipped_gene_count_statistics_{{batch}}_{DATE_STAMP}.txt"
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 {input.join_zoom_script_f:q} clip \
      --h5ad       {input.h5ad_f:q} \
      --cells      {input.cells_f:q} \
      --chunk-size {params.chunk_size} \
      --trend      {input.trend_f:q} \
      --output     {output.statistics_f:q}
    """


rule join_zoom_rank_and_select_highly_variable_genes:
  input:
    statistics_fs      = expand(zoom_clip_pattern, batch = ZOOM_BATCHES),
    trend_f            = zoom_trend_f,
    join_zoom_script_f = JOIN_ZOOM_SCRIPT_F,
    hvgs_script_f      = HVGS_SCRIPT_F,
    exclude_fs        = (
      [JOIN_ZOOM_PARAMS['hvg']['hvg_exclude_from_file']]
      if JOIN_ZOOM_PARAMS['hvg'].get('hvg_exclude_from_file') else [])
  output:
    hvgs_f = HVG_FS[0]
  params:
    n_hvgs       = N_HVGS,
    exclude_flag = '--exclude' if JOIN_ZOOM_PARAMS['hvg'].get('hvg_exclude_from_file') else ''
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'join_zoom_rank_and_select_highly_variable_genes', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'join_zoom_rank_and_select_highly_variable_genes', 'time', attempt)
  log:
    f"{logs_dir}/join_zoom_rank_and_select_highly_variable_genes_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_zoom_rank_and_select_highly_variable_genes_{DATE_STAMP}.txt"
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 {input.join_zoom_script_f:q} rank \
      --inputs {input.statistics_fs:q} \
      --trend  {input.trend_f:q} \
      --n-hvgs {params.n_hvgs} \
      --output {output.hvgs_f:q} \
      {params.exclude_flag} {input.exclude_fs:q}
    """
