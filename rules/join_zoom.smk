# Joined zooms prepare source tables, then reuse the existing join rules unchanged
# for matrix assembly, native PCA, integration, pseudobulks, markers and reports.
# load modules
import os
import pathlib
import sys

scprocess_dir = pathlib.Path(config.pop('scprocess_dir'))
sys.path.insert(0, str(scprocess_dir / 'scripts'))
from scprocess_utils import (get_join_zoom_specs, check_join_zoom_config,
  get_join_zoom_parameters, get_join_zoom_source_selection_parameters)

# check configuration
scdata_dir    = pathlib.Path(os.environ['SCPROCESS_DATA_DIR'])
ZOOM_NAME     = config.pop('join_zoom_name')
parent_config = config
ZOOM_SPECS    = get_join_zoom_specs(parent_config, scprocess_dir)
if ZOOM_NAME not in ZOOM_SPECS:
  raise ValueError(f'Unknown joined zoom {ZOOM_NAME!r}; available: {list(ZOOM_SPECS)}')
ZOOM_SPEC = ZOOM_SPECS[ZOOM_NAME]
config    = check_join_zoom_config(parent_config, ZOOM_SPEC, scprocess_dir, scdata_dir)

# get parameter sets
JOIN_ZOOM_PARAMS             = get_join_zoom_parameters(parent_config, config, ZOOM_SPEC)
ZOOM_SOURCE_SELECTION_PARAMS = get_join_zoom_source_selection_parameters(config, JOIN_ZOOM_PARAMS)
config['scprocess_dir'] = str(scprocess_dir)

# shared join rules supply matrix assembly, PCA, integration, markers and reports
include: 'join.smk'

# unpack zoom parameters
ZOOM_BATCHES = list(JOIN_ZOOM_PARAMS['batches'])

# specify locations
JOIN_ZOOM_SCRIPT_F       = str(scprocess_dir / 'scripts/join_zoom.py')
HVGS_SCRIPT_F            = str(scprocess_dir / 'scripts/hvgs.py')
zoom_stats_pattern       = f"{join_int_dir}/tmp_raw_gene_count_statistics_{{batch}}_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
zoom_clip_pattern        = f"{join_int_dir}/tmp_clipped_gene_count_statistics_{{batch}}_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
zoom_trend_f             = f"{join_int_dir}/estimated_variances_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
zoom_batch_trend_pattern = f"{join_int_dir}/tmp_estimated_variances_{{batch}}_{JOIN_TAG}_{DATE_STAMP}.csv.gz"

rule join_zoom_prepare:
  input:
    integration_f      = JOIN_ZOOM_PARAMS['paths']['integration'],
    coldata_f          = JOIN_ZOOM_PARAMS['paths']['coldata'],
    metadata_f         = JOIN_ZOOM_PARAMS['paths']['metadata'],
    qc_f               = JOIN_ZOOM_PARAMS['paths']['qc'],
    labels_f           = JOIN_ZOOM_PARAMS['paths']['labels'],
    manifest_f         = JOIN_ZOOM_PARAMS['paths']['manifest'],
    h5ads              = [b['path'] for b in JOIN_ZOOM_PARAMS['batches'].values()],
    JOIN_ZOOM_SCRIPT_F = JOIN_ZOOM_SCRIPT_F
  output:
    cells        = JOIN_ZOOM_PARAMS['selected_f'],
    summary      = JOIN_ZOOM_PARAMS['summary_f'],
    batch_cells  = [b['cells_f'] for b in JOIN_ZOOM_PARAMS['batches'].values()],
    integrations = INTEGRATED_FS,
    metadata     = SAMPLE_META_FS,
    qc           = QC_ALL_FS,
    manifests    = H5ADS_YAML_FS
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
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_prepare', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_prepare', 'time', attempt)
  log:
    f"{logs_dir}/join_zoom_prepare_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_zoom_prepare_{DATE_STAMP}.txt"
  conda:
    '../envs/hvgs.yaml'
  # A script transports the manifest through Snakemake's job object, avoiding
  # shell argument limits for joins with thousands of H5AD files.
  script: '../scripts/join_zoom.py'

rule join_zoom_calculate_raw_gene_count_statistics:
  input:
    cells              = lambda wildcards: JOIN_ZOOM_PARAMS['batches'][wildcards.batch]['cells_f'],
    h5ad               = lambda wildcards: JOIN_ZOOM_PARAMS['batches'][wildcards.batch]['path'],
    JOIN_ZOOM_SCRIPT_F = JOIN_ZOOM_SCRIPT_F
  output: zoom_stats_pattern
  params: chunk_size = JOIN_ZOOM_PARAMS['hvg']['hvg_chunk_size']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_calculate_raw_gene_count_statistics', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_calculate_raw_gene_count_statistics', 'time', attempt)
  log:
    f"{logs_dir}/join_zoom_calculate_raw_gene_count_statistics_{{batch}}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_zoom_calculate_raw_gene_count_statistics_{{batch}}_{DATE_STAMP}.txt"
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 {input.JOIN_ZOOM_SCRIPT_F:q} stats --h5ad {input.h5ad:q} --cells {input.cells:q} \
      --chunk-size {params.chunk_size} --output {output:q}
    """

rule join_zoom_fit_expected_gene_variances:
  input:
    stats              = expand(zoom_stats_pattern, batch=ZOOM_BATCHES),
    JOIN_ZOOM_SCRIPT_F = JOIN_ZOOM_SCRIPT_F,
    HVGS_SCRIPT_F      = HVGS_SCRIPT_F
  output:
    trend   = zoom_trend_f,
    batches = expand(zoom_batch_trend_pattern, batch=ZOOM_BATCHES)
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_fit_expected_gene_variances', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_fit_expected_gene_variances', 'time', attempt)
  log:
    f"{logs_dir}/join_zoom_fit_expected_gene_variances_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_zoom_fit_expected_gene_variances_{DATE_STAMP}.txt"
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 {input.JOIN_ZOOM_SCRIPT_F:q} trend --inputs {input.stats:q} --output {output.trend:q} --batch-outputs {output.batches:q}
    """

rule join_zoom_calculate_clipped_gene_count_statistics:
  input:
    cells              = lambda wildcards: JOIN_ZOOM_PARAMS['batches'][wildcards.batch]['cells_f'],
    h5ad               = lambda wildcards: JOIN_ZOOM_PARAMS['batches'][wildcards.batch]['path'],
    trend              = zoom_batch_trend_pattern,
    JOIN_ZOOM_SCRIPT_F = JOIN_ZOOM_SCRIPT_F
  output: zoom_clip_pattern
  params: chunk_size = JOIN_ZOOM_PARAMS['hvg']['hvg_chunk_size']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_calculate_clipped_gene_count_statistics', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_calculate_clipped_gene_count_statistics', 'time', attempt)
  log:
    f"{logs_dir}/join_zoom_calculate_clipped_gene_count_statistics_{{batch}}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_zoom_calculate_clipped_gene_count_statistics_{{batch}}_{DATE_STAMP}.txt"
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 {input.JOIN_ZOOM_SCRIPT_F:q} clip --h5ad {input.h5ad:q} --cells {input.cells:q} \
      --chunk-size {params.chunk_size} --trend {input.trend:q} --output {output:q}
    """

rule join_zoom_rank_and_select_highly_variable_genes:
  input:
    stats              = expand(zoom_clip_pattern, batch=ZOOM_BATCHES),
    trend              = zoom_trend_f,
    JOIN_ZOOM_SCRIPT_F = JOIN_ZOOM_SCRIPT_F,
    HVGS_SCRIPT_F      = HVGS_SCRIPT_F,
    exclude            = [JOIN_ZOOM_PARAMS['hvg']['hvg_exclude_from_file']] if JOIN_ZOOM_PARAMS['hvg'].get('hvg_exclude_from_file') else []
  output: HVG_FS[0]
  params:
    n_hvgs       = N_HVGS,
    exclude_flag = '--exclude' if JOIN_ZOOM_PARAMS['hvg'].get('hvg_exclude_from_file') else ''
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_rank_and_select_highly_variable_genes', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_rank_and_select_highly_variable_genes', 'time', attempt)
  log:
    f"{logs_dir}/join_zoom_rank_and_select_highly_variable_genes_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_zoom_rank_and_select_highly_variable_genes_{DATE_STAMP}.txt"
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 {input.JOIN_ZOOM_SCRIPT_F:q} rank --inputs {input.stats:q} --trend {input.trend:q} \
      --n-hvgs {params.n_hvgs} --output {output:q} {params.exclude_flag} {input.exclude:q}
    """

# Matrix assembly, PCA, integration, pseudobulks, markers and reports use the
# rules from join.smk, included above.
