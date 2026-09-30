# Joined zooms prepare source tables, then reuse the existing join rules unchanged
# for matrix assembly, native PCA, integration, pseudobulks, markers and reports.
import json
import os
import pathlib
import sys

_scprocess_dir = pathlib.Path(config['scprocess_dir'])
sys.path.insert(0, str(_scprocess_dir / 'scripts'))
from scprocess_utils import get_join_zoom_context, get_join_zoom_prepare_params

_zoom_name = config.pop('join_zoom_name')
config.pop('scprocess_dir')
config, JOIN_ZOOM_CONTEXT = get_join_zoom_context(
    config, _zoom_name, _scprocess_dir, pathlib.Path(os.environ['SCPROCESS_DATA_DIR']))
_zoom_effective_config = json.dumps(config, default=str, sort_keys=True)
config['scprocess_dir'] = str(_scprocess_dir)

include: 'join.smk'

_z = JOIN_ZOOM_CONTEXT
_zoom_context_json = json.dumps(_z, sort_keys=True)
_zoom_prepare_params = get_join_zoom_prepare_params(config, _z)
_zoom_script = str(scprocess_dir / 'scripts/join_zoom.py')
_stats_pattern = f"{join_int_dir}/hvg_stats/{{batch}}.csv.gz"
_clip_pattern = f"{join_int_dir}/hvg_clipped/{{batch}}.csv.gz"
_trend_f = f"{join_int_dir}/hvg_trend.csv.gz"
_batch_trend_pattern = f"{join_int_dir}/hvg_trends/{{batch}}.csv.gz"
_zoom_batches = list(_z['batches'])

rule join_zoom_prepare:
  input:
    parent_files = list(dict.fromkeys(_z['paths'].values())),
    h5ads = [b['path'] for b in _z['batches'].values()],
    script_f = _zoom_script
  output:
    cells = _z['selected_f'],
    summary = _z['summary_f'],
    batch_cells = [b['cells_f'] for b in _z['batches'].values()],
    integrations = INTEGRATED_FS,
    metadata = SAMPLE_META_FS,
    qc = QC_ALL_FS,
    manifests = H5ADS_YAML_FS
  params:
    context = _zoom_prepare_params['context'],
    effective = _zoom_prepare_params['effective']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_prepare', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_prepare', 'time', attempt)
  log: f"{logs_dir}/join_zoom_prepare_{DATE_STAMP}.log"
  benchmark: f"{benchmark_dir}/join_zoom_prepare_{DATE_STAMP}.txt"
  conda: '../envs/hvgs.yaml'
  # A script transports the manifest through Snakemake's job object, avoiding
  # shell argument limits for joins with thousands of H5AD files.
  script: '../scripts/join_zoom.py'

localrules: join_zoom_record_config

rule join_zoom_record_config:
  """Record requested settings without invalidating analytical outputs."""
  input:
    cells = _z['selected_f']
  output:
    effective = _z['effective_f']
  params:
    context = _zoom_context_json,
    effective = _zoom_effective_config
  log: f"{logs_dir}/join_zoom_record_config_{DATE_STAMP}.log"
  script: '../scripts/join_zoom_record_config.py'


rule join_zoom_hvg_stats:
  input:
    cells = lambda wc: _z['batches'][wc.batch]['cells_f'],
    h5ad = lambda wc: _z['batches'][wc.batch]['path'],
    script_f = _zoom_script
  output: _stats_pattern
  params: chunk_size = _z['hvg']['hvg_chunk_size']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_hvg_stats', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_hvg_stats', 'time', attempt)
  log: f"{logs_dir}/join_zoom_hvg_stats_{{batch}}_{DATE_STAMP}.log"
  benchmark: f"{benchmark_dir}/join_zoom_hvg_stats_{{batch}}_{DATE_STAMP}.txt"
  conda: '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 {input.script_f:q} stats --h5ad {input.h5ad:q} --cells {input.cells:q} \
      --chunk-size {params.chunk_size} --output {output:q}
    """

rule join_zoom_hvg_trend:
  input:
    stats = expand(_stats_pattern, batch=_zoom_batches),
    script_f = _zoom_script,
    hvg_script = str(scprocess_dir / 'scripts/hvgs.py')
  output:
    trend = _trend_f,
    batches = expand(_batch_trend_pattern, batch=_zoom_batches)
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_hvg_trend', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_hvg_trend', 'time', attempt)
  log: f"{logs_dir}/join_zoom_hvg_trend_{DATE_STAMP}.log"
  benchmark: f"{benchmark_dir}/join_zoom_hvg_trend_{DATE_STAMP}.txt"
  conda: '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 {input.script_f:q} trend --inputs {input.stats:q} --output {output.trend:q} --batch-outputs {output.batches:q}
    """

rule join_zoom_hvg_clip:
  input:
    cells = lambda wc: _z['batches'][wc.batch]['cells_f'],
    h5ad = lambda wc: _z['batches'][wc.batch]['path'],
    trend = _batch_trend_pattern,
    script_f = _zoom_script
  output: _clip_pattern
  params: chunk_size = _z['hvg']['hvg_chunk_size']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_hvg_clip', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_hvg_clip', 'time', attempt)
  log: f"{logs_dir}/join_zoom_hvg_clip_{{batch}}_{DATE_STAMP}.log"
  benchmark: f"{benchmark_dir}/join_zoom_hvg_clip_{{batch}}_{DATE_STAMP}.txt"
  conda: '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 {input.script_f:q} clip --h5ad {input.h5ad:q} --cells {input.cells:q} \
      --chunk-size {params.chunk_size} --trend {input.trend:q} --output {output:q}
    """

rule join_zoom_hvgs:
  input:
    stats = expand(_clip_pattern, batch=_zoom_batches),
    trend = _trend_f,
    script_f = _zoom_script,
    hvg_script = str(scprocess_dir / 'scripts/hvgs.py'),
    exclude = [_z['hvg']['hvg_exclude_from_file']] if _z['hvg'].get('hvg_exclude_from_file') else []
  output: HVG_FS[0]
  params:
    n_hvgs = N_HVGS,
    exclude_flag = '--exclude' if _z['hvg'].get('hvg_exclude_from_file') else ''
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_hvgs', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_zoom_hvgs', 'time', attempt)
  log: f"{logs_dir}/join_zoom_hvgs_{DATE_STAMP}.log"
  benchmark: f"{benchmark_dir}/join_zoom_hvgs_{DATE_STAMP}.txt"
  conda: '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 {input.script_f:q} rank --inputs {input.stats:q} --trend {input.trend:q} \
      --n-hvgs {params.n_hvgs} --output {output:q} {params.exclude_flag} {input.exclude:q}
    """
