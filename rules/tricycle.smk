# Fixed-reference tricycle projection and pooled phase estimation.
# Preserve the standalone default target without overriding a parent workflow.
TRICYCLE_DEFAULT_TARGET      = workflow.default_target is None

CELL_CYCLE_ENABLED            = 'cell_cycle' in config
TRICYCLE_SCORES_F             = f'{cell_cycle_dir}/tricycle_scores_{FULL_TAG}_{DATE_STAMP}.csv.gz'
TRICYCLE_MARKER_EXPRESSION_F  = f'{cell_cycle_dir}/tricycle_marker_expression_{FULL_TAG}_{DATE_STAMP}.csv.gz'
TRICYCLE_ORIGIN_F             = f'{cell_cycle_dir}/tricycle_origin_{FULL_TAG}_{DATE_STAMP}.csv'
TRICYCLE_ORIGIN_SENSITIVITY_F = f'{cell_cycle_dir}/tricycle_origin_sensitivity_{FULL_TAG}_{DATE_STAMP}.csv'
TRICYCLE_DIAGNOSTICS_F        = f'{cell_cycle_dir}/tricycle_origin_diagnostics_{FULL_TAG}_{DATE_STAMP}.csv.gz'

TRICYCLE_ELLIPSE_F            = f'{cell_cycle_dir}/tricycle_blob_{FULL_TAG}_{DATE_STAMP}.csv'

CELL_CYCLE_SCRIPT_F           = f'{scprocess_dir}/scripts/cell_cycle.R'
TRICYCLE_SCRIPT_F             = f'{scprocess_dir}/scripts/tricycle.py'


rule calculate_tricycle_projection_one_run:
  input:
    counts_h5_f = lambda wildcards: get_filtered_counts_file(
      wildcards.run, config['ambient']['ambient_method'], amb_dir, DATE_STAMP),
    coldata_f   = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    rowdata_f   = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    scores_f            = f'{cell_cycle_dir}/run_{{run}}_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    marker_expression_f = f'{cell_cycle_dir}/run_{{run}}_marker_expression_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    summary_f           = f'{cell_cycle_dir}/run_{{run}}_summary_{FULL_TAG}_{DATE_STAMP}.csv'
  params:
    species = lambda wildcards: config['cell_cycle']['_species'],
    run_var = RUN_VAR
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'calculate_tricycle_projection_one_run', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'calculate_tricycle_projection_one_run', 'time', attempt)
  conda:
    '../envs/tricycle.yaml'
  benchmark:
    f'{benchmark_dir}/cell_cycle/tricycle_run_{{run}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/cell_cycle/tricycle_run_{{run}}_{DATE_STAMP}.log'
  shell: """
    exec &>> {log}
    Rscript --vanilla -e 'source("{scprocess_dir}/scripts/tricycle.R"); estimate_run_tricycle(
      counts_h5_f = "{input.counts_h5_f}",
      coldata_f = "{input.coldata_f}", rowdata_f = "{input.rowdata_f}",
      run_column = "{params.run_var}", run_value = "{wildcards.run}",
      species = "{params.species}", out_scores_f = "{output.scores_f}",
      out_marker_expression_f = "{output.marker_expression_f}",
      out_summary_f = "{output.summary_f}")'
    """


rule aggregate_tricycle:
  input:
    score_fs             = expand(
      f'{cell_cycle_dir}/run_{{run}}_{FULL_TAG}_{DATE_STAMP}.csv.gz', run=RUNS),
    marker_expression_fs = expand(
      f'{cell_cycle_dir}/run_{{run}}_marker_expression_{FULL_TAG}_{DATE_STAMP}.csv.gz',
      run=RUNS),
    summary_fs           = expand(
      f'{cell_cycle_dir}/run_{{run}}_summary_{FULL_TAG}_{DATE_STAMP}.csv', run=RUNS),
    coldata_f            = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    hvg_mat_f            = f'{hvg_dir}/top_hvgs_counts_{FULL_TAG}_{DATE_STAMP}.h5',
    dbl_hvg_mat_f        = f'{hvg_dir}/top_hvgs_doublet_counts_{FULL_TAG}_{DATE_STAMP}.h5',
    tricycle_script_f    = TRICYCLE_SCRIPT_F
  output:
    scores_f            = TRICYCLE_SCORES_F,
    marker_expression_f = TRICYCLE_MARKER_EXPRESSION_F,
    origin_f            = TRICYCLE_ORIGIN_F,
    sensitivity_f       = TRICYCLE_ORIGIN_SENSITIVITY_F,
    diagnostics_f       = TRICYCLE_DIAGNOSTICS_F
  params:
    score_args                   = lambda wildcards, input: ' '.join(f'--score_f {path}' for path in input.score_fs),
    marker_expression_args       = lambda wildcards, input: ' '.join(
      f'--marker_expression_f {path}' for path in input.marker_expression_fs),
    summary_args                 = lambda wildcards, input: ' '.join(f'--summary_f {path}' for path in input.summary_fs),
    bandwidth_multiplier         = lambda wildcards: config['cell_cycle']['cyc_bandwidth_multiplier'],
    kde_grid_size                = lambda wildcards: config['cell_cycle']['cyc_kde_grid_size'],
    density_percentile           = lambda wildcards: config['cell_cycle']['cyc_density_percentile'],
    density_percentile_grid_args = lambda wildcards: ' '.join(
      f"--density_percentile_grid {value}"
      for value in config['cell_cycle']['cyc_density_percentile_grid']),
    origin_args                  = lambda wildcards: (
      '--origin_manual ' + ' '.join(str(value) for value in config['cell_cycle']['cyc_origin_manual'])
      if 'cyc_origin_manual' in config['cell_cycle'] else ''),
    seed                         = lambda wildcards: config['cell_cycle']['cyc_seed']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'aggregate_tricycle', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'aggregate_tricycle', 'time', attempt)
  conda:
    '../envs/hvgs.yaml'
  benchmark:
    f'{benchmark_dir}/cell_cycle/tricycle_aggregate_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/cell_cycle/tricycle_aggregate_{DATE_STAMP}.log'
  shell: """
    exec &>> {log}
    python3 {input.tricycle_script_f} \
      {params.score_args} \
      {params.marker_expression_args} \
      {params.summary_args} \
      --coldata_f {input.coldata_f} \
      --required_h5_f {input.hvg_mat_f} \
      --required_h5_f {input.dbl_hvg_mat_f} \
      --bandwidth_multiplier {params.bandwidth_multiplier} \
      --kde_grid_size {params.kde_grid_size} \
      --density_percentile {params.density_percentile} \
      {params.density_percentile_grid_args} \
      {params.origin_args} \
      --seed {params.seed} \
      --out_scores_f {output.scores_f} \
      --out_marker_expression_f {output.marker_expression_f} \
      --out_origin_f {output.origin_f} \
      --out_sensitivity_f {output.sensitivity_f} \
      --out_diagnostics_f {output.diagnostics_f}
    """


rule estimate_tricycle_reference_ellipse:
  input:
    diagnostics_f       = TRICYCLE_DIAGNOSTICS_F,
    cell_cycle_script_f = CELL_CYCLE_SCRIPT_F
  output:
    ellipse_f = TRICYCLE_ELLIPSE_F
  default_target: TRICYCLE_DEFAULT_TARGET
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'estimate_tricycle_reference_ellipse', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(
      RESOURCE_PARAMS, rules, input, 'estimate_tricycle_reference_ellipse', 'time', attempt)
  log:
    f'{logs_dir}/cell_cycle/tricycle_blob_{DATE_STAMP}.log'
  benchmark:
    f'{benchmark_dir}/cell_cycle/tricycle_blob_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/rlibs.yaml'
  shell: """
    exec &>> {log}
    Rscript --vanilla -e "source('{input.cell_cycle_script_f}'); setDTthreads({threads}); \
      estimate_cell_cycle_blob('{input.diagnostics_f}', '{output.ellipse_f}')"
    """
