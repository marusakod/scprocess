"""Rules for fixed-reference tricycle projection and pooled phase estimation."""

CELL_CYCLE_ENABLED = 'cell_cycle' in config
TRICYCLE_SCORES_F = f'{int_dir}/tricycle_scores_{FULL_TAG}_{DATE_STAMP}.csv.gz'
TRICYCLE_ORIGIN_F = f'{int_dir}/tricycle_origin_{FULL_TAG}_{DATE_STAMP}.csv'
TRICYCLE_DIAGNOSTICS_F = f'{int_dir}/tricycle_origin_diagnostics_{FULL_TAG}_{DATE_STAMP}.csv.gz'


def _tricycle_runs_for_sample(sample):
  return [run for run in RUNS if sample in RUNS_TO_SAMPLES[run]]


def _tricycle_count_paths(runs):
  return [get_filtered_counts_file(run, config['ambient']['ambient_method'], amb_dir, DATE_STAMP)
          for run in runs]


def _tricycle_named_r_vector(runs):
  entries = [f'"{run}"="{path}"' for run, path in zip(runs, _tricycle_count_paths(runs))]
  return 'c(' + ','.join(entries) + ')'


TRICYCLE_FALLBACK_RUNS = RUNS if (
  config['multiplexing']['demux_type'] in ('hto', 'custom')
) else []


if CELL_CYCLE_ENABLED:
  rule estimate_sample_tricycle:
    input:
      counts_h5_fs = lambda wildcards: _tricycle_count_paths(_tricycle_runs_for_sample(wildcards.sample)),
      coldata_f = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz',
      rowdata_f = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
    output:
      scores_f = f'{int_dir}/tricycle/sample_{{sample}}_{FULL_TAG}_{DATE_STAMP}.csv.gz',
      summary_f = f'{int_dir}/tricycle/sample_{{sample}}_summary_{FULL_TAG}_{DATE_STAMP}.csv'
    params:
      counts_r = lambda wildcards: _tricycle_named_r_vector(_tricycle_runs_for_sample(wildcards.sample)),
      species = config['cell_cycle']['species']
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb = lambda wildcards, attempt, input: get_resources(
        RESOURCE_PARAMS, rules, input, 'run_integration', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(
        RESOURCE_PARAMS, rules, input, 'run_integration', 'time', attempt)
    conda:
      '../envs/tricycle.yaml'
    benchmark:
      f'{benchmark_dir}/integration/tricycle_sample_{{sample}}_{DATE_STAMP}.benchmark.txt'
    log:
      f'{logs_dir}/integration/tricycle_sample_{{sample}}_{DATE_STAMP}.log'
    shell: """
      exec &>> {log}
      Rscript -e "source('{scprocess_dir}/scripts/tricycle.R'); estimate_sample_tricycle(
        counts_h5_fs = {params.counts_r},
        coldata_f = '{input.coldata_f}', rowdata_f = '{input.rowdata_f}',
        sample_id = '{wildcards.sample}',
        species = '{params.species}', out_scores_f = '{output.scores_f}',
        out_summary_f = '{output.summary_f}')"
      """


  if TRICYCLE_FALLBACK_RUNS:
    rule estimate_unassigned_doublet_tricycle:
      input:
        counts_h5_f = lambda wildcards: _tricycle_count_paths([wildcards.run])[0],
        coldata_f = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz',
        rowdata_f = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
      output:
        scores_f = f'{int_dir}/tricycle/unassigned_{{run}}_{FULL_TAG}_{DATE_STAMP}.csv.gz',
        summary_f = f'{int_dir}/tricycle/unassigned_{{run}}_summary_{FULL_TAG}_{DATE_STAMP}.csv'
      params:
        counts_r = lambda wildcards: _tricycle_named_r_vector([wildcards.run]),
        species = config['cell_cycle']['species'],
        run_var = RUN_VAR
      threads: 1
      retries: config['resources']['retries']
      resources:
        mem_mb = lambda wildcards, attempt, input: get_resources(
          RESOURCE_PARAMS, rules, input, 'run_integration', 'memory', attempt),
        runtime = lambda wildcards, attempt, input: get_resources(
          RESOURCE_PARAMS, rules, input, 'run_integration', 'time', attempt)
      conda:
        '../envs/tricycle.yaml'
      benchmark:
        f'{benchmark_dir}/integration/tricycle_unassigned_{{run}}_{DATE_STAMP}.benchmark.txt'
      log:
        f'{logs_dir}/integration/tricycle_unassigned_{{run}}_{DATE_STAMP}.log'
      shell: """
        exec &>> {log}
        Rscript -e "source('{scprocess_dir}/scripts/tricycle.R'); estimate_tricycle_group(
          counts_h5_fs = {params.counts_r}, coldata_f = '{input.coldata_f}',
          rowdata_f = '{input.rowdata_f}',
          group_column = '{params.run_var}', group_value = '{wildcards.run}',
          species = '{params.species}', out_scores_f = '{output.scores_f}',
          out_summary_f = '{output.summary_f}', output_unassigned_only = TRUE)"
        """


  rule aggregate_tricycle:
    input:
      score_fs = (
        expand(f'{int_dir}/tricycle/sample_{{sample}}_{FULL_TAG}_{DATE_STAMP}.csv.gz', sample=SAMPLES) +
        expand(f'{int_dir}/tricycle/unassigned_{{run}}_{FULL_TAG}_{DATE_STAMP}.csv.gz', run=TRICYCLE_FALLBACK_RUNS)
      ),
      summary_fs = (
        expand(f'{int_dir}/tricycle/sample_{{sample}}_summary_{FULL_TAG}_{DATE_STAMP}.csv', sample=SAMPLES) +
        expand(f'{int_dir}/tricycle/unassigned_{{run}}_summary_{FULL_TAG}_{DATE_STAMP}.csv', run=TRICYCLE_FALLBACK_RUNS)
      ),
      coldata_f = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz',
      hvg_mat_f = f'{hvg_dir}/top_hvgs_counts_{FULL_TAG}_{DATE_STAMP}.h5',
      dbl_hvg_mat_f = f'{hvg_dir}/top_hvgs_doublet_counts_{FULL_TAG}_{DATE_STAMP}.h5'
    output:
      scores_f = TRICYCLE_SCORES_F,
      origin_f = TRICYCLE_ORIGIN_F,
      diagnostics_f = TRICYCLE_DIAGNOSTICS_F
    params:
      score_args = lambda wildcards, input: ' '.join(f'--score_f {path}' for path in input.score_fs),
      summary_args = lambda wildcards, input: ' '.join(f'--summary_f {path}' for path in input.summary_fs),
      bandwidth_multiplier = config['cell_cycle']['bandwidth_multiplier'],
      kde_grid_size = config['cell_cycle']['kde_grid_size'],
      target_cells = config['cell_cycle']['target_cells'],
      seed = config['cell_cycle']['seed']
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb = 8000,
      runtime = 30
    conda:
      '../envs/hvgs.yaml'
    benchmark:
      f'{benchmark_dir}/integration/tricycle_aggregate_{DATE_STAMP}.benchmark.txt'
    log:
      f'{logs_dir}/integration/tricycle_aggregate_{DATE_STAMP}.log'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/tricycle.py \
        {params.score_args} \
        {params.summary_args} \
        --coldata_f {input.coldata_f} \
        --required_h5_f {input.hvg_mat_f} \
        --required_h5_f {input.dbl_hvg_mat_f} \
        --bandwidth_multiplier {params.bandwidth_multiplier} \
        --kde_grid_size {params.kde_grid_size} \
        --target_cells {params.target_cells} \
        --seed {params.seed} \
        --out_scores_f {output.scores_f} \
        --out_origin_f {output.origin_f} \
        --out_diagnostics_f {output.diagnostics_f}
      """
