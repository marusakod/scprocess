# QC and doublet identification for each run
from contextlib import redirect_stdout, redirect_stderr
import polars as pl
import pandas as pd
import gzip
import os
import numpy as np

localrules: make_qc_thresholds_csv, check_qc_quality

# mini wrapper function
def get_all_batches_str(wildcards):
  return ','.join(RUNS_TO_BATCHES[wildcards.run])


rule make_qc_thresholds_csv:
  input:
    metadata_f = config['project']['sample_metadata'],
    custom_f   = config['project']['custom_sample_params'] if 'custom_sample_params' in config['project'] else []
  output:
    cuts_f = f'{qc_dir}/qc_thresholds_by_{BATCH_VAR}_{FULL_TAG}_{DATE_STAMP}.csv'
  params:
    qc_by_batch = {batch: values['qc'] for batch, values in BATCH_PARAMS.items()},
    batch_var   = BATCH_VAR
  log:
    f'{logs_dir}/qc/make_qc_thresholds_csv_{DATE_STAMP}.log'
  run:
    with open(str(log), "a") as f, redirect_stdout(f), redirect_stderr(f):
      write_qc_thresholds(params.qc_by_batch, params.batch_var, output.cuts_f)


rule check_demux_ids:
  input:
    metadata_f = config['project']['sample_metadata'],
    demux_f    = config['multiplexing']['demux_output'] if config['multiplexing']['demux_type'] == 'custom' else []
  output:
    check_f = f'{qc_dir}/demux_ids_check_{FULL_TAG}_{DATE_STAMP}.ok'
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'check_demux_ids', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'check_demux_ids', 'time', attempt)
  log:
    f'{logs_dir}/qc/check_demux_ids_{DATE_STAMP}.log'
  conda:
    '../envs/py_env.yaml'
  shell: """
    exec &>> {log}

    PYTHONPATH={scprocess_dir}/scripts python -c "from scprocess_utils import check_demux_ids; \
      check_demux_ids('{input.demux_f}', '{input.metadata_f}', '{output.check_f}')"
    """


rule run_qc_one_run:
  input:
    af_gtf_dt_f = config['mapping_af']['gene_info_f'],
    metadata_f  = config['project']['sample_metadata'],
    af_h5_f     = f'{af_dir}/af_{{run}}/{af_rna_dir}af_counts_mat.h5',
    cuts_f      = f'{qc_dir}/qc_thresholds_by_{BATCH_VAR}_{FULL_TAG}_{DATE_STAMP}.csv',
    run_stats_f = f'{amb_dir}/ambient_run_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    amb_yaml_f  = f'{amb_dir}/ambient_{{run}}/ambient_{{run}}_{DATE_STAMP}_output_paths.yaml',
    demux_check = f'{qc_dir}/demux_ids_check_{FULL_TAG}_{DATE_STAMP}.ok' \
      if config['multiplexing']['demux_type'] == 'custom' else [],
    demux_f     = (f'{demux_dir}/sce_cells_htos_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds') \
      if config['multiplexing']['demux_type'] == 'hto' else \
      config['multiplexing']['demux_output'] if config['multiplexing']['demux_type'] == 'custom' else []
  output:
    qc_f             = temp(f'{qc_dir}/tmp_qc_dt_{{run}}_{FULL_TAG}_{DATE_STAMP}.csv.gz'),
    coldata_f        = temp(f'{qc_dir}/tmp_coldata_dt_{{run}}_{FULL_TAG}_{DATE_STAMP}.csv.gz'),
    rowdata_f        = temp(f'{qc_dir}/tmp_rowdata_dt_{{run}}_{FULL_TAG}_{DATE_STAMP}.csv.gz'),
    dbl_f            = f'{dbl_dir}/dbl_{{run}}/scDblFinder_{{run}}_outputs_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    dbl_artificial_f = f'{dbl_dir}/dbl_{{run}}/scDblFinder_{{run}}_artificial_doublets_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    dbl_summary_f    = f'{dbl_dir}/dbl_{{run}}/scDblFinder_{{run}}_summary_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  params:
    all_batches_str = get_all_batches_str,
    mito_str        = config['mapping_af']['af_mito_str'],
    ambient_method  = config['ambient']['ambient_method'],
    exclude_mito    = config['qc']['exclude_mito'],
    hard_min_counts = config['qc']['qc_hard_min_counts'],
    hard_min_feats  = config['qc']['qc_hard_min_feats'],
    hard_max_mito   = config['qc']['qc_hard_max_mito'],
    run_var         = RUN_VAR,
    demux_type      = config['multiplexing']['demux_type'],
    batch_var       = BATCH_VAR,
    dbl_min_feats   = config['qc']['dbl_min_feats'],
    dbl_dbr_sd      = lambda wildcards: 'NULL' if RUN_PARAMS[wildcards.run]['qc']['dbl_dbr_sd'] is None else RUN_PARAMS[wildcards.run]['qc']['dbl_dbr_sd'],
    dbl_clusters    = 'NULL' if config['qc']['dbl_clusters'] is None else str(config['qc']['dbl_clusters']).upper()
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_qc_one_run', 'memory', attempt, wildcards.run),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_qc_one_run', 'time', attempt, wildcards.run)
  benchmark:
    f'{benchmark_dir}/qc/run_qc_{{run}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/qc/run_qc_one_run_{{run}}_{DATE_STAMP}.log'
  conda:
    '../envs/rlibs.yaml'
  shell: """
    exec &>> {log}

    Rscript --vanilla -e "source('{scprocess_dir}/scripts/SampleQC.R'); source('{scprocess_dir}/scripts/utils.R'); \
      main_qc( \
        run_name         = '{wildcards.run}', \
        metadata_f       = '{input.metadata_f}', \
        cuts_f           = '{input.cuts_f}', \
        amb_yaml_f       = '{input.amb_yaml_f}', \
        run_stats_f      = '{input.run_stats_f}', \
        demux_f          = '{input.demux_f}', \
        gtf_dt_f         = '{input.af_gtf_dt_f}', \
        ambient_method   = '{params.ambient_method}', \
        all_batches_str  = '{params.all_batches_str}', \
        rowdata_f        = '{output.rowdata_f}', \
        qc_f             = '{output.qc_f}', \
        coldata_f        = '{output.coldata_f}', \
        dbl_f            = '{output.dbl_f}', \
        dbl_artificial_f = '{output.dbl_artificial_f}', \
        dbl_summary_f    = '{output.dbl_summary_f}', \
        mito_str         = '{params.mito_str}', \
        exclude_mito     = '{params.exclude_mito}', \
        hard_min_counts  = {params.hard_min_counts}, \
        hard_min_feats   = {params.hard_min_feats}, \
        hard_max_mito    = {params.hard_max_mito}, \
        run_var          = '{params.run_var}', \
        demux_type       = '{params.demux_type}', \
        batch_var        = '{params.batch_var}', \
        dbl_min_feats    = {params.dbl_min_feats}, \
        dbl_dbr_sd       = {params.dbl_dbr_sd}, \
        dbl_clusters     = {params.dbl_clusters} \
      )"
    """


rule merge_doublet_diagnostics:
  input:
    real_fs       = expand(
      f'{dbl_dir}/dbl_{{run}}/scDblFinder_{{run}}_outputs_{FULL_TAG}_{DATE_STAMP}.csv.gz',
      run = RUNS),
    artificial_fs = expand(
      f'{dbl_dir}/dbl_{{run}}/scDblFinder_{{run}}_artificial_doublets_{FULL_TAG}_{DATE_STAMP}.csv.gz',
      run = RUNS),
    summary_fs    = expand(
      f'{dbl_dir}/dbl_{{run}}/scDblFinder_{{run}}_summary_{FULL_TAG}_{DATE_STAMP}.csv.gz',
      run = RUNS)
  output:
    real_f       = f'{dbl_dir}/scDblFinder_real_all_runs_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    artificial_f = f'{dbl_dir}/scDblFinder_artificial_doublets_all_runs_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    summary_f    = f'{dbl_dir}/scDblFinder_summary_all_runs_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'merge_doublet_diagnostics', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'merge_doublet_diagnostics', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/qc/merge_doublet_diagnostics_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/qc/merge_doublet_diagnostics_{DATE_STAMP}.log'
  conda:
    '../envs/py_env.yaml'
  shell:
    r"""
    exec &> {log}

    python {scprocess_dir}/scripts/merge_doublet_diagnostics.py \
      --real-output {output.real_f:q} \
      --artificial-output {output.artificial_f:q} \
      --summary-output {output.summary_f:q} \
      --real-inputs {input.real_fs:q} \
      --artificial-inputs {input.artificial_fs:q} \
      --summary-inputs {input.summary_fs:q}
    """


rule merge_qc:
  input:
    qc_fs      = expand(f'{qc_dir}/tmp_qc_dt_{{run}}_{FULL_TAG}_{DATE_STAMP}.csv.gz', run = RUNS),
    coldata_fs = expand(f'{qc_dir}/tmp_coldata_dt_{{run}}_{FULL_TAG}_{DATE_STAMP}.csv.gz', run = RUNS)
  output:
    qc_merged_f      = f'{qc_dir}/qc_all_samples_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    coldata_merged_f = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  params:
    metadata_f = config['project']['sample_metadata']
  threads: 1
  retries: config['resources']['retries']
  benchmark:
    f'{benchmark_dir}/qc/merge_qc_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/qc/merge_qc_{DATE_STAMP}.log'
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'merge_qc', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'merge_qc', 'time', attempt)
  conda:
    '../envs/rlibs.yaml'
  shell: """
    exec &>> {log}

    QC_LIST="{input.qc_fs}" COL_LIST="{input.coldata_fs}" \
      Rscript --vanilla -e "source('{scprocess_dir}/scripts/SampleQC.R'); \
        qc_files  = trimws(strsplit(Sys.getenv('QC_LIST'), ' ')[[1]]); \
        col_files = trimws(strsplit(Sys.getenv('COL_LIST'), ' ')[[1]]); \
        merge_qc_files( \
          qc_files  = qc_files, \
          col_files = col_files, \
          qc_out    = '{output.qc_merged_f}', \
          col_out   = '{output.coldata_merged_f}' \
        )"
    """


rule merge_rowdata:
  input:
    rowdata_fs = expand(f'{qc_dir}/tmp_rowdata_dt_{{run}}_{FULL_TAG}_{DATE_STAMP}.csv.gz', run = RUNS)
  output:
    rowdata_merged_f = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'merge_rowdata', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'merge_rowdata', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/qc/merge_rowdata_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/qc/merge_rowdata_{DATE_STAMP}.log'
  run:
    with open(str(log), "a") as f, redirect_stdout(f), redirect_stderr(f):

      # read all nonempty rowdata files
      rows_df_ls  = [pl.read_csv(f) for f in input.rowdata_fs if os.path.getsize(f) > 0 ]

      # check if identical
      first_df    = rows_df_ls[0]
      all_ident   = all(first_df.equals(df) for df in rows_df_ls[1:])

      # save only one df
      if all_ident:
        with gzip.open(output.rowdata_merged_f, 'wb') as f:
          first_df.write_csv(f)
      else:
        raise ValueError("error: rowdata for all sce objects not identical.")


checkpoint get_qc_sample_statistics:
  input:
    run_stats_f = f'{amb_dir}/ambient_run_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    qc_merged_f = f'{qc_dir}/qc_all_samples_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    cuts_f      = f'{qc_dir}/qc_thresholds_by_{BATCH_VAR}_{FULL_TAG}_{DATE_STAMP}.csv'
  output:
    qc_stats_f = f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv'
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'get_qc_sample_statistics', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'get_qc_sample_statistics', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/qc/get_qc_sample_statistics_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/qc/get_qc_sample_statistics_{DATE_STAMP}.log'
  run:
    with open(str(log), "a") as f, redirect_stdout(f), redirect_stderr(f):

      sample_stats_df = extract_qc_sample_statistics(input.run_stats_f, input.qc_merged_f, input.cuts_f,
        config, BATCHES, RUNS_TO_BATCHES, BATCH_VAR, RUN_VAR)
      sample_stats_df.write_csv(output.qc_stats_f)


rule check_qc_quality:
  input:
    qc_stats_f = f"{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv"
  output:
    qc_flag_f = touch(f"{qc_dir}/qc_passed_{FULL_TAG}_{DATE_STAMP}.flag")
  params:
    batch_var = BATCH_VAR
  log:
    f'{logs_dir}/qc/check_qc_quality_{DATE_STAMP}.log'
  run:
    with open(str(log), "a") as f, redirect_stdout(f), redirect_stderr(f):

      import polars as pl
      qc_df     = pl.read_csv(input.qc_stats_f)
      num_ok    = qc_df[f'bad_{params.batch_var}'].not_().sum()

      if num_ok < 2:
        # This will crash the job and prevent downstream rules from starting
        raise Exception(f"CRITICAL: Only {num_ok} samples passed QC. Stopping.")
