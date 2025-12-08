# snakemake rule for doing QC on sce object
import polars as pl
import pandas as pd
import gzip
import os
import numpy as np

localrules: make_qc_thresholds_csv, make_tmp_sce_paths_yaml

# get output file paths as string
def _get_qc_files_str(run, RUNS_TO_BATCHES, qc_dir, FULL_TAG, DATE_STAMP):
  # make lists
  sce_fs_ls   = []
  batches_ls  = []
  for b in RUNS_TO_BATCHES[run]:
    sce_fs_ls.append(f"{qc_dir}/sce_cells_tmp_{b}_{FULL_TAG}_{DATE_STAMP}.rds")
    batches_ls.append(b)

  # concatenate them
  sce_str   = ','.join(sce_fs_ls)
  batch_str  = ','.join(batches_ls)

  # make out dictionary
  out_dc = {
    "batch_str":  batch_str, 
    "sce_str":    sce_str
  }
  return out_dc


# mini wrapper functions
def get_all_batches_str(wildcards):
  return _get_qc_files_str(wildcards.run, RUNS_TO_BATCHES, qc_dir, FULL_TAG, DATE_STAMP)['batch_str']


def get_sce_fs_str(wildcards):
  return _get_qc_files_str(wildcards.run, RUNS_TO_BATCHES, qc_dir, FULL_TAG, DATE_STAMP)['sce_str']


def extract_qc_sample_statistics(run_stats_f, qc_merged_f, cuts_f, config, BATCHES, RUNS_TO_BATCHES, BATCH_VAR, RUN_VAR):
  # load the merged qc file, also thresholds
  qc_df     = pl.read_csv(qc_merged_f)
  cuts_df   = pl.read_csv(cuts_f)

  # count the number of cells that passed qc per sample
  qc_df     = qc_df.filter(pl.col("keep") == True)
  sample_df = qc_df.group_by(BATCH_VAR).agg( pl.col(BATCH_VAR).count().alias('n_cells') )

  cuts_sample_df = cuts_df.select(BATCH_VAR)
  sample_df = cuts_sample_df.join(sample_df, on=BATCH_VAR, how='left') 
  sample_df = sample_df.fill_null(0)

  # check all samples present in sample_df
  if set(sample_df[BATCH_VAR].to_list()) != set(cuts_df[BATCH_VAR].to_list()):
    raise ValueError("sample_df and cuts_df have different sets of values for BATCH_VAR")

  # load cuts, merge and use to identify samples that do not meet the minimum cell threshold
  sample_df = sample_df.join( cuts_df.select([BATCH_VAR, "min_cells"]), on = BATCH_VAR, how = "left")
  sample_df = sample_df.with_columns( (pl.col('n_cells') < pl.col('min_cells')).alias('bad_qc') )

  # handle samples excluded after cellbender
  if config['ambient']['ambient_method'] == 'cellbender':
    # load ambient sample stats
    amb_stats   = pl.read_csv(run_stats_f)

    # get bad pools or samples
    bad_bender  = amb_stats.filter( pl.col('bad_run') )[ RUN_VAR ].to_list()

    # check for these
    if config['multiplexing']['demux_type'] != "none":
      bad_bender_batches = []
      for bad_run in bad_bender:
        if bad_run in RUNS_TO_BATCHES:
          bad_bender_batches.extend(RUNS_TO_BATCHES[bad_run])
      if not all(b in BATCHES for b in bad_bender_batches):
        raise ValueError("Some bad bender samples are not in the BATCHES list.")
    else:
      bad_bender_batches = bad_bender

    # record samples where cellbender went wrong
    sample_df = sample_df.with_columns( pl.lit(False).alias('bad_bender'))
    if len(bad_bender_batches) != 0: 
       # add bad_bender column to sample_df
      bad_bender_df = pl.DataFrame({
        BATCH_VAR:    bad_bender_batches,
        'n_cells':    None, 
        'min_cells':  None, 
        'bad_qc':     False, 
        'bad_bender': True
      }).with_columns(
        pl.col("n_cells").cast(sample_df["n_cells"].dtype),
        pl.col("min_cells").cast(sample_df["min_cells"].dtype),
      )
      sample_df = sample_df.vstack(bad_bender_df)

    # check we didn't miss anything
    if not all([b in sample_df[BATCH_VAR].to_list() for b in BATCHES]):
      raise KeyError("some samples missing from sample_df")

    # label as bad if bad_bender or bad_qc
    sample_df   = sample_df.with_columns( (pl.col('bad_bender') | pl.col('bad_qc')).alias(f'bad_{BATCH_VAR}') )

  else:
    sample_df   = sample_df.with_columns( pl.col('bad_qc').alias(f'bad_{BATCH_VAR}') )

  return sample_df


rule make_qc_thresholds_csv:
  output:
    cuts_f      = f'{qc_dir}/qc_thresholds_by_{BATCH_VAR}_{FULL_TAG}_{DATE_STAMP}.csv'
  params:
    batch_var   = BATCH_VAR
  run:
    # make polars dataframe from dictionary of parameters
    rows_data   = []
    for batch, param_ls in BATCH_PARAMS.items():
      qc_params   = param_ls['qc']
      row_data    = {params.batch_var: batch}
      row_data.update(qc_params)
      rows_data.append(row_data)
    cuts_df     = pl.DataFrame(rows_data)
    
    # rename the columns
    def rename_fn(col):
      if col.startswith("qc_"):
        return col.replace("qc_", "")
      else:
        return col
    cuts_df     = cuts_df.rename( lambda col_name: rename_fn(col_name) )

    # save to csv
    cuts_df.write_csv(output.cuts_f)


rule run_qc_one_run:
  input:
    af_h5_f     = af_dir  + '/af_{run}/' + af_rna_dir + 'af_counts_mat.h5', 
    cuts_f      = f'{qc_dir}/qc_thresholds_by_{BATCH_VAR}_{FULL_TAG}_{DATE_STAMP}.csv',
    run_stats_f = amb_dir + '/ambient_run_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    amb_yaml_f  = amb_dir + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml',
    demux_f     = (demux_dir + '/sce_cells_htos_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.rds') \
      if config['multiplexing']['demux_type'] == 'hto' else \
      config['multiplexing']['demux_output'] if config['multiplexing']['demux_type'] == 'custom' else []
  output:
    qc_f         = temp(qc_dir + '/tmp_qc_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'), 
    coldata_f    = temp(qc_dir + '/tmp_coldata_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'),
    rowdata_f    = temp(qc_dir + '/tmp_rowdata_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'), 
    # dimred_f     = dbl_dir + '/dbl_{run}/scDblFinder_{run}_dimreds_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    dbl_f        = f'{dbl_dir}/dbl_{{run}}/scDblFinder_{{run}}_outputs_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  params:
    metadata_f      = config['project']['sample_metadata'],
    af_gtf_dt_f     = config['mapping']['af_gtf_dt_f'],
    all_batches_str = get_all_batches_str,
    sce_fs_str      = get_sce_fs_str,
    mito_str        = config['mapping']['af_mito_str'],
    ambient_method  = config['ambient']['ambient_method'],
    exclude_mito    = config['qc']['exclude_mito'],
    hard_min_counts = config['qc']['qc_hard_min_counts'],
    hard_min_feats  = config['qc']['qc_hard_min_feats'],
    hard_max_mito   = config['qc']['qc_hard_max_mito'],
    run_var         = RUN_VAR,
    demux_type      = config['multiplexing']['demux_type'],
    batch_var       = BATCH_VAR,
    dbl_min_feats   = config['qc']['dbl_min_feats']
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('run_qc', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run)*(1.5**(attempt-1)),
    runtime = lambda wildcards, input: get_resources('run_qc', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_qc/run_qc_{run}_' + DATE_STAMP + '.benchmark.txt'
  conda:
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/SampleQC.R'); source('scripts/utils.R'); \
      main_qc( \
        run_name        = '{wildcards.run}', \
        metadata_f      = '{params.metadata_f}', \
        cuts_f          = '{input.cuts_f}', \
        amb_yaml_f      = '{input.amb_yaml_f}', \
        run_stats_f     = '{input.run_stats_f}', \
        demux_f         = '{input.demux_f}', \
        gtf_dt_f        = '{params.af_gtf_dt_f}', \
        ambient_method  = '{params.ambient_method}', \
        sce_fs_str      = '{params.sce_fs_str}', \
        all_batches_str = '{params.all_batches_str}', \
        rowdata_f       = '{output.rowdata_f}', \
        qc_f            = '{output.qc_f}', \
        coldata_f       = '{output.coldata_f}', \
        dbl_f           = '{output.dbl_f}', \
        mito_str        = '{params.mito_str}', \
        exclude_mito    = '{params.exclude_mito}', \
        hard_min_counts =  {params.hard_min_counts}, \
        hard_min_feats  =  {params.hard_min_feats}, \
        hard_max_mito   =  {params.hard_max_mito}, \
        run_var         = '{params.run_var}', \
        demux_type      = '{params.demux_type}', \
        batch_var       = '{params.batch_var}', \
        dbl_min_feats   =  {params.dbl_min_feats} \
      )"
    """


rule merge_qc:
  input:
    qc_fs      = expand(f'{qc_dir}/tmp_qc_dt_{{run}}_{FULL_TAG}_{DATE_STAMP}.csv.gz', run = RUNS),
    coldata_fs = expand(f'{qc_dir}/tmp_coldata_dt_{{run}}_{FULL_TAG}_{DATE_STAMP}.csv.gz', run = RUNS)
  output:
    qc_merged_f      = f'{qc_dir}/qc_all_samples_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    coldata_merged_f = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  threads: 1
  retries: config['resources']['retries']
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_qc/merge_qc_' + DATE_STAMP + '.benchmark.txt'
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('merge_qc', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('merge_qc', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  run:
    # read all nonempty input files and concatenate them
    qc_df_ls    = [ pl.read_csv(f, schema_overrides = {"log_N": pl.Float64}) for f in input.qc_fs if os.path.getsize(f) > 0 ]
    qc_df_all   = pl.concat(qc_df_ls)

    cols_df_ls  = [ pl.read_csv(f) for f in input.coldata_fs if os.path.getsize(f) > 0 ]
    cols_df_all = pl.concat(cols_df_ls)

    # save merged dataframes to output files
    with gzip.open(output.qc_merged_f, 'wb') as f:
      qc_df_all.write_csv(f)
    with gzip.open(output.coldata_merged_f, 'wb') as f:
      cols_df_all.write_csv(f)


rule merge_rowdata:
  input:
    rowdata_fs = expand(qc_dir + '/tmp_rowdata_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', run = RUNS)
  output:
    rowdata_merged_f = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('merge_rowdata', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('merge_rowdata', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_qc/merge_rowdata_' + DATE_STAMP + '.benchmark.txt'
  run:
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


rule get_qc_sample_statistics:
  input:
    run_stats_f   = f'{amb_dir}/ambient_run_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    qc_merged_f   = f'{qc_dir}/qc_all_samples_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    cuts_f        = f'{qc_dir}/qc_thresholds_by_{BATCH_VAR}_{FULL_TAG}_{DATE_STAMP}.csv'
  output:
    qc_stats_f    = f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv'
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('get_qc_sample_statistics', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('get_qc_sample_statistics', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_qc/get_qc_sample_statistics_' + DATE_STAMP + '.benchmark.txt'
  run:
    sample_stats_df = extract_qc_sample_statistics(input.run_stats_f, input.qc_merged_f, input.cuts_f,
      config, BATCHES, RUNS_TO_BATCHES, BATCH_VAR, RUN_VAR)
    sample_stats_df.write_csv(output.qc_stats_f)

