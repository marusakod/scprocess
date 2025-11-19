# snakemake rule for doing QC on sce object
import polars as pl
import pandas as pd
import gzip
import os
import numpy as np

localrules: make_qc_thresholds_csv, make_tmp_sce_paths_yaml

# get output file paths as string
def get_qc_files_str(run, SAMPLES_TO_RUNS, qc_dir, FULL_TAG, DATE_STAMP):
  # make lists
  sce_fs_ls = []
  smpls_ls  = []
  for s in SAMPLES_TO_RUNS[run]:
    sce_fs_ls.append(f"{qc_dir}/sce_cells_tmp_{s}_{FULL_TAG}_{DATE_STAMP}.rds")
    smpls_ls.append(s)

  # concatenate them
  sce_str   = ','.join(sce_fs_ls)
  smpl_str  = ','.join(smpls_ls)

  # make out dictionary
  out_dc = {
    "smpl_str": smpl_str, 
    "sce_str":  sce_str
  }
  return out_dc


# mini wrapper functions
def get_all_smpls_str(wildcards):
  return get_qc_files_str(wildcards.run, SAMPLES_TO_RUNS, qc_dir, FULL_TAG, DATE_STAMP)['smpl_str']

def get_sce_fs_str(wildcards):
  return get_qc_files_str(wildcards.run, SAMPLES_TO_RUNS, qc_dir, FULL_TAG, DATE_STAMP)['sce_str']


def extract_qc_sample_statistics(run_stats_f, qc_merged_f, cuts_f, config, SAMPLES, RUN_VAR, SAMPLES_TO_RUNS):
  # load the merged qc file, also thresholds
  qc_df     = pl.read_csv(qc_merged_f)
  cuts_df   = pl.read_csv(cuts_f)

  # count the number of cells that passed qc per sample
  qc_df     = qc_df.filter(pl.col("keep") == True)
  sample_df = qc_df.group_by('sample_id').agg( pl.col("sample_id").count().alias('n_cells') )

  # check all samples present in sample_df
  if set(sample_df['sample_id'].to_list()) != set(cuts_df['sample_id'].to_list()):
    raise ValueError("sample_df and cuts_df have different sets of values for 'sample_id'")

  # load cuts, merge and use to identify samples that do not meet the minimum cell threshold
  sample_df = sample_df.join( cuts_df.select(["sample_id", "min_cells"]), on = "sample_id", how = "left")
  sample_df = sample_df.with_columns( (pl.col('n_cells') < pl.col('min_cells')).alias('bad_qc') )

  # handle samples excluded after cellbender
  if config['ambient']['ambient_method'] == 'cellbender':
    # load ambient sample stats
    amb_stats   = pl.read_csv(run_stats_f)

    # get bad pools or samples
    bad_bender  = amb_stats.filter( pl.col('bad_run') )[ RUN_VAR ].to_list()

    # check for these
    if config['multiplexing']['demux_type'] != "none":
      bad_bender_samples = []
      for bad_run in bad_bender:
        if bad_run in SAMPLES_TO_RUNS:
          bad_bender_samples.extend(SAMPLES_TO_RUNS[bad_run])
      if not all(s in SAMPLES for s in bad_bender_samples):
        raise ValueError("Some bad bender samples are not in the SAMPLES list.")
    else:
      bad_bender_samples = bad_bender

    # record samples where cellbender went wrong
    sample_df = sample_df.with_columns( pl.lit(False).alias('bad_bender'))
    if len(bad_bender_samples) != 0: 
       # add bad_bender column to sample_df
      bad_bender_df = pl.DataFrame({
        'sample_id':  bad_bender_samples, 
        'n_cells':    np.nan, 
        'bad_qc':     False, 
        'bad_bender': True
      })
      sample_df = sample_df.vstack(bad_bender_df)

    # check we didn't miss anything
    if not all([s in sample_df['sample_id'].to_list() for s in SAMPLES]):
      raise KeyError("some samples missing from sample_df")

    # label as bad if bad_bender or bad_qc
    sample_df   = sample_df.with_columns( (pl.col('bad_bender') | pl.col('bad_qc')).alias('bad_sample') )

  else:
    sample_df   = sample_df.with_columns( pl.col('bad_qc').alias('bad_sample') )

  return sample_df


rule make_qc_thresholds_csv:
  output:
    cuts_f      = qc_dir  + '/qc_thresholds_by_sample_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  run:
    # make polars dataframe from dictionary of parameters
    rows_data   = []
    for sample_id, param_ls in SAMPLE_PARAMS.items():
      qc_params   = param_ls['qc']
      row_data    = {'sample_id': sample_id}
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


rule run_qc:
  input:
    af_h5_f     = af_dir  + '/af_{run}/' + af_rna_dir + 'af_counts_mat.h5', 
    cuts_f      = qc_dir  + '/qc_thresholds_by_sample_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
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
    dbl_f        = dbl_dir + '/dbl_{run}/scDblFinder_{run}_outputs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  params:
    metadata_f      = config['project']['sample_metadata'],
    af_gtf_dt_f     = config['mapping']['af_gtf_dt_f'],
    all_smpls_str   = get_all_smpls_str,
    sce_fs_str      = get_sce_fs_str,
    mito_str        = config['mapping']['af_mito_str'],
    ambient_method  = config['ambient']['ambient_method'],
    exclude_mito    = config['qc']['exclude_mito'],
    hard_min_counts = config['qc']['qc_hard_min_counts'],
    hard_min_feats  = config['qc']['qc_hard_min_feats'],
    hard_max_mito   = config['qc']['qc_hard_max_mito'],
    run_var         = RUN_VAR,
    demux_type      = config['multiplexing']['demux_type'],
    dbl_min_feats   = config['qc']['dbl_min_feats']
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt, input: (
      attempt * (
      config['resources'].get('gb_run_qc', None) * MB_PER_GB
      if config['resources'].get('gb_run_qc') is not None
      else (1522.029 + 34.180 * (os.path.getsize(input.af_h5_f)//MB_PER_GB**2)) + 2*MB_PER_GB # lm + buffer
      )
    ), 
    runtime = lambda wildcards, input:
      config['resources'].get('min_run_qc', None)
      if config['resources'].get('min_run_qc') is not None
      else (35.03 + 1.687*(os.path.getsize(input.af_h5_f)//MB_PER_GB**2))/60 + 5 # lm + 5 min buffer
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_qc/run_qc_{run}_' + DATE_STAMP + '.benchmark.txt'
  conda:
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/SampleQC.R'); source('scripts/ambient.R'); \
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
        all_samples_str = '{params.all_smpls_str}', \
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
        dbl_min_feats   =  {params.dbl_min_feats} \
      )"
    """


rule merge_qc:
  input:
    qc_fs      = expand(qc_dir + '/tmp_qc_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', run = RUNS),
    coldata_fs = expand(qc_dir + '/tmp_coldata_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', run = RUNS)
  output:
    qc_merged_f      = qc_dir  + '/qc_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz',
    coldata_merged_f = qc_dir  + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  threads: 1
  retries: config['resources']['retries']
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_qc/merge_qc_' + DATE_STAMP + '.benchmark.txt'
  resources:
    mem_mb = lambda wildcards, attempt: (
      attempt * (
      config['resources'].get('gb_merge_qc', None) * MB_PER_GB
      if config['resources'].get('gb_merge_qc') is not None
      else (145.174 + 0.162 * len(SAMPLES)) + 2*MB_PER_GB # lm + buffer
      )
    ), 
    runtime = lambda wildcards:
      config['resources'].get('min_merge_qc', None) 
      if config['resources'].get('min_merge_qc') is not None
      else (-1.151 + 0.603*len(SAMPLES))/60 + 5 # lm + 5 min buffer
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
    mem_mb = lambda wildcards, attempt: (
      attempt * (
      config['resources'].get('gb_merge_rowdata', None) * MB_PER_GB
      if config['resources'].get('gb_merge_rowdata') is not None
      else (111.267 + 7.465 * len(SAMPLES)) + 2*MB_PER_GB # lm + buffer
      )
    ), 
    runtime = lambda wildcards:
      config['resources'].get('min_merge_rowdata', None)
      if config['resources'].get('min_merge_rowdata') is not None
      else (0.593 + 0.056*len(SAMPLES))/60 + 5 # lm + 5 min buffer
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
    run_stats_f   = amb_dir + '/ambient_run_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    qc_merged_f   = qc_dir  + '/qc_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz',
    cuts_f        = qc_dir  + '/qc_thresholds_by_sample_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    qc_stats_f    = qc_dir + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt: (
      attempt * (
      config['resources'].get('gb_get_qc_sample_statistics', None) * MB_PER_GB
      if config['resources'].get('gb_get_qc_sample_statistics') is not None
      else (78.85 + 3.092 * len(SAMPLES)) + 2*MB_PER_GB # lm + buffer
      )
    ), 
    runtime = lambda wildcards:
      config['resources'].get('min_get_qc_sample_statistics', None)
      if config['resources'].get('min_get_qc_sample_statistics') is not None
      else (-0.005 + 0.033*len(SAMPLES))/60 + 5 # lm + 5 min buffer
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_qc/get_qc_sample_statistics_' + DATE_STAMP + '.benchmark.txt'
  run:
    sample_stats_df = extract_qc_sample_statistics(input.run_stats_f, input.qc_merged_f, input.cuts_f,
      config, SAMPLES, RUN_VAR, SAMPLES_TO_RUNS)
    sample_stats_df.write_csv(output.qc_stats_f)


# write sce objects paths to a yaml file
rule make_tmp_sce_paths_yaml:
  input:
    qc_stats_f  = qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv' # so that this runs after get_qc_sample_statistics
  output:
    sces_yaml_f = temp(qc_dir  + '/sce_tmp_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml')
  threads: 1
  retries: config['resources']['retries']
  run:
    # split paths and sample names
    fs = [f"{qc_dir}/sce_cells_tmp_{s}_{FULL_TAG}_{DATE_STAMP}.rds" for s in SAMPLES]
    
    # check that all files exist
    for f in fs:
      if not os.path.isfile(f):
        raise FileNotFoundError(f"file {f} doesn't exist")

    # create a dictionary
    fs_dict = dict(zip(SAMPLES, fs))

    # write to yaml
    with open(output.sces_yaml_f, 'w') as f:
      yaml.dump(fs_dict, f, default_flow_style=True)

