# rule to aggregate single cell data into pseudobulk matrices

import polars as pl


localrules: make_empty_pb_input_df

# for empty pseudobulks
rule make_empty_pb_input_df:
  input:
    af_mat_ls   = expand( [af_dir + '/af_{run}/' + af_rna_dir + 'af_counts_mat.h5'], run = RUNS), 
    af_knee_ls  = expand( [af_dir + '/af_{run}/' + af_rna_dir + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz'], run = RUNS),
    run_stats_f = amb_dir + '/ambient_run_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    af_paths_f  = pb_dir + '/af_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  run:
    # make dataframe with alevin outputs
    df          = pl.DataFrame({
      RUN_VAR:      RUNS,
      'af_mat_f':   input.af_mat_ls,
      'af_knee_f':  input.af_knee_ls
    })
    
    # add bad sample labels if cellbender
    if config['ambient']['ambient_method'] == 'cellbender':
      run_stats_df  = pl.read_csv(input.run_stats_f).select([RUN_VAR, 'bad_run'])
      df            = df.join(run_stats_df, on = RUN_VAR)

    # add output file paths
    pb_tmp_pat  = f"{pb_dir}/tmp_pb_empties_{"{}"}_{FULL_TAG}_{DATE_STAMP}.rds"
    df          = df.with_columns(
      pl.format(pb_tmp_pat, pl.col(RUN_VAR)).alias('pb_tmp_f')
    )

    # save dataframe
    df.write_csv(output.af_paths_f)


# make empties per sample then combine
rule make_one_pb_empty:
  input:
    af_h5_f       = af_dir  + '/af_{run}/' + af_rna_dir + 'af_counts_mat.h5', 
    af_paths_f    = pb_dir +  '/af_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    pb_empty_f    = temp(pb_dir + '/tmp_pb_empties_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.rds')
  params:
    ambient_method  = config['ambient']['ambient_method']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt, input: (
      attempt * (
      config['resources'].get('gb_make_one_pb_empty', None) * MB_PER_GB
      if config['resources'].get('gb_make_one_pb_empty') is not None
      else (446.932 + 26.863 * (os.path.getsize(input.af_h5_f)//MB_PER_GB**2)) + 2*MB_PER_GB # lm + buffer
      )
    ), 
    runtime = lambda wildcards, input:
      config['resources'].get('min_make_one_pb_empty', None)
      if config['resources'].get('min_make_one_pb_empty') is not None
      else (26.849 + 0.427*(os.path.getsize(input.af_h5_f)//MB_PER_GB**2))/60 + 5 # lm + 5 min buffer
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_pb_empties/make_one_pb_empty_{run}_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/ambient.R'); source('scripts/pseudobulk_and_empties.R'); \
    make_pb_empty(
      sel_run         = '{wildcards.run}', 
      af_paths_f      = '{input.af_paths_f}', 
      pb_empty_f      = '{output.pb_empty_f}', 
      ambient_method  = '{params.ambient_method}',
      run_var         = '{RUN_VAR}')"
    """

rule merge_pb_empty:
  input:
    pb_empty_fs   = expand(pb_dir + '/tmp_pb_empties_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.rds', run = RUNS), 
    rowdata_f     = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', 
    af_paths_f    = pb_dir +  '/af_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    pb_empty_f    = pb_dir + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  params:
    ambient_method  = config['ambient']['ambient_method']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * config['resources']['gb_merge_pb_empty'] * MB_PER_GB, 
    runtime     = config['resources']['min_merge_pb_empty']
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_pb_empties/merge_pb_empty_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/ambient.R'); source('scripts/pseudobulk_and_empties.R'); \
    merge_empty_pbs( \
      af_paths_f      = '{input.af_paths_f}', 
      rowdata_f       = '{input.rowdata_f}',
      empty_pbs_f     = '{output.pb_empty_f}', 
      ambient_method  = '{params.ambient_method}'
    )"
    """


# make pseudobulk with all cells that passed qc
rule make_pb_all:
  input:
    sces_yaml_f = qc_dir  + '/sce_tmp_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml', 
    pb_empty_f  = pb_dir  + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
    qc_stats_f  = qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    pb_all_f    = pb_dir  + '/pb_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt: (
      attempt * (
      config['resources'].get('gb_make_pb_all', None) * MB_PER_GB
      if config['resources'].get('gb_make_pb_all') is not None
      else (5379.419 + 213.689*len(SAMPLES)) + 2*MB_PER_GB # lm + buffer
      )
    ), 
    runtime = lambda wildcards: 
      config['resources'].get('min_make_pb_all', None)
      if config['resources'].get('min_make_pb_all') is not None
      else (19.151 + 4.049*len(SAMPLES))/60 + 5 # lm + 5 min buffer
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_pb_empties/make_pb_all_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/ambient.R'); source('scripts/pseudobulk_and_empties.R'); \
    make_pb_cells( \
      sce_fs_yaml = '{input.sces_yaml_f}',
      qc_stats_f  = '{input.qc_stats_f}',
      pb_f        = '{output.pb_all_f}',
      n_cores     =  {threads})"
    """


# calculate ambient genes across all samples or per group
rule calculate_ambient_genes:
  input:
    pb_empty_f  = pb_dir + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
    pb_all_f    = pb_dir + '/pb_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  output:
    empty_gs_f  = empty_dir + '/edger_empty_genes_' + FULL_TAG + '_all_' + DATE_STAMP + '.txt.gz' # one file per group
  params:
    fdr_thr     = config['pb_empties']['ambient_genes_fdr_thr'],
    logfc_thr   = config['pb_empties']['ambient_genes_logfc_thr']
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt: (
      attempt * (
      config['resources'].get('gb_calculate_ambient_genes', None) * MB_PER_GB
      if config['resources'].get('gb_calculate_ambient_genes') is not None
      else (900.397 + 4.057*len(SAMPLES)) + 2*MB_PER_GB # lm + buffer
      )
    ), 
    runtime = lambda wildcards: 
      config['resources'].get('min_calculate_ambient_genes', None)
      if config['resources'].get('min_calculate_ambient_genes') is not None
      else (19.290 + 0.384*len(SAMPLES))/60 + 5 # lm + 5 min buffer
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_pb_empties/calculate_ambient_genes_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    calc_empty_genes(
      pb_cells_f = '{input.pb_all_f}',
      pb_empty_f = '{input.pb_empty_f}',
      fdr_thr    =  {params.fdr_thr},
      logfc_thr  =  {params.logfc_thr},
      empty_gs_f = '{output.empty_gs_f}')"
    """
