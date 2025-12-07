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
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('make_one_pb_empty', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run),
    runtime = lambda wildcards, input: get_resources('make_one_pb_empty', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_pb_empties/make_one_pb_empty_{run}_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
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
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('merge_pb_empty', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('merge_pb_empty', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_pb_empties/merge_pb_empty_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
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
    qc_stats_f  = f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv'
  output:
    pb_all_f    = pb_dir  + '/pb_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  params:
    batch_var   = BATCH_VAR
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('make_pb_all', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('make_pb_all', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_pb_empties/make_pb_all_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    make_pb_cells( \
      sce_fs_yaml = '{input.sces_yaml_f}',
      qc_stats_f  = '{input.qc_stats_f}',
      pb_f        = '{output.pb_all_f}',
      batch_var   = '{params.batch_var}',
      n_cores     =  {threads})"
    """


# calculate ambient genes across all samples or per group
rule calculate_ambient_genes:
  input:
    pb_empty_f  = pb_dir + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
    pb_all_f    = pb_dir + '/pb_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  output:
    empty_gs_f  = empty_dir + '/edger_empty_genes_' + FULL_TAG + '_all_' + DATE_STAMP + '.csv.gz' # one file per group
  params:
    fdr_thr     = config['pb_empties']['ambient_genes_fdr_thr'],
    logfc_thr   = config['pb_empties']['ambient_genes_logfc_thr']
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('calculate_ambient_genes', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('calculate_ambient_genes', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
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
