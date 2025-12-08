# rule to aggregate single cell data into pseudobulk matrices

import polars as pl


localrules: make_runs_to_batches_df, make_tmp_pb_cells_df, make_empty_pb_input_df

# for empty pseudobulks
rule make_empty_pb_input_df:
  input:
    af_mat_ls   = expand( [f'{af_dir}/af_{{run}}/{af_rna_dir}af_counts_mat.h5'], run = RUNS), 
    af_knee_ls  = expand( [f'{af_dir}/af_{{run}}/{af_rna_dir}knee_plot_data_{{run}}_{DATE_STAMP}.txt.gz'], run = RUNS),
    run_stats_f = f'{amb_dir}/ambient_run_statistics_{FULL_TAG}_{DATE_STAMP}.csv'
  output:
    af_paths_f  = f'{pb_dir}/af_paths_{FULL_TAG}_{DATE_STAMP}.csv'
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
    af_h5_f       = f'{af_dir}/af_{{run}}/{af_rna_dir}af_counts_mat.h5', 
    af_paths_f    = f'{pb_dir}/af_paths_{FULL_TAG}_{DATE_STAMP}.csv'
  output:
    pb_empty_f    = temp(f'{pb_dir}/tmp_pb_empties_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds')
  params:
    ambient_method  = config['ambient']['ambient_method']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('make_one_pb_empty', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run),
    runtime = lambda wildcards, input: get_resources('make_one_pb_empty', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_pb_empties/make_one_pb_empty_{{run}}_{DATE_STAMP}.benchmark.txt'
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
    pb_empty_fs   = expand(f'{pb_dir}/tmp_pb_empties_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds', run = RUNS),
    af_paths_f    = f'{pb_dir}/af_paths_{FULL_TAG}_{DATE_STAMP}.csv',
    rowdata_f     = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    pb_empty_f    = f'{pb_dir}/pb_empties_{FULL_TAG}_{DATE_STAMP}.rds'
  params:
    ambient_method  = config['ambient']['ambient_method']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('merge_pb_empty', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('merge_pb_empty', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_pb_empties/merge_pb_empty_{DATE_STAMP}.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    merge_pbs_empty( \
      af_paths_f      = '{input.af_paths_f}', 
      rowdata_f       = '{input.rowdata_f}',
      pb_empty_f      = '{output.pb_empty_f}', 
      ambient_method  = '{params.ambient_method}'
    )"
    """


rule make_runs_to_batches_df:
  output:
    batch_lu_f  = f'{pb_dir}/runs_to_batches_{FULL_TAG}_{DATE_STAMP}.csv'
  params:
    batches     = RUNS_TO_BATCHES
  run:
    # make df
    lu_ls     = []
    for r, bs in RUNS_TO_BATCHES.items():
      lu_tmp    = pl.DataFrame({
        "run_var":    r,
        "batch_var":  bs
      })
      lu_ls.append(lu_tmp)
    lu_df     = pl.concat(lu_ls)

    # save
    lu_df.write_csv(output.batch_lu_f)


rule make_one_pb_cells:
  input:
    batch_lu_f  = f'{pb_dir}/runs_to_batches_{FULL_TAG}_{DATE_STAMP}.csv',
    af_h5_f     = f'{af_dir}/af_{{ run }}/{af_rna_dir}af_counts_mat.h5',
    qc_stats_f  = f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    h5_paths_f  = f'{amb_dir}/paths_h5_filtered_{FULL_TAG}_{DATE_STAMP}.csv',
    qc_all_f    = f'{qc_dir}/qc_all_samples_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    pb_cells_f  = temp(f'{pb_dir}/tmp_pb_cells_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds')
  params:
    run_var     = RUN_VAR,
    batch_var   = BATCH_VAR
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('make_one_pb_cells', 'memory', 
      lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run),
    runtime = lambda wildcards, input: get_resources('make_one_pb_cells', 'time', 
      lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_pb_empties/make_one_pb_cells_{{run}}_{DATE_STAMP}.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    make_pb_cells(
      sel_run     = '{wildcards.run}',
      batch_lu_f  = '{input.batch_lu_f}',
      qc_stats_f  = '{input.qc_stats_f}',
      h5_paths_f  = '{input.h5_paths_f}', 
      qc_f        = '{input.qc_all_f}',
      run_var     = '{params.run_var}',
      batch_var   = '{params.batch_var}',
      pb_cells_f  = '{output.pb_cells_f}'
    )"
    """


rule make_tmp_pb_cells_df:
  input:
    pb_cells_fs   = expand(f'{pb_dir}/tmp_pb_cells_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds', run = RUNS)
  output:
    cells_paths_f = temp(f'{pb_dir}/tmp_pb_cells_paths_{FULL_TAG}_{DATE_STAMP}.csv')
  params:
    run_var     = RUN_VAR,
    runs        = RUNS
  run:
    # make df
    paths_df    = pl.DataFrame({
      params.run_var: params.runs,
      "pb_path":      input.pb_cells_fs
    })
    paths_df    = paths_df.filter( pl.col("pb_path").map_elements(os.path.getsize, return_dtype=pl.Int64) > 0 )

    # save
    paths_df.write_csv(output.cells_paths_f)


rule merge_pb_cells:
  input:
    cells_paths_f = f'{pb_dir}/tmp_pb_cells_paths_{FULL_TAG}_{DATE_STAMP}.csv',
    pb_cells_fs   = expand(f'{pb_dir}/tmp_pb_cells_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds', run = RUNS), 
    rowdata_f     = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    pb_cells_f    = f'{pb_dir}/pb_cells_all_{FULL_TAG}_{DATE_STAMP}.rds'
  params:
    batch_var     = BATCH_VAR
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('merge_pb_empty', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('merge_pb_empty', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_pb_empties/merge_pb_cells_{DATE_STAMP}.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    merge_pbs_cells( \
      cells_paths_f = '{input.cells_paths_f}', 
      rowdata_f     = '{input.rowdata_f}',
      batch_var     = '{params.batch_var}',
      pb_cells_f    = '{output.pb_cells_f}'
    )"
    """


# calculate ambient genes across all samples or per group
rule calculate_ambient_genes:
  input:
    pb_empty_f  = f'{pb_dir}/pb_empties_{FULL_TAG}_{DATE_STAMP}.rds', 
    pb_cells_f  = f'{pb_dir}/pb_cells_all_{FULL_TAG}_{DATE_STAMP}.rds'
  output:
    empty_gs_f  = f'{empty_dir}/edger_empty_genes_all_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  params:
    fdr_thr     = config['pb_empties']['ambient_genes_fdr_thr'],
    logfc_thr   = config['pb_empties']['ambient_genes_logfc_thr']
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('calculate_ambient_genes', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('calculate_ambient_genes', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_pb_empties/calculate_ambient_genes_{DATE_STAMP}.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    calc_empty_genes(
      pb_cells_f = '{input.pb_cells_f}',
      pb_empty_f = '{input.pb_empty_f}',
      fdr_thr    =  {params.fdr_thr},
      logfc_thr  =  {params.logfc_thr},
      empty_gs_f = '{output.empty_gs_f}')"
    """
