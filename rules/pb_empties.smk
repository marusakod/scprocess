# rule to aggregate single cell data into pseudobulk matrices


localrules: make_pb_input_df

rule make_pb_input_df: # for empty pseudobulks
  input:
    af_mat_ls   = expand( [af_dir + '/af_{run}/' + af_rna_dir + 'af_counts_mat.h5'], run = runs), 
    af_knee_ls  = expand( [af_dir + '/af_{run}/' + af_rna_dir + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz'], run = runs), 
    amb_stats_f = amb_dir + '/ambient_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    af_paths_f  = pb_dir + '/af_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  run:

    # make dataframe with alevin outputs
    df          = pd.DataFrame({
      SAMPLE_VAR:   runs,
      'af_mat_f':   input.af_mat_ls,
      'af_knee_f':  input.af_knee_ls
    })
    
    # add bad sample labels if cellbender
    if AMBIENT_METHOD == 'cellbender':
      amb_stats_df = pd.read_csv(input.amb_stats_f)
      bad_samples  = amb_stats_df.loc[ amb_stats_df['bad_sample'] == True, SAMPLE_VAR].tolist()
      df['bad_sample'] = df[SAMPLE_VAR].apply(lambda run: run in bad_samples)

    # add output file paths
    df['pb_tmp_f'] = df[SAMPLE_VAR].apply(lambda r: f"{pb_dir}/pb_empties_{r}_{FULL_TAG}_{DATE_STAMP}.rds")
    
    # save dataframe
    df.to_csv(output.af_paths_f, index = False)


# make empties per sample then combine
rule make_one_pb_empty:
  input:
    af_paths_f    = pb_dir +  '/af_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    pb_empty_f    = temp(pb_dir + '/pb_empties_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.rds')
  threads: 1
  retries: RETRIES 
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * MB_PB_MAKE_PBS
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_pb_empties/make_pb_empty_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell:
    """
    Rscript -e "source('scripts/utils.R'); source('scripts/ambient.R'); source('scripts/pseudobulk_and_empties.R'); \
    make_pb_empty( \
      sel_s       = '{wildcards.run}', 
      af_paths_f  = '{input.af_paths_f}', 
      pb_empty_f  = '{output.pb_empty_f}', 
      ambient_method = '{AMBIENT_METHOD}',
      sample_var  = '{SAMPLE_VAR}')"
    """

rule merge_pb_empty:
  input:
    pb_empty_fs   = expand(pb_dir + '/pb_empties_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.rds', run = runs), 
    rowdata_f     = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    af_paths_f    = pb_dir +  '/af_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    pb_empty_f    = pb_dir + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  threads: 1
  retries: RETRIES 
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * MB_PB_MAKE_PBS
  conda: 
    '../envs/rlibs.yaml'
  shell:
    """
    Rscript -e "source('scripts/utils.R'); source('scripts/ambient.R'); source('scripts/pseudobulk_and_empties.R'); \
    merge_empty_pbs( \
      af_paths_f  = '{input.af_paths_f}', 
      rowdata_f   = '{input.rowdata_f}',
      empty_pbs_f = '{output.pb_empty_f}', 
      ambient_method = '{AMBIENT_METHOD}')"
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
  retries: RETRIES 
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * MB_PB_MAKE_PBS
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_pb_empties/make_pb_all_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell:
    """
    Rscript -e "source('scripts/utils.R'); source('scripts/ambient.R'); source('scripts/pseudobulk_and_empties.R'); \
    make_pb_cells( \
      sce_fs_yaml = '{input.sces_yaml_f}',
      qc_stats_f  = '{input.qc_stats_f}',
      pb_f        = '{output.pb_all_f}',
      n_cores     = {threads})"
    """


# calculate ambient genes across all samples or per group
rule calculate_ambient_genes:
  input:
    pb_empty_f  = pb_dir + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
    pb_all_f    = pb_dir + '/pb_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  output:
    empty_gs_f  = empty_dir + '/edger_empty_genes_' + FULL_TAG + '_all_' + DATE_STAMP + '.txt.gz' # one file per group
  threads: 4
  retries: RETRIES 
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * MB_PB_MAKE_PBS
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_pb_empties/calculate_ambient_genes_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell:
    """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    calc_empty_genes(
      pb_cells_f = '{input.pb_all_f}',
      pb_empty_f = '{input.pb_empty_f}',
      fdr_thr    = {AMBIENT_GENES_FDR_THR},
      logfc_thr  = {AMBIENT_GENES_LOGFC_THR},
      empty_gs_f = '{output.empty_gs_f}')"
    """
