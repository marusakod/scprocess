# rule to aggregate single cell data into pseudobulk matrices


localrules: make_pb_input_df

rule make_pb_input_df: # for empty pseudobulks
  input:
    af_mat_ls   = expand( [af_dir + '/af_{run}/' + af_rna_dir + 'af_counts_mat.h5'], run = runs), 
    af_knee_ls  = expand( [af_dir + '/af_{run}/' + af_rna_dir + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz'], run = runs) 
  output:
    af_paths_f  = pb_dir + '/af_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  run:

    # make pandas dataframe of cellbender outputs
    df          = pd.DataFrame({
      SAMPLE_VAR:   runs,
      'af_mat_f':   input.af_mat_ls,
      'af_knee_f':  input.af_knee_ls
    })
    
    # save dataframe
    df.to_csv(output.af_paths_f, index = False)



rule make_pb_empty:
  input:
    amb_stats_f   = amb_dir + '/ambient_sample_statistics_' + DATE_STAMP + '.txt', 
    af_paths_f    = pb_dir +  '/af_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    rowdata_f     = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    pb_empty_f    = pb_dir + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  threads: 8
  retries: RETRIES 
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * MB_PB_MAKE_PBS
  conda: 
    '../envs/rlibs.yml'
  shell:
    """
    Rscript -e "source('scripts/utils.R'); source('scripts/ambient.R'); source('scripts/pseudobulk_and_empties.R'); \
    make_pb_empty( \
      af_paths_f  = '{input.af_paths_f}', \
      rowdata_f   = '{input.rowdata_f}', \
      amb_stats_f = '{input.amb_stats_f}', \
      pb_empty_f  = '{output.pb_empty_f}', \
      ambient_method = '{AMBIENT_METHOD}', \
      sample_var  = '{SAMPLE_VAR}', \
      n_cores     =  {threads})"
    """


# make pseudobulk with all cells that passed qc
rule make_pb_all:
  input:
    sces_yaml_f = qc_dir  + '/sce_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml', 
    pb_empty_f  = pb_dir  + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
    qc_stats_f  = qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.txt'
  output:
    pb_all_f    = pb_dir  + '/pb_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  threads: 4
  retries: RETRIES 
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * MB_PB_MAKE_PBS
  conda: 
    '../envs/rlibs.yml'
  shell:
    """
    Rscript -e "source('scripts/utils.R'); source('scripts/ambient.R'); source('scripts/pseudobulk_and_empties.R'); \
    make_pb_cells( \
      sce_fs_yaml = '{input.sces_yaml_f}', \
      qc_stats_f  = '{input.qc_stats_f}', \
      pb_f        = '{output.pb_all_f}', \
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
  conda: 
    '../envs/rlibs.yml'
  shell:
    """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    calc_empty_genes(
      pb_cells_f = '{input.pb_all_f}', \
      pb_empty_f = '{input.pb_empty_f}', \
      fdr_thr    = {AMBIENT_GENES_FDR_THR}, \
      logfc_thr  = {AMBIENT_GENES_LOGFC_THR}, \
      empty_gs_f = '{output.empty_gs_f}')"
    """
