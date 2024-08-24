# rule to aggregate single cell data into pseudobulk matrices

# make_pb_input_df
localrules: make_pb_input_df
rule make_pb_input_df:
  input:
    samples_f   = sce_dir + '/sce_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    af_mat_ls   = expand( [af_dir + '/af_{sample}/af_counts_mat.h5'], sample = SAMPLES),
    af_knee_ls  = expand( [af_dir + '/af_{sample}/knee_plot_data_{sample}_' + DATE_STAMP + '.txt.gz'], sample = SAMPLES)
  output:
    af_paths_f  = pb_dir + '/af_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  run:
    # get list of ok samples
    samples_df  = pd.read_csv(input.samples_f)
    ok_samples  = samples_df.sample_id.tolist()

    # make pandas dataframe of cellbender outputs
    df          = pd.DataFrame({
      'sample_id':  SAMPLES,
      'af_mat_f':   input.af_mat_ls,
      'af_knee_f':  input.af_knee_ls
    })
    
    # restrict to ok samples
    df          = df[ df['sample_id'].isin(ok_samples) ]

    # save dataframe
    df.to_csv(output.af_paths_f, index = False)


# make_pb_subset
rule make_pb_subset:
  params:
    subset      = '{subset}'
  input:
    sce_sub_f   = lbl_dir + '/sce_subset_' + FULL_TAG + '_{subset}_' + DATE_STAMP + '.rds',
    af_paths_f  = pb_dir + '/af_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    pb_subset_f = pb_dir + '/pb_subset_' + FULL_TAG + '_{subset}_' + DATE_STAMP + '.rds'
  threads: 4
  resources:
    mem_mb      = 8192
  conda: 
    '../envs/rlibs.yml'
  shell:
    """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    make_pb_subset(sce_f = '{input.sce_sub_f}', af_paths_f = '{input.af_paths_f}', \
      pb_f = '{output.pb_subset_f}', n_cores = {threads})"
    """


# make_pb_empty
rule make_pb_empty:
  input:
    sce_all_f     = int_dir + '/sce_clean_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    sce_ls        = expand( [lbl_dir +'/sce_subset_' + FULL_TAG + '_{subset}_' + DATE_STAMP + '.rds'], 
      subset = None if PB_SUBSETS is None else [*PB_SUBSETS] ),
    af_paths_f    = pb_dir + '/af_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    pb_empty_f    = pb_dir + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    empty_locs_f  = pb_dir + '/empty_plateau_locations_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  threads: 4
  resources:
    mem_mb      = 8192
  conda: 
    '../envs/rlibs.yml'
  shell:
    """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    make_pb_empty( \
      af_paths_f = '{input.af_paths_f}', gtf_dt_f = '{AF_GTF_DT_F}', \
      custom_empties_f = '{PB_CUSTOM_EMPTIES_F}', empty_locs_f = '{output.empty_locs_f}', \
      pb_empty_f = '{output.pb_empty_f}', n_cores = {threads})"
    """

# calc_empty_genes
rule calc_empty_genes:
  params:
    subset      = '{subset}'
  input:
    pb_subset_f = pb_dir + '/pb_subset_' + FULL_TAG + '_{subset}_' + DATE_STAMP + '.rds',
    pb_empty_f  = pb_dir + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  output:
    empty_gs_f  = empty_dir + '/edger_empty_genes_' + FULL_TAG + '_{subset}_' + DATE_STAMP + '.txt.gz'
  resources:
    mem_mb      = 8192
  shell:
    """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    calc_empty_genes(pb_cells_f = '{input.pb_subset_f}', \
      pb_empty_f = '{input.pb_empty_f}', empty_gs_f = '{output.empty_gs_f}')"
    """


# make_pb_all
rule make_pb_all:
  input:
    sce_all_f   = int_dir + '/sce_clean_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    af_paths_f  = pb_dir + '/af_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    pb_empty_f  = pb_dir + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  output:
    pb_all_f    = pb_dir + '/pb_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    empty_gs_f  = empty_dir + '/edger_empty_genes_' + FULL_TAG + '_all_' + DATE_STAMP + '.txt.gz'
  threads: 4
  resources:
    mem_mb      = 8192
  conda: 
    '../envs/rlibs.yml'
  shell:
    """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    make_pb_subset(sce_f = '{input.sce_all_f}', af_paths_f = '{input.af_paths_f}', \
      pb_f = '{output.pb_all_f}', n_cores = {threads})"
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    calc_empty_genes(pb_cells_f = '{output.pb_all_f}', \
      pb_empty_f = '{input.pb_empty_f}', empty_gs_f = '{output.empty_gs_f}')"
    """
