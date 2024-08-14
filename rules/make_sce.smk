# snakemake rule for joining ambient outputs into sce

localrules: make_sce_input_df

# function for make_sce_input_df that takes different ambient methods as input
def make_sce_input_df(AMBIENT_METHOD, smpl_stats_f, samples_ls, ambient_outs_yamls):
    # Load sample stats (only needed for cellbender)
    if AMBIENT_METHOD == 'cellbender':
        sample_stats = pd.read_csv(smpl_stats_f, sep=',')
        ok_samples = sample_stats[~sample_stats['bad_sample']]['sample_id'].tolist()

    # Initialize an empty list to collect dataframes
    df_list = []

    for sample, ambient_outs_yaml in zip(samples_ls, ambient_outs_yamls):
        # Load ambient outs YAML file
        with open(ambient_outs_yaml) as f:
            amb_outs = yaml.load(f, Loader=yaml.FullLoader)

        # Determine the correct matrix file paths based on ambient method
        if AMBIENT_METHOD == 'cellbender':
            sce_df = pd.DataFrame({
                'sample_id': [sample],
                'cb_full': [amb_outs['cb_full_f']],
                'cb_filt': [amb_outs['cb_filt_f']]
            })
            # Restrict to ok samples
            sce_df = sce_df[sce_df['sample_id'].isin(ok_samples)]
            df_list.append(sce_df)
        elif AMBIENT_METHOD == 'decontx':
            sce_df = pd.DataFrame({
                'sample_id': [sample],
                'dcx_filt': [amb_outs['dcx_filt_f']]
            })
            df_list.append(sce_df)
        else:
            sce_df = pd.DataFrame({
                'sample_id': [sample],
                'bcs_filt': [amb_outs['cell_filt_f']]
            })
            df_list.append(sce_df)

    # Concatenate all collected dataframes into one
    final_sce_df = pd.concat(df_list, ignore_index=True)
    return final_sce_df





# input strings to make_sce_object were suuuuuper long, so we use a df instead
rule make_sce_input_df:
  input:
     smpl_stats_f    = amb_dir + '/ambient_sample_statistics_' + DATE_STAMP + '.txt',
     amb_yaml_fs = expand( amb_dir + '/ambient_{sample}/ambient_{sample}_' + DATE_STAMP + '_output_paths.yaml', sample = SAMPLES)
  output:
    sce_df      = sce_dir + '/sce_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  run:
    # anything else needed for AMBIENT METHOD??
    df =  make_sce_input_df(AMBIENT_METHOD, input.smpl_stats_f, SAMPLES, input.amb_yaml_fs)
    # save dataframe
    df.to_csv(output.sce_df, index = False)



if AMBIENT_METHOD == 'cellbender':
  rule make_sce_object:
    input:
      sce_df      = sce_dir + '/sce_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
    output:
      sce_all_f   = sce_dir + '/sce_cells_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
    threads: 16
    resources:
      mem_mb    = 8192
    conda:
      '../envs/rlibs.yml'
    shell:
    """
        # save sce object
        Rscript -e "source('scripts/make_sce.R'); \
          save_cellbender_as_sce( \
            sce_df      = '{input.sce_df}', \
            metadata_f  = '{METADATA_F}', \
            gtf_dt_f    = '{AF_GTF_DT_F}', \
            mito_str    = '{AF_MITO_STR}', \
            sce_f       = '{output.sce_all_f}', \
            bender_prob = {SCE_BENDER_PROB}, \
            n_cores     = {threads})"
    """
else:
  localrules: make_sce_object
  rule mke_sce_object:
    input:
      sce_df      = sce_dir + '/sce_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
    output:
      sce_all_f   = sce_dir + '/sce_cells_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
    threads: 16
    resources:
      mem_mb    = 8192
    conda:
      '../envs/rlibs.yml'
    shell:
    """
        Rscript -e "source('scripts/make_sce.R'); \
          save_noncb_as_sce ( \
            sce_df              = '{input.sce_df}', \
            ambient_method      = '{AMBIENT_METHOD}', \
            metadata_f          = '{METADATA_F}', \
            gtf_dt_f            = '{AF_GTF_DT_F}', \
            mito_str            = '{AF_MITO_STR}', \
            sce_f               = '{output.sce_all_f}', \
            min_counts          = {QC_HARD_MIN_COUNTS}, \
            n_cores             = {threads})"
    """

