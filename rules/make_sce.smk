# snakemake rule for joining ambient outputs into sce

localrules: make_sce_input_df


def make_sce_input_df(AMBIENT_METHOD, DEMUX_TYPE, SCPROCESS_DATA_DIR, SAMPLE_CHEMISTRY,
                     smpl_stats_f, samples_ls, ambient_outs_yamls, hto_h5_fs):

        # get the right run variable
    if DEMUX_TYPE is not None:
      run = 'pool_id'
    else:
      run = 'sample_id'

    # Load sample stats (only needed for cellbender)
    if AMBIENT_METHOD == 'cellbender':
        sample_stats = pd.read_csv(smpl_stats_f, delimiter='\t')
        ok_samples = sample_stats[~sample_stats['bad_sample']][run].tolist()

    # Initialize an empty list to collect dataframes
    df_list = []

    # Adjust zip input based on DEMUX_TYPE
    if DEMUX_TYPE == 'af':
        zipped_data = zip(samples_ls, ambient_outs_yamls, hto_h5_fs)
    else:
        zipped_data = zip(samples_ls, ambient_outs_yamls)

    for data in zipped_data:
        if DEMUX_TYPE == 'af':
            sample, ambient_outs_yaml, hto_h5_f = data
        else:
            sample, ambient_outs_yaml = data

        # Load ambient outs YAML file
        with open(ambient_outs_yaml) as f:
            amb_outs = yaml.load(f, Loader=yaml.FullLoader)

        # Determine the correct matrix file paths based on ambient method
        if AMBIENT_METHOD == 'cellbender':
            sce_df = pd.DataFrame({
                 run : [sample],
                'cb_full': [amb_outs['cb_full_f']],
                'cb_filt': [amb_outs['cb_filt_f']], 
                'bcs_csv': [amb_outs['cb_bcs_f']]
            })
            # Restrict to ok samples
            sce_df = sce_df[sce_df[run].isin(ok_samples)]
        elif AMBIENT_METHOD == 'decontx':
            sce_df = pd.DataFrame({
                run : [sample],
                'dcx_filt': [amb_outs['dcx_filt_f']],
                'bcs_csv' : [amb_outs['dcx_bcs_f']]
            })
        else:
            sce_df = pd.DataFrame({
                 run : [sample],
                'bcs_filt': [amb_outs['cell_filt_f']], 
                'bcs_csv' : [amb_outs['cell_bcs_f']]
            })

        # Add path to af hto file and barcode tranlation file if DEMUX_TYPE is "af"
        if DEMUX_TYPE == 'af':
            sce_df['hto_f'] = hto_h5_f

            wl_df_f = os.path.join(SCPROCESS_DATA_DIR, 'cellranger_ref/cellranger_whitelists.csv')
            wl_df = pd.read_csv(wl_df_f)
            wl_trans_f = wl_df.loc[wl_df['chemistry'] == SAMPLE_CHEMISTRY, 'translation_f'].values[0]
            wl_trans_f = os.path.join(SCPROCESS_DATA_DIR, 'cellranger_ref', wl_trans_f)
            sce_df['wl_trans_f'] = wl_trans_f
        

        df_list.append(sce_df)

    # Concatenate all collected dataframes into one
    final_sce_df = pd.concat(df_list, ignore_index=True)
    return final_sce_df


# input strings to make_sce_object were suuuuuper long, so we use a df instead
rule make_sce_input_df:
  input:
     smpl_stats_f    = amb_dir + '/ambient_sample_statistics_' + DATE_STAMP + '.txt',
     amb_yaml_fs = expand( amb_dir + '/ambient_{sample}/ambient_{sample}_' + DATE_STAMP + '_output_paths.yaml', sample = runs)
     hto_h5_fs   = expand(af_dir + '/af_{sample}/hto/af_hto_counts_mat.h5', sample = runs) if DEMUX_TYPE == "af" else []
  output:
    sce_df      = sce_dir + '/sce_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  run:
    # anything else needed for AMBIENT METHOD??
    df =  make_sce_input_df(AMBIENT_METHOD, DEMUX_TYPE, input.smpl_stats_f, runs, input.amb_yaml_fs)
    # save dataframe
    df.to_csv(output.sce_df, index = False)


if DEMUX_TYPE == 'af':
# save sce object with hto counts and demultiplex
  rule make_hto_sce_object: 
    input: 
      sce_df      = sce_dir + '/sce_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
    output:
      sce_hto_f   = sce_dir + '/sce_cells_htos_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
    threads: 4
    retries: RETRIES
      mem_mb = lambda wildcards, attempt: attempt * MB_MAKE_SCE_OBJECT
    conda:
     '../envs/rlibs.yml'
    shell: 
     """
       # save hto sce with demultiplexing info
        Rscript -e "source('scripts/make_sce.R'); \
          save_hto_sce( \
            sce_df_f    = '{input.sce_df}', \
            sce_f       = '{output.sce_hto_f}' \
            n_cores     = {threads})"
     """


if AMBIENT_METHOD == 'cellbender':
  rule make_sce_object:
    input:
      sce_df      = sce_dir + '/sce_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
      demux_f     = hto_sce_f if DEMUX_TYPE == 'af' else DEMUX_F
    output:
      sce_all_f   = sce_dir + '/sce_cells_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
    threads: 4
    retries: RETRIES
    resources:
      mem_mb    =  lambda wildcards, attempt: attempt * MB_MAKE_SCE_OBJECT
    conda:
      '../envs/rlibs.yml'
    params:
      demux_type = "" if DEMUX_TYPE is None else DEMUX_TYPE
    shell:
      """
        # save sce object
        Rscript -e "source('scripts/make_sce.R'); \
          save_cellbender_as_sce( \
            sce_df_f       = '{input.sce_df}', \
            metadata_f     = '{METADATA_F}', \
            gtf_dt_f       = '{AF_GTF_DT_F}', \
            mito_str       = '{AF_MITO_STR}', \
            sce_f          = '{output.sce_all_f}', \
            bender_prob    = {SCE_BENDER_PROB}, \
            n_cores        = {threads}, \
            demux_type     = '{params.demux_type}', \
            keep_smpls_str = '{SAMPLE_STR}')"
      """
else:
  rule make_sce_object:
    input:
      sce_df      = sce_dir + '/sce_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
      demux_f     = hto_sce_f if DEMUX_TYPE == 'af' else DEMUX_F
    output:
      sce_all_f   = sce_dir + '/sce_cells_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
    threads: 4
    retries: RETRIES 
    resources:
      mem_mb    = lambda wildcards, attempt: attempt * MB_MAKE_SCE_OBJECT
    conda:
      '../envs/rlibs.yml'
    params:
      demux_type = "" if DEMUX_TYPE is None else DEMUX_TYPE
    shell:
      """
        Rscript -e "source('scripts/make_sce.R'); \
          save_noncb_as_sce ( \
            sce_df_f            = '{input.sce_df}', \
            ambient_method      = '{AMBIENT_METHOD}', \
            metadata_f          = '{METADATA_F}', \
            gtf_dt_f            = '{AF_GTF_DT_F}', \
            mito_str            = '{AF_MITO_STR}', \
            sce_f               = '{output.sce_all_f}', \
            min_counts          = {QC_HARD_MIN_COUNTS}, \
            n_cores             = {threads}, \
            demux_type          = '{params.demux_type}', \
            keep_smpls_str      = '{SAMPLE_STR}')"
      """

