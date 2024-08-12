# snakemake rule for joining cellbender outputs into sce

localrules: make_sce_input_df

# input strings to make_sce_object were suuuuuper long, so we use a df instead
rule make_sce_input_df:
  input:
    cb_bad_f    = cb_dir + '/bender_bad_samples_' + DATE_STAMP + '.txt',
    cb_full_ls  = expand(cb_dir + '/bender_{sample}/bender_{sample}_' + DATE_STAMP + '.h5', sample = SAMPLES),
    cb_filt_ls  = expand(cb_dir + '/bender_{sample}/bender_{sample}_' + DATE_STAMP + '_filtered.h5', sample = SAMPLES),
    cb_qc_f     = expand(cb_dir + '/bender_{sample}/bender_qc_metrics_{sample}_' + DATE_STAMP + '.txt.gz', sample = SAMPLES) ## added here to have it running before the sce and to clean the dependencies chain
  output:
    sce_df      = sce_dir + '/sce_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  run:
    # get list of "bad" samples, i.e. too many barcodes called as cells by cellbender
    cb_bad_df   = pd.read_csv(input.cb_bad_f, sep = '\t')
    ok_samples  = cb_bad_df[ ~cb_bad_df['bad_sample'] ].sample_id.tolist()

    # make pandas dataframe of cellbender outputs
    df          = pd.DataFrame({
      'sample_id':  SAMPLES,
      'cb_full':    input.cb_full_ls,
      'cb_filt':    input.cb_filt_ls
    })
    # restrict to ok samples
    df          = df[ df['sample_id'].isin(ok_samples) ]

    # save dataframe
    df.to_csv(output.sce_df, index = False)


if DO_CELLBENDER:
  rule make_sce_object:
    input:
      sce_df      = sce_dir + '/sce_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
    output:
      sce_all_f   = sce_dir + '/sce_' + ('bender' if DO_CELLBENDER else 'alevin') + '_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
    threads: 16
      retries: 5
    resources:
      mem_mb    = 8192
    conda:
      '../envs/rlibs.yml'
    shell:
      """
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
  rule make_sce_object:
    input:
      sce_df      = sce_dir + '/sce_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
    output:
      sce_all_f   = sce_dir + '/sce_' + ('bender' if DO_CELLBENDER else 'alevin') + '_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
    threads: 16
      retries: 5
    resources:
      mem_mb    = 8192
    conda:
      '../envs/rlibs.yml'
    shell:
      """
        Rscript -e "source('scripts/make_sce.R'); \
          save_alevin_h5_as_sce( \
            sce_df      = '{input.sce_df}', \
            af_dir      = '{af_dir}', \
            metadata_f  = '{METADATA_F}', \
            gtf_dt_f    = '{AF_GTF_DT_F}', \
            mito_str    = '{AF_MITO_STR}', \
            sce_f       = '{output.sce_all_f}', \
            min_counts  = {QC_HARD_MIN_COUNTS}, \
            n_cores     = {threads})"
      """

