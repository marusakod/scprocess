# snakemake rule for running scDblFinder on sce object

localrules: make_dbl_files_df

rule run_scDblFinder:
  input:
    smpl_stats_f  = amb_dir + '/ambient_sample_statistics_' + DATE_STAMP + '.txt',
    sce_all_f     = sce_dir + '/sce_cells_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  output:
    dbl_f       = dbl_dir + '/dbl_{sample}/scDblFinder_{sample}_outputs_' + FULL_TAG +'_' + DATE_STAMP + '.txt.gz',
    dimred_f    = dbl_dir + '/dbl_{sample}/scDblFinder_{sample}_dimreds_' + FULL_TAG +'_' + DATE_STAMP + '.txt.gz'
  threads: 1
  retries: RETRIES 
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * MB_RUN_SCDBLFINDER
  conda: 
   '../envs/rlibs.yml'
  shell:
   """
    # run scDblFinder
    Rscript -e "source('scripts/doublet_id.R'); \
      main_doublet_id( \
        sel_sample      = '{wildcards.sample}', \
        sce_f           = '{input.sce_all_f}', \
        sample_stats_f  = '{input.smpl_stats_f}' \
        ambient_method  = '{AMBIENT_METHOD}', \
        dbl_f           = '{output.dbl_f}', \
        dimred_f        = '{output.dimred_f}', \
        min_feats       = {DBL_MIN_FEATS})"
   """


# input strings to make_sce_object were suuuuuper long, so we use a df instead
rule make_dbl_files_df:
  input:
    scdbl_ls    = expand(dbl_dir + '/dbl_{sample}/scDblFinder_{sample}_outputs_' + FULL_TAG +'_' + DATE_STAMP + '.txt.gz', sample = runs),
    dimred_ls   = expand(dbl_dir + '/dbl_{sample}/scDblFinder_{sample}_dimreds_' + FULL_TAG +'_' + DATE_STAMP + '.txt.gz', sample = runs)
  output:
    dbl_fs_f    = dbl_dir + '/doublet_id_files_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  run:
    # make pandas dataframe of cellbender outputs
    if DEMUX_TYPE is not None:
     sample_var = 'pool_id'
    else:
     sample_var = 'sample_id'

    df = pd.DataFrame({
      sample_var:   runs,
      'dbl_f':      input.scdbl_ls,
      'dimred_f':   input.dimred_ls
    })

    # save dataframe
    df.to_csv(output.dbl_fs_f, index = False)


rule combine_scDblFinder_outputs:
  input:
    dbl_fs_f        = dbl_dir + '/doublet_id_files_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    # doublet_id
    combn_dbl_f     = dbl_dir + '/scDblFinder_combined_outputs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    combn_dimred_f  = dbl_dir + '/scDblFinder_combined_dimreds_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  threads: 8
  retries: RETRIES 
  resources:
    mem_mb    = lambda wildcards, attempt: attempt * MB_COMBINE_SCDBLFINDER_OUTPUTS
  params:
      demux_type = "" if DEMUX_TYPE is None else DEMUX_TYPE
  conda: 
    '../envs/rlibs.yml'
  shell:
    """
    # combine scDblFinder outputs
    Rscript -e "source('scripts/doublet_id.R'); \
      combine_scDblFinder_outputs( \
        dbl_fs_f        = '{input.dbl_fs_f}',
        combn_dbl_f     = '{output.combn_dbl_f}', \
        combn_dimred_f  = '{output.combn_dimred_f}', \
        demux_type      = '{params.demux_type}', \
        n_cores = {threads})"
    """
