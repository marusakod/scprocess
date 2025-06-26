# snakemake rule for integrating samples with harmony

rule run_integration:
  input:
    hvg_mat_f     = hvg_dir + '/top_hvgs_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
    dbl_hvg_mat_f = hvg_dir + '/top_hvgs_doublet_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
    sample_qc_f   = qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    coldata_f     = qc_dir  + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    integration_f = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  threads: 8
  retries: RETRIES 
  resources:
    mem_mb   = lambda wildcards, attempt: attempt * MB_RUN_INTEGRATION
  conda: 
    '../envs/rlibs.yml'
  shell:
    """
    # run harmony
    Rscript -e "source('scripts/integration.R'); source('scripts/ambient.R');
      run_integration( 
        hvg_mat_f        = '{input.hvg_mat_f}', 
        dbl_hvg_mat_f    = '{input.dbl_hvg_mat_f}', 
        sample_qc_f      = '{input.sample_qc_f}', 
        coldata_f        = '{input.coldata_f}', 
        demux_type       = '{DEMUX_TYPE}', 
        exclude_mito     = '{EXCLUDE_MITO}', 
        reduction        = '{INT_REDUCTION}',
        n_dims           = {INT_N_DIMS}, 
        cl_method        = '{INT_CL_METHOD}', 
        dbl_res          = {INT_DBL_RES}, 
        dbl_cl_prop      = {INT_DBL_CL_PROP}, 
        theta            = {INT_THETA}, 
        res_ls_concat    = '{INT_RES_LS}', 
        integration_f    = '{output.integration_f}', 
        batch_var        = '{BATCH_VAR}', 
        n_cores          = {threads})"
    """
