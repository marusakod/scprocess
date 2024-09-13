# snakemake rule for integrating samples with harmony

rule run_harmony:
  input:
    sce_all_f   = sce_dir + '/sce_cells_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    keep_f      = qc_dir  + '/keep_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    dbl_f       = dbl_dir + '/scDblFinder_combined_outputs_' + FULL_TAG +'_' + DATE_STAMP + '.txt.gz'
  output:
    harmony_f   = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    hvgs_f      = int_dir + '/harmony_hvgs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    sce_clean_f = int_dir + '/sce_clean_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  threads: 8
  retries: RETRIES 
  resources:
    mem_mb   = lambda wildcards, attempt: attempt * MB_RUN_HARMONY
  conda: 
    '../envs/rlibs.yml'
  shell:
    """
    # run harmony
    Rscript -e "source('scripts/integration.R'); \
      run_harmony(sce_all_f = '{input.sce_all_f}', \
        keep_f = '{input.keep_f}', dbl_f = '{input.dbl_f}', \
        exc_regex = '{INT_EXC_REGEX}', n_hvgs = {INT_N_HVGS}, n_dims = {INT_N_DIMS}, \
        dbl_res = {INT_DBL_RES}, dbl_cl_prop = {INT_DBL_CL_PROP}, \
        theta = {INT_THETA}, res_ls_concat = '{INT_RES_LS}', \
        harmony_f = '{output.harmony_f}', hvgs_f = '{output.hvgs_f}', \
        sce_clean_f = '{output.sce_clean_f}', n_cores = {threads})"
    """
