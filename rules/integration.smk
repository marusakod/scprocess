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
    '../envs/rlibs.yaml'
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_integration/run_integration_' + DATE_STAMP + '.benchmark.txt'
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


# rule to create sce objects without any doublets (and delete temporary sce objects in the qc directory)
rule make_clean_sces: 
  input:
    sces_yaml_f   = qc_dir  + '/sce_tmp_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml', 
    integration_f = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    clean_sce_f   = int_dir + '/sce_cells_clean_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  threads: 1
  retries: RETRIES
  resources:
    mem_mb = lambda wildcards, attempt: attempt * MB_MAKE_CLEAN_SCES
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_integration/make_clean_sces_{sample}_' + DATE_STAMP + '.benchmark.txt'
  conda:
    '../envs/rlibs.yaml'
  shell:
    """

    Rscript -e "source('scripts/integration.R');
      make_clean_sces(
        sel_s         = '{wildcards.sample}', 
        integration_f = '{input.integration_f}', 
        sces_yaml_f   = '{input.sces_yaml_f}', 
        clean_sce_f   = '{output.clean_sce_f}')"
    
    """

# make a yaml with all clean sce file paths
rule make_clean_sce_paths_yaml:
   input:
    clean_sce_f = expand(int_dir + '/sce_cells_clean_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.rds', sample = SAMPLES) # not used
   output:
    sces_yaml_f = int_dir + '/sce_clean_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml'
   benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_integration/make_clean_sce_paths_yaml_' + DATE_STAMP + '.benchmark.txt'
   run:
    # split paths and sample names
    fs = [f"{int_dir}/sce_cells_clean_{s}_{FULL_TAG}_{DATE_STAMP}.rds" for s in SAMPLES]
    
    # check that all files exist
    for f in fs:
     assert os.path.isfile(f), \
      f"File {f} doesn't exist"

    # create a dictionary
    fs_dict = dict(zip(SAMPLES, fs))

    # write to yaml
    with open(output.sces_yaml_f, 'w') as f:
     yaml.dump(fs_dict, f, default_flow_style=False)
