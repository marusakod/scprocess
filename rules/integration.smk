# snakemake rule for integrating samples with harmony

rule run_integration:
  input:
    hvg_mat_f     = hvg_dir + '/top_hvgs_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
    dbl_hvg_mat_f = hvg_dir + '/top_hvgs_doublet_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
    sample_qc_f   = qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    coldata_f     = qc_dir  + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  output:
    integration_f = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  params:
    demux_type      = config['multiplexing']['demux_type'],
    exclude_mito    = config['qc']['exclude_mito'],
    int_embedding   = config['integration']['int_embedding'],
    int_theta       = config['integration']['int_theta'],
    int_batch_var   = config['integration']['int_batch_var'],
    int_n_dims      = config['integration']['int_n_dims'],
    int_dbl_res     = config['integration']['int_dbl_res'],
    int_dbl_cl_prop = config['integration']['int_dbl_cl_prop'],
    int_cl_method   = config['integration']['int_cl_method'],
    int_res_ls      = config['integration']['int_res_ls']
  threads: 1
  retries: config['resources']['retries'] 
  resources:
    mem_mb   = lambda wildcards, attempt: attempt * config['resources']['gb_run_integration'] * MB_PER_GB
  conda: 
    '../envs/integration.yaml'
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_integration/run_integration_' + DATE_STAMP + '.benchmark.txt'
  shell: """
    # check whether gpu available
    if [ -n "$CUDA_VISIBLE_DEVICES" ]; then
      use_gpu=1
    else
      use_gpu=0
    fi

    python3 scripts/integration.py \
      --hvg_mat_f     {input.hvg_mat_f} \
      --dbl_hvg_mat_f {input.dbl_hvg_mat_f} \
      --sample_qc_f   {input.sample_qc_f} \
      --coldata_f     {input.coldata_f} \
      --demux_type    {params.demux_type} \
      --exclude_mito  "{params.exclude_mito}" \
      --embedding     {params.int_embedding} \
      --n_dims        {params.int_n_dims} \
      --cl_method     {params.int_cl_method} \
      --dbl_res       {params.int_dbl_res} \
      --dbl_cl_prop   {params.int_dbl_cl_prop} \
      --theta         {params.int_theta} \
      --res_ls_concat "{params.int_res_ls}" \
      --integration_f {output.integration_f} \
      --batch_var     {params.int_batch_var} \
      --use_gpu       $use_gpu
    """

# if INT_USE_GPU:
#   rule run_gpu_integration:
#     input:
#       hvg_mat_f     = hvg_dir + '/top_hvgs_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
#       dbl_hvg_mat_f = hvg_dir + '/top_hvgs_doublet_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
#       sample_qc_f   = qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
#       coldata_f     = qc_dir  + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
#     output:
#       integration_f = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
#     params:
#       demux_type      = config['multiplexing']['demux_type'],
#       exclude_mito    = config['qc']['exclude_mito'],
#       int_embedding   = config['integration']['int_embedding'],
#       int_theta       = config['integration']['int_theta'],
#       int_batch_var   = config['integration']['int_batch_var'],
#       int_n_dims      = config['integration']['int_n_dims'],
#       int_dbl_res     = config['integration']['int_dbl_res'],
#       int_dbl_cl_prop = config['integration']['int_dbl_cl_prop'],
#       int_cl_method   = config['integration']['int_cl_method'],
#       int_res_ls      = config['integration']['int_res_ls']
#     threads: 1
#     retries: RETRIES 
#     resources:
#       mem_mb   = lambda wildcards, attempt: attempt * MB_RUN_INTEGRATION
#     conda: 
#       '../envs/integration.yaml'
#     benchmark:
#       benchmark_dir + '/' + SHORT_TAG + '_integration/run_gpu_integration_' + DATE_STAMP + '.benchmark.txt'
#     shell: """
#       python3 scripts/integration.py \
#         --hvg_mat_f       {input.hvg_mat_f} \
#         --dbl_hvg_mat_f   {input.dbl_hvg_mat_f} \
#         --sample_qc_f     {input.sample_qc_f} \
#         --coldata_f       {input.coldata_f} \
#         --demux_type      {params.demux_type} \
#         --exclude_mito    "{params.exclude_mito}" \
#         --embedding       {params.int_embedding} \
#         --n_dims          {params.int_n_dims} \
#         --cl_method       {params.int_cl_method} \
#         --dbl_res         {params.int_dbl_res} \
#         --dbl_cl_prop     {params.int_dbl_cl_prop} \
#         --theta           {params.int_theta} \
#         --res_ls_concat   "{params.int_res_ls}" \
#         --integration_f   {output.integration_f} \
#         --batch_var       {params.batch_var} \
#         --gpu
#       """
# else:
#   rule run_integration:
#     input:
#       hvg_mat_f     = hvg_dir + '/top_hvgs_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
#       dbl_hvg_mat_f = hvg_dir + '/top_hvgs_doublet_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
#       sample_qc_f   = qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
#       coldata_f     = qc_dir  + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
#     output:
#       integration_f = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
#     params:
#       demux_type      = config['multiplexing']['demux_type'],
#       exclude_mito    = config['qc']['exclude_mito'],
#       int_embedding   = config['integration']['int_embedding'],
#       int_theta       = config['integration']['int_theta'],
#       int_batch_var   = config['integration']['int_batch_var'],
#       int_n_dims      = config['integration']['int_n_dims'],
#       int_dbl_res     = config['integration']['int_dbl_res'],
#       int_dbl_cl_prop = config['integration']['int_dbl_cl_prop'],
#       int_cl_method   = config['integration']['int_cl_method'],
#       int_res_ls      = config['integration']['int_res_ls']
#     threads: 1
#     retries: RETRIES 
#     resources:
#       mem_mb   = lambda wildcards, attempt: attempt * MB_RUN_INTEGRATION
#     conda: 
#       '../envs/integration.yaml'
#     benchmark:
#       benchmark_dir + '/' + SHORT_TAG + '_integration/run_integration_' + DATE_STAMP + '.benchmark.txt'
#     shell: """
#       python3 scripts/integration.py \
#         --hvg_mat_f       {input.hvg_mat_f} \
#         --dbl_hvg_mat_f   {input.dbl_hvg_mat_f} \
#         --sample_qc_f     {input.sample_qc_f} \
#         --coldata_f       {input.coldata_f} \
#         --demux_type      {params.demux_type} \
#         --exclude_mito    "{params.exclude_mito}" \
#         --embedding       {params.int_embedding} \
#         --n_dims          {params.int_n_dims} \
#         --cl_method       {params.int_cl_method} \
#         --dbl_res         {params.int_dbl_res} \
#         --dbl_cl_prop     {params.int_dbl_cl_prop} \
#         --theta           {params.int_theta} \
#         --res_ls_concat   "{params.int_res_ls}" \
#         --integration_f   {output.integration_f} \
#         --batch_var       {params.batch_var}
#       """



# rule to create sce objects without any doublets (and delete temporary sce objects in the qc directory)
rule make_clean_sces: 
  input:
    sces_yaml_f   = qc_dir  + '/sce_tmp_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml', 
    integration_f = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  output:
    clean_sce_f   = int_dir + '/sce_cells_clean_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_make_clean_sces'] * MB_PER_GB
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_integration/make_clean_sces_{sample}_' + DATE_STAMP + '.benchmark.txt'
  conda:
    '../envs/rlibs.yaml'
  shell: """
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
   run:
    # split paths and sample names
    fs = [f"{int_dir}/sce_cells_clean_{s}_{FULL_TAG}_{DATE_STAMP}.rds" for s in SAMPLES]
    
    # check that all files exist
    for f in fs:
      if not os.path.isfile(f):
        raise FileNotFoundError(f"File {f} doesn't exist")

    # create a dictionary
    fs_dict = dict(zip(SAMPLES, fs))

    # write to yaml
    with open(output.sces_yaml_f, 'w') as f:
      yaml.dump(fs_dict, f, default_flow_style=False)

