# snakemake rule for integrating samples with harmony

def get_run_for_one_batch(batch, RUNS_TO_BATCHES):
  sel_run   = [ run for run, run_batches in RUNS_TO_BATCHES.items() if batch in run_batches ]
  return sel_run

localrules: make_clean_sce_paths_yaml

rule run_integration:
  input:
    hvg_mat_f     = f'{hvg_dir}/top_hvgs_counts_{FULL_TAG}_{DATE_STAMP}.h5', 
    dbl_hvg_mat_f = f'{hvg_dir}/top_hvgs_doublet_counts_{FULL_TAG}_{DATE_STAMP}.h5', 
    qc_stats_f    = f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    coldata_f     = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    integration_f = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  params:
    demux_type      = config['multiplexing']['demux_type'],
    exclude_mito    = config['qc']['exclude_mito'],
    int_embedding   = config['integration']['int_embedding'],
    int_cl_method   = config['integration']['int_cl_method'],
    int_theta       = config['integration']['int_theta'],
    int_batch_var   = config['integration']['int_batch_var'],
    int_n_dims      = config['integration']['int_n_dims'],
    int_dbl_res     = config['integration']['int_dbl_res'],
    int_dbl_cl_prop = config['integration']['int_dbl_cl_prop'],
    int_res_ls      = config['integration']['int_res_ls']
  threads: 1
  retries: config['resources']['retries'] 
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_integration', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_integration', 'time', attempt)
  conda: 
    '../envs/integration.yaml'
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_integration/run_integration_{DATE_STAMP}.benchmark.txt'
  shell: """
    set +u
    # if GPU is available, use it
    USE_GPU_FLAG=""
    if [ -n "$CUDA_VISIBLE_DEVICES" ]; then
       USE_GPU_FLAG="--use_gpu"
    fi
    set -u
    
    python3 scripts/integration.py run_integration \
      --hvg_mat_f     {input.hvg_mat_f} \
      --dbl_hvg_mat_f {input.dbl_hvg_mat_f} \
      --sample_qc_f   {input.qc_stats_f} \
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
      $USE_GPU_FLAG
    """


# rule to create sce objects without any doublets (and delete temporary sce objects in the qc directory)
rule make_clean_h5ads: 
  input:
    integration_f = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    h5_paths_f    = f'{hvg_dir}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv',
    coldata_f     = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    rowdata_f     = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
  output:
    clean_h5ad_f  = f'{int_dir}/anndata_cells_clean_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5ad'
  params:
    sel_run       = lambda wildcards: get_run_for_one_batch(wildcards.batch, RUNS_TO_BATCHES),
    run_var       = RUN_VAR,
    batch_var     = BATCH_VAR
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'make_clean_h5ads', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'make_clean_h5ads', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_integration/make_clean_h5ads_{{batch}}_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/integration.yaml'
  shell: """
    python3 scripts/make_clean_h5ad.py \
      {wildcards.batch} \
      {params.sel_run} \
      {input.integration_f} \
      {input.h5_paths_f} \
      {input.coldata_f} \
      {input.rowdata_f} \
      {params.run_var} \
      {params.batch_var} \
      {output.clean_h5ad_f}
    """


# make a yaml with all clean sce file paths
rule make_clean_h5ad_paths_yaml:
  input:
    clean_h5ad_fs  = expand(f'{int_dir}/anndata_cells_clean_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5ad', batch = BATCHES)
  output:
    h5ads_yaml_f   = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml'
  run:
    # split paths and batch names
    fs = [f"{int_dir}/anndata_cells_clean_{batch}_{FULL_TAG}_{DATE_STAMP}.h5ad" for batch in BATCHES]
    fs_dict = dict(zip(BATCHES, fs))

    # write to yaml
    with open(output.h5ads_yaml_f, 'w') as f:
      yaml.dump(fs_dict, f, default_flow_style=False)


rule convert_h5ad_to_sce: 
  input:
    clean_h5ad_f  = f'{int_dir}/anndata_cells_clean_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5ad'
  output:
    clean_sce_f   = f'{int_dir}/sce_cells_clean_{{batch}}_{FULL_TAG}_{DATE_STAMP}.rds'
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'convert_h5ad_to_sce', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'convert_h5ad_to_sce', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_integration/convert_h5ad_to_sce_{{batch}}_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/rlibs.yaml'
  shell:"""
    Rscript -e "source('scripts/integration.R');
    make_clean_sce_from_h5ad(
      sel_batch  = '{wildcards.batch}', 
      adata_f    = '{input.clean_h5ad_f}',
      sce_f      = '{output.clean_sce_f}'
    )"
    """