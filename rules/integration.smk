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


def int_get_filt_counts_f(run):
  ambient_method = config['ambient']['ambient_method']
  if ambient_method == "cellbender":
    return f'{amb_dir}/ambient_{run}/bender_{run}_{DATE_STAMP}_filtered.h5'
  elif ambient_method == "decontx":
    return f'{amb_dir}/ambient_{run}/decontx_{run}_{DATE_STAMP}_filtered.h5'
  elif ambient_method == "none":
    return f'{amb_dir}/ambient_{run}/uncorrected_{run}_{DATE_STAMP}_filtered.h5'

# rule to create sce objects without any doublets (and delete temporary sce objects in the qc directory)
rule make_clean_sces: 
  input:
    integration_f = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    h5_paths_f    = f'{hvg_dir}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv',
    h5_filt_fs    = [int_get_filt_counts_f(run) for run in RUNS],
    coldata_f     = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz',
  output:
    clean_sce_f   = f'{int_dir}/sce_cells_clean_{{batch}}_{FULL_TAG}_{DATE_STAMP}.rds'
  params:
    sel_run       = lambda wildcards: get_run_for_one_batch(wildcards.batch, RUNS_TO_BATCHES),
    gtf_dt_f      = config['mapping']['af_gtf_dt_f'],
    run_var       = RUN_VAR,
    batch_var     = BATCH_VAR,
    mito_str      = config['mapping']['af_mito_str'],
    exclude_mito  = config['qc']['exclude_mito']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'make_clean_sces', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'make_clean_sces', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_integration/make_clean_sces_{{batch}}_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/utils.R'); source('scripts/SampleQC.R');
      make_clean_sces(
        sel_b         = '{wildcards.batch}',
        sel_run       = '{params.sel_run}',
        integration_f = '{input.integration_f}',
        h5_paths_f    = '{input.h5_paths_f}',
        coldata_f     = '{input.coldata_f}',
        gtf_dt_f      = '{params.gtf_dt_f}',
        run_var       = '{params.run_var}',
        batch_var     = '{params.batch_var}',
        mito_str      = '{params.mito_str}',
        exclude_mito  = '{params.exclude_mito}',
        clean_sce_f   = '{output.clean_sce_f}'
      )"
    """


# make a yaml with all clean sce file paths
rule make_clean_sce_paths_yaml:
  input:
    clean_sce_fs  = expand(f'{int_dir}/sce_cells_clean_{{batch}}_{FULL_TAG}_{DATE_STAMP}.rds', batch = BATCHES)
  output:
    sces_yaml_f   = f'{int_dir}/sce_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml'
  run:
    # split paths and batch names
    fs = [f"{int_dir}/sce_cells_clean_{batch}_{FULL_TAG}_{DATE_STAMP}.rds" for batch in BATCHES]
    
    # check that all files exist
    for f in fs:
      if not os.path.isfile(f):
        raise FileNotFoundError(f"File {f} doesn't exist")

    # create a dictionary
    fs_dict = dict(zip(BATCHES, fs))

    # write to yaml
    with open(output.sces_yaml_f, 'w') as f:
      yaml.dump(fs_dict, f, default_flow_style=False)

