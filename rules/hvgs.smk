
import yaml  

localrules: make_hvg_df

rule make_hvg_df:
  input:
    ambient_yml_out = expand([f'{amb_dir}/ambient_{{run}}/ambient_{{run}}_{DATE_STAMP}_output_paths.yaml'], run = RUNS)
  output:
    hvg_paths_f     = f'{hvg_dir}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv'
  params:
    run_var         = RUN_VAR,
    batch_var       = BATCH_VAR,
    demux_type      = config['multiplexing']['demux_type']
  run:
    hvg_df = make_hvgs_input_df(RUNS, input.ambient_yml_out, params.run_var, params.batch_var, 
      RUNS_TO_BATCHES, params.demux_type, FULL_TAG, DATE_STAMP, hvg_dir)
    hvg_df.write_csv(output.hvg_paths_f)



# create temporary csr h5 files
rule make_tmp_csr_matrix:
  input:
    hvg_paths_f     = f'{hvg_dir}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv', 
    cell_filter_f   = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    qc_stats_f      = f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    rowdata_f       = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    clean_h5_f      = temp(expand(f'{hvg_dir}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5', batch = BATCHES))
  params:
    run_var         = RUN_VAR,
    batch_var       = BATCH_VAR,
    demux_type      = config['multiplexing']['demux_type'],
    hvg_chunk_size  = config['hvg']['hvg_chunk_size']
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'make_tmp_csr_matrix', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'make_tmp_csr_matrix', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_hvgs/make_tmp_csr_matrix_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    python3 scripts/hvgs.py get_csr_counts \
      {input.hvg_paths_f} \
      {input.cell_filter_f} \
      "keep" \
      {input.qc_stats_f} \
      {input.rowdata_f} \
      {params.run_var} \
      {params.batch_var} \
      {params.demux_type} \
      --chunksize {params.hvg_chunk_size} \
      --ncores {threads}
    """


if config['hvg']['hvg_method'] == 'sample': 
  localrules: merge_sample_std_var_stats
  # calculate stats for each sample separately  
  rule get_stats_for_std_variance_for_sample:
    input:
      clean_h5_f      = f'{hvg_dir}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5', 
      qc_stats_f      = f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
      rowdata_f       = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
    output:
      std_var_stats_f = temp(f'{hvg_dir}/tmp_std_var_stats_{{batch}}_{FULL_TAG}_{DATE_STAMP}.csv.gz')
    params:
      batch_var       = BATCH_VAR
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'get_stats_for_std_variance_for_sample', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'get_stats_for_std_variance_for_sample', 'time', attempt)
    benchmark:
      f'{benchmark_dir}/{SHORT_TAG}_hvgs/get_stats_for_std_variance_for_sample_{{batch}}_{DATE_STAMP}.benchmark.txt'
    conda:
      '../envs/hvgs.yaml'
    shell: """
      python3 scripts/hvgs.py calculate_std_var_stats_for_sample \
        {wildcards.batch} \
        {params.batch_var} \
        {input.qc_stats_f} \
        {input.clean_h5_f} \
        {input.rowdata_f} \
        {output.std_var_stats_f}
      """

  rule merge_sample_std_var_stats:
    input:                 
      std_var_stats_f = expand(f'{hvg_dir}/tmp_std_var_stats_{{batch}}_{FULL_TAG}_{DATE_STAMP}.csv.gz', batch = BATCHES)
    output:
      std_var_stats_merged_f = f'{hvg_dir}/standardized_variance_stats_{FULL_TAG}_{DATE_STAMP}.csv.gz'
    threads: 1
    benchmark:
      f'{benchmark_dir}/{SHORT_TAG}_hvgs/merge_sample_std_var_stats_{DATE_STAMP}.benchmark.txt'
    run:
      merge_tmp_files(input.std_var_stats_f, output.std_var_stats_merged_f)


else:
  # define some rules that don't need the cluster
  localrules: merge_group_mean_var, merge_group_std_var_stats

  rule get_mean_var_for_group:
    input:
      clean_h5_f      = expand(f'{hvg_dir}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5', batch = BATCHES),
      hvg_paths_f     = f'{hvg_dir}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv',
      rowdata_f       = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
      qc_stats_f      = f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv'
    output: 
      mean_var_f      = temp(f'{hvg_dir}/tmp_mean_var_{{group}}_chunk_{{chunk}}_{FULL_TAG}_{DATE_STAMP}.csv.gz')
    params:
      metadata_f      = config['project']['sample_metadata'],
      batch_var       = BATCH_VAR,
      hvg_method      = config['hvg']['hvg_method'],
      hvg_chunk_size  = config['hvg']['hvg_chunk_size'],
      hvg_group_var   = config['hvg']['hvg_metadata_split_var']
    threads: 8
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'get_mean_var_for_group', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'get_mean_var_for_group', 'time', attempt)
    benchmark:
      f'{benchmark_dir}/{SHORT_TAG}_hvgs/get_mean_var_for_group_{{group}}_chunk_{{chunk}}_{DATE_STAMP}.benchmark.txt'
    conda:
      '../envs/hvgs.yaml'
    shell: """
      GROUPVAR_FLAG=""
      if [ "{params.hvg_method}" = "groups" ]; then
        GROUPVAR_FLAG="--groupvar {params.hvg_group_var}"
      fi

      python3 scripts/hvgs.py calculate_mean_var_for_chunk \
        {input.hvg_paths_f} \
        {input.rowdata_f} \
        {params.metadata_f} \
        {input.qc_stats_f} \
        {output.mean_var_f} \
        {wildcards.chunk} \
        {params.hvg_method} \
        {params.batch_var} \
        --chunksize {params.hvg_chunk_size} \
        --group {wildcards.group} \
        --ncores {threads} \
        $GROUPVAR_FLAG
      """


  rule merge_group_mean_var:
    input:
      mean_var_f  = expand(
        f'{hvg_dir}/tmp_mean_var_{{group}}_chunk_{{chunk}}_{FULL_TAG}_{DATE_STAMP}.csv.gz',
        group = config['hvg']['hvg_group_names'], 
        chunk = range(config['hvg']['hvg_num_chunks'])
        )
    output:
      mean_var_merged_f = temp(f'{hvg_dir}/means_variances_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz')
    threads: 1
    run:
      merge_tmp_files(input.mean_var_f, output.mean_var_merged_f)


  rule get_estimated_variances:
    input:
      mean_var_merged_f = f'{hvg_dir}/means_variances_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
    output:
      estim_vars_f      = temp(f'{hvg_dir}/estimated_variances_{FULL_TAG}_{DATE_STAMP}.csv.gz')
    params:
      hvg_method        = config['hvg']['hvg_method'],
      batch_var         = BATCH_VAR
    threads: 1
    retries: config['resources']['retries']
    conda:
      '../envs/hvgs.yaml'
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'get_estimated_variances', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'get_estimated_variances', 'time', attempt)
    benchmark:
      f'{benchmark_dir}/{SHORT_TAG}_hvgs/get_estimated_variances_{DATE_STAMP}.benchmark.txt'
    shell: """
      python3 scripts/hvgs.py calculate_estimated_vars \
        {output.estim_vars_f} \
        {params.hvg_method} \
        {params.batch_var} \
        {input.mean_var_merged_f}
      """


  rule get_stats_for_std_variance_for_group:
    input: 
      clean_h5_fs     = expand(f'{hvg_dir}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5', batch = BATCHES),
      estim_vars_f    = f'{hvg_dir}/estimated_variances_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
      hvg_paths_f     = f'{hvg_dir}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv',
      rowdata_f       = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
      qc_stats_f      = f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    output:
      std_var_stats_f = temp(f'{hvg_dir}/tmp_std_var_stats_{{group}}_chunk_{{chunk}}_{FULL_TAG}_{DATE_STAMP}.csv.gz')
    params:
      metadata_f      = config['project']['sample_metadata'],
      batch_var       = BATCH_VAR,
      hvg_method      = config['hvg']['hvg_method'],
      hvg_chunk_size  = config['hvg']['hvg_chunk_size'],
      hvg_group_var   = config['hvg']['hvg_metadata_split_var']
    threads: 8
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'get_stats_for_std_variance_for_group', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'get_stats_for_std_variance_for_group', 'time', attempt)
    benchmark:
      f'{benchmark_dir}/{SHORT_TAG}_hvgs/get_stats_for_std_variance_for_group_{{group}}_chunk_{{chunk}}_{DATE_STAMP}.benchmark.txt'
    conda:
      '../envs/hvgs.yaml'
    shell: """
      python3 scripts/hvgs.py calculate_std_var_stats_for_chunk \
        {input.hvg_paths_f} \
        {input.rowdata_f} \
        {params.metadata_f} \
        {input.qc_stats_f} \
        {output.std_var_stats_f} \
        {input.estim_vars_f} \
        {wildcards.chunk} \
        {params.hvg_method} \
        {params.batch_var} \
        --chunksize {params.hvg_chunk_size} \
        --group {wildcards.group} \
        --groupvar {params.hvg_group_var} \
        --ncores {threads}
      """


  rule merge_group_std_var_stats:
    input:                 
      std_var_stats_f = expand(
        f'{hvg_dir}/tmp_std_var_stats_{{group}}_chunk_{{chunk}}_{FULL_TAG}_{DATE_STAMP}.csv.gz',
        group = config['hvg']['hvg_group_names'], 
        chunk = range(config['hvg']['hvg_num_chunks'])
        )
    output:
      std_var_stats_merged_f = temp(f'{hvg_dir}/standardized_variance_stats_{FULL_TAG}_{DATE_STAMP}.csv.gz')
    benchmark:
      f'{benchmark_dir}/{SHORT_TAG}_hvgs/merge_group_std_var_stats_{DATE_STAMP}.benchmark.txt'
    run:
      merge_tmp_files(input.std_var_stats_f, output.std_var_stats_merged_f)


rule get_highly_variable_genes:
  input:
    std_var_stats_f = f'{hvg_dir}/standardized_variance_stats_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    empty_gs_fs     = f'{empty_dir}/edger_empty_genes_all_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    hvg_f           = f'{hvg_dir}/hvg_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  threads: 1
  retries: config['resources']['retries']
  params:
    hvg_method  = config['hvg']['hvg_method'],
    batch_var   = BATCH_VAR,
    n_hvgs      = config['hvg']['hvg_n_hvgs'],
    no_ambient  = config['hvg']['hvg_exclude_ambient_genes']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'get_highly_variable_genes', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'get_highly_variable_genes', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_hvgs/get_highly_variable_genes_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    NOAMBIENT_FLAG=""
    if [ "{params.no_ambient}" = "True" ]; then
      NOAMBIENT_FLAG="--noambient"
    fi

    python3 scripts/hvgs.py calculate_hvgs \
      {input.std_var_stats_f} \
      {output.hvg_f} \
      {input.empty_gs_fs} \
      {params.hvg_method} \
      {params.batch_var} \
      {params.n_hvgs} \
      $NOAMBIENT_FLAG
     """


rule create_hvg_matrix:
  input: 
    clean_h5_f  = expand(f'{hvg_dir}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5', batch = BATCHES),
    qc_stats_f  = f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    hvg_paths_f = f'{hvg_dir}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv', 
    hvg_f       = f'{hvg_dir}/hvg_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
  output:
    hvg_mat_f   = f'{hvg_dir}/top_hvgs_counts_{FULL_TAG}_{DATE_STAMP}.h5'
  params:
    demux_type  = config['multiplexing']['demux_type'],
    batch_var   = BATCH_VAR
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'create_hvg_matrix', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'create_hvg_matrix', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_hvgs/create_hvg_matrix_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    python3 scripts/hvgs.py create_hvg_matrix \
      {input.qc_stats_f} \
      {input.hvg_paths_f} \
      {input.hvg_f} \
      {output.hvg_mat_f} \
      {params.demux_type} \
      {params.batch_var}
    """


rule create_doublets_hvg_matrix:
  input: 
    hvg_paths_f   = f'{hvg_dir}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv', 
    hvg_f         = f'{hvg_dir}/hvg_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    qc_f          = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    qc_stats_f    = f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv'
  output: 
    dbl_hvg_mat_f = f'{hvg_dir}/top_hvgs_doublet_counts_{FULL_TAG}_{DATE_STAMP}.h5'
  params:
    run_var       = RUN_VAR,
    demux_type    = config['multiplexing']['demux_type'],
    batch_var     = BATCH_VAR
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'create_doublets_hvg_matrix', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'create_doublets_hvg_matrix', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_hvgs/create_doublets_hvg_matrix_{DATE_STAMP}.benchmark.txt'
  conda: 
    '../envs/hvgs.yaml'
  shell: """
    python3 scripts/hvgs.py create_doublets_matrix \
      {input.hvg_paths_f} \
      {input.hvg_f} \
      {input.qc_f} \
      {input.qc_stats_f} \
      {output.dbl_hvg_mat_f} \
      {params.run_var} \
      {params.demux_type} \
      {params.batch_var}
    """

