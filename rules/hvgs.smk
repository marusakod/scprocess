
import yaml  


localrules: make_hvg_df


rule make_hvg_df:
  input:
    ambient_yaml_out=expand([amb_dir + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml'], run=RUNS)
  output:
    hvg_paths_f     = hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv' 
  params:
    run_var         = RUN_VAR,
    demux_type      = config['multiplexing']['demux_type']
  run:
    hvg_df = make_hvgs_input_df(params.demux_type, params.run_var, RUNS, input.ambient_yaml_out,
      SAMPLES_TO_RUNS, FULL_TAG, DATE_STAMP, hvg_dir)
    hvg_df.to_csv(output.hvg_paths_f, index=False)


# create temporary csr h5 files
rule make_tmp_csr_matrix:
  input:
    hvg_paths_f       = hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    cell_filter_f     = qc_dir  + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz',
    qc_sample_stats_f = qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    rowdata_f         = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  output:
    clean_h5_f        = temp(expand(hvg_dir + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5', sample = SAMPLES))
  params:
    run_var           = RUN_VAR,
    demux_type        = config['multiplexing']['demux_type'],
    hvg_chunk_size    = config['hvg']['hvg_chunk_size']
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_run_hvgs'] * MB_PER_GB
  conda:
    '../envs/hvgs.yaml'
  shell: """
    python3 scripts/hvgs.py get_csr_counts \
      {input.hvg_paths_f} \
      {input.cell_filter_f} \
      "keep" \
      "True" \
      {input.qc_sample_stats_f} \
      {input.rowdata_f} \
      {params.run_var} \
      {params.demux_type} \
      --chunksize {params.hvg_chunk_size} \
      --ncores {threads}
    """


if config['hvg']['hvg_method'] == 'sample': 
  localrules: merge_sample_std_var_stats
  # calculate stats for each sample separately  
  rule get_stats_for_std_variance_for_sample:
    input:
      clean_h5_f      = hvg_dir + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
      qc_smpl_stats_f = qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
      rowdata_f       = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
    output:
      std_var_stats_f = temp(hvg_dir + '/tmp_std_var_stats_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_run_hvgs'] * MB_PER_GB
    conda:
      '../envs/hvgs.yaml'
    shell: """
      python3 scripts/hvgs.py calculate_std_var_stats_for_sample \
        {wildcards.sample} \
        {input.qc_smpl_stats_f} \
        {input.clean_h5_f} \
        {input.rowdata_f} \
        {output.std_var_stats_f}
      """

  rule merge_sample_std_var_stats:
    input:                 
      std_var_stats_f = expand(hvg_dir + '/tmp_std_var_stats_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', sample = SAMPLES)
    output:
      std_var_stats_merged_f = hvg_dir + '/standardized_variance_stats_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_run_hvgs'] * MB_PER_GB
    run:
      merge_tmp_files(input.std_var_stats_f, output.std_var_stats_merged_f)


else:
  # define some rules that don't need the cluster
  localrules: merge_group_mean_var, merge_group_std_var_stats


  rule get_mean_var_for_group:
    input:
      clean_h5_f      = expand(hvg_dir + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5', sample = SAMPLES),
      hvg_paths_f     = hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
      rowdata_f       = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', 
      qc_smpl_stats_f = qc_dir + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
    output: 
      mean_var_f      = temp(hvg_dir + '/tmp_mean_var_{group}_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
    params:
      hvg_method      = config['hvg']['hvg_method'],
      hvg_chunk_size  = config['hvg']['hvg_chunk_size'],
      hvg_group_var   = config['hvg']['hvg_metadata_split_var']
    threads: 8
    resources:
      mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_run_hvgs'] * MB_PER_GB
    conda:
      '../envs/hvgs.yaml'
    shell: """
      python3 scripts/hvgs.py calculate_mean_var_for_chunk \
        {input.hvg_paths_f} \
        {input.rowdata_f} \
        {METADATA_F} \
        {input.qc_smpl_stats_f} \
        {output.mean_var_f} \
        {wildcards.chunk} \
        {params.hvg_method} \
        {params.hvg_chunk_size} \
        --group {wildcards.group} \
        --groupvar {params.hvg_group_var} \
        --ncores {threads} 
    
      """


  rule merge_group_mean_var:
    input:
      mean_var_f  = expand(
        hvg_dir + '/tmp_mean_var_{group}_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
        group = config['hvg']['hvg_group_names'], 
        chunk = range(config['hvg']['hvg_num_chunks'])
        )
    output:
      mean_var_merged_f = temp(hvg_dir + '/means_variances_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_run_hvgs'] * MB_PER_GB
    run:
      merge_tmp_files(input.mean_var_f, output.mean_var_merged_f)


  rule get_estimated_variances:
    input:
      mean_var_merged_f = hvg_dir + '/means_variances_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
    output:
      estim_vars_f      = temp(hvg_dir + '/estimated_variances_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
    params:
      hvg_method        = config['hvg']['hvg_method']
    threads: 1
    retries: config['resources']['retries']
    conda:
      '../envs/hvgs.yaml'
    resources:
      mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_run_hvgs'] * MB_PER_GB
    shell: """
      python3 scripts/hvgs.py calculate_estimated_vars \
        {output.estim_vars_f} \
        {params.hvg_method} \
        {input.mean_var_merged_f}
      """


  rule get_stats_for_std_variance_for_group:
    input: 
      clean_h5_fs     = expand(hvg_dir + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5', sample = SAMPLES),
      estim_vars_f    = hvg_dir + '/estimated_variances_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
      hvg_paths_f     = hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
      rowdata_f       = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz', 
      qc_smpl_stats_f = qc_dir + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
    output:
      std_var_stats_f = temp(hvg_dir + '/tmp_std_var_stats_{group}_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
    params:
      hvg_method      = config['hvg']['hvg_method'],
      hvg_chunk_size  = config['hvg']['hvg_chunk_size'],
      hvg_group_var   = config['hvg']['hvg_metadata_split_var']
    threads: 8
    resources:
      mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_run_hvgs'] * MB_PER_GB
    conda:
      '../envs/hvgs.yaml'
    shell: """
      python3 scripts/hvgs.py calculate_std_var_stats_for_chunk \
        {input.hvg_paths_f} \
        {input.rowdata_f} \
        {METADATA_F} \
        {input.qc_smpl_stats_f} \
        {output.std_var_stats_f} \
        {input.estim_vars_f} \
        {wildcards.chunk} \
        {params.hvg_method} \
        --chunksize {params.hvg_chunk_size} \
        --group {wildcards.group} \
        --groupvar {params.hvg_group_var} \
        --ncores {threads}
      """


  rule merge_group_std_var_stats:
    input:                 
      std_var_stats_f = expand(
        hvg_dir + '/tmp_std_var_stats_{group}_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
        group = config['hvg']['hvg_group_names'], 
        chunk = range(config['hvg']['hvg_num_chunks'])
        )
    output:
      std_var_stats_merged_f = temp(hvg_dir + '/standardized_variance_stats_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_run_hvgs'] * MB_PER_GB
    run:
      merge_tmp_files(input.std_var_stats_f, output.std_var_stats_merged_f)


rule get_highly_variable_genes:
  input:
    std_var_stats_f = hvg_dir + '/standardized_variance_stats_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    empty_gs_fs     = empty_dir + '/edger_empty_genes_' + FULL_TAG + '_all_' + DATE_STAMP + '.txt.gz' 
  output:
    hvg_f = hvg_dir + '/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  threads: 1
  retries: config['resources']['retries']
  params:
    hvg_method = config['hvg']['hvg_method'],
    n_hvgs     = config['hvg']['hvg_n_hvgs'],
    no_ambient = config['hvg']['hvg_exclude_ambient_genes']
  resources:
     mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_run_hvgs'] * MB_PER_GB
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
      {params.n_hvgs} \
      $NOAMBIENT_FLAG
     """


rule create_hvg_matrix:
  input: 
    clean_h5_f      = expand(hvg_dir + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5', sample = SAMPLES),
    qc_smpl_stats_f = qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    hvg_paths_f     = hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    hvg_f           = hvg_dir + '/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
  output:
    hvg_mat_f       = hvg_dir + '/top_hvgs_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5'
  params:
    demux_type      = config['multiplexing']['demux_type']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_run_hvgs'] * MB_PER_GB
  conda:
    '../envs/hvgs.yaml'
  shell: """
    python3 scripts/hvgs.py read_top_genes \
      {input.qc_smpl_stats_f} \
      {input.hvg_paths_f} \
      {input.hvg_f} \
      {output.hvg_mat_f} \
      {params.demux_type}
    """


rule create_doublets_hvg_matrix:
  input: 
    hvg_paths_f       = hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    hvg_f             = hvg_dir + '/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    qc_f              = qc_dir  + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz',
    qc_sample_stats_f = qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output: 
    dbl_hvg_mat_f     = hvg_dir + '/top_hvgs_doublet_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5'
  params:
    run_var           = RUN_VAR,
    demux_type        = config['multiplexing']['demux_type']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_run_hvgs'] * MB_PER_GB
  conda: 
    '../envs/hvgs.yaml'
  shell: """
    python3 scripts/hvgs.py create_doublets_matrix \
    {input.hvg_paths_f} \
    {input.hvg_f} \
    {input.qc_f} \
    {input.qc_sample_stats_f} \
    {output.dbl_hvg_mat_f} \
    {params.run_var} \
    {params.demux_type}
    """

