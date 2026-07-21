# rule for flex mapping - runs per pool/library
rule run_mapping_flex:
  params:
    arv_instance  = config['project'].get('arv_instance', ""),
    af_home_dir   = config['mapping_af']['alevin_fry_home'],
    af_index_dir  = config['mapping_af']['af_index_dir'],
    probeset_f    = config['mapping_af']['probeset_f'],
    probe_bcs_f   = config['mapping_af']['probe_bcs_f'],
    af_chemistry  = lambda wildcards: LIB_PARAMS[wildcards.lib]["mapping_af"]["af_chemistry"],
    whitelist_f   = lambda wildcards: LIB_PARAMS[wildcards.lib]["mapping_af"]["gex_whitelist_f"],
    where         = lambda wildcards: LIB_PARAMS[wildcards.lib]["mapping_af"]["where"],
    R1_fs         = lambda wildcards: LIB_PARAMS[wildcards.lib]["mapping_af"]["R1_fs"],
    R2_fs         = lambda wildcards: LIB_PARAMS[wildcards.lib]["mapping_af"]["R2_fs"],
    lib_pool_dir  = lib_pool_dir
  output:
    rad_f         = temp(f'{af_dir}/{lib_pool_dir}af_{{lib}}/flex/af_map/map.rad'),
    collate_rad_f = temp(f'{af_dir}/{lib_pool_dir}af_{{lib}}/flex/af_quant/map.collated.rad'),
    fry_dir       = directory(f'{af_dir}/{lib_pool_dir}af_{{lib}}/flex/af_quant/'),
    mtx_f         = f'{af_dir}/{lib_pool_dir}af_{{lib}}/flex/af_quant/alevin/quants_mat.mtx',
    cols_f        = f'{af_dir}/{lib_pool_dir}af_{{lib}}/flex/af_quant/alevin/quants_mat_cols.txt',
    rows_f        = f'{af_dir}/{lib_pool_dir}af_{{lib}}/flex/af_quant/alevin/quants_mat_rows.txt'
  benchmark:
    f'{benchmark_dir}/mapping/run_mapping_flex_{{lib}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/mapping/run_mapping_flex_{{lib}}_{DATE_STAMP}.log'
  threads: config['resources']['n_run_mapping_flex']
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_mapping_flex', 'memory', attempt, wildcards.lib),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_mapping_flex', 'time', attempt, wildcards.lib)
  conda:
    '../envs/alevin_fry.yaml'
  shell:"""
    exec &>> {log}

    ARV_ARG=""
    if [[ "{params.arv_instance}" != "" ]]; then
      ARV_ARG="--arv_instance {params.arv_instance}"
    fi

    python3 {scprocess_dir}/scripts/mapping.py map_flex_fastqs_to_counts {wildcards.lib} \
      --af_dir          "{af_dir}" \
      --lib_pool_dir    "{params.lib_pool_dir}" \
      --af_home_dir     "{params.af_home_dir}" \
      --where           "{params.where}" \
      --R1_fs           {params.R1_fs} \
      --R2_fs           {params.R2_fs} \
      --threads         {threads} \
      --af_index_dir    "{params.af_index_dir}" \
      --af_chemistry    "{params.af_chemistry}" \
      --gex_whitelist_f "{params.whitelist_f}" \
      --probeset_f      "{params.probeset_f}" \
      --probe_bc_f      "{params.probe_bcs_f}" \
      $ARV_ARG
    """
