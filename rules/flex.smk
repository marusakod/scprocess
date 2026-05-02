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

    python3 scripts/mapping.py map_flex_fastqs_to_counts {wildcards.lib} \
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


# rule for saving per-sample h5, knee data, and ambient params - runs per sample/run
rule save_alevin_flex_to_h5:
  input:
    fry_dir     = lambda wildcards: f'{af_dir}/{lib_pool_dir}af_{RUNS_TO_LIBS[wildcards.run]}/flex/af_quant/'
  output:
    af_h5_f     = f'{af_dir}/af_{{run}}/flex/af_counts_mat.h5',
    amb_yaml_f  = f'{af_dir}/af_{{run}}/flex/ambient_params_{{run}}_{DATE_STAMP}.yaml',
    knee_data_f = f'{af_dir}/af_{{run}}/flex/knee_plot_data_{{run}}_{DATE_STAMP}.csv.gz'
  params:
    probe_id      = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["probe_id"],
    knee1         = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["knee1"],
    shin1         = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["shin1"],
    knee2         = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["knee2"],
    shin2         = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["shin2"],
    exp_cells     = lambda wildcards: RUN_PARAMS[wildcards.run]["ambient"]["cb_expected_cells"],
    total_inc     = lambda wildcards: RUN_PARAMS[wildcards.run]["ambient"]["cb_total_droplets_included"],
    low_count_thr = lambda wildcards: RUN_PARAMS[wildcards.run]["ambient"]["cb_low_count_threshold"]
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'save_alevin_flex_to_h5', 'memory', attempt, wildcards.run),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'save_alevin_flex_to_h5', 'time', attempt, wildcards.run)
  benchmark:
    f'{benchmark_dir}/mapping/save_alevin_flex_to_h5_{{run}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/mapping/save_alevin_flex_to_h5_{{run}}_{DATE_STAMP}.log'
  conda:
    '../envs/rlibs.yaml'
  shell: """
    exec &>> {log}

    Rscript -e "source('scripts/mapping.R');
      save_flex_alevin_h5_ambient_params(
        run           = '{wildcards.run}',
        fry_dir       = '{input.fry_dir}',
        probe_id      = '{params.probe_id}',
        h5_f          = '{output.af_h5_f}',
        cb_yaml_f     = '{output.amb_yaml_f}',
        knee_data_f   = '{output.knee_data_f}',
        run_var       = '{RUN_VAR}',
        knee1         = '{params.knee1}',
        shin1         = '{params.shin1}',
        knee2         = '{params.knee2}',
        shin2         = '{params.shin2}',
        exp_cells     = '{params.exp_cells}',
        total_included= '{params.total_inc}',
        low_count_thr = '{params.low_count_thr}'
      )"
    """
