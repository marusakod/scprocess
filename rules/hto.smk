
# build index for HTO mapping
rule build_hto_index:
  input:
    feature_ref_f = config["multiplexing"]["feature_ref"]
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'build_hto_index', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'build_hto_index', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_hto/build_hto_index_{DATE_STAMP}.benchmark.txt'
  output: 
    hto_f       = f'{af_dir}/hto.tsv',
    t2g_f       = f'{af_dir}/t2g_hto.tsv',
    idx_log_f   = f'{af_dir}/hto_index/ref_indexing.log',
    hto_idx_dir = directory(f'{af_dir}/hto_index')
  conda:
    '../envs/alevin_fry.yaml'
  shell: """
    cd {af_dir}
    # create a tsv file from feature ref file
    awk -F "," '
      NR == 1 {{
        for (i = 1; i <= NF; i++) {{
          if ($i == "hto_id") hto_id_col = i
          if ($i == "sequence") seq_col = i
        }}
      }}
      NR > 1 {{
        print $hto_id_col "\t" $seq_col
      }}' {input.feature_ref_f} > {output.hto_f}

    # 
    salmon index -t {output.hto_f} -i hto_index --features -k7
    
    # Make gene-to-transcript mapping file
    awk '{{print $1"\t"$1;}}' {output.hto_f} > {output.t2g_f}
    """

# run mapping for HTO files
rule run_mapping_hto:
  input: 
    hto_idx_dir   = f'{af_dir}/hto_index', 
    chem_stats_f  = f'{af_dir}/af_{{run}}/{af_rna_dir}chemistry_statistics.yaml'
  output:
    fry_dir       = directory(f'{af_dir}/af_{{run}}/hto/af_quant/'),
    rad_f         = temp(f'{af_dir}/af_{{run}}/hto/af_map/map.rad'),
    collate_rad_f = temp(f'{af_dir}/af_{{run}}/hto/af_quant/map.collated.rad'),
    mtx_f         = f'{af_dir}/af_{{run}}/hto/af_quant/alevin/quants_mat.mtx',
    cols_f        = f'{af_dir}/af_{{run}}/hto/af_quant/alevin/quants_mat_cols.txt',
    rows_f        = f'{af_dir}/af_{{run}}/hto/af_quant/alevin/quants_mat_rows.txt'
  params:
    demux_type    = config['multiplexing']['demux_type'],
    af_home_dir   = config['mapping']['alevin_fry_home'],
    wl_lu_f       = config['mapping']['wl_lu_f'], 
    where         = lambda wildcards: RUN_PARAMS[wildcards.run]["multiplexing"]["where"],
    R1_fs         = lambda wildcards: RUN_PARAMS[wildcards.run]["multiplexing"]["R1_fs"],
    R2_fs         = lambda wildcards: RUN_PARAMS[wildcards.run]["multiplexing"]["R2_fs"],
    af_chemistry  = lambda wildcards: RUN_PARAMS[wildcards.run]["multiplexing"]["af_chemistry"],
    expected_ori  = "fw",
    whitelist_f   = lambda wildcards: RUN_PARAMS[wildcards.run]["multiplexing"]["whitelist_f"]
  conda:
    '../envs/alevin_fry.yaml'
  threads: config['resources']['n_run_mapping']
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_mapping_hto', 'memory', attempt, wildcards.run),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_mapping_hto', 'time', attempt, wildcards.run)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_hto/run_mapping_hto_{{run}}_{DATE_STAMP}.benchmark.txt'
  shell:"""
    # check whether doing arvados
    ARV_REGEX="^arkau-[0-9a-z]{{5}}-[0-9a-z]{{15}}$"
    if [[ "{params.where}" =~ $ARV_REGEX ]]; then
      ml arvados
      arv-env arkau
    fi

    # get chemistry parameters from input yaml file
    AF_CHEMISTRY=$(grep "selected_af_chemisty:" {input.chem_stats_f} | sed 's/selected_af_chemisty: //')
    EXP_ORI=$(grep "selected_ori:" {input.chem_stats_f} | sed 's/selected_ori: //')
    WHITELIST_F=$(grep "selected_whitelist:" {input.chem_stats_f} | sed 's/selected_whitelist: //')

    # run mapping
    python3 scripts/mapping.py {wildcards.run} \
      --af_dir          {af_dir}\
      --demux_type      {params.demux_type} \
      --what            "hto" \
      --af_home_dir     {params.af_home_dir} \
      --where           {params.where} \
      --R1_fs           {params.R1_fs} \
      --R2_fs           {params.R2_fs} \
      --threads         {threads} \
      --af_index_dir    {input.hto_idx_dir} \
      --wl_lu_f         {params.wl_lu_f} \
      --tenx_chemistry  $AF_CHEMISTRY \
      --exp_ori         $EXP_ORI \
      --whitelist_f     $WHITELIST_F
    """


rule save_alevin_hto_to_h5:
  input: 
    fry_dir     = f'{af_dir}/af_{{run}}/hto/af_quant/'
  output: 
    h5_f        = f'{af_dir}/af_{{run}}/hto/af_hto_counts_mat.h5',
    knee_data_f = f'{af_dir}/af_{{run}}/hto/knee_plot_data_{{run}}_{DATE_STAMP}.csv.gz'
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'save_alevin_hto_to_h5', 'memory', attempt, wildcards.run),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'save_alevin_hto_to_h5', 'time', attempt, wildcards.run)
  benchmark: 
    f'{benchmark_dir}/{SHORT_TAG}_hto/save_alevin_hto_to_h5_{{run}}_{DATE_STAMP}.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/mapping.R');
      save_alevin_h5_knee_params_df(
        run         = '{wildcards.run}', 
        fry_dir     = '{input.fry_dir}', 
        hto_mat     = 1, 
        run_var     = '{RUN_VAR}',
        h5_f        = '{output.h5_f}',
        knee_data_f = '{output.knee_data_f}'
      )"
    """


# save sce object with hto counts and demultiplex
rule make_hto_sce_objects: 
  input: 
    smpl_stats_f = f'{amb_dir}/ambient_run_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    amb_yaml_f   = f'{amb_dir}/ambient_{{run}}/ambient_{{run}}_{DATE_STAMP}_output_paths.yaml',
    hto_h5_f     = f'{af_dir}/af_{{run}}/hto/af_hto_counts_mat.h5', 
    chem_stats_f = f'{af_dir}/af_{{run}}/{af_rna_dir}chemistry_statistics.yaml'
  params:
    ambient_method    = config['ambient']['ambient_method'],
    seurat_quantile   = config['multiplexing']['seurat_quantile']
  output:
    sce_hto_f   = f'{demux_dir}/sce_cells_htos_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds'
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'make_hto_sce_objects', 'memory', attempt, wildcards.run),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'make_hto_sce_objects', 'time', attempt, wildcards.run)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_hto/make_hto_sce_objects_{{run}}_{DATE_STAMP}.benchmark.txt'
  conda:
   '../envs/rlibs.yaml'
  shell: """

  # get translation file from chemistry stats
  WHITELIST_TRANS_F=$(grep "selected_translation_f:" {input.chem_stats_f} | sed 's/selected_translation_f: //')
  
  # save hto sce with demultiplexing info
  Rscript -e "source('scripts/multiplexing.R'); source('scripts/utils.R'); 
    get_one_hto_sce( 
      sel_pool        = '{wildcards.run}', 
      sample_stats_f  = '{input.smpl_stats_f}', 
      amb_yaml_f      = '{input.amb_yaml_f}', 
      hto_mat_f       = '{input.hto_h5_f}', 
      trans_f         = '$WHITELIST_TRANS_F', 
      hto_sce_f       = '{output.sce_hto_f}', 
      ambient_method  = '{params.ambient_method}',
      seurat_quantile =  {params.seurat_quantile}
    )"
  """

