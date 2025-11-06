
# build index for HTO mapping
rule build_hto_index:
  input:
    feature_ref_f = config["multiplexing"]["feature_ref"]
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb = 4 * MB_PER_GB
  output: 
    hto_f       = af_dir + '/hto.tsv',
    t2g_f       = af_dir + '/t2g_hto.tsv',
    idx_log_f   = af_dir + '/hto_index/ref_indexing.log',
    hto_idx_dir = directory(af_dir + '/hto_index')
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
    hto_idx_dir   = af_dir + '/hto_index'
  output:
    fry_dir       = directory(af_dir + '/af_{run}/hto/af_quant/'),
    rad_f         = temp(af_dir + '/af_{run}/hto/af_map/map.rad'),
    mtx_f         = af_dir + '/af_{run}/hto/af_quant/alevin/quants_mat.mtx',
    cols_f        = af_dir + '/af_{run}/hto/af_quant/alevin/quants_mat_cols.txt',
    rows_f        = af_dir + '/af_{run}/hto/af_quant/alevin/quants_mat_rows.txt'
  params:
    demux_type    = config['multiplexing']['demux_type'],
    af_home_dir   = config['mapping']['alevin_fry_home'],
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
    mem_mb        = lambda wildcards, attempt: attempt * config['resources']['gb_run_mapping'] * MB_PER_GB
  shell:"""
    # check whether doing arvados
    ARV_REGEX="^arkau-[0-9a-z]{{5}}-[0-9a-z]{{15}}$"
    if [[ "{params.where}" =~ $ARV_REGEX ]]; then
      ml arvados
      arv-env arkau
    fi
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
      --tenx_chemistry  {params.af_chemistry} \
      --exp_ori         {params.expected_ori} \
      --whitelist_f     {params.whitelist_f}
    """


rule save_alevin_hto_to_h5:
  input: 
    fry_dir     = af_dir + '/af_{run}/hto/af_quant/'
  output: 
    h5_f        = af_dir + '/af_{run}/hto/af_hto_counts_mat.h5',
    knee_data_f = af_dir + '/af_{run}/hto/knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz'
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_save_alevin_to_h5'] * MB_PER_GB
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
    smpl_stats_f = amb_dir + '/ambient_run_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    amb_yaml_f   = amb_dir + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml',
    hto_h5_f     = af_dir + '/af_{run}/hto/af_hto_counts_mat.h5'
  params:
    whitelist_trans_f = lambda wildcards: RUN_PARAMS[wildcards.run]["multiplexing"]["whitelist_trans_f"],
    ambient_method    = config['ambient']['ambient_method'],
    seurat_quantile   = config['multiplexing']['seurat_quantile']
  output:
    sce_hto_f   = demux_dir + '/sce_cells_htos_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_make_hto_sce_objects'] * MB_PER_GB
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_make_hto_sce/make_hto_sce_objects_{run}_' + DATE_STAMP + '.benchmark.txt'
  conda:
   '../envs/rlibs.yaml'
  shell: """
  # save hto sce with demultiplexing info
  Rscript -e "source('scripts/multiplexing.R'); source('scripts/ambient.R'); 
    get_one_hto_sce( 
      sel_pool        = '{wildcards.run}', 
      sample_stats_f  = '{input.smpl_stats_f}', 
      amb_yaml_f      = '{input.amb_yaml_f}', 
      hto_mat_f       = '{input.hto_h5_f}', 
      trans_f         = '{params.whitelist_trans_f}', 
      hto_sce_f       = '{output.sce_hto_f}', 
      ambient_method  = '{params.ambient_method}',
      seurat_quantile =  {params.seurat_quantile}
    )"
  """

