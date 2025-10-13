# snakemake rule for making sce objects with hto information (for multiplexed samples)

if DEMUX_TYPE == "hto":
# save sce object with hto counts and demultiplex
  rule make_hto_sce_objects: 
    input: 
      smpl_stats_f = amb_dir + '/ambient_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
      amb_yaml_f   = amb_dir + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml',
      hto_h5_f     = af_dir + '/af_{run}/hto/af_hto_counts_mat.h5'
    params:
      translation_f = lambda wildcards: parse_alevin_params(CUSTOM_SAMPLE_PARAMS_F, CHEMISTRY, SCPROCESS_DATA_DIR, wildcards.run)[3]
    output:
      sce_hto_f   = demux_dir + '/sce_cells_htos_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_make_hto_sce_objects'] * MB_PER_GB
    conda:
     '../envs/rlibs.yaml'
    shell: """
    # save hto sce with demultiplexing info
    Rscript -e "source('scripts/multiplexing.R'); source('scripts/ambient.R'); 
      get_one_hto_sce( 
        sel_sample      = '{wildcards.run}', 
        sample_var      = '{SAMPLE_VAR}', 
        sample_stats_f  = '{input.smpl_stats_f}', 
        amb_yaml_f      = '{input.amb_yaml_f}', 
        hto_mat_f       = '{input.hto_h5_f}', 
        trans_f         = '{params.translation_f}', 
        hto_sce_f       = '{output.sce_hto_f}', 
        ambient_method  = '{AMBIENT_METHOD}'
      )"
    """
