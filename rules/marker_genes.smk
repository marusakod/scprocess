rule run_marker_genes:
  input:
    sces_yaml_f    = int_dir  + '/sce_clean_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml',
    integration_f  = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  output:
    pb_f      = mkr_dir + '/pb_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.rds',
    mkrs_f    = mkr_dir + '/pb_marker_genes_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.txt.gz',
    pb_hvgs_f = mkr_dir + '/pb_hvgs_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.txt.gz'
  params:
    af_gtf_dt_f     = config['mapping']['af_gtf_dt_f'],
    mkr_gsea_dir    = config['marker_genes']['mkr_gsea_dir'],
    mkr_sel_res     = config['marker_genes']['mkr_sel_res'],
    mkr_min_cl_size = config['marker_genes']['mkr_min_cl_size'],
    mkr_min_cells   = config['marker_genes']['mkr_min_cells']
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('run_marker_genes', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('run_marker_genes', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_marker_genes/run_marker_genes_' + DATE_STAMP + '.benchmark.txt'
  conda: '../envs/rlibs.yaml'
  shell:"""
    Rscript -e "source('scripts/utils.R'); source('scripts/marker_genes.R'); calculate_marker_genes(
      integration_f = '{input.integration_f}', 
      sces_yaml_f   = '{input.sces_yaml_f}',
      pb_f          = '{output.pb_f}',
      mkrs_f        = '{output.mkrs_f}',
      pb_hvgs_f     = '{output.pb_hvgs_f}',
      gtf_dt_f      = '{params.af_gtf_dt_f}',
      sel_res       = '{params.mkr_sel_res}',
      min_cl_size   =  {params.mkr_min_cl_size},
      min_cells     =  {params.mkr_min_cells},
      n_cores       =  {threads})"
    """

rule run_fgsea:
  input:
    mkrs_f        = mkr_dir + '/pb_marker_genes_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.txt.gz'
  output:
    fgsea_go_bp_f = mkr_dir + '/fgsea_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_go_bp_' + DATE_STAMP + '.txt.gz', 
    fgsea_go_cc_f = mkr_dir + '/fgsea_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_go_cc_' + DATE_STAMP + '.txt.gz',
    fgsea_go_mf_f = mkr_dir + '/fgsea_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_go_mf_' + DATE_STAMP + '.txt.gz'
  params:
    species         = config['project']['species'],
    mkr_gsea_dir    = config['marker_genes']['mkr_gsea_dir'],
    mkr_min_cpm_go  = config['marker_genes']['mkr_min_cpm_go'],
    mkr_max_zero_p  = config['marker_genes']['mkr_max_zero_p'],
    mkr_gsea_cut    = config['marker_genes']['mkr_gsea_cut'], 
    mkr_not_ok_re   = config['marker_genes']['mkr_not_ok_re']
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('run_fgsea', 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('run_fgsea', 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_marker_genes/run_fgsea_' + DATE_STAMP + '.benchmark.txt'
  conda: '../envs/rlibs.yaml'
  shell:"""
    Rscript -e "source('scripts/utils.R'); source('scripts/fgsea.R'); run_fgsea(
      mkrs_f        = '{input.mkrs_f}', 
      fgsea_go_bp_f = '{output.fgsea_go_bp_f}', 
      fgsea_go_cc_f = '{output.fgsea_go_cc_f}', 
      fgsea_go_mf_f = '{output.fgsea_go_mf_f}', 
      species       = '{params.species}', 
      gsea_dir      = '{params.mkr_gsea_dir}', 
      min_cpm_go    = {params.mkr_min_cpm_go}, 
      max_zero_p    = {params.mkr_max_zero_p},
      gsea_cut      = {params.mkr_gsea_cut},
      not_ok_re     = '{params.mkr_not_ok_re}',
      n_cores       =  {threads})"
    """
