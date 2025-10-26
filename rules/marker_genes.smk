# snakemake rule for calculating marker genes

def get_conditional_outputs(species):
  if species in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020']:
    return {
      'fgsea_go_bp_f': mkr_dir + '/fgsea_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_go_bp_' + DATE_STAMP + '.txt.gz',
      'fgsea_go_cc_f': mkr_dir + '/fgsea_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_go_cc_' + DATE_STAMP + '.txt.gz',
      'fgsea_go_mf_f': mkr_dir + '/fgsea_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_go_mf_' + DATE_STAMP + '.txt.gz',
      'fgsea_paths_f': mkr_dir + '/fgsea_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_paths_' + DATE_STAMP + '.txt.gz',
      'fgsea_hlmk_f':  mkr_dir + '/fgsea_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_hlmk_' + DATE_STAMP + '.txt.gz'
    }
  else:
    return {}


rule run_marker_genes:
  input:
    sces_yaml_f    = int_dir  + '/sce_clean_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml',
    integration_f  = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    pb_f      = mkr_dir + '/pb_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.rds',
    mkrs_f    = mkr_dir + '/pb_marker_genes_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.txt.gz',
    pb_hvgs_f = mkr_dir + '/pb_hvgs_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.txt.gz',
    **get_conditional_outputs(config['project']['species'])
  params:
    species         = config['project']['species'],
    af_gtf_dt_f     = config['mapping']['af_gtf_dt_f'],
    mkr_gsea_dir    = config['marker_genes']['mkr_gsea_dir'],
    mkr_sel_res     = config['marker_genes']['mkr_sel_res'],
    mkr_min_cl_size = config['marker_genes']['mkr_min_cl_size'],
    mkr_min_cells   = config['marker_genes']['mkr_min_cells'],
    mkr_not_ok_re   = config['marker_genes']['mkr_not_ok_re'],
    mkr_min_cpm_go  = config['marker_genes']['mkr_min_cpm_go'],
    mkr_max_zero_p  = config['marker_genes']['mkr_max_zero_p'],
    mkr_gsea_cut    = config['marker_genes']['mkr_gsea_cut'],
    fgsea_args = lambda wildcards, output: ", ".join([
        f"fgsea_go_bp_f = '{output.get('fgsea_go_bp_f', '')}'",
        f"fgsea_go_cc_f = '{output.get('fgsea_go_cc_f', '')}'",
        f"fgsea_go_mf_f = '{output.get('fgsea_go_mf_f', '')}'",
        f"fgsea_paths_f = '{output.get('fgsea_paths_f', '')}'",
        f"fgsea_hlmk_f = '{output.get('fgsea_hlmk_f', '')}',"
    ])
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_run_marker_genes'] * MB_PER_GB
  conda: '../envs/rlibs.yaml'
  shell:"""
    Rscript -e "source('scripts/utils.R'); source('scripts/marker_genes.R'); calculate_marker_genes(
      integration_f = '{input.integration_f}', 
      sces_yaml_f   = '{input.sces_yaml_f}',
      pb_f          = '{output.pb_f}',
      mkrs_f        = '{output.mkrs_f}',
      pb_hvgs_f     = '{output.pb_hvgs_f}',
      {params.fgsea_args}
      species       = '{params.species}',
      gtf_dt_f      = '{params.af_gtf_dt_f}',
      gsea_dir      = '{params.mkr_gsea_dir}',
      sel_res       = '{params.mkr_sel_res}',
      min_cl_size   =  {params.mkr_min_cl_size},
      min_cells     =  {params.mkr_min_cells},
      not_ok_re     = '{params.mkr_not_ok_re}',
      min_cpm_go    =  {params.mkr_min_cpm_go},
      max_zero_p    =  {params.mkr_max_zero_p},
      gsea_cut      =  {params.mkr_gsea_cut},
      n_cores       =  {threads})"
    """
