# rules to render html files

# # render_html_index --> render from index.Rmd file that is created when workflow R project is started ?
# rule render_html_index:
#   output:
#     rmd_f       = f'{rmd_dir}/index.Rmd',
#     html_f      = f'{docs_dir}/index.html'
#   threads: 1
#   retries: RETRIES 
#   resources:
#     mem_mb      =  lambda wildcards, attempt: attempt * 4096
#   run:
#     # define what we will substitute in
#     print('setting up template')
#     sub_dict    = {
#       'YOUR_NAME':        YOUR_NAME,
#       'AFFILIATION':      AFFILIATION,
#       'SHORT_TAG':        SHORT_TAG
#       }
#     # make and render Rmd file
#     template_f  = 'resources/rmd_templates/marker_genes.Rmd.template'
#     print('rendering template')
#     render_html(PROJ_DIR, template_f, sub_dict, output.rmd_f)


# rule render_html_mapping
rule render_html_mapping:
  input:
    knee_fs   = expand(af_dir + '/af_{run}/' + af_rna_dir + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz', run=RUNS)
  output:
    r_utils_f = f"{code_dir}/utils.R",
    r_map_f   = f"{code_dir}/mapping.R",
    r_amb_f   = f"{code_dir}/ambient.R",
    rmd_f     = f"{rmd_dir}/{SHORT_TAG}_mapping.Rmd",
    html_f    = f"{docs_dir}/{SHORT_TAG}_mapping.html"
  params:
    your_name       = config['project']['your_name'],
    affiliation     = config['project']['affiliation'],
    short_tag       = config['project']['short_tag'],
    date_stamp      = config['project']['date_stamp'],
    proj_dir        = config['project']['proj_dir'],
    ambient_method  = config['ambient']['ambient_method'],
    run_var         = RUN_VAR,
    runs_str        = ','.join(RUNS)
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb    =  lambda wildcards, attempt: attempt * config['resources']['gb_render_htmls'] * MB_PER_GB
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_render_htmls/render_html_mapping_' + DATE_STAMP + '.benchmark.txt'
  conda:
    '../envs/rlibs.yaml'
  shell: """
    # copy R code over
    echo "copying relevant R files over"
    cp scripts/utils.R {output.r_utils_f}
    cp scripts/mapping.R {output.r_map_f}
    cp scripts/ambient.R {output.r_amb_f}

    # define rule and template
    template_f=$(realpath resources/rmd_templates/mapping.Rmd.template)
    rule="mapping"

    # rendering html
    Rscript --vanilla -e "source('scripts/render_htmls.R'); \
      render_html(
        rule_name       = '$rule', 
        proj_dir        = '{params.proj_dir}', 
        temp_f          =  '$template_f', 
        rmd_f           = '{output.rmd_f}', 
        your_name       = '{params.your_name}', 
        affiliation     = '{params.affiliation}', 
        short_tag       = '{params.short_tag}', 
        date_stamp      = '{params.date_stamp}', 
        threads         =  {threads},
        runs_str        = '{params.runs_str}', 
        ambient_method  = '{params.ambient_method}', 
        run_var         = '{params.run_var}', 
        af_dir          = '{af_dir}', 
        af_rna_dir      = '{af_rna_dir}'
      )"
    """

if config['multiplexing']['demux_type'] == "hto":
  # rule render_html_multiplexing
  rule render_html_multiplexing:
    input:
      r_utils_f   = code_dir + '/utils.R', 
      r_amb_f     = code_dir + '/ambient.R',
      hto_knee_fs = expand(af_dir + '/af_{run}/hto/' + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz', run=RUNS), 
      sce_hto_fs  = expand(demux_dir + '/sce_cells_htos_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.rds', run=RUNS)
    output:
      r_demux_f   = f"{code_dir}/multiplexing.R",
      rmd_f       = f"{rmd_dir}/{SHORT_TAG}_demultiplexing.Rmd",
      html_f      = f"{docs_dir}/{SHORT_TAG}_demultiplexing.html"
    params:
      your_name       = config['project']['your_name'],
      affiliation     = config['project']['affiliation'],
      short_tag       = config['project']['short_tag'],
      date_stamp      = config['project']['date_stamp'],
      proj_dir        = config['project']['proj_dir'],
      metadata_f      = config['project']['sample_metadata'],
      ambient_method  = config['ambient']['ambient_method'],
      run_var         = RUN_VAR,
      runs_str        = ','.join(RUNS)
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb      =  lambda wildcards, attempt: attempt * config['resources']['gb_render_htmls'] * MB_PER_GB
    benchmark:
      benchmark_dir + '/' + SHORT_TAG + '_render_htmls/render_html_multiplexing_' + DATE_STAMP + '.benchmark.txt'
    conda:
      '../envs/rlibs.yaml'
    shell: """
      # copy R code over
      echo "copying relevant R files over"
      cp scripts/multiplexing.R {output.r_demux_f}
    
      # make and render Rmd file
      template_f=$(realpath resources/rmd_templates/multiplexing.Rmd.template)
      rule="multiplexing"

      # rendering html
      Rscript --vanilla -e "source('scripts/render_htmls.R'); \
        render_html(
          rule_name       = '$rule', 
          proj_dir        = '{params.proj_dir}', 
          temp_f          =  '$template_f', 
          rmd_f           = '{output.rmd_f}', 
          your_name       = '{params.your_name}', 
          affiliation     = '{params.affiliation}', 
          short_tag       = '{params.short_tag}', 
          date_stamp      = '{params.date_stamp}', 
          threads         =  {threads},
          runs_str        = '{params.runs_str}', 
          metadata_f      = '{params.metadata_f}', 
          ambient_method  = '{params.ambient_method}', 
          run_var         = '{params.run_var}', 
          af_dir          = '{af_dir}', 
          demux_dir       = '{demux_dir}'
        )"
      """

# render_html_ambient
rule render_html_ambient:
  input:
    r_amb_f       = f"{code_dir}/ambient.R", 
    r_utils_f     = f"{code_dir}/utils.R",
    run_stats_f  = amb_dir + '/ambient_run_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    rmd_f         = f"{rmd_dir}/{SHORT_TAG}_ambient.Rmd",
    html_f        = f"{docs_dir}/{SHORT_TAG}_ambient.html"
  params:
    your_name         = config['project']['your_name'],
    affiliation       = config['project']['affiliation'],
    short_tag         = config['project']['short_tag'],
    date_stamp        = config['project']['date_stamp'],
    proj_dir          = config['project']['proj_dir'],
    metadata_f        = config['project']['sample_metadata'],
    run_var           = RUN_VAR,
    runs_str          = ','.join(RUNS),
    ambient_method    = config['ambient']['ambient_method'],
    cb_max_prop_kept  = config['ambient']['cb_max_prop_kept']
  threads: 4
  retries: config['resources']['retries'] 
  conda:
    '../envs/rlibs.yaml'
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * config['resources']['gb_render_htmls'] * MB_PER_GB
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_render_htmls/render_html_ambient_' + DATE_STAMP + '.benchmark.txt'
  shell: """
    # define rule and template
    template_f=$(realpath resources/rmd_templates/ambient.Rmd.template)
    rule="ambient"

    # rendering html
    Rscript --vanilla -e "source('scripts/render_htmls.R'); \
      render_html(
        rule_name         = '$rule', 
        temp_f            = '$template_f', 
        rmd_f             = '{output.rmd_f}', 
        your_name         = '{params.your_name}', 
        affiliation       = '{params.affiliation}', 
        proj_dir          = '{params.proj_dir}', 
        short_tag         = '{params.short_tag}', 
        date_stamp        = '{params.date_stamp}', 
        threads           =  {threads}, 
        run_stats_f       = '{input.run_stats_f}', 
        run_var           = '{params.run_var}',
        runs_str          = '{params.runs_str}',
        ambient_method    = '{params.ambient_method}', 
        cb_prop_max_kept  =  {params.cb_max_prop_kept}
      )"    
    """

# render_html_qc
rule render_html_qc:
  input:
    r_utils_f   = f"{code_dir}/utils.R",
    qc_dt_f     = qc_dir  + '/qc_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz',
    cuts_f      = qc_dir  + '/qc_thresholds_by_sample_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  params:
    metadata_f          = config['project']['sample_metadata'],
    your_name           = config['project']['your_name'],
    affiliation         = config['project']['affiliation'],
    short_tag           = config['project']['short_tag'],
    date_stamp          = config['project']['date_stamp'],
    proj_dir            = config['project']['proj_dir'],
    min_cells           = config['qc']['qc_min_cells'],
    qc_hard_min_counts  = config['qc']['qc_hard_min_counts'],
    qc_hard_min_feats   = config['qc']['qc_hard_min_feats'],
    qc_hard_max_mito    = config['qc']['qc_hard_max_mito']
  output:
    r_qc_f      = f"{code_dir}/qc.R",
    rmd_f       = f"{rmd_dir}/{SHORT_TAG}_qc.Rmd",
    html_f      = f"{docs_dir}/{SHORT_TAG}_qc.html"
  threads: 1
  retries: config['resources']['retries'] 
  conda:
    '../envs/rlibs.yaml'
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * config['resources']['gb_render_htmls'] * MB_PER_GB
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_render_htmls/render_html_qc_' + DATE_STAMP + '.benchmark.txt'
  shell: """
    # copy R code over
    echo "copying relevant R files over"
    cp scripts/SampleQC.R {output.r_qc_f}

    # define rule and template
    template_f=$(realpath resources/rmd_templates/SampleQC.Rmd.template)
    rule="qc"

    # rendering html
    Rscript --vanilla -e "source('scripts/render_htmls.R'); \
    render_html(
      rule_name           = '$rule',
      proj_dir            = '{params.proj_dir}', 
      temp_f              = '$template_f',
      rmd_f               = '{output.rmd_f}',
      your_name           = '{params.your_name}', 
      affiliation         = '{params.affiliation}', 
      short_tag           = '{params.short_tag}', 
      date_stamp          = '{params.date_stamp}', 
      threads             =  {threads}, 
      metadata_f          = '{params.metadata_f}',
      qc_dt_f             = '{input.qc_dt_f}',
      cuts_f              = '{input.cuts_f}',
      min_cells           =  {params.min_cells},
      qc_hard_min_counts  =  {params.qc_hard_min_counts},
      qc_hard_min_feats   =  {params.qc_hard_min_feats},
      qc_hard_max_mito    =  {params.qc_hard_max_mito}
      )"
    """

# render_html_hvgs
rule render_html_hvgs:
  input:
    r_utils_f   = f"{code_dir}/utils.R",
    hvgs_f      = hvg_dir   + '/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    empty_gs_f  = empty_dir + '/edger_empty_genes_' + FULL_TAG + '_all_' + DATE_STAMP + '.txt.gz',
    pb_empty_f  = pb_dir  + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  output:
    r_hvgs_f    = f"{code_dir}/hvgs.R",
    rmd_f       = f"{rmd_dir}/{SHORT_TAG}_hvgs.Rmd",
    html_f      = f"{docs_dir}/{SHORT_TAG}_hvgs.html"
  params:
    your_name   = config['project']['your_name'],
    affiliation = config['project']['affiliation'],
    short_tag   = config['project']['short_tag'],
    date_stamp  = config['project']['date_stamp'],
    proj_dir    = config['project']['proj_dir']
  threads: 1
  retries: config['resources']['retries'] 
  conda:
    '../envs/rlibs.yaml'
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * config['resources']['gb_render_htmls'] * MB_PER_GB
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_render_htmls/render_html_hvgs_' + DATE_STAMP + '.benchmark.txt'
  shell: """
    # copy R code over
    echo "copying relevant R files over"
    cp scripts/hvgs.R {output.r_hvgs_f}

    # define rule and template
    template_f=$(realpath resources/rmd_templates/hvgs.Rmd.template)
    rule="hvg"

    # rendering html
    Rscript --vanilla -e "source('scripts/render_htmls.R'); \
    render_html(
      rule_name   = '$rule',
      temp_f      = '$template_f',
      rmd_f       = '{output.rmd_f}',
      your_name   = '{params.your_name}',
      affiliation = '{params.affiliation}',
      proj_dir    = '{params.proj_dir}',
      short_tag   = '{params.short_tag}',
      date_stamp  = '{params.date_stamp}',
      threads     =  {threads},
      hvgs_f      = '{input.hvgs_f}',
      empty_gs_f  = '{input.empty_gs_f}',
      pb_empty_f  = '{input.pb_empty_f}'
    )"   
    """

# render_html_integration
rule render_html_integration:
  input:
    r_utils_f     = f"{code_dir}/utils.R",
    r_amb_f       = f"{code_dir}/ambient.R",
    qc_dt_f       = qc_dir  + '/qc_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz',
    integration_f = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    r_int_f     = f"{code_dir}/integration.R",
    rmd_f       = f"{rmd_dir}/{SHORT_TAG}_integration.Rmd",
    html_f      = f"{docs_dir}/{SHORT_TAG}_integration.html"
  params:
    your_name       = config['project']['your_name'],
    affiliation     = config['project']['affiliation'],
    short_tag       = config['project']['short_tag'],
    date_stamp      = config['project']['date_stamp'],
    proj_dir        = config['project']['proj_dir'],
    int_res_ls_str  = ','.join(map(str, config['integration']['int_res_ls'])),
    int_dbl_cl_prop = config['integration']['int_dbl_cl_prop'],
    int_embedding   = config['integration']['int_embedding'],
    demux_type      = config['multiplexing']['demux_type']
  threads: 1
  retries: config['resources']['retries']
  conda:
    '../envs/rlibs.yaml'
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * config['resources']['gb_render_htmls'] * MB_PER_GB
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_render_htmls/render_html_integration_' + DATE_STAMP + '.benchmark.txt'
  shell: """
    # copy R code over
    echo "copying relevant R files over"
    cp scripts/integration.R {output.r_int_f}

    # define rule and template
    template_f=$(realpath resources/rmd_templates/integration.Rmd.template)
    rule="integration"
   
    # rendering html
    Rscript --vanilla -e "source('scripts/render_htmls.R'); \
    render_html(
      rule_name       = '$rule',
      temp_f          = '$template_f',
      rmd_f           = '{output.rmd_f}',
      your_name       = '{params.your_name}',
      affiliation     = '{params.affiliation}',
      proj_dir        = '{params.proj_dir}',
      short_tag       = '{params.short_tag}',
      date_stamp      = '{params.date_stamp}',
      qc_dt_f         = '{input.qc_dt_f}',
      threads         =  {threads},
      integration_f   = '{input.integration_f}',
      int_res_ls_str  = '{params.int_res_ls_str}',
      int_dbl_cl_prop =  {params.int_dbl_cl_prop},
      int_embedding   = '{params.int_embedding}',
      demux_type      = '{params.demux_type}'
    )"
    """

# render_html_marker_genes
rule render_html_marker_genes:
  input:
    r_utils_f     = f"{code_dir}/utils.R",
    r_int_f       = f'{code_dir}/integration.R',
    pb_f          = mkr_dir + '/pb_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.rds',
    mkrs_f        = mkr_dir + '/pb_marker_genes_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.txt.gz',
    integration_f = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    hvgs_f        = mkr_dir + '/pb_hvgs_' + FULL_TAG + f'_{config['marker_genes']['mkr_sel_res']}_' + DATE_STAMP + '.txt.gz',
    ambient_f     = empty_dir + '/edger_empty_genes_' + FULL_TAG + '_all_' + DATE_STAMP + '.txt.gz',
    **get_conditional_outputs(config['project']['species'])
  output:
    r_mkr_f       = f"{code_dir}/marker_genes.R",
    rmd_f         = f'{rmd_dir}/{SHORT_TAG}_marker_genes_{config['marker_genes']['mkr_sel_res']}.Rmd',
    html_f        = f'{docs_dir}/{SHORT_TAG}_marker_genes_{config['marker_genes']['mkr_sel_res']}.html'
  threads: 8
  retries: config['resources']['retries']
  params:
    your_name         = config['project']['your_name'],
    affiliation       = config['project']['affiliation'],
    short_tag         = config['project']['short_tag'],
    date_stamp        = config['project']['date_stamp'],
    proj_dir          = config['project']['proj_dir'],
    metadata_f        = config['project']['sample_metadata'],
    species           = config['project']['species'],
    meta_vars         = ','.join(config['project']['metadata_vars']),
    af_gtf_dt_f       = config['mapping']['af_gtf_dt_f'],
    custom_mkr_names  = config['marker_genes']['custom_mkr_names'],
    custom_mkr_paths  = config['marker_genes']['custom_mkr_paths'],
    mkr_sel_res       = config['marker_genes']['mkr_sel_res'],
    mkr_not_ok_re     = config['marker_genes']['mkr_not_ok_re'],
    mkr_min_cpm_mkr   = config['marker_genes']['mkr_min_cpm_mkr'],
    mkr_min_cells     = config['marker_genes']['mkr_min_cells'],
    mkr_gsea_cut      = config['marker_genes']['mkr_gsea_cut'],
    fgsea_args = lambda wildcards, input: " ".join(
      [
        f"fgsea_go_bp_f = '{input.get('fgsea_go_bp_f', '')}',",
        f"fgsea_go_cc_f = '{input.get('fgsea_go_cc_f', '')}',",
        f"fgsea_go_mf_f = '{input.get('fgsea_go_mf_f', '')}',",
        f"fgsea_paths_f = '{input.get('fgsea_paths_f', '')}',",
        f"fgsea_hlmk_f  = '{input.get('fgsea_hlmk_f', '')}',"
      ]
    ).strip()
  conda: '../envs/rlibs.yaml'
  resources:
    mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_render_htmls'] * MB_PER_GB
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_render_htmls/render_html_marker_genes_' + DATE_STAMP + '.benchmark.txt'
  shell: """
    # copy R code over
    echo "copying relevant R files over"
    cp scripts/marker_genes.R {output.r_mkr_f}

    # define rule and template
    template_f=$(realpath resources/rmd_templates/marker_genes.Rmd.template)
    rule="markers"

    # rendering html
    Rscript --vanilla -e "source('scripts/render_htmls.R'); \
    render_html(
      rule_name         = '$rule',
      temp_f            = '$template_f',
      rmd_f             = '{output.rmd_f}',
      your_name         = '{params.your_name}',
      affiliation       = '{params.affiliation}',
      proj_dir          = '{params.proj_dir}',
      short_tag         = '{params.short_tag}',
      date_stamp        = '{params.date_stamp}',
      threads           =  {threads},
      metadata_f        = '{params.metadata_f}',
      meta_vars_ls      = '{params.meta_vars}',
      gtf_dt_f          = '{params.af_gtf_dt_f}',
      ambient_f         = '{input.ambient_f}',
      integration_f     = '{input.integration_f}',
      pb_f              = '{input.pb_f}',
      hvgs_f            = '{input.hvgs_f}',
      mkrs_f            = '{input.mkrs_f}',
      custom_mkr_names  = '{params.custom_mkr_names}',
      custom_mkr_paths  = '{params.custom_mkr_paths}',
      mkr_sel_res       =  {params.mkr_sel_res},
      mkr_not_ok_re     = '{params.mkr_not_ok_re}',
      mkr_min_cpm_mkr   =  {params.mkr_min_cpm_mkr},
      mkr_min_cells     =  {params.mkr_min_cells},
      mkr_gsea_cut      =  {params.mkr_gsea_cut},
      {params.fgsea_args}
      species           = '{params.species}'
    )"
    """

if "label_celltypes" in config:
  # render_html_label_celltypes
  rule render_html_label_celltypes:
    input:
      r_utils_f   = f"{code_dir}/utils.R",
      r_mkr_f     = f'{code_dir}/marker_genes.R',
      guesses_f   = f'{lbl_dir}/cell_annotations_{FULL_TAG}_{DATE_STAMP}.txt.gz'
    output:
      r_lbl_f     = f'{code_dir}/label_celltypes.R',
      rmd_f       = f'{rmd_dir}/{SHORT_TAG}_label_celltypes.Rmd',
      html_f      = f'{docs_dir}/{SHORT_TAG}_label_celltypes.html'
    params:
      your_name       = config['project']['your_name'],
      affiliation     = config['project']['affiliation'],
      short_tag       = config['project']['short_tag'],
      date_stamp      = config['project']['date_stamp'],
      proj_dir        = config['project']['proj_dir'],
      lbl_xgb_f       = config['label_celltypes']['lbl_xgb_f'], 
      lbl_xgb_cls_f   = config['label_celltypes']['lbl_xgb_cls_f'], 
      lbl_tissue      = config['label_celltypes']['lbl_tissue'], 
      lbl_sel_res_cl  = config['label_celltypes']['lbl_sel_res_cl'], 
      lbl_min_pred    = config['label_celltypes']['lbl_min_pred'], 
      lbl_min_cl_prop = config['label_celltypes']['lbl_min_cl_prop'], 
      lbl_min_cl_size = config['label_celltypes']['lbl_min_cl_size']
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb      =  lambda wildcards, attempt: attempt * config['resources']['gb_render_htmls'] * MB_PER_GB
    conda:
      '../envs/rlibs.yaml'
    benchmark:
      benchmark_dir + '/' + SHORT_TAG + '_render_htmls/render_html_label_celltypes_' + DATE_STAMP + '.benchmark.txt'
    shell: """
      # copy R code over
      echo "copying relevant R files over"
      cp scripts/label_celltypes.R {output.r_lbl_f}

      template_f=$(realpath resources/rmd_templates/label_celltypes.Rmd.template)
      rule="cell_labels"

      Rscript --vanilla -e "source('scripts/render_htmls.R'); \
      render_html(
        rule_name       = '$rule', 
        temp_f          =  '$template_f', 
        rmd_f           = '{output.rmd_f}', 
        your_name       = '{params.your_name}',
        affiliation     = '{params.affiliation}',
        proj_dir        = '{params.proj_dir}',
        short_tag       = '{params.short_tag}',
        date_stamp      = '{params.date_stamp}',
        threads         =  {threads}, 
        guesses_f       = '{input.guesses_f}', 
        lbl_xgb_f       = '{params.lbl_xgb_f}', 
        lbl_xgb_cls_f   = '{params.lbl_xgb_cls_f}', 
        lbl_tissue      = '{params.lbl_tissue}', 
        lbl_sel_res_cl  = '{params.lbl_sel_res_cl}', 
        lbl_min_pred    =  {params.lbl_min_pred}, 
        lbl_min_cl_prop =  {params.lbl_min_cl_prop}, 
        lbl_min_cl_size =  {params.lbl_min_cl_size}
      )"
      """
