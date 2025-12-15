
def get_conditional_fgsea_files(species, do_gsea):
  if (species in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020']) & do_gsea:
    return {
      'fgsea_go_bp_f': f'{mkr_dir}/fgsea_{FULL_TAG}_{ config['marker_genes']['mkr_sel_res'] }_go_bp_{DATE_STAMP}.csv.gz',
      'fgsea_go_cc_f': f'{mkr_dir}/fgsea_{FULL_TAG}_{ config['marker_genes']['mkr_sel_res'] }_go_cc_{DATE_STAMP}.csv.gz',
      'fgsea_go_mf_f': f'{mkr_dir}/fgsea_{FULL_TAG}_{ config['marker_genes']['mkr_sel_res'] }_go_mf_{DATE_STAMP}.csv.gz'
    }
  else:
    return {}



rule render_html_index:
  input:
    html_reports  = glob.glob(f'{docs_dir}/{SHORT_TAG}*.html')
  output:
    rmd_f  = f"{rmd_dir}/index.Rmd",
    html_f = f"{docs_dir}/index.html"
  threads: 1
  resources: 
    mem_mb  = 2* MB_PER_GB, 
    runtime = 2
  params: 
    your_name       = config['project']['your_name'],
    affiliation     = config['project']['affiliation'],
    proj_dir        = config['project']['proj_dir'],
    short_tag       = SHORT_TAG, 
    full_tag        = config['project']['full_tag'],
    date_stamp      = config['project']['date_stamp'], 
    mkr_sel_res     = config['marker_genes']['mkr_sel_res']
  conda:
    '../envs/rlibs.yaml'
  shell:"""
    # Template for the RMarkdown file
    template_f=$(realpath resources/rmd_templates/index.Rmd.template)
    rule="index"
    rmd_f={rmd_dir}/index.Rmd

    Rscript --vanilla -e "source('scripts/render_htmls.R'); \
    render_html(
      rule_name       = '$rule', 
      proj_dir        = '{params.proj_dir}', 
      your_name       = '{params.your_name}', 
      affiliation     = '{params.affiliation}',
      docs_dir        = '{docs_dir}', 
      short_tag       = '{params.short_tag}', 
      full_tag        = '{params.full_tag}',  
      date_stamp      = '{params.date_stamp}',
      mkr_sel_res     = '{params.mkr_sel_res}', 
      temp_f          = '$template_f',
      rmd_f           = '{output.rmd_f}'
    )"
    """


# rule render_html_mapping
rule render_html_mapping:
  input:
    knee_fs         = expand(f'{af_dir}/af_{{run}}/{af_rna_dir}knee_plot_data_{{run}}_{DATE_STAMP}.csv.gz', run=RUNS)
  output:
    r_utils_f       = f"{code_dir}/utils.R",
    r_map_f         = f"{code_dir}/mapping.R",
    rmd_f           = f"{rmd_dir}/{SHORT_TAG}_mapping.Rmd",
    html_f          = f"{docs_dir}/{SHORT_TAG}_mapping.html"
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
    mem_mb    =  lambda wildcards, attempt, input: attempt * get_resources('render_html_mapping', rules, 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS), 
    runtime   =  lambda wildcards, input: get_resources('render_html_mapping', rules, 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_render_htmls/render_html_mapping_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/rlibs.yaml'
  shell: """
    # copy R code over
    echo "copying relevant R files over"
    cp scripts/utils.R {output.r_utils_f}
    cp scripts/mapping.R {output.r_map_f}

    # define rule and template
    template_f=$(realpath resources/rmd_templates/mapping.Rmd.template)
    rule="mapping"

    # rendering html
    Rscript --vanilla -e "source('scripts/render_htmls.R'); \
      render_html(
        rule_name       = '$rule', 
        proj_dir        = '{params.proj_dir}', 
        temp_f          = '$template_f', 
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
      r_utils_f   = f'{code_dir}/utils.R', 
      hto_knee_fs = expand(f'{af_dir}/af_{{run}}/hto/knee_plot_data_{{run}}_{DATE_STAMP}.csv.gz', run=RUNS), 
      sce_hto_fs  = expand(f'{demux_dir}/sce_cells_htos_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds', run=RUNS)
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
      batch_var       = BATCH_VAR,
      runs_str        = ','.join(RUNS),
      r_map_f         = f'{code_dir}/mapping.R'
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('render_html_multiplexing', rules, 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
      runtime = lambda wildcards, input: get_resources('render_html_multiplexing', rules, 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
    benchmark:
      f'{benchmark_dir}/{SHORT_TAG}_render_htmls/render_html_multiplexing_{DATE_STAMP}.benchmark.txt'
    conda:
      '../envs/rlibs.yaml'
    shell: """
      # copy R code over
      echo "copying relevant R files over"
      cp scripts/mapping.R {params.r_map_f}
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
          batch_var       = '{params.batch_var}',
          af_dir          = '{af_dir}', 
          demux_dir       = '{demux_dir}'
        )"
      """

# render_html_ambient
rule render_html_ambient:
  input: 
    r_utils_f   = f"{code_dir}/utils.R",
    run_stats_f = f'{amb_dir}/ambient_run_statistics_{FULL_TAG}_{DATE_STAMP}.csv'
  output:
    r_amb_f     = f"{code_dir}/ambient.R",
    rmd_f       = f"{rmd_dir}/{SHORT_TAG}_ambient.Rmd",
    html_f      = f"{docs_dir}/{SHORT_TAG}_ambient.html"
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
    '../envs/ambientr.yaml'
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('render_html_ambient', rules, 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('render_html_ambient', rules, 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_render_htmls/render_html_ambient_{DATE_STAMP}.benchmark.txt'
  shell: """
    # copy R code over
    echo "copying relevant R files over"
    cp scripts/ambient.R {output.r_amb_f}
    
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
    qc_dt_f     = f'{qc_dir}/qc_all_samples_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    cuts_f      = f'{qc_dir}/qc_thresholds_by_{BATCH_VAR}_{FULL_TAG}_{DATE_STAMP}.csv'
  params:
    metadata_f          = config['project']['sample_metadata'],
    your_name           = config['project']['your_name'],
    affiliation         = config['project']['affiliation'],
    short_tag           = config['project']['short_tag'],
    date_stamp          = config['project']['date_stamp'],
    proj_dir            = config['project']['proj_dir'],
    min_cells           = config['qc']['qc_min_cells'],
    batch_var           = BATCH_VAR,
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
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('render_html_qc', rules, 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('render_html_qc', rules, 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_render_htmls/render_html_qc_{DATE_STAMP}.benchmark.txt'
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
      batch_var           = '{params.batch_var}',
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
    hvgs_f      = f'{hvg_dir}/hvg_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    empty_gs_f  = f'{empty_dir}/edger_empty_genes_all_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    pb_empty_f  = f'{pb_dir}/pb_empties_{FULL_TAG}_{DATE_STAMP}.rds'
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
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('render_html_hvgs', rules, 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('render_html_hvgs', rules, 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_render_htmls/render_html_hvgs_{DATE_STAMP}.benchmark.txt'
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
    qc_dt_f       = f'{qc_dir}/qc_all_samples_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    integration_f = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
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
    demux_type      = config['multiplexing']['demux_type'],
    batch_var       = BATCH_VAR
  threads: 1
  retries: config['resources']['retries']
  conda:
    '../envs/rlibs.yaml'
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('render_html_integration', rules, 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('render_html_integration', rules, 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_render_htmls/render_html_integration_{DATE_STAMP}.benchmark.txt'
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
      demux_type      = '{params.demux_type}',
      batch_var       = '{params.batch_var}'
    )"
    """

# render_html_marker_genes
rule render_html_marker_genes:
  input:
    r_utils_f     = f"{code_dir}/utils.R",
    r_int_f       = f'{code_dir}/integration.R',
    pb_f          = f'{mkr_dir}/pb_{FULL_TAG}_{ config['marker_genes']['mkr_sel_res'] }_{DATE_STAMP}.rds',
    mkrs_f        = f'{mkr_dir}/pb_marker_genes_{FULL_TAG}_{ config['marker_genes']['mkr_sel_res'] }_{DATE_STAMP}.csv.gz',
    integration_f = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    hvgs_f        = f'{mkr_dir}/pb_hvgs_{FULL_TAG}_{ config['marker_genes']['mkr_sel_res'] }_{DATE_STAMP}.csv.gz',
    empty_gs_f    = f'{empty_dir}/edger_empty_genes_all_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    **get_conditional_fgsea_files(config['project']['species'], config['marker_genes']['mkr_do_gsea'])
  output:
    r_mkr_f       = f"{code_dir}/marker_genes.R",
    r_fgsea_f     = f"{code_dir}/fgsea.R",
    rmd_f         = f'{rmd_dir}/{SHORT_TAG}_marker_genes_{ config['marker_genes']['mkr_sel_res'] }.Rmd',
    html_f        = f'{docs_dir}/{SHORT_TAG}_marker_genes_{ config['marker_genes']['mkr_sel_res'] }.html'
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
    mkr_do_gsea       = config['marker_genes']['mkr_do_gsea'], 
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
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('render_html_marker_genes', rules, 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
    runtime = lambda wildcards, input: get_resources('render_html_marker_genes', rules,'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_render_htmls/render_html_marker_genes_{DATE_STAMP}.benchmark.txt'
  shell: """
    # copy R code over
    echo "copying relevant R files over"
    cp scripts/marker_genes.R {output.r_mkr_f}
    cp scripts/fgsea.R {output.r_fgsea_f}

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
      ambient_f         = '{input.empty_gs_f}',
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
      species           = '{params.species}', 
      do_gsea           = '{params.mkr_do_gsea}'
    )"
    """

if "label_celltypes" in config:
  # render_html_label_celltypes
  rule render_html_label_celltypes:
    input:
      r_utils_f   = f"{code_dir}/utils.R",
      r_int_f     = f'{code_dir}/integration.R',
      int_f       = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
      guess_f_ls  = expand(f'{lbl_dir}/labels_{{labeller}}_model_{{model}}_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
        zip, 
          labeller  = [ entry['labeller'] for entry in LABELLER_PARAMS],
          model     = [ entry['model']    for entry in LABELLER_PARAMS]
        )
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
      labeller_ls     = [ entry['labeller']     for entry in LABELLER_PARAMS],
      model_ls        = [ entry['model']        for entry in LABELLER_PARAMS],
      hi_res_cl_ls    = [ entry['hi_res_cl']    for entry in LABELLER_PARAMS], 
      min_cl_prop_ls  = [ entry['min_cl_prop']  for entry in LABELLER_PARAMS], 
      min_cl_size_ls  = [ entry['min_cl_size']  for entry in LABELLER_PARAMS]
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('render_html_label_celltypes', rules, 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS),
      runtime = lambda wildcards, input: get_resources('render_html_label_celltypes', rules, 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS)
    conda:
      '../envs/rlibs.yaml'
    benchmark:
      f'{benchmark_dir}/{SHORT_TAG}_render_htmls/render_html_label_celltypes_{DATE_STAMP}.benchmark.txt'
    shell: """
      # copy R code over
      echo "copying relevant R files over"
      cp scripts/label_celltypes.R {output.r_lbl_f}

      template_f=$(realpath resources/rmd_templates/label_celltypes.Rmd.template)
      rule="label_celltypes"

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
        threads         =  {threads},
        int_f           = '{input.int_f}',
        guess_f_ls      = '{input.guess_f_ls}',
        labeller_ls     = '{params.labeller_ls}',
        model_ls        = '{params.model_ls}',
        hi_res_cl_ls    = '{params.hi_res_cl_ls}',
        min_cl_prop_ls  = '{params.min_cl_prop_ls}',
        min_cl_size_ls  = '{params.min_cl_size_ls}'
      )"
      """
