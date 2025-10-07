# rules to render html files

localrules: copy_r_code

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
    knee_fs   = expand(af_dir + '/af_{run}/' + af_rna_dir + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz', run=runs)
  output:
    r_utils_f = f"{code_dir}/utils.R",
    r_map_f   = f"{code_dir}/mapping.R",
    r_amb_f   = f"{code_dir}/ambient.R",
    rmd_f     = f"{rmd_dir}/{SHORT_TAG}_mapping.Rmd",
    html_f    = f"{docs_dir}/{SHORT_TAG}_mapping.html"
  threads: 4
  retries: RETRIES
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * MB_RENDER_HTMLS
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
    rule="af"

    # rendering html
    Rscript --vanilla -e "source('scripts/render_htmls.R'); \
      render_html(
        rule_name      = '$rule', 
        proj_dir       = '{PROJ_DIR}', 
        temp_f         =  '$template_f', 
        rmd_f          = '{output.rmd_f}', 
        threads        = {threads}, 
        YOUR_NAME      = '{YOUR_NAME}', 
        AFFILIATION    = '{AFFILIATION}', 
        SHORT_TAG      = '{SHORT_TAG}', 
        DATE_STAMP     = '{DATE_STAMP}', 
        RUNS_STR       = '{RUNS_STR}', 
        PROJ_DIR       = '{PROJ_DIR}', 
        AMBIENT_METHOD = '{AMBIENT_METHOD}', 
        SAMPLE_VAR     = '{SAMPLE_VAR}', 
        af_dir         = '{af_dir}', 
        af_rna_dir     = '{af_rna_dir}'
      )"
    """

if DEMUX_TYPE == 'af':
  # rule render_html_multiplexing
  rule render_html_multiplexing:
    input:
      r_utils_f   = code_dir + '/utils.R', 
      r_amb_f     = code_dir + '/ambient.R',
      hto_knee_fs = expand(af_dir + '/af_{run}/hto/' + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz', run=runs), 
      sce_hto_fs  = expand(demux_dir + '/sce_cells_htos_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.rds', run = runs)
    output:
      r_demux_f   = f"{code_dir}/multiplexing.R",
      rmd_f       = f"{rmd_dir}/{SHORT_TAG}_demultiplexing.Rmd",
      html_f      = f"{docs_dir}/{SHORT_TAG}_demultiplexing.html"
    threads: 1
    retries: RETRIES
    resources:
      mem_mb      =  lambda wildcards, attempt: attempt * MB_RENDER_HTMLS
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
          rule_name      = '$rule', 
          proj_dir       = '{PROJ_DIR}', 
          temp_f         =  '$template_f', 
          rmd_f          = '{output.rmd_f}', 
          YOUR_NAME      = '{YOUR_NAME}', 
          AFFILIATION    = '{AFFILIATION}', 
          PROJ_DIR       = '{PROJ_DIR}', 
          SHORT_TAG      = '{SHORT_TAG}', 
          DATE_STAMP     = '{DATE_STAMP}', 
          RUNS_STR       = '{RUNS_STR}', 
          METADATA_F     = '{METADATA_F}', 
          AMBIENT_METHOD = '{AMBIENT_METHOD}', 
          SAMPLE_VAR     = '{SAMPLE_VAR}', 
          af_dir         = '{af_dir}', 
          demux_dir      = '{demux_dir}')"

      """

# render_html_ambient
rule render_html_ambient:
  input:
    r_amb_f       = f"{code_dir}/ambient.R", 
    r_utils_f     = f"{code_dir}/utils.R",
    smpl_stats_f  = amb_dir + '/ambient_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    rmd_f         = f"{rmd_dir}/{SHORT_TAG}_ambient.Rmd",
    html_f        = f"{docs_dir}/{SHORT_TAG}_ambient.html"
  threads: 4
  retries: RETRIES 
  conda:
    '../envs/rlibs.yaml'
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * MB_RENDER_HTMLS
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_render_htmls/render_html_ambient_' + DATE_STAMP + '.benchmark.txt'
  shell: 
    """
    # define rule and template
    template_f=$(realpath resources/rmd_templates/ambient.Rmd.template)
    rule="ambient"

    # rendering html
    Rscript --vanilla -e "source('scripts/render_htmls.R'); \
    render_html(
      rule_name     = '$rule', 
      proj_dir      = '{PROJ_DIR}', 
      temp_f        = '$template_f', 
      rmd_f         = '{output.rmd_f}', 
      YOUR_NAME     = '{YOUR_NAME}', 
      AFFILIATION   = '{AFFILIATION}', 
      PROJ_DIR      = '{PROJ_DIR}', 
      SHORT_TAG     = '{SHORT_TAG}', 
      DATE_STAMP    = '{DATE_STAMP}', 
      threads       = {threads}, 
      smpl_stats_f  = '{input.smpl_stats_f}', 
      SAMPLE_VAR    = '{SAMPLE_VAR}',
      RUNS_STR        = '{RUNS_STR}',
      AMBIENT_METHOD  = '{AMBIENT_METHOD}', 
      CELLBENDER_PROP_MAX_KEPT = {CELLBENDER_PROP_MAX_KEPT}
      )"    
    """

# render_html_qc
rule render_html_qc:
  input:
    r_utils_f   = f"{code_dir}/utils.R",
    qc_dt_f     = qc_dir  + '/qc_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    r_qc_f      = f"{code_dir}/qc.R",
    rmd_f       = f"{rmd_dir}/{SHORT_TAG}_qc.Rmd",
    html_f      = f"{docs_dir}/{SHORT_TAG}_qc.html"
  threads: 1
  retries: RETRIES 
  conda:
    '../envs/rlibs.yaml'
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * MB_RENDER_HTMLS
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
      rule_name          = '$rule', 
      proj_dir           = '{PROJ_DIR}', 
      temp_f             =  '$template_f', 
      rmd_f              = '{output.rmd_f}', 
      YOUR_NAME          = '{YOUR_NAME}', 
      AFFILIATION        = '{AFFILIATION}', 
      PROJ_DIR           = '{PROJ_DIR}', 
      SHORT_TAG          = '{SHORT_TAG}', 
      DATE_STAMP         = '{DATE_STAMP}', 
      threads            = {threads}, 
      meta_f             = '{METADATA_F}', 
      qc_dt_f            = '{input.qc_dt_f}', 
      QC_HARD_MIN_COUNTS = {QC_HARD_MIN_COUNTS}, 
      QC_HARD_MIN_FEATS  = {QC_HARD_MIN_FEATS}, 
      QC_HARD_MAX_MITO   = {QC_HARD_MAX_MITO}, 
      QC_MIN_COUNTS      = {QC_MIN_COUNTS}, 
      QC_MIN_FEATS       = {QC_MIN_FEATS}, 
      QC_MIN_MITO        = {QC_MIN_MITO}, 
      QC_MAX_MITO        = {QC_MAX_MITO}, 
      QC_MIN_SPLICE      = {QC_MIN_SPLICE}, 
      QC_MAX_SPLICE      = {QC_MAX_SPLICE}, 
      QC_MIN_CELLS       = {QC_MIN_CELLS}
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
  threads: 1
  retries: RETRIES 
  conda:
    '../envs/rlibs.yaml'
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * MB_RENDER_HTMLS
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
      proj_dir    = '{PROJ_DIR}', 
      temp_f      =  '$template_f', 
      rmd_f       = '{output.rmd_f}', 
      YOUR_NAME   = '{YOUR_NAME}', 
      AFFILIATION = '{AFFILIATION}', 
      PROJ_DIR    = '{PROJ_DIR}', 
      SHORT_TAG   = '{SHORT_TAG}', 
      DATE_STAMP  = '{DATE_STAMP}', 
      threads     = {threads}, 
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
    qc_dt_f       = qc_dir  + '/qc_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    integration_f = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    r_int_f     = f"{code_dir}/integration.R",
    rmd_f       = f"{rmd_dir}/{SHORT_TAG}_integration.Rmd",
    html_f      = f"{docs_dir}/{SHORT_TAG}_integration.html"
  params: 
    int_res_ls = ','.join(map(str, INT_RES_LS))
  threads: 1
  retries: RETRIES 
  conda:
    '../envs/rlibs.yaml'
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * MB_RENDER_HTMLS
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
      rule_name     = '$rule', 
      proj_dir      = '{PROJ_DIR}', 
      temp_f        = '$template_f', 
      rmd_f         = '{output.rmd_f}', 
      YOUR_NAME     = '{YOUR_NAME}', 
      AFFILIATION   = '{AFFILIATION}', 
      SHORT_TAG     = '{SHORT_TAG}', 
      PROJ_DIR      = '{PROJ_DIR}',
      DATE_STAMP    = '{DATE_STAMP}', 
      threads       = {threads}, 
      qc_dt_f       = '{input.qc_dt_f}', 
      integration_f = '{input.integration_f}', 
      INT_RES_LS      = '{params.int_res_ls}', 
      INT_DBL_CL_PROP = {INT_DBL_CL_PROP}, 
      INT_REDUCTION   = '{INT_REDUCTION}', 
      DEMUX_TYPE      = '{DEMUX_TYPE}'
    )"
    """

# render_html_marker_genes
rule render_html_marker_genes:
  input:
    r_utils_f     = f"{code_dir}/utils.R",
    r_int_f       = f'{code_dir}/integration.R',
    pb_f          = mkr_dir + '/pb_' + FULL_TAG + f'_{MKR_SEL_RES}_' + DATE_STAMP + '.rds',
    mkrs_f        = mkr_dir + '/pb_marker_genes_' + FULL_TAG + f'_{MKR_SEL_RES}_' + DATE_STAMP + '.txt.gz',
    integration_f = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    hvgs_f        = mkr_dir + '/pb_hvgs_' + FULL_TAG + f'_{MKR_SEL_RES}_' + DATE_STAMP + '.txt.gz',
    ambient_f     = empty_dir + '/edger_empty_genes_' + FULL_TAG + '_all_' + DATE_STAMP + '.txt.gz',
    **get_conditional_outputs(SPECIES)
  output:
    r_mkr_f       = f"{code_dir}/marker_genes.R",
    rmd_f         = f'{rmd_dir}/{SHORT_TAG}_marker_genes_{MKR_SEL_RES}.Rmd',
    html_f        = f'{docs_dir}/{SHORT_TAG}_marker_genes_{MKR_SEL_RES}.html'
  threads: 8
  retries: RETRIES
  params:
    meta_vars  = ','.join(METADATA_VARS),
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
    mem_mb = lambda wildcards, attempt: attempt * MB_RENDER_HTMLS
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
      proj_dir          = '{PROJ_DIR}',
      temp_f            = '$template_f',
      rmd_f             = '{output.rmd_f}',
      YOUR_NAME         = '{YOUR_NAME}',
      AFFILIATION       = '{AFFILIATION}',
      PROJ_DIR          = '{PROJ_DIR}',
      SHORT_TAG         = '{SHORT_TAG}',
      DATE_STAMP        = '{DATE_STAMP}',
      threads           = {threads},
      meta_f            = '{METADATA_F}',
      meta_vars_ls      = '{params.meta_vars}',
      gtf_dt_f          = '{AF_GTF_DT_F}',
      ambient_f         = '{input.ambient_f}',
      integration_f     = '{input.integration_f}',
      pb_f              = '{input.pb_f}',
      mkrs_f            = '{input.mkrs_f}',
      CUSTOM_MKR_NAMES  = '{CUSTOM_MKR_NAMES}',
      CUSTOM_MKR_PATHS  = '{CUSTOM_MKR_PATHS}',
      hvgs_f            = '{input.hvgs_f}',
      {params.fgsea_args}
      MKR_SEL_RES       = {MKR_SEL_RES},
      MKR_NOT_OK_RE     = '{MKR_NOT_OK_RE}',
      MKR_MIN_CPM_MKR   = {MKR_MIN_CPM_MKR},
      MKR_MIN_CELLS     = {MKR_MIN_CELLS},
      MKR_GSEA_CUT      = {MKR_GSEA_CUT},
      SPECIES           = '{SPECIES}'
    )"
    """

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
  threads: 1
  retries: RETRIES
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * MB_RENDER_HTMLS
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_render_htmls/render_html_label_celltypes_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    # copy R code over
    echo "copying relevant R files over"
    cp scripts/label_celltypes.R {output.r_lbl_f}

    template_f=$(realpath resources/rmd_templates/label_celltypes.Rmd.template)
    rule="cell_labels"

    Rscript --vanilla -e "source('scripts/render_htmls.R'); \
    render_html(
      rule_name       = '$rule', 
      proj_dir        = '{PROJ_DIR}', 
      temp_f          =  '$template_f', 
      rmd_f           = '{output.rmd_f}', 
      YOUR_NAME       = '{YOUR_NAME}', 
      AFFILIATION     = '{AFFILIATION}', 
      SHORT_TAG       = '{SHORT_TAG}', 
      PROJ_DIR        = '{PROJ_DIR}', 
      DATE_STAMP      = '{DATE_STAMP}',
      threads         = {threads}, 
      guesses_f       = '{input.guesses_f}', 
      LBL_XGB_F       = '{LBL_XGB_F}', 
      LBL_XGB_CLS_F   = '{LBL_XGB_CLS_F}', 
      LBL_TISSUE      = '{LBL_TISSUE}', 
      LBL_SEL_RES_CL  = '{LBL_SEL_RES_CL}', 
      LBL_MIN_PRED    = {LBL_MIN_PRED}, 
      LBL_MIN_CL_PROP = {LBL_MIN_CL_PROP}, 
      LBL_MIN_CL_SIZE = {LBL_MIN_CL_SIZE}
    )"
    """
