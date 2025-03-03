
# rules to render html files

# localrules: render_html_index

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
#     template_f  = 'templates/marker_genes.Rmd.template'
#     print('rendering template')
#     render_html(PROJ_DIR, template_f, sub_dict, output.rmd_f)


# copy R scripts to code directory (not for rules that only run if specifically called)rule copy_r_code:
rule copy_r_code:
  output: 
    r_utils_f   = f"{code_dir}/utils.R",
    r_amb_f     = f"{code_dir}/ambient.R",
    r_qc_f      = f"{code_dir}/qc.R",
    r_sce_f     = f"{code_dir}/make_sce.R",
    r_dbl_f     = f"{code_dir}/doublet_id.R",
    r_int_f     = f"{code_dir}/integration.R",
    r_mkr_f     = f"{code_dir}/marker_genes.R",
    r_zoom_f	  = f"{code_dir}/zoom.R", 
    r_demux_f   = f"{code_dir}/multiplexing.R"
  shell: 
    """
    echo "copying relevant R files over"
    
    cp scripts/utils.R {output.r_utils_f}
    cp scripts/ambient.R {output.r_amb_f}
    cp scripts/SampleQC.R {output.r_qc_f}
    cp scripts/make_sce.R {output.r_sce_f}
    cp scripts/doublet_id.R {output.r_dbl_f}
    cp scripts/integration.R {output.r_int_f}
    cp scripts/marker_genes.R {output.r_mkr_f}
    cp scripts/zoom.R {output.r_zoom_f}    
    cp scripts/multiplexing.R {output.r_demux_f}    
    """

 


rule render_html_alevin_fry:
  input:
    expand(af_dir + '/af_{sample}/' + af_rna_dir + 'knee_plot_data_{sample}_' + DATE_STAMP + '.txt.gz', sample=runs)
  output:
    rmd_f       = f"{rmd_dir}/{SHORT_TAG}_alevin_fry.Rmd",
    html_f      = f"{docs_dir}/{SHORT_TAG}_alevin_fry.html"
  threads: 1
  retries: RETRIES
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * 4096
  conda:
    '../envs/rlibs.yml'
  shell:
    """
    echo "copying relevant R files over"
    
    # make and render Rmd file
    template_f=$(realpath templates/alevin_fry.Rmd.template)
    rule="af"

    # rendering html
    Rscript --vanilla -e "source('scripts/render_reports.R'); \
        render_reports(
        rule_name      = '$rule', \
        proj_dir       = '{PROJ_DIR}', \
        temp_f         =  '$template_f', \
        rmd_f          = '{output.rmd_f}', \
        YOUR_NAME      = '{YOUR_NAME}', \
        AFFILIATION    = '{AFFILIATION}', \
        SHORT_TAG      = '{SHORT_TAG}', \
        DATE_STAMP     = '{DATE_STAMP}', \
        RUNS_STR       = '{RUNS_STR}', \
        AMBIENT_METHOD = '{AMBIENT_METHOD}', \
        SAMPLE_VAR     = '{SAMPLE_VAR}', \
        af_dir         = '{af_dir}', \
        af_rna_dir     = '{af_rna_dir}')"

    """

if DEMUX_TYPE == 'af':
  rule render_html_multiplexing:
    input:
      expand(af_dir + '/af_{sample}/hto/' + 'knee_plot_data_{sample}_' + DATE_STAMP + '.txt.gz', sample=runs), 
      sce_hto_f = sce_dir + '/sce_cells_htos_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
    output:
      rmd_f       = f"{rmd_dir}/{SHORT_TAG}_multiplexing.Rmd",
      html_f      = f"{docs_dir}/{SHORT_TAG}_multiplexing.html"
    threads: 1
    retries: RETRIES
    resources:
      mem_mb      =  lambda wildcards, attempt: attempt * 4096
    conda:
      '../envs/rlibs.yml'
    shell:
      """
      echo "copying relevant R files over"
    
      # make and render Rmd file
      template_f=$(realpath templates/multiplexing.Rmd.template)
      rule="multiplexing"

      # rendering html
      Rscript --vanilla -e "source('scripts/render_reports.R'); \
          render_reports(
          rule_name      = '$rule', \
          proj_dir       = '{PROJ_DIR}', \
          temp_f         =  '$template_f', \
          rmd_f          = '{output.rmd_f}', \
          sce_hto_f      = '{input.sce_hto_f}', \
          YOUR_NAME      = '{YOUR_NAME}', \
          AFFILIATION    = '{AFFILIATION}', \
          SHORT_TAG      = '{SHORT_TAG}', \
          DATE_STAMP     = '{DATE_STAMP}', \
          RUNS_STR       = '{RUNS_STR}', \
          METADATA_F     = '{METADATA_F}', \
          AMBIENT_METHOD = '{AMBIENT_METHOD}', \
          SAMPLE_VAR     = '{SAMPLE_VAR}', \
          af_dir         = '{af_dir}')"

      """


# render_html_ambient
rule render_html_ambient: # some outputs are the same as outputs in render_html_alevin_fry
  input:
    expand( amb_dir + '/ambient_{sample}/barcodes_qc_metrics_{sample}_' + DATE_STAMP + '.txt.gz', sample = runs )
  output:
    rmd_f       = f"{rmd_dir}/{SHORT_TAG}_ambient.Rmd",
    html_f      = f"{docs_dir}/{SHORT_TAG}_ambient.html"
  threads: 1
  retries: RETRIES 
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * 4096
  conda:
    '../envs/rlibs.yml'
  shell:
    """
        template_f=$(realpath templates/ambient.Rmd.template)
        rule="ambient"
        
        Rscript --vanilla -e "source('scripts/render_reports.R'); \
        render_reports(
        rule_name = '$rule', \
        proj_dir = '{PROJ_DIR}', \
        temp_f =  '$template_f', \
        rmd_f = '{output.rmd_f}', \
        YOUR_NAME = '{YOUR_NAME}', \
        AFFILIATION = '{AFFILIATION}', \
        SHORT_TAG = '{SHORT_TAG}', \
        DATE_STAMP = '{DATE_STAMP}', \
        SAMPLE_STR = '{RUNS_STR}', \
        AMBIENT_METHOD = '{AMBIENT_METHOD}', \
        af_dir = '{af_dir}')"
    """

rule render_html_qc:
  input:
    qc_f        = qc_dir  + '/qc_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    keep_f      = qc_dir  + '/keep_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    rmd_f       = f"{rmd_dir}/{SHORT_TAG}_qc.Rmd",
    html_f      = f"{docs_dir}/{SHORT_TAG}_qc.html"
  threads: 1
  retries: RETRIES 
  conda:
    '../envs/rlibs.yml'
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * 4096
  shell:
    """
  
    #define rule and template
    template_f=$(realpath templates/SampleQC.Rmd.template)
    rule="qc"

    # rendering html
    Rscript --vanilla -e "source('scripts/render_reports.R'); \
    render_reports(
    rule_name = '$rule', \
    proj_dir = '{PROJ_DIR}', \
    temp_f =  '$template_f', \
    rmd_f = '{output.rmd_f}', \
    YOUR_NAME = '{YOUR_NAME}', \
    AFFILIATION = '{AFFILIATION}', \
    SHORT_TAG = '{SHORT_TAG}', \
    DATE_STAMP = '{DATE_STAMP}', \
    threads = {threads}, \
    meta_f = '{METADATA_F}', \
    qc_dt_f = '{input.qc_f}', \
    qc_keep_f = '{input.keep_f}', \
    AMBIENT_METHOD = '{AMBIENT_METHOD}', \
    QC_HARD_MIN_COUNTS = {QC_HARD_MIN_COUNTS}, \
    QC_HARD_MIN_FEATS = {QC_HARD_MIN_FEATS}, \
    QC_HARD_MAX_MITO = {QC_HARD_MAX_MITO}, \
    QC_MIN_COUNTS = {QC_MIN_COUNTS}, \
    QC_MIN_FEATS = {QC_MIN_FEATS}, \
    QC_MIN_MITO = {QC_MIN_MITO}, \
    QC_MAX_MITO = {QC_MAX_MITO}, \
    QC_MIN_SPLICE = {QC_MIN_SPLICE}, \
    QC_MAX_SPLICE = {QC_MAX_SPLICE}, \
    QC_MIN_CELLS = {QC_MIN_CELLS}, \
    QC_FILTER_BENDER = '{QC_FILTER_BENDER}')"
    
    """



rule render_html_integration:
  input:
    sce_f       = sce_dir + '/sce_cells_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    qc_f        = qc_dir  + '/qc_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    keep_f      = qc_dir  + '/keep_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    harmony_f   = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    hvgs_f      = int_dir   + '/harmony_hvgs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    rmd_f       = f"{rmd_dir}/{SHORT_TAG}_integration.Rmd",
    html_f      = f"{docs_dir}/{SHORT_TAG}_integration.html"
  params: 
    int_res_ls = ','.join(map(str, INT_RES_LS))
  threads: 1
  retries: RETRIES 
  conda:
    '../envs/rlibs.yml'
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * 4096
  shell:
    """
    template_f=$(realpath templates/integration.Rmd.template)
    rule="integration"
    
    Rscript --vanilla -e "source('scripts/render_reports.R'); \
    render_reports(
    rule_name = '$rule', \
    proj_dir = '{PROJ_DIR}', \
    temp_f =  '$template_f', \
    rmd_f = '{output.rmd_f}', \
    YOUR_NAME = '{YOUR_NAME}', \
    AFFILIATION = '{AFFILIATION}', \
    SHORT_TAG = '{SHORT_TAG}', \
    DATE_STAMP = '{DATE_STAMP}', \
    threads = {threads}, \
    sce_all_f = '{input.sce_f}', \
    qc_dt_f = '{input.qc_f}', \
    qc_keep_f = '{input.keep_f}', \
    hmny_f = '{input.harmony_f}', \
    hmny_hvgs_f = '{input.hvgs_f}', \
    INT_RES_LS = '{params.int_res_ls}', \
    INT_SEL_RES = '{INT_SEL_RES}', \
    INT_DBL_CL_PROP = {INT_DBL_CL_PROP})"
    """


# # render_html_label_celltypes
rule render_html_label_celltypes:
  input:
    r_mkr_f     = f'{code_dir}/marker_genes.R',
    guesses_f   = f'{lbl_dir}/cell_annotations_{FULL_TAG}_{DATE_STAMP}.txt.gz'
  output:
    r_lbl_f     = f'{code_dir}/label_celltypes.R',
    rmd_f       = f'{rmd_dir}/{SHORT_TAG}_label_celltypes.Rmd',
    html_f      = f'{docs_dir}/{SHORT_TAG}_label_celltypes.html'
  threads: 1
  retries: RETRIES
  resources:
    mem_mb      =  lambda wildcards, attempt: attempt * 4096
  conda: 
    '../envs/rlibs.yml'
  shell:
    """
    # copy R code over
    echo "copying relevant R files over"
    cp scripts/label_celltypes.R {output.r_lbl_f}

    template_f=$(realpath templates/label_celltypes.Rmd.template)
    rule="cell_labels"

    Rscript --vanilla -e "source('scripts/render_reports.R'); \
    render_reports(
    rule_name = '$rule', \
    proj_dir = '{PROJ_DIR}', \
    temp_f =  '$template_f', \
    rmd_f = '{output.rmd_f}', \
    YOUR_NAME = '{YOUR_NAME}', \
    AFFILIATION = '{AFFILIATION}', \
    SHORT_TAG = '{SHORT_TAG}', \
    DATE_STAMP = '{DATE_STAMP}', \
    threads = {threads}, \
    guesses_f = '{input.guesses_f}', \
    LBL_XGB_F = '{LBL_XGB_F}', \
    CUSTOM_LABELS_F = '{CUSTOM_LABELS_F}', \
    INT_SEL_RES = '{INT_SEL_RES}', \
    LBL_TISSUE = '{LBL_TISSUE}', \
    LBL_SEL_RES_CL = '{LBL_SEL_RES_CL}', \
    LBL_MIN_PRED = {LBL_MIN_PRED}, \
    LBL_MIN_CL_PROP = {LBL_MIN_CL_PROP}, \
    LBL_MIN_CL_SIZE = {LBL_MIN_CL_SIZE})"
    """


# # render_html_marker_genes
rule render_html_marker_genes:
  input:
    r_int_f     = f'{code_dir}/integration.R',
    pb_f        = mkr_dir + '/pb_'              + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.rds', 
    mkrs_f      = mkr_dir   + '/pb_marker_genes_' + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.txt.gz', 
    harmony_f   = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    hvgs_f      = mkr_dir   + '/pb_hvgs_'         + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.txt.gz',
    **(
        {  # Include FGSEA outputs **only if** SPECIES is in the allowed list
            'fgsea_go_bp_f':   mkr_dir + '/fgsea_' + FULL_TAG + f'_{INT_SEL_RES}_go_bp_' + DATE_STAMP + '.txt.gz',
            'fgsea_go_cc_f':   mkr_dir + '/fgsea_' + FULL_TAG + f'_{INT_SEL_RES}_go_cc_' + DATE_STAMP + '.txt.gz',
            'fgsea_go_mf_f':   mkr_dir + '/fgsea_' + FULL_TAG + f'_{INT_SEL_RES}_go_mf_' + DATE_STAMP + '.txt.gz',
            'fgsea_paths_f':   mkr_dir + '/fgsea_' + FULL_TAG + f'_{INT_SEL_RES}_paths_' + DATE_STAMP + '.txt.gz',
            'fgsea_hlmk_f':    mkr_dir + '/fgsea_' + FULL_TAG + f'_{INT_SEL_RES}_hlmk_' + DATE_STAMP + '.txt.gz'
        } if SPECIES in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020'] else {}
    )
  output:
    rmd_f       = f'{rmd_dir}/{SHORT_TAG}_marker_genes_{INT_SEL_RES}.Rmd',
    html_f      = f'{docs_dir}/{SHORT_TAG}_marker_genes_{INT_SEL_RES}.html'
  threads: 8
  retries: RETRIES
  params: 
    meta_vars = ','.join(METADATA_VARS)
  conda: 
    '../envs/rlibs.yml'
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * MB_HTML_MARKER_GENES
  shell:
    """

    template_f=$(realpath templates/marker_genes.Rmd.template)
    rule="markers"

    # define what we will substitute in

    Rscript --vanilla -e "source('scripts/render_reports.R'); \
    render_reports(
    rule_name = '$rule', \
    proj_dir = '{PROJ_DIR}', \
    temp_f =  '$template_f', \
    rmd_f = '{output.rmd_f}', \
    YOUR_NAME = '{YOUR_NAME}', \
    AFFILIATION = '{AFFILIATION}', \
    SHORT_TAG = '{SHORT_TAG}', \
    DATE_STAMP = '{DATE_STAMP}', \
    threads = {threads}, \
    meta_f = '{METADATA_F}', \
    meta_vars_ls = '{params.meta_vars}', \
    gtf_dt_f = '{AF_GTF_DT_F}', \
    hmny_f = '{input.harmony_f}', \
    pb_f = '{input.pb_f}', \
    mkrs_f = '{input.mkrs_f}', \
    CUSTOM_MKR_NAMES = '{CUSTOM_MKR_NAMES}', \
    CUSTOM_MKR_PATHS = '{CUSTOM_MKR_PATHS}', \
    hvgs_f = '{input.hvgs_f}', \
    {f'fgsea_go_bp_f = {output.fgsea_go_bp_f},' if 'fgsea_go_bp_f' in output else ''} \
    {f'fgsea_go_cc_f = {output.fgsea_go_cc_f},' if 'fgsea_go_cc_f' in output else ''} \
    {f'fgsea_go_mf_f = {output.fgsea_go_mf_f},' if 'fgsea_go_mf_f' in output else ''} \
    {f'fgsea_paths_f = {output.fgsea_paths_f},' if 'fgsea_paths_f' in output else ''} \
    {f'fgsea_hlmk_f  = {output.fgsea_hlmk_f},' if 'fgsea_hlmk_f' in output else ''} \
    INT_EXC_REGEX = '{INT_EXC_REGEX}', \
    INT_SEL_RES = {INT_SEL_RES}, \
    MKR_NOT_OK_RE = '{MKR_NOT_OK_RE}', \
    MKR_MIN_CPM_MKR = {MKR_MIN_CPM_MKR}, \
    MKR_MIN_CELLS = {MKR_MIN_CELLS}, \
    MKR_GSEA_CUT = {MKR_GSEA_CUT}, \
    SPECIES = '{SPECIES}')"

    """


# # render_html_zoom
rule render_html_zoom:
  input:
    zoom_sce_sub_f      = zoom_dir + '/{zoom_name}/' + 'zoom_sce_clean_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.rds',
    zoom_hmny_f         = zoom_dir + '/{zoom_name}/' + 'zoom_integrated_dt_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.txt.gz',
    zoom_pb_f           = zoom_dir + '/{zoom_name}/' + 'zoom_pb_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.rds',
    zoom_mkrs_f         = zoom_dir + '/{zoom_name}/' + 'zoom_pb_marker_genes_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.txt.gz',
    zoom_hvgs_f         = zoom_dir + '/{zoom_name}/' + 'zoom_pb_hvgs_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.txt.gz',
    zoom_fgsea_go_bp_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_go_bp_' + DATE_STAMP +'.txt.gz',
    zoom_fgsea_go_cc_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_go_cc_' + DATE_STAMP +'.txt.gz',
    zoom_fgsea_go_mf_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_go_mf_' + DATE_STAMP +'.txt.gz',
    zoom_fgsea_paths_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_paths_' + DATE_STAMP +'.txt.gz',
    zoom_fgsea_hlmk_f   = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_hlmk_' + DATE_STAMP +'.txt.gz'
  output:
    rmd_f       = rmd_dir + '/' + SHORT_TAG + '_zoom' + '_{zoom_name}_{zoom_res}.Rmd',
    html_f      = docs_dir + '/' + SHORT_TAG + '_zoom' + '_{zoom_name}_{zoom_res}.html'
  params:
    zoom_name   = '{zoom_name}',
    zoom_res    = '{zoom_res}', 
    meta_vars   = ','.join(METADATA_VARS)
  threads: 1
  retries: RETRIES 
  conda: 
    '../envs/rlibs.yml'
  resources:
    mem_mb =  lambda wildcards, attempt: attempt * 8192
  shell:
    """

    template_f=$(realpath templates/zoom.Rmd.template)
    rule="zoom"

    Rscript --vanilla -e "source('scripts/render_reports.R'); \
    render_reports(
    rule_name = '$rule', \
    proj_dir = '{PROJ_DIR}', \
    temp_f =  '$template_f', \
    rmd_f = '{output.rmd_f}', \
    YOUR_NAME = '{YOUR_NAME}', \
    AFFILIATION = '{AFFILIATION}', \
    SHORT_TAG = '{SHORT_TAG}', \
    DATE_STAMP = '{DATE_STAMP}', \
    threads = {threads}, \
    zoom_name = '{params.zoom_name}', \
    zoom_res = {params.zoom_res}, \
    meta_f = '{METADATA_F}', \
    meta_vars_ls = '{params.meta_vars}', \
    gtf_dt_f = '{AF_GTF_DT_F}', \
    sce_sub_f = '{input.zoom_sce_sub_f}', \
    hmny_f = '{input.zoom_hmny_f}', \
    pb_f = '{input.zoom_pb_f}', \
    mkrs_f = '{input.zoom_mkrs_f}', \
    hvgs_f = '{input.zoom_hvgs_f}', \
    canon_f = '{MKR_CANON_F}', \
    fgsea_go_bp_f = '{input.zoom_fgsea_go_bp_f}', \
    fgsea_go_cc_f = '{input.zoom_fgsea_go_cc_f}', \
    fgsea_go_mf_f = '{input.zoom_fgsea_go_mf_f}', \
    fgsea_paths_f = '{input.zoom_fgsea_paths_f}', \
    fgsea_hlmk_f = '{input.zoom_fgsea_hlmk_f}', \
    INT_DBL_CL_PROP = {INT_DBL_CL_PROP}, \
    INT_EXC_REGEX = '{INT_EXC_REGEX}', \
    MKR_NOT_OK_RE = '{MKR_NOT_OK_RE}', \
    MKR_MIN_CPM_MKR = {MKR_MIN_CPM_MKR}, \
    MKR_MIN_CELLS = {MKR_MIN_CELLS}, \
    MKR_GSEA_CUT = {MKR_GSEA_CUT}, \
    SPECIES = '{SPECIES}')"

    """


# render_html_empties
rule render_html_empties:
  input:
    r_int_f     = f'{code_dir}/{SHORT_TAG}06_integration.R',
    guesses_f   = lbl_dir + '/xgboost_guesses_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    empty_csv_f = f'{lbl_dir}/empties_{FULL_TAG}_{DATE_STAMP}.csv',
    empty_gs_ls = expand( [lbl_dir + '/empty_genes' + '/empty_genes_' + FULL_TAG + '_{subset}_' + DATE_STAMP + '.txt.gz'], \
    subset = None if LBL_SCE_SUBSETS is None else [*LBL_SCE_SUBSETS] )
  output:
    r_pb_f      = f'{code_dir}/pseudobulk_and_empties.R',
    rmd_f       = f'{rmd_dir}/{SHORT_TAG}_empties.Rmd',
    html_f      = f'{docs_dir}/{SHORT_TAG}_empties.html'
  threads: 1
  retries: RETRIES 
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * 4096
  conda: 
    '../envs/rlibs.yml'
  shell:
    """
    # copy R code over
    echo "copying relevant R files over"
    cp scripts/pseudobulk_and_empties.R {output.r_pb_f}

    template_f=$(realpath templates/empties.Rmd.template)
    rule="pb_empties"

    
    Rscript --vanilla -e "source('scripts/render_reports.R'); \
    render_reports(
    rule_name = '$rule', \
    proj_dir = '{PROJ_DIR}', \
    temp_f =  '$template_f', \
    rmd_f = '{output.rmd_f}', \
    YOUR_NAME = '{YOUR_NAME}', \
    AFFILIATION = '{AFFILIATION}', \
    SHORT_TAG = '{SHORT_TAG}', \
    DATE_STAMP = '{DATE_STAMP}', \
    threads = {threads}, \
    guesses_f = '{input.guesses_f}', \
    empty_csv_f = '{input.empty_csv_f}', \
    LBL_XGB_F = '{LBL_XGB_F}', \
    LBL_SEL_RES_CL = {LBL_SEL_RES_CL}, \
    LBL_MIN_PRED = {LBL_MIN_PRED}, \
    LBL_MIN_CL_PROP = {LBL_MIN_CL_PROP}, \
    LBL_MIN_CL_SIZE = {LBL_MIN_CL_SIZE})"
     
    """

