# rules to render html files

# localrules: render_html_index

# # render_html_index --> render from index.Rmd file that is created when workflow R project is started ?
# rule render_html_index:
#   output:
#     rmd_f       = f'{rmd_dir}/index.Rmd',
#     html_f      = f'{docs_dir}/index.html'
#   threads: 1
#   resources:
#     mem_mb      = 4096
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





# # render_html_alevin_fry
# rule render_html_alevin_fry:
#   output:
#     r_utils_f   = f"{code_dir}/utils.R",
#     r_amb_f      = f"{code_dir}/ambient.R",
#     rmd_f       = f"{rmd_dir}/{SHORT_TAG}_alevin_fry.Rmd",
#     html_f      = f"{docs_dir}/{SHORT_TAG}_alevin_fry.html"
#   threads: 1
#   resources:
#     mem_mb      = 4096
#   conda:
#     '../envs/rlibs.yml'
#   shell:
#     """
#     # copy R code over
#     echo "copying relevant R files over"
    
#     cp scripts/utils.R {output.r_utils_f}
#     cp scripts/ambient.R {output.r_amb_f}

#     if [ "{EXC_SAMPLES}" == "None" ]; then
#         exc_samples_ls='""'
#     else
#         exc_samples_ls='c("'$(IFS='","'; echo "${EXC_SAMPLES[*]}")'")'
#     fi


#     meta_f=$(realpath {METADATA_F})
#     af_dir=$(realpath {af_dir})
    

#     # define what we will substitute in
#     echo "setting up template"

#     sub_ls=$(jq -n \
#             --arg YOUR_NAME "{YOUR_NAME}" \
#             --arg AFFILIATION "{AFFILIATION}" \
#             --arg SHORT_TAG "{SHORT_TAG}" \
#             --arg DATE_STAMP "{DATE_STAMP}" \
#             --arg meta_f "$meta_f" \
#             --arg exc_samples_ls "$exc_samples_ls" \
#             --arg af_dir "$af_dir" \
#             '{
#                 YOUR_NAME: $YOUR_NAME,
#                 AFFILIATION: $AFFILIATION,
#                 SHORT_TAG: $SHORT_TAG,
#                 DATE_STAMP: $DATE_STAMP,
#                 meta_f: $meta_f,
#                 exc_samples_ls: $exc_samples_ls,
#                 af_dir: $af_dir
#             }')


#     # make and render Rmd file
#     template_f="templates/alevin_fry.Rmd.template"
#     echo "rendering html"


#     Rscript -e "source('scripts/render_reports.R'); \
#     render_reports(
#     proj_dir = '{PROJ_DIR}', \
#     temp_f = '{template_f}', \
#     temp_ls = $sub_ls, \
#     rmd_f = '{output.rmd_f}')"

#     """


# render_html_ambient
rule render_html_ambient: # some outputs are the same as outputs in render_html_alevin_fry
  input:
    expand( amb_dir + '/ambient_{sample}/barcodes_qc_metrics_{sample}_' + DATE_STAMP + '.txt.gz', sample = SAMPLES )
  output:
    rmd_f       = f"{rmd_dir}/{SHORT_TAG}_ambient.Rmd"
    #html_f      = f"{docs_dir}/{SHORT_TAG}_ambient.html"
  threads: 1
  resources:
    mem_mb      = 4096
  conda:
    '../envs/rlibs.yml'
  shell:
    """
        template_f=$(realpath templates/ambient.Rmd.template)
        rule="ambient"
        
        echo "Rendering html"
        Rscript -e "source('scripts/render_reports.R'); \
        render_reports(
        rule_name = '$rule', \
        proj_dir = '{proj_dir}', \
        temp_f =  '$template_f', \
        rmd_f = '{output.rmd_f}', \
        YOUR_NAME = '{YOUR_NAME}', \
        AFFILIATION = '{AFFILIATION}', \
        SHORT_TAG = '{SHORT_TAG}', \
        DATE_STAMP = '{DATE_STAMP}', \
        SAMPLE_STR = '{SAMPLE_STR}', \
        AMBIENT_METHOD = '{AMBIENT_METHOD}', \
        CUSTOM_PARAMS_F = '{CUSTOM_PARAMS_F}', \
        af_dir = '{af_dir}')"
    """

# rule render_html_qc:
#   input:
#     qc_f        = qc_dir  + '/qc_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
#     keep_f      = qc_dir  + '/keep_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
#   output:
#     r_qc_f      = f'{code_dir}/qc.R',
#     rmd_f       = f"{rmd_dir}/{SHORT_TAG}_qc.Rmd",
#     html_f      = f"{docs_dir}/{SHORT_TAG}_qc.html"
#   threads: 8
#   retries: 5
#   resources:
#     mem_mb      = lambda wildcards, attempt: attempt * 4096
#   run:
#     # copy R code over
#     print('copying relevant R files over')
#     import shutil
#     shutil.copyfile('scripts/SampleQC.R', output.r_qc_f)

#     # define what we will substitute in
#     print('setting up template')
#     sub_dict     = {
#       'YOUR_NAME':          YOUR_NAME,
#       'AFFILIATION':        AFFILIATION,
#       'SHORT_TAG':          SHORT_TAG,
#       'DATE_STAMP':         DATE_STAMP,
#       'threads':            threads,
#       'meta_f':             os.path.relpath(METADATA_F, PROJ_DIR),
#       'qc_dt_f':            os.path.relpath(qc_dir  + '/qc_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', PROJ_DIR),
#       'qc_keep_f':          os.path.relpath(qc_dir  + '/keep_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', PROJ_DIR),
#       'QC_HARD_MIN_COUNTS': QC_HARD_MIN_COUNTS,
#       'QC_HARD_MIN_FEATS':  QC_HARD_MIN_FEATS,
#       'QC_HARD_MAX_MITO':   QC_HARD_MAX_MITO,
#       'QC_MIN_COUNTS':      QC_MIN_COUNTS,
#       'QC_MIN_FEATS':       QC_MIN_FEATS,
#       'QC_MIN_MITO':        QC_MIN_MITO,
#       'QC_MAX_MITO':        QC_MAX_MITO,
#       'QC_MIN_SPLICE':      QC_MIN_SPLICE,
#       'QC_MAX_SPLICE':      QC_MAX_SPLICE,
#       'QC_MIN_CELLS':       QC_MIN_CELLS,
#       'QC_FILTER_BENDER':   QC_FILTER_BENDER
#       }
#     # make and render Rmd file
#     template_f  = 'templates/SampleQC.Rmd.template'
#     print('rendering template')
#     render_html(PROJ_DIR, RLIBS_DIR, template_f, sub_dict, output.rmd_f)


# # render_html_integration
# rule render_html_integration:
#   input:
#     qc_f        = qc_dir  + '/qc_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
#     keep_f      = qc_dir  + '/keep_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
#     harmony_f   = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
#   output:
#     r_sce_f     = f"{code_dir}/make_sce.R",
#     r_dbl_f     = f"{code_dir}/doublet_id.R",
#     r_int_f     = f"{code_dir}/integration.R",
#     rmd_f       = f"{rmd_dir}/{SHORT_TAG}_integration.Rmd",
#     html_f      = f"{docs_dir}/{SHORT_TAG}_integration.html"
#   threads: 8
#   retries: 5
#   resources:
#     mem_mb      = lambda wildcards, attempt: attempt * 4096
#   run:
#     # copy R code over
#     print('copying relevant R files over')
#     import shutil
#     shutil.copyfile('scripts/make_sce.R', output.r_sce_f)
#     shutil.copyfile('scripts/doublet_id.R', output.r_dbl_f)
#     shutil.copyfile('scripts/integration.R', output.r_int_f)

#     # define what we will substitute in
#     print('setting up template')
#     sub_dict    = {
#       'YOUR_NAME':        YOUR_NAME,
#       'AFFILIATION':      AFFILIATION,
#       'SHORT_TAG':        SHORT_TAG,
#       'DATE_STAMP':       DATE_STAMP,
#       'threads':          threads,
#       'sce_all_f':        os.path.relpath(sce_dir + '/sce_' + ('bender' if DO_CELLBENDER else 'alevin') + \
#         '_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds', PROJ_DIR),
#       'qc_dt_f':          os.path.relpath(qc_dir  + '/qc_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', PROJ_DIR),
#       'qc_keep_f':        os.path.relpath(qc_dir  + '/keep_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', PROJ_DIR),
#       'hmny_f':           os.path.relpath(int_dir + '/integrated_dt_'       + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', PROJ_DIR),
#       'hmny_hvgs_f':      os.path.relpath(int_dir + '/harmony_hvgs_'        + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', PROJ_DIR),
#       'INT_RES_LS':       'c(' + ', '.join(map(str, INT_RES_LS)) + ')',
#       'INT_SEL_RES':      INT_SEL_RES,
#       'INT_DBL_CL_PROP':  INT_DBL_CL_PROP
#       }
#     # make and render Rmd file
#     template_f  = 'templates/integration.Rmd.template'
#     print('rendering template')
#     render_html(PROJ_DIR, RLIBS_DIR, template_f, sub_dict, output.rmd_f)


# # render_html_label_celltypes
# rule render_html_label_celltypes:
#   input:
#     r_mkr_f     = f'{code_dir}/marker_genes.R',
#     guesses_f   = f'{lbl_dir}/xgboost_guesses_{FULL_TAG}_{DATE_STAMP}.txt.gz'
#   output:
#     r_lbl_f     = f'{code_dir}/label_celltypes.R',
#     rmd_f       = f'{rmd_dir}/{SHORT_TAG}_label_celltypes.Rmd',
#     html_f      = f'{docs_dir}/{SHORT_TAG}_label_celltypes.html'
#   threads: 1
#   retries: 5
#   resources:
#     mem_mb      = lambda wildcards, attempt: attempt * 4096
#   run:
#     # copy R code over
#     print('copying relevant R files over')
#     import shutil
#     shutil.copyfile('scripts/label_celltypes.R', output.r_lbl_f)

#     # define what we will substitute in
#     print('setting up template')
#     sub_dict    = {
#       'YOUR_NAME':        YOUR_NAME,
#       'AFFILIATION':      AFFILIATION,
#       'SHORT_TAG':        SHORT_TAG,
#       'DATE_STAMP':       DATE_STAMP,
#       'threads':          threads,
#       'guesses_f':        os.path.relpath(input.guesses_f, PROJ_DIR),
#       'LBL_XGB_F':        LBL_XGB_F,
#       'LBL_SEL_RES_CL':   LBL_SEL_RES_CL,
#       'LBL_MIN_PRED':     LBL_MIN_PRED,
#       'LBL_MIN_CL_PROP':  LBL_MIN_CL_PROP,
#       'LBL_MIN_CL_SIZE':  LBL_MIN_CL_SIZE
#       }

#     # make and render Rmd file
#     template_f  = 'templates/label_celltypes.Rmd.template'
#     print('rendering template')
#     render_html(PROJ_DIR, RLIBS_DIR, template_f, sub_dict, output.rmd_f)


# # render_html_marker_genes
# rule render_html_marker_genes:
#   input:
#     r_int_f     = f'{code_dir}/integration.R',
#     pb_f        = mkr_dir + '/pb_'              + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.rds'
#   output:
#     r_mkr_f     = f'{code_dir}/marker_genes.R',
#     rmd_f       = f'{rmd_dir}/{SHORT_TAG}_marker_genes_{INT_SEL_RES}.Rmd',
#     html_f      = f'{docs_dir}/{SHORT_TAG}_marker_genes_{INT_SEL_RES}.html'
#   threads: 8
#   retries: 5
#   resources:
#     mem_mb      = lambda wildcards, attempt: attempt * MB_HTML_MARKER_GENES
#   run:
#     # copy R code over
#     print('copying relevant R files over')
#     import shutil
#     shutil.copyfile('scripts/marker_genes.R', output.r_mkr_f)

#     # define what we will substitute in
#     print('setting up template')
#     sub_dict    = {
#       'YOUR_NAME':        YOUR_NAME,
#       'AFFILIATION':      AFFILIATION,
#       'SHORT_TAG':        SHORT_TAG,
#       'DATE_STAMP':       DATE_STAMP,
#       'threads':          threads,
#       'meta_f':           os.path.relpath(METADATA_F, PROJ_DIR),
#       'meta_vars_ls':     'c("' + '", "'.join(map(str, METADATA_VARS)) + '")',
#       'gtf_dt_f':         AF_GTF_DT_F,
#       'hmny_f':           os.path.relpath(int_dir + '/integrated_dt_'   + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', PROJ_DIR),
#       'pb_f':             os.path.relpath(mkr_dir + '/pb_'              + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.rds', PROJ_DIR),
#       'mkrs_f':           os.path.relpath(mkr_dir + '/pb_marker_genes_' + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.txt.gz', PROJ_DIR),
#       'canon_f':          MKR_CANON_F,
#       'hvgs_f':           os.path.relpath(mkr_dir + '/pb_hvgs_'         + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.txt.gz', PROJ_DIR),
#       'fgsea_go_bp_f':    os.path.relpath(mkr_dir + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'go_bp_' + DATE_STAMP + '.txt.gz', PROJ_DIR),
#       'fgsea_go_cc_f':    os.path.relpath(mkr_dir + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'go_cc_' + DATE_STAMP + '.txt.gz', PROJ_DIR),
#       'fgsea_go_mf_f':    os.path.relpath(mkr_dir + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'go_mf_' + DATE_STAMP + '.txt.gz', PROJ_DIR),
#       'fgsea_paths_f':    os.path.relpath(mkr_dir + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'paths_' + DATE_STAMP + '.txt.gz', PROJ_DIR),
#       'fgsea_hlmk_f':     os.path.relpath(mkr_dir + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'hlmk_' + DATE_STAMP + '.txt.gz', PROJ_DIR),
#       'INT_EXC_REGEX':    INT_EXC_REGEX,
#       'INT_SEL_RES':      INT_SEL_RES,
#       'MKR_NOT_OK_RE':    MKR_NOT_OK_RE,
#       'MKR_MIN_CPM_MKR':  MKR_MIN_CPM_MKR,
#       'MKR_MIN_CELLS':    MKR_MIN_CELLS,
#       'MKR_GSEA_CUT':     MKR_GSEA_CUT
#       }
#     # make and render Rmd file
#     template_f  = 'templates/marker_genes.Rmd.template'
#     print('rendering template')
#     render_html(PROJ_DIR, RLIBS_DIR, template_f, sub_dict, output.rmd_f)


# # render_html_zoom
# rule render_html_zoom:
#   input:
#     zoom_sce_sub_f      = zoom_dir + '/{zoom_name}/' + 'zoom_sce_clean_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.rds',
#     zoom_hmny_f         = zoom_dir + '/{zoom_name}/' + 'zoom_integrated_dt_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.txt.gz',
#     zoom_pb_f           = zoom_dir + '/{zoom_name}/' + 'zoom_pb_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.rds',
#     zoom_mkrs_f         = zoom_dir + '/{zoom_name}/' + 'zoom_pb_marker_genes_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.txt.gz',
#     zoom_hvgs_f         = zoom_dir + '/{zoom_name}/' + 'zoom_pb_hvgs_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.txt.gz',
#     zoom_fgsea_go_bp_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_go_bp_' + DATE_STAMP +'.txt.gz',
#     zoom_fgsea_go_cc_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_go_cc_' + DATE_STAMP +'.txt.gz',
#     zoom_fgsea_go_mf_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_go_mf_' + DATE_STAMP +'.txt.gz',
#     zoom_fgsea_paths_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_paths_' + DATE_STAMP +'.txt.gz',
#     zoom_fgsea_hlmk_f   = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_hlmk_' + DATE_STAMP +'.txt.gz'
#   output:
#     # r_zoom_f    = f"{code_dir}/zoom.R",
#     rmd_f       = rmd_dir + '/' + SHORT_TAG + '_zoom' + '_{zoom_name}_{zoom_res}.Rmd',
#     html_f      = docs_dir + '/' + SHORT_TAG + '_zoom' + '_{zoom_name}_{zoom_res}.html'
#   params:
#     zoom_name   = '{zoom_name}',
#     zoom_res    = '{zoom_res}'
#   threads: 1
#   retries: 5
#   resources:
#     mem_mb      = lambda wildcards, attempt: attempt * 8192
#   run:
#     # copy R code over
#     print('copying relevant R files over')
#     import shutil
#     shutil.copyfile('scripts/zoom.R', f"{code_dir}/zoom.R")

#     # define what we will substitute in
#     print('setting up template')
#     sub_dict    = {
#       'YOUR_NAME':        YOUR_NAME,
#       'AFFILIATION':      AFFILIATION,
#       'SHORT_TAG':        SHORT_TAG,
#       'DATE_STAMP':       DATE_STAMP,
#       'threads':          threads,
#       'zoom_name':        params.zoom_name,
#       'zoom_res':         params.zoom_res,
#       'SHORT_TAG':        SHORT_TAG,
#       'meta_f':           os.path.relpath(METADATA_F, PROJ_DIR),
#       'meta_vars_ls':     'c("' + '", "'.join(map(str, METADATA_VARS)) + '")',
#       'gtf_dt_f':         AF_GTF_DT_F,
#       'sce_sub_f':        os.path.relpath(input.zoom_sce_sub_f, PROJ_DIR),
#       'hmny_f':           os.path.relpath(input.zoom_hmny_f, PROJ_DIR),
#       'pb_f':             os.path.relpath(input.zoom_pb_f, PROJ_DIR),
#       'mkrs_f':           os.path.relpath(input.zoom_mkrs_f, PROJ_DIR),
#       'hvgs_f':           os.path.relpath(input.zoom_hvgs_f, PROJ_DIR),
#       'canon_f':          MKR_CANON_F,
#       'fgsea_go_bp_f':    os.path.relpath(input.zoom_fgsea_go_bp_f, PROJ_DIR),
#       'fgsea_go_cc_f':    os.path.relpath(input.zoom_fgsea_go_cc_f, PROJ_DIR),
#       'fgsea_go_mf_f':    os.path.relpath(input.zoom_fgsea_go_mf_f, PROJ_DIR),
#       'fgsea_paths_f':    os.path.relpath(input.zoom_fgsea_paths_f, PROJ_DIR),
#       'fgsea_hlmk_f':     os.path.relpath(input.zoom_fgsea_hlmk_f, PROJ_DIR),
#       'INT_DBL_CL_PROP':  INT_DBL_CL_PROP,
#       'INT_EXC_REGEX':    INT_EXC_REGEX,
#       'MKR_NOT_OK_RE':    MKR_NOT_OK_RE,
#       'MKR_MIN_CPM_MKR':  MKR_MIN_CPM_MKR,
#       'MKR_MIN_CELLS':    MKR_MIN_CELLS,
#       'MKR_GSEA_CUT':     MKR_GSEA_CUT
#       }

#     # make and render Rmd file
#     template_f  = 'templates/zoom.Rmd.template'
#     print('rendering template')
#     render_html(PROJ_DIR, RLIBS_DIR, template_f, sub_dict, output.rmd_f)


# # render_html_empties
# rule render_html_empties:
#   input:
#     r_int_f     = f'{code_dir}/{SHORT_TAG}06_integration.R',
#     guesses_f   = lbl_dir + '/xgboost_guesses_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
#     empty_csv_f = f'{lbl_dir}/empties_{FULL_TAG}_{DATE_STAMP}.csv',
#     empty_gs_ls = expand( [lbl_dir + '/empty_genes' + '/empty_genes_' + FULL_TAG + '_{subset}_' + DATE_STAMP + '.txt.gz'], \
#       subset = None if LBL_SCE_SUBSETS is None else [*LBL_SCE_SUBSETS] )
#   output:
#     r_pb_f      = f'{code_dir}/pseudobulk_and_empties.R',
#     rmd_f       = f'{rmd_dir}/{SHORT_TAG}_empties.Rmd',
#     html_f      = f'{docs_dir}/{SHORT_TAG}_empties.html'
#   threads: 1
#   retries: 5
#   resources:
#     mem_mb      = lambda wildcards, attempt: attempt * 4096
#   run:
#     # copy R code over
#     print('copying relevant R files over')
#     import shutil
#     shutil.copyfile('scripts/pseudobulk_and_empties.R', output.r_pb_f)

#     # define what we will substitute in
#     print('setting up template')
#     sub_dict    = {
#       'YOUR_NAME':      YOUR_NAME,
#       'AFFILIATION':    AFFILIATION,
#       'SHORT_TAG':      SHORT_TAG,
#       'DATE_STAMP':     DATE_STAMP,
#       'threads':        threads,
#       'guesses_f':      os.path.relpath(input.guesses_f, PROJ_DIR),
#       'empty_csv_f':    os.path.relpath(input.empty_csv_f, PROJ_DIR),
#       'LBL_XGB_F':      LBL_XGB_F,
#       'LBL_SEL_RES_CL': LBL_SEL_RES_CL,
#       'LBL_MIN_PRED':   LBL_MIN_PRED,
#       'LBL_MIN_CL_PROP':LBL_MIN_CL_PROP,
#       'LBL_MIN_CL_SIZE':LBL_MIN_CL_SIZE
#       }

#     # make and render Rmd file
#     template_f  = 'templates/empties.Rmd.template'
#     print('rendering template')
#     render_html(PROJ_DIR, RLIBS_DIR, template_f, sub_dict, output.rmd_f)

