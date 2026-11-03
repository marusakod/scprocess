#!/usr/bin/env Rscript
# Classifier results have their own report: zoom and join templates substitute fully
# without classifier placeholders, and the standalone template accepts zoom inputs.
# Usage: Rscript --vanilla tests/test_xgboost_report_separation.R
source('scripts/render_htmls.R')

# substitute one template with placeholder values for its required arguments
substitute_template <- function(rule, template_f, args) {
  params = do.call(get_sub_ls, c(list(rule = rule, proj_dir = tempdir()), args))
  rmd_f  = tempfile(fileext = '.Rmd')
  make_rmd_from_temp(rule, template_f, params, rmd_f)
  readLines(rmd_f)
}

# every required name gets a placeholder; a few need meaningful values
dummy_args <- function(arg_names, overrides) {
  args = setNames(as.list(rep('x', length(arg_names))), arg_names)
  args[names(overrides)] = overrides
  args
}

zoom_names = c('your_name', 'affiliation', 'short_tag', 'date_stamp',
  'threads', 'zoom_dir', 'zoom_name', 'metadata_f', 'meta_vars_ls',
  'gtf_dt_f', 'qc_f', 'r_label_celltypes_f', 'zoom_lbls_f', 'zoom_report_lbls_f',
  'zoom_original_lbls_f', 'zoom_lbls_col', 'zoom_label_hi_res_cl', 'zoom_sel_labels',
  'zoom_qc_min_counts', 'zoom_qc_max_counts', 'zoom_qc_min_feats', 'zoom_qc_max_feats',
  'zoom_qc_min_mito', 'zoom_qc_max_mito', 'zoom_qc_min_splice', 'zoom_qc_max_splice',
  'cell_hvgs_f', 'int_f', 'pb_plot_data_f', 'pb_hvgs_f', 'mkrs_f', 'empty_gs_f', 'pb_empty_f',
  'fgsea_go_bp_f', 'fgsea_go_cc_f', 'fgsea_go_mf_f', 'int_res_ls',
  'custom_mkr_names', 'custom_mkr_paths', 'mkr_not_ok_re', 'mkr_min_cpm_mkr', 'mkr_sel_res',
  'mkr_min_cells', 'mkr_gsea_var', 'mkr_gsea_cut', 'ref_txome', 'batch_var', 'do_gsea')
zoom_lines = substitute_template('zoom', 'resources/rmd_templates/zoom.Rmd.template',
  dummy_args(zoom_names, list(threads = 1L, meta_vars_ls = '', mkr_sel_res = '0.5',
    do_gsea = 'FALSE', ref_txome = 'none')))

join_names = c('your_name', 'affiliation', 'join_name', 'join_tag',
  'join_int_dir', 'join_mkr_dir', 'ref_txome', 'mkr_sel_res',
  'int_res_ls', 'scprocess_dir', 'date_stamp', 'metadata_vars',
  'custom_mkr_names', 'custom_mkr_paths',
  'label_f_ls', 'labeller_ls', 'model_ls', 'hi_res_cl_ls', 'min_cl_prop_ls',
  'mkr_min_cpm_mkr', 'mkr_min_cells', 'mkr_gsea_cut',
  'integration_f', 'qc_f', 'sample_meta_f', 'mkrs_f', 'pb_hvgs_f', 'pb_plot_data_f',
  'fgsea_go_bp_f', 'fgsea_go_cc_f', 'fgsea_go_mf_f')
join_lines = substitute_template('join', 'resources/rmd_templates/join.Rmd.template',
  dummy_args(join_names, list()))

# zoom and join reports: fully substituted, with no classifier content
for (lines in list(zoom_lines, join_lines)) {
  stopifnot(!any(grepl('${', lines, fixed = TRUE)),
    !any(grepl('xgb|xgboost', lines, ignore.case = TRUE)))
}

# standalone classifier report with the arguments passed for a project zoom
xgb_args = list(your_name = 'test', affiliation = 'test', short_tag = 'proj_zoom_neurons',
  ref_tag = 'neurons_v1', predictions_f = 'p.csv.gz', importance_f = 'i.csv',
  pseudobulk_f = 'pb.h5ad', model_card_f = 'card.json', label_counts_f = 'counts.csv',
  fulldata_predictions_f = 'full.csv.gz', integration_f = 'int.csv.gz',
  min_cells = '10', n_cores = '1')
xgb_lines = substitute_template('train_xgboost',
  'resources/rmd_templates/train_xgboost.Rmd.template', xgb_args)
stopifnot(!any(grepl('${', xgb_lines, fixed = TRUE)),
  any(grepl('full.csv.gz', xgb_lines, fixed = TRUE)))

message('Zoom and join reports omit classifier results; standalone report substitutes zoom inputs.')
