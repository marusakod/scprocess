
suppressPackageStartupMessages({
    library('stringr')
    library('assertthat')
    library('workflowr')
    library('glue')
    library('stringi')
})

# function that replaces placeholder strings in template Rmds and renders html reports
render_html <- function(rule_name, proj_dir, temp_f, rmd_f, ...) {
  setwd(proj_dir)

  # get list with all values that need to be replaced in the template
  temp_ls = get_sub_ls(rule_name, ...)

  # make Rmd file
  message('Creating Rmd file from template ', temp_f)
  make_rmd_from_temp(temp_f, temp_ls, rmd_f)

  # convert to html
  message('Rendering html')
  workflowr::wflow_build(
    files = rmd_f,
    view = FALSE,
    verbose = TRUE,
    delete_cache = TRUE
    )
}

make_rmd_from_temp <- function(temp_f, temp_ls, rmd_f) {
  if (!file.exists(rmd_f)) {
    # read remplate file
    temp_str = readLines(temp_f, warn = FALSE)

    # convert into one long string
    temp_str = paste(temp_str, collapse = "\n")

    # substitute placeholders
    rmd_str = glue(temp_str, .envir = as.environment(temp_ls), .open = "${", .close = "}")

    # write to rmd file
    writeLines(rmd_str, rmd_f)
  }
}

get_sub_ls <- function(rule = c('af', 'multiplexing', 'ambient', 'qc', 'hvg', 'integration', 
  'markers', 'cell_labels', 'zoom', 'pb_empties'), ...) {

  # get arguments
  sel_rule = match.arg(rule)
  add_args = list(...)
  add_args_names = names(add_args)

  # check if all extra args for a specific rule are present
  if (sel_rule == 'ambient') {
    req_names = c('YOUR_NAME','AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 'smpl_stats_f',
      'threads','SAMPLE_VAR', 'AMBIENT_METHOD', 'DATE_STAMP', 'RUNS_STR', 'CELLBENDER_PROP_MAX_KEPT')

    assert_that(all(req_names %in% add_args_names))
    
    if (add_args[['AMBIENT_METHOD']] != "none") {
      plot_removed_title = "## How many reads were removed as ambient?"
      plot_removed_txt = "Plots show what proportion of reads were removed from all barcodes called as cells."
    } else{
      plot_removed_title = ""
      plot_removed_txt = ""
    }

    if (add_args[['AMBIENT_METHOD']] == "cellbender") {
      tbl_removed_title = "## Samples excluded by ambient removal step"
      tbl_removed_txt = paste0("CellBender includes all barcodes in the analysis up to the `total_droplets` threshold", 
        " covering both cell-containing droplets and empty droplets.", " If CellBender calls the majority of these included droplets as cells,", 
        " it may indicate an underlying issue. This typically occurs in low-quality samples where cell-containing barcodes and empty droplets cannot", 
        " be clearly distinguished in the barcode rank plot. The table below shows the proportion of included droplets that were classified as cells by CellBender.",
        " Samples where this proportion exceeded ", add_args[['CELLBENDER_PROP_MAX_KEPT']]*100,  "% were excluded from further analysis.")
    } else{
      tbl_removed_title = ""
      tbl_removed_txt = ""
    }

    params_ls = c(
      add_args[setdiff(req_names, 'CELLBENDER_PROP_MAX_KEPT')],
      list(plot_removed_title = plot_removed_title, 
           plot_removed_txt   = plot_removed_txt, 
           tbl_removed_title  = tbl_removed_title, 
           tbl_removed_txt    = tbl_removed_txt))

  } else if (sel_rule == 'af') {
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'RUNS_STR','AMBIENT_METHOD','SAMPLE_VAR',
      'af_dir', 'af_rna_dir')

    assert_that(all(req_names %in% add_args_names))

    params_ls = add_args

  } else if (sel_rule == 'multiplexing') {
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'RUNS_STR','AMBIENT_METHOD','METADATA_F',
      'SAMPLE_VAR', 'af_dir', 'demux_dir')

    assert_that(all(req_names %in% add_args_names))

    params_ls = add_args

  } else if (sel_rule == 'qc') {
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'threads', 'meta_f', 'qc_dt_f',
      'QC_HARD_MIN_COUNTS', 'QC_HARD_MIN_FEATS', 'QC_HARD_MAX_MITO',
      'QC_MIN_COUNTS', 'QC_MIN_FEATS', 'QC_MIN_MITO', 'QC_MAX_MITO',
      'QC_MIN_SPLICE', 'QC_MAX_SPLICE', 'QC_MIN_CELLS')

    assert_that(all(req_names %in% add_args_names))

    params_ls = add_args

  } else if (sel_rule == 'hvg') {
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'threads', 'hvgs_f', 'empty_gs_f', 'pb_empty_f')
    
    assert_that(all(req_names %in% add_args_names))

    params_ls = add_args

  } else if (sel_rule == 'integration') {
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'threads', 'qc_dt_f', 'integration_f', 'INT_REDUCTION', 'DEMUX_TYPE', 
      'INT_RES_LS', 'INT_DBL_CL_PROP')

    assert_that(all(req_names %in% add_args_names))

    params_ls = add_args[req_names]

  } else if (sel_rule == 'markers') {
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'threads', 'meta_f','meta_vars_ls',
      'gtf_dt_f', 'integration_f', 'pb_f', 'mkrs_f', 'hvgs_f', 'ambient_f',
      'fgsea_go_bp_f', 'fgsea_go_cc_f', 'fgsea_go_mf_f','fgsea_paths_f', 'fgsea_hlmk_f',
      'MKR_SEL_RES', 'CUSTOM_MKR_NAMES', 'CUSTOM_MKR_PATHS',
      'MKR_NOT_OK_RE', 'MKR_MIN_CPM_MKR', 'MKR_MIN_CELLS', 'MKR_GSEA_CUT', 'SPECIES')

    assert_that(all(req_names %in% add_args_names))

    metadata_vars = add_args[['meta_vars_ls']] %>% 
    str_split(pattern = ",") %>% unlist()
    if (length(metadata_vars) > 0){
      meta_bars_title = "## Cluster splits by metadata variables"
      meta_bars_txt   = paste0("For each cluster the proportion of cells coming from samples associated with", 
       " specific values of ", paste(metadata_vars, collapse = ', ') %>% stri_replace_last_fixed(",", " and"), '.')
      meta_umap_title = "## Metadata variables over UMAP{.tabset}"
      meta_umap_txt   = paste0("The plot shows a binned UMAP with facets corresponding to specific values of ", 
       (if (length(metadata_vars) == 1) print(metadata_vars) else print("different metadata variables")),
       " which allows the evaluation of whether cells sharing certain annotations are particularly abundant in some clusters.")
   }else{
      meta_bars_title = ""
      meta_bars_txt   = ""
      meta_umap_title = ""
      meta_umap_txt   = ""
   }
  
    if(add_args[['SPECIES']] %in% c('human_2024', 'human_2020', 'mouse_2024', 'mouse_2020')){
      fgsea_title = "## GSEA characterisation of clusters{.tabset}"
      fgsea_txt   = paste0("GSEA was performed on marker genes for each cluster, using log fold change as the ranking variable.", 
       " Top x pathways grouped into 5 categories with some threshold are shown for each cluster.")
    }else{
      fgsea_title = ""
      fgsea_txt   = ""
    }
  
    params_ls = c(
      add_args[req_names],
      list(meta_bars_title = meta_bars_title, 
           meta_bars_txt   = meta_bars_txt, 
           meta_umap_title = meta_umap_title, 
           meta_umap_txt   = meta_umap_txt, 
           fgsea_title     = fgsea_title,
           fgsea_txt       = fgsea_txt))

  } else if (sel_rule == 'cell_labels') {
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'threads' ,'guesses_f',
      'LBL_TISSUE', 'LBL_XGB_F', 'LBL_XGB_CLS_F', 'LBL_SEL_RES_CL', 'LBL_MIN_PRED',
      'LBL_MIN_CL_PROP', 'LBL_MIN_CL_SIZE')
      
     assert_that(all(req_names %in% add_args_names))

      if (add_args[["LBL_TISSUE"]] == 'human_cns') {
        train_data_str = "whole brain human single nuclei atlas (Siletti et al. 2023)"
      } else if (add_args[["LBL_TISSUE"]] == 'mouse_cns') {
        train_data_str = "whole brain mouse single nuclei atlas (Langlieb et al. 2023)"
      } else {
        train_data_str = "insert name of study here"
      }

      params_ls = add_args
      params_ls = c(params_ls, train_data_str = train_data_str)

  } else if (sel_rule == 'zoom') {
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 'DATE_STAMP', 
      'threads', 'zoom_dir', 'zoom_name', 'meta_f', 'meta_vars_ls',
      'gtf_dt_f', 'qc_f', 'int_f', 'pb_f', 'mkrs_f', 'mkrs_hvgs_f', 'hvgs_f', 'empty_gs_f', 'pb_empty_f', 
      'fgsea_go_bp_f','fgsea_go_cc_f', 'fgsea_go_mf_f', 'fgsea_paths_f', 'fgsea_hlmk_f',
      'CUSTOM_MKR_NAMES', 'CUSTOM_MKR_PATHS', 'MKR_NOT_OK_RE', 'MKR_MIN_CPM_MKR', 'MKR_SEL_RES',
      'MKR_MIN_CELLS', 'MKR_GSEA_CUT', 'SPECIES')
    
    assert_that(all(req_names %in% add_args_names))

    metadata_vars = add_args[['meta_vars_ls']] %>% 
    str_split(pattern = ",") %>% unlist()
    if (length(metadata_vars) > 0){
      meta_bars_title = "### Cluster splits by metadata variables"
      meta_bars_txt   = paste0("For each cluster the proportion of cells coming from samples associated with", 
       " specific values of ", paste(metadata_vars, collapse = ', ') %>% stri_replace_last_fixed(",", " and"), '.')
      meta_umap_title = "### Metadata variables over UMAP{.tabset}"
      meta_umap_txt   = paste0("The plot shows a binned UMAP with facets corresponding to specific values of ", 
       (if (length(metadata_vars) == 1) print(metadata_vars) else print("different metadata variables")),
       " which allows the evaluation of whether cells sharing certain annotations are particularly abundant in some clusters.")
   }else{
      meta_bars_title = ""
      meta_bars_txt   = ""
      meta_umap_title = ""
      meta_umap_txt   = ""
   }

    if(add_args[['SPECIES']] %in% c('human_2024', 'human_2020', 'mouse_2024', 'mouse_2020')){
      fgsea_title = "### GSEA characterisation of clusters{.tabset}"
      fgsea_txt   = paste0("GSEA was performed on marker genes for each cluster, using log fold change as the ranking variable.", 
       " Top x pathways grouped into 5 categories with some threshold are shown for each cluster.")
    }else{
      fgsea_title = ""
      fgsea_txt   = ""
    }

    params_ls = c(
      add_args[req_names],
      list(meta_bars_title = meta_bars_title, 
           meta_bars_txt   = meta_bars_txt, 
           meta_umap_title = meta_umap_title, 
           meta_umap_txt   = meta_umap_txt, 
           fgsea_title     = fgsea_title,
           fgsea_txt       = fgsea_txt))

  } else if (sel_rule == 'pb_empties') {
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'threads', 'guesses_f', 'empty_csv_f',
      'LBL_XGB_F', 'LBL_SEL_RES_CL', 'LBL_MIN_PRED', 'LBL_MIN_CL_PROP',
      'LBL_MIN_CL_SIZE', 'LBL_MIN_CL_SIZE')
    
    assert_that(all(req_names %in% add_args_names))
    params_ls = add_args[req_names]
  }

  return(params_ls)
}
