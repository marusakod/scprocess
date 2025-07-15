
suppressPackageStartupMessages({
    library('stringr')
    library('assertthat')
    library('workflowr')
    library('glue')
})

# function that replaces placeholder strings in template Rmds and renders html reports
render_reports <- function(rule_name, proj_dir, temp_f, rmd_f, ...){
  setwd(proj_dir)

  # get list with all values that need to be replaced in the template
  temp_ls = get_sub_ls(rule_name, ...)

  # make Rmd file
  message('Creating Rmd file from template ', temp_f)
  make_rmd_from_temp(temp_f, temp_ls, rmd_f)

  message('Rendering html')

  workflowr::wflow_build(
    files = rmd_f,
    view = FALSE,
    verbose = TRUE,
    delete_cache = TRUE
    )
}

make_rmd_from_temp <- function(temp_f, temp_ls, rmd_f){
  if(!file.exists(rmd_f)){
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

get_sub_ls <- function(rule = c('af', 'multiplexing', 'cellbender', 'qc', 'integration', 
  'markers', 'cell_labels', 'zoom', 'pb_empties'), ...){

  # get arguments
  sel_rule = match.arg(rule)
  add_args = list(...)
  add_args_names = names(add_args)

  # check if all extra args for a specific rule are present
  if(sel_rule == 'cellbender'){
    req_names = c('YOUR_NAME','AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 'stats_f',
      'CELLBENDER_PROP_MAX_KEPT','DATE_STAMP', 'RUNS_STR', 'af_dir')
    
    assert_that(all(req_names %in% add_args_names))
    
    params_ls = add_args
  }else if(sel_rule == 'af'){
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'RUNS_STR','AMBIENT_METHOD','SAMPLE_VAR',
      'af_dir', 'af_rna_dir')

    assert_that(all(req_names %in% add_args_names))

    params_ls = add_args

  }else if(sel_rule == 'multiplexing'){
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'RUNS_STR','AMBIENT_METHOD','METADATA_F',
      'SAMPLE_VAR', 'af_dir', 'demux_dir')

    assert_that(all(req_names %in% add_args_names))

    params_ls = add_args

  }else if(sel_rule == 'qc'){
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'threads', 'meta_f', 'qc_dt_f',
      'QC_HARD_MIN_COUNTS', 'QC_HARD_MIN_FEATS', 'QC_HARD_MAX_MITO',
      'QC_MIN_COUNTS', 'QC_MIN_FEATS', 'QC_MIN_MITO', 'QC_MAX_MITO',
      'QC_MIN_SPLICE', 'QC_MAX_SPLICE', 'QC_MIN_CELLS')

    assert_that(all(req_names %in% add_args_names))

    params_ls = add_args

  }else if(sel_rule == 'integration'){
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'threads', 'qc_dt_f', 'integration_f', 
      'INT_RES_LS', 'INT_DBL_CL_PROP')

    assert_that(all(req_names %in% add_args_names))

    params_ls = add_args[req_names]

  }else if(sel_rule == 'markers'){
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'threads', 'meta_f','meta_vars_ls',
      'gtf_dt_f', 'integration_f', 'pb_f', 'mkrs_f', 'hvgs_f', 'fgsea_go_bp_f',
      'fgsea_go_cc_f', 'fgsea_go_mf_f','fgsea_paths_f', 'fgsea_hlmk_f',
      'MKR_SEL_RES', 'CUSTOM_MKR_NAMES', 'CUSTOM_MKR_PATHS',
      'MKR_NOT_OK_RE', 'MKR_MIN_CPM_MKR', 'MKR_MIN_CELLS', 'MKR_GSEA_CUT', 'SPECIES')

    assert_that(all(req_names %in% add_args_names))

    # based on species determine whether code chunks with gsea results shoudl be eval or not
    if(add_args[['SPECIES']] %in% c('human_2024', 'human_2020', 'mouse_2024', 'mouse_2020')){
      eval_fgsea = TRUE
    }else{
      eval_fgsea = FALSE
    }

    params_ls = c(add_args[setdiff(req_names, 'SPECIES')], list(eval_fgsea = eval_fgsea))

  }else if(sel_rule == 'cell_labels'){
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'threads' ,'guesses_f', 'CUSTOM_LABELS_F', 'MKR_SEL_RES',
      'LBL_TISSUE', 'LBL_XGB_F', 'LBL_SEL_RES_CL', 'LBL_MIN_PRED',
      'LBL_MIN_CL_PROP', 'LBL_MIN_CL_SIZE')
      
    assert_that(all(req_names %in% add_args_names))

    if(add_args[['CUSTOM_LABELS_F']] == ""){
      eval_boost = TRUE
      if(add_args[["LBL_TISSUE"]] == 'human_cns'){
        train_data_str = "whole brain human single nuclei atlas (Siletti et al. 2023)"
      }else if(add_args[["LBL_TISSUE"]] == 'mouse_cns'){
        train_data_str = "whole brain mouse singlce nuclei atlas (Langlieb et al. 2023)"
      }else if(add_args[["LBL_TISSUE"]] == 'human pbmc'){
        train_data_str = "insert name of study here"
      }else{
        train_data_str = "insert name of study here"
      }

      main_lbl_txt = paste0("Annotation of cell types was performed using",
        " an XGBoost classifier trained on annotated data from the ", train_data_str)
      hmap_title = "`harmony`clusters vs `XGBoost`-predicted cell types{.tabset}"
      hmap_txt = "heatmap with xgboost annotations vs harmony clusters"
      umap_title = "XGBoost`-predicted cell types over UMAP{.tabset}"
      umap_txt = paste0("For each cluster the most frequent prediction labels are shown")

    }else{
      eval_boost = FALSE
      main_lbl_txt = paste0("Annotation of cell types is based on a user defined table")
      hmap_title = "`harmony` clusters vs user defined cell type annotations"
      hmap_txt = "harmony clusters vs user defined cell type annotations"
      umap_title = "user defined cell types over UMAP"
      umap_txt = paste0("User defined annotations are displayed ",
                        "over a UMAP, together with clusters at resolution ",
                        add_args[["MKR_SEL_RES"]])
    }

      params_ls = add_args
      params_ls = c(
      params_ls,
      list(eval_boost   = eval_boost,
            main_lbl_txt = main_lbl_txt,
            hmap_title   = hmap_title,
            hmap_txt     = hmap_txt,
            umap_title   = umap_title,
            umap_txt     = umap_txt)
    )

  }else if(sel_rule == 'zoom'){
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 'DATE_STAMP', 
      'threads', 'zoom_name', 'zoom_res', 'meta_f', 'meta_vars_ls', # meta_vars_ls should be made one string
      'gtf_dt_f', 'sce_sub_f', 'hmny_f', 'pb_f', 'mkrs_f', 'hvgs_f', 'canon_f', 
      'fgsea_go_bp_f','fgsea_go_cc_f', 'fgsea_go_mf_f', 'fgsea_paths_f', 'fgsea_hlmk_f',
      'INT_DBL_CL_PROP', 'INT_EXC_REGEX', 'MKR_NOT_OK_RE', 'MKR_MIN_CPM_MKR', 
      'MKR_MIN_CELLS', 'MKR_GSEA_CUT', 'SPECIES')
    assert_that(all(req_names %in% add_args_names))

    if(add_args[['SPECIES']] %in% c('human_2024', 'human_2020', 'mouse_2024', 'mouse_2020')){
      eval_fgsea = TRUE
    }else{
      eval_fgsea = FALSE
    }
    params_ls = c(add_args[setdiff(req_names, 'SPECIES')], list(eval_fgsea = eval_fgsea))

  }else if(sel_rule == 'pb_empties'){
    req_names = c('YOUR_NAME', 'AFFILIATION', 'SHORT_TAG', 'PROJ_DIR', 
      'DATE_STAMP', 'threads', 'guesses_f', 'empty_csv_f',
      'LBL_XGB_F', 'LBL_SEL_RES_CL', 'LBL_MIN_PRED', 'LBL_MIN_CL_PROP',
      'LBL_MIN_CL_SIZE', 'LBL_MIN_CL_SIZE')
    
    assert_that(all(req_names %in% add_args_names))
    params_ls = add_args[req_names]
  }

  return(params_ls)
}
