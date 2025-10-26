
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
  temp_ls = get_sub_ls(rule_name, proj_dir, ...)

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

get_sub_ls <- function(rule = c('mapping', 'multiplexing', 'ambient', 'qc', 'hvg', 'integration', 
  'markers', 'cell_labels', 'zoom', 'pb_empties'), proj_dir, ...) {
  # get arguments
  sel_rule = match.arg(rule)
  add_args = list(...)
  add_args_names = names(add_args)

  # check if all extra args for a specific rule are present
  if (sel_rule == 'ambient') {
    req_names = c('your_name','affiliation', 'short_tag', 'run_stats_f',
      'threads','run_var', 'ambient_method', 'date_stamp', 'runs_str', 'cb_prop_max_kept')

    assert_that(all(req_names %in% add_args_names))
    
    if (add_args[['ambient_method']] != "none") {
      plot_removed_title = "## How many reads were removed as ambient?"
      plot_removed_txt = "Plots show what proportion of reads were removed from all barcodes called as cells."
      spl_umis_pl_txt = paste0("The plots show the relationship between the number of UMIs and the percentage of spliced reads, ", 
        "separated by sample, both before and after ambient RNA removal. Ambient RNA is typically characterized by a high proportion of spliced reads, ", 
        "making the spliced percentage a useful metric for evaluating the effectiveness of an ambient RNA removal method. ", 
        "By examining changes in spliced percentages, we can assess how well the method performed in reducing ambient RNA contamination.")
    } else{
      spl_umis_pl_txt = paste0("The plot shows the relationship between the number of UMIs and the percentage of spliced reads, separated by sample. ", 
        "It can help diagnose potential ambient RNA contamination, as ambient RNA is typically characterized by a high proportion of spliced reads.")
      plot_removed_title = ""
      plot_removed_txt = ""
    }

    if (add_args[['ambient_method']] == "cellbender") {
      tbl_removed_title = "## Samples excluded by ambient removal step"
      tbl_removed_txt = paste0("CellBender includes all barcodes in the analysis up to the `total_droplets` threshold", 
        " covering both cell-containing droplets and empty droplets.", " If CellBender calls the majority of these included droplets as cells,", 
        " it may indicate an underlying issue. This typically occurs in low-quality samples where cell-containing barcodes and empty droplets cannot", 
        " be clearly distinguished in the barcode rank plot. The table below shows the proportion of included droplets that were classified as cells by CellBender.",
        " Samples where this proportion exceeded ", add_args[['CB_PROP_MAX_KEPT']]*100,  "% were excluded from further analysis.")
    } else{
      tbl_removed_title = ""
      tbl_removed_txt = ""
    }

    params_ls = c(
      add_args[setdiff(req_names, 'cb_prop_max_kept')],
      list(spl_umis_pl_txt    = spl_umis_pl_txt, 
           plot_removed_title = plot_removed_title, 
           plot_removed_txt   = plot_removed_txt, 
           tbl_removed_title  = tbl_removed_title, 
           tbl_removed_txt    = tbl_removed_txt))

  } else if (sel_rule == 'mapping') {
    req_names = c('your_name', 'affiliation', 'short_tag', 
      'date_stamp', 'runs_str','ambient_method','run_var',
      'af_dir', 'af_rna_dir')

    assert_that(all(req_names %in% add_args_names))

    params_ls = add_args

  } else if (sel_rule == 'multiplexing') {
    req_names = c('your_name', 'affiliation', 'short_tag', 
      'date_stamp', 'runs_str','ambient_method','metadata_f',
      'run_var', 'af_dir', 'demux_dir')

    assert_that(all(req_names %in% add_args_names))

    params_ls = add_args

  } else if (sel_rule == 'qc') {
    req_names = c('your_name', 'affiliation', 'short_tag', 
      'date_stamp', 'threads', 'metadata_f', 'qc_dt_f', 'cuts_f',
      'min_cells', 'qc_hard_min_counts', 'qc_hard_min_feats', 'qc_hard_max_mito')
    setdiff(req_names, add_args_names)
    assert_that(all(req_names %in% add_args_names))

    params_ls = add_args

  } else if (sel_rule == 'hvg') {
    req_names = c('your_name', 'affiliation', 'short_tag', 
      'date_stamp', 'threads', 'hvgs_f', 'empty_gs_f', 'pb_empty_f')
    
    assert_that(all(req_names %in% add_args_names))

    params_ls = add_args

  } else if (sel_rule == 'integration') {
    req_names = c('your_name', 'affiliation', 'short_tag', 
      'date_stamp', 'threads', 'qc_dt_f', 'integration_f', 'int_embedding', 'demux_type', 
      'int_res_ls_str', 'int_dbl_cl_prop')

    assert_that(all(req_names %in% add_args_names))

    params_ls = add_args[req_names]

  } else if (sel_rule == 'markers') {
    req_names = c('your_name', 'affiliation', 'short_tag', 
      'date_stamp', 'threads', 'metadata_f','meta_vars_ls',
      'gtf_dt_f', 'integration_f', 'pb_f', 'mkrs_f', 'hvgs_f', 'ambient_f',
      'fgsea_go_bp_f', 'fgsea_go_cc_f', 'fgsea_go_mf_f','fgsea_paths_f', 'fgsea_hlmk_f',
      'mkr_sel_res', 'custom_mkr_names', 'custom_mkr_paths',
      'mkr_not_ok_re', 'mkr_min_cpm_mkr', 'mkr_min_cells', 'mkr_gsea_cut', 'species')

    assert_that(all(req_names %in% add_args_names))

    metadata_vars = add_args[['meta_vars_ls']] %>% 
    str_split(pattern = ",") %>% unlist()
    if (length(metadata_vars) > 0){
      meta_bars_title = "## Cluster splits by metadata variables"
      meta_bars_txt   = paste0("For each cluster the proportion of cells coming from samples associated with", 
       " specific values of ", paste(metadata_vars, collapse = ', ') %>% stri_replace_last_fixed(",", " and"), ' is shown.')
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
  
    if(add_args[['species']] %in% c('human_2024', 'human_2020', 'mouse_2024', 'mouse_2020')){
      fgsea_title = "## GSEA characterisation of clusters{.tabset}"
      fgsea_txt   = paste0("Gene Set Enrichment Analysis (GSEA) was performed on marker genes for each cluster, using log fold change as the ranking variable.", 
      " The top 10 pathways, grouped into five categories and selected based on a significance threshold of 0.05, are displayed for each cluster.")
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
    req_names = c('your_name', 'affiliation', 'short_tag', 
      'date_stamp', 'threads' ,'guesses_f',
      'lbl_tissue', 'lbl_xgb_f', 'lbl_xgb_cls_f', 'lbl_sel_res_cl', 'lbl_min_pred',
      'lbl_min_cl_prop', 'lbl_min_cl_size')
      
     assert_that(all(req_names %in% add_args_names))

      if (add_args[["lbl_tissue"]] == 'human_cns') {
        train_data_str = "whole brain human single nuclei atlas (siletti et al. 2023)"
      } else if (add_args[["lbl_tissue"]] == 'mouse_cns') {
        train_data_str = "whole brain mouse single nuclei atlas (Langlieb et al. 2023)"
      } else {
        train_data_str = "insert name of study here"
      }

      params_ls = add_args
      params_ls = c(params_ls, train_data_str = train_data_str)

  } else if (sel_rule == 'zoom') {
    req_names = c('your_name', 'affiliation', 'short_tag', 'date_stamp', 
      'threads', 'zoom_dir', 'zoom_name', 'metadata_f', 'meta_vars_ls',
      'gtf_dt_f', 'qc_f', 'cell_hvgs_f', 'int_f', 'pb_f', 'pb_hvgs_f', 'mkrs_f', 'empty_gs_f', 'pb_empty_f', 
      'fgsea_go_bp_f','fgsea_go_cc_f', 'fgsea_go_mf_f', 'fgsea_paths_f', 'fgsea_hlmk_f', 'int_res_ls',
      'custom_mkr_names', 'custom_mkr_paths', 'mkr_not_ok_re', 'mkr_min_cpm_mkr', 'mkr_sel_res',
      'mkr_min_cells', 'mkr_gsea_cut', 'species')
    
    assert_that(all(req_names %in% add_args_names))

    metadata_vars = add_args[['meta_vars_ls']] %>% 
    str_split(pattern = ",") %>% unlist()
    if (length(metadata_vars) > 0){
      meta_bars_title = "### Cluster splits by metadata variables"
      meta_bars_txt   = paste0("For each cluster (resolution: ", add_args[['mkr_sel_res']], ")", " the proportion of cells coming from samples associated with", 
       " specific values of ", paste(metadata_vars, collapse = ', ') %>% stri_replace_last_fixed(",", " and"), ' is shown.')
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

    if(add_args[['species']] %in% c('human_2024', 'human_2020', 'mouse_2024', 'mouse_2020')){
      fgsea_title = "### GSEA characterisation of clusters{.tabset}"
      fgsea_txt   = paste0("Gene Set Enrichment Analysis (GSEA) was performed on marker genes for each cluster, using log fold change as the ranking variable.", 
      " The top 10 pathways, grouped into five categories and selected based on a significance threshold of 0.05, are displayed for each cluster.")
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
    req_names = c('your_name', 'affiliation', 'short_tag', 
      'date_stamp', 'threads', 'guesses_f', 'empty_csv_f',
      'lbl_xgb_f', 'lbl_sel_res_cl', 'lbl_min_pred', 'lbl_min_cl_prop',
      'lbl_min_cl_size', 'lbl_min_cl_size')
    
    assert_that(all(req_names %in% add_args_names))
  }
  params_ls$proj_dir = proj_dir

  return(params_ls)
}
