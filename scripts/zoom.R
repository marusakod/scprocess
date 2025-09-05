
suppressPackageStartupMessages({
  library("data.table")
  library("magrittr")
  library("assertthat")
  library("SingleCellExperiment")
  library("stringr")
})

run_zoom_integration <- function(hvg_mat_f, coldata_f, smpl_stats_f, demux_type, exclude_mito,
   reduction, n_dims, cl_method, theta, res_ls_concat, integration_f, batch_var, n_cores = 4) {
  
  exclude_mito = as.logical(exclude_mito)
  
  # unpack inputs
  res_ls      = res_ls_concat %>% str_split(" ") %>% unlist %>% as.numeric
  
  # change clustering method to integer
  if(cl_method == 'louvain'){
    cl_method = 1
  } else if (cl_method == 'leiden') { 
    cl_method = 4
    stop("sorry leiden doesn't work yet :(")
  }else{
    message('Using louvain clustering as default')
    cl_method = 1
  }
  
  message('running integration')
  
  message('  setting up cluster')
  plan("multicore", workers = n_cores)
  options( future.globals.maxSize = 2^35 )
  
  message('  loading relevant cell ids')
  smpl_stats_dt = fread(smpl_stats_f)
  assert_that("sample_id" %in% colnames(smpl_stats_dt))
  ok_samples  = smpl_stats_dt[ bad_sample == FALSE ]$sample_id
  
  all_coldata = fread(coldata_f) %>%
    setkey(cell_id)
  assert_that("sample_id" %in% colnames(all_coldata))

  message('  loading hvg matrix')
  hvg_mat     = .get_alevin_mx(hvg_mat_f, sel_s = '')
  
  # subset coldata
  assert_that( all(colnames(hvg_mat) %in% all_coldata$cell_id) )
  subset_coldata = all_coldata[colnames(hvg_mat), ]
  
  message('  normalizing hvg matrix')
  hvg_mat_norm =  normalize_hvg_mat(
    hvg_mat, subset_coldata, exclude_mito, scale_f = 10000
  )
  
  # turn into seurat object
  message('  prepping Seurat object')
  suppressWarnings({
    meta = data.frame(subset_coldata)
    rownames(meta) = meta$cell_id
    seu     = Seurat::CreateSeuratObject(
      counts      = hvg_mat,
      meta.data   = meta,
      project     = "dummy"
    )
    # add normalized counts 
    seu[['RNA']]$data = hvg_mat_norm
  })
  
  
  # run harmony
  int_dt = .run_one_integration(
    seu, batch_var, cl_method, n_dims, 
    theta = theta, res_ls, reduction)        
  
  # save outputs
  fwrite(int_dt, file = integration_f)
  
  message('done!')
}

make_subset_sces <- function(sel_s, clean_sce_f, integration_f, smpl_stats_f,
  sces_yaml_f, subset_f, subset_col, subset_str){
  
  # get all good samples
  zoom_stats_dt = fread(smpl_stats_f)
  keep_samples  = zoom_stats_dt[bad_sample == FALSE, sample_id]
  
  # get list of input files
  all_sce_paths = yaml::read_yaml(sces_yaml_f) %>% unlist()
  
  if(!sel_s %in% keep_samples){
    # write empty file if sample was excluded
    message('Sample ', sel_s, ' was excluded. Creating empty sce file')
    file.create(clean_sce_f)
  }else{
    message('Creating sce file for sample ', sel_s)
    # get input sce file
    in_sce_f = all_sce_paths[[sel_s]]
    assert_that(file.exists(in_sce_f))
    in_sce   = readRDS(in_sce_f)
    
    # get cell ids to extract
    subset_vals = str_split(subset_str, pattern = ',') %>% unlist
    subset_dt   = fread(subset_f)
    assert_that(subset_col %in% colnames(subset_dt))
    assert_that(all(c("cell_id", "sample_id") %in% colnames(subset_dt)))
    
    subset_dt = subset_dt %>%
      .[sample_id == sel_s] %>%
      .[get(subset_col) %in% subset_vals] %>%
      .[, c('sample_id', 'cell_id', subset_col), with = FALSE]
    
    assert_that(all(subset_dt$cell_id %in% colnames(in_sce)))

    sce_zoom = in_sce[, subset_dt$cell_id]
    
    # get integration results
    int_dt    = fread(integration_f) 
    smpl_int  = int_dt %>% .[sample_id == sel_s]
  
    # remove umap and clustering cols from before and add new ones
    rm_cols = c('UMAP1', 'UMAP2', str_subset(names(colData(sce_zoom)), "RNA_snn_res"))
    new_coldata = colData(sce_zoom) %>% as.data.table %>%
      .[ , (rm_cols) := NULL] %>%
    # add label
      merge(subset_dt, by = c('cell_id', 'sample_id')) %>%
      setnames(new = 'label', old = subset_col) %>%
      setkey(cell_id)
    
    sce_zoom = sce_zoom[, smpl_int$cell_id]
    new_coldata = new_coldata[smpl_int$cell_id, ]

    colData(sce_zoom) = DataFrame(as.data.frame(new_coldata))
    
    # get useful integration variables
    int_vs   = c('UMAP1', 'UMAP2', str_subset(names(smpl_int), "RNA_snn_res"))
  
    # add these to sce object
    for (v in int_vs) {
      if (str_detect(v, "RNA_snn_res")) {
        colData(sce_zoom)[[ v ]] = smpl_int[[ v ]] %>% factor
      } else {
        colData(sce_zoom)[[ v ]] = smpl_int[[ v ]]
      }
    }
    
    saveRDS(sce_zoom, clean_sce_f, compress = FALSE)
  }
  
  message('done!')
}

# 
# zoom_integrate_within_group <- function(full_tag, date_stamp, zoom_dir, 
#   hmny_f, sce_all_f, dbl_f, species, gtf_dt_f, 
#   sel_res, exc_regex, dbl_res, dbl_cl_prop, theta, 
#   not_ok_re, gsea_dir, min_cpm_go, max_zero_p, gsea_cut, min_cells, 
#   zoom_name, zoom_sel_cls_concat, zoom_res, zoom_n_hvgs, zoom_n_dims, 
#   zoom_min_n_sample, zoom_min_n_cl, zoom_n_train, 
#   n_cores = 4, seed = 123, overwrite = FALSE) {
#   # check some inputs
#   assert_that(
#     is.character(zoom_res),
#     is.character(full_tag),
#     is.character(date_stamp),
#     is.character(zoom_dir),
#     is.character(hmny_f),
#     is.character(sce_all_f),
#     is.character(dbl_f),
#     is.character(species),
#     is.character(gtf_dt_f),
#     is.character(exc_regex),
#     is.character(not_ok_re),
#     is.character(gsea_dir),
#     is.character(zoom_name)
#   )
#   assert_that(
#     dir.exists(zoom_dir),
#     file.exists(hmny_f),
#     file.exists(sce_all_f),
#     file.exists(dbl_f),
#     file.exists(gtf_dt_f),
#     dir.exists(gsea_dir)
#   )
#   assert_that(
#     is.numeric(sel_res), 
#     is.numeric(dbl_res), 
#     is.numeric(dbl_cl_prop), 
#     is.numeric(theta), 
#     is.numeric(min_cpm_go), 
#     is.numeric(max_zero_p), 
#     is.numeric(gsea_cut), 
#     is.numeric(min_cells), 
#     is.numeric(zoom_n_hvgs), 
#     is.numeric(zoom_n_dims), 
#     is.numeric(zoom_min_n_sample), 
#     is.numeric(zoom_min_n_cl), 
#     is.numeric(zoom_n_train), 
#     is.numeric(n_cores), 
#     is.numeric(seed)
#   )
#   # get nice clusters
#   hmny_dt   = fread(hmny_f) %>% .[ !is.na(UMAP1) ] %>% 
#     .[, cluster := get(paste0("RNA_snn_res.", sel_res))] %>% 
#     .[, .(sample_id, cell_id, integration, UMAP1, UMAP2, cluster) ]
#   zoom_sel_cls  = zoom_sel_cls_concat %>% str_split(" ") %>% unlist
#   assert_that( all(zoom_sel_cls %in% unique(hmny_dt$cluster)) )
# 
#   # define files to save
#   this_dir     = sprintf("%s/%s", zoom_dir, zoom_name)
#   if (!dir.exists(this_dir))
#     dir.create(this_dir)
# 
#   # define files for harmony
#   keep_sub_f  = sprintf("%s/zoom_cluster_subset_%s_%s_%s.txt.gz", 
#     this_dir, full_tag, zoom_name, date_stamp)
#   hmny_sub_f  = sprintf("%s/zoom_integrated_dt_%s_%s_%s_%s.txt.gz", 
#     this_dir, full_tag, zoom_name, zoom_res, date_stamp)
#   hvgs_sub_f  = sprintf("%s/zoom_harmony_hvgs_%s_%s_%s.txt.gz", 
#     this_dir, full_tag, zoom_name, date_stamp)
#   sce_sub_f   = sprintf("%s/zoom_sce_clean_%s_%s_%s_%s.rds", 
#     this_dir, full_tag, zoom_name, zoom_res, date_stamp)
# 
#   # define files for marker genes
#   pb_f          = sprintf("%s/zoom_pb_%s_%s_%s_%s.rds",
#     this_dir, full_tag, zoom_name, zoom_res, date_stamp)
#   mkrs_f        = sprintf("%s/zoom_pb_marker_genes_%s_%s_%s_%s.txt.gz",
#     this_dir, full_tag, zoom_name, zoom_res, date_stamp)
#   hvgs_f        = sprintf("%s/zoom_pb_hvgs_%s_%s_%s_%s.txt.gz",
#     this_dir, full_tag, zoom_name, zoom_res, date_stamp)
#   fgsea_go_bp_f = sprintf("%s/zoom_fgsea_%s_%s_%s_go_bp_%s.txt.gz",
#     this_dir, full_tag, zoom_name, zoom_res, date_stamp)
#   fgsea_go_cc_f = sprintf("%s/zoom_fgsea_%s_%s_%s_go_cc_%s.txt.gz",
#     this_dir, full_tag, zoom_name, zoom_res, date_stamp)
#   fgsea_go_mf_f = sprintf("%s/zoom_fgsea_%s_%s_%s_go_mf_%s.txt.gz",
#     this_dir, full_tag, zoom_name, zoom_res, date_stamp)
#   fgsea_paths_f = sprintf("%s/zoom_fgsea_%s_%s_%s_paths_%s.txt.gz",
#     this_dir, full_tag, zoom_name, zoom_res, date_stamp)
#   fgsea_hlmk_f  = sprintf("%s/zoom_fgsea_%s_%s_%s_hlmk_%s.txt.gz",
#     this_dir, full_tag, zoom_name, zoom_res, date_stamp)
# 
#   # xgboost parameters
#   hvg_pcs_f   = sprintf("%s/zoom_hvg_pcs_%s_%s_%s.rds", 
#     this_dir, full_tag, zoom_name, date_stamp)
#   boost_f     = sprintf("%s/zoom_xgboost_obj_%s_%s_%s_%s.rds", 
#     this_dir, full_tag, zoom_name, zoom_res, date_stamp)
#   pairwise_f  = sprintf("%s/zoom_pairwise_dt_%s_%s_%s_%s.txt.gz", 
#     this_dir, full_tag, zoom_name, zoom_res, date_stamp)
#   imputed_f   = sprintf("%s/zoom_imputed_dt_%s_%s_%s_%s.txt.gz", 
#     this_dir, full_tag, zoom_name, zoom_res, date_stamp)
# 
#   # split into groups, save into separate txt files
#   message("  saving just these")
#   hmny_grp    = hmny_dt[ cluster %in% zoom_sel_cls ] %>% 
#     .[, N_sample := .N, by = sample_id ] %>% 
#     .[ N_sample >= zoom_min_n_sample ]
#   fwrite(hmny_grp, file = keep_sub_f)
# 
#   # run harmony on this subset
#   res_ls_concat = zoom_res %>% paste(collapse = " ")
#   if ( !all(file.exists(c(keep_sub_f, hmny_sub_f, sce_sub_f))) | overwrite ) {
#     message("  running harmony on them")
#     set.seed(seed)
#     run_harmony(sce_all_f, keep_sub_f, dbl_f, 
#       exc_regex, zoom_n_hvgs, zoom_n_dims, dbl_res, dbl_cl_prop, theta, res_ls_concat,
#       hmny_sub_f, hvgs_sub_f, sce_sub_f, n_cores)      
#   }
# 
#   # infer labels for missing cells
#   if ( !all(file.exists(c(hvg_pcs_f, boost_f, pairwise_f, imputed_f))) | overwrite ) {
#     # get parameters for imputing cells
#     message("  inferring labels for missing cells")
#     set.seed(seed)
#     impute_clusters_for_small_samples(hmny_dt, zoom_res, zoom_sel_cls, hmny_sub_f, sce_all_f, 
#       hvg_pcs_f, boost_f, pairwise_f, imputed_f,
#       exc_regex, zoom_n_hvgs, zoom_n_dims, zoom_min_n_cl, zoom_n_train, n_cores)
#   }
# 
#   # find marker genes for subset clusters
#   if ( !all(file.exists(c(pb_f, mkrs_f, hvgs_f, 
#     fgsea_go_bp_f, fgsea_go_cc_f, fgsea_go_mf_f, fgsea_paths_f, fgsea_hlmk_f))) | overwrite ) {
#     message("  finding marker genes")
#     calculate_marker_genes(sce_sub_f, pb_f, mkrs_f, hvgs_f,
#       fgsea_go_bp_f, fgsea_go_cc_f, fgsea_go_mf_f, fgsea_paths_f, fgsea_hlmk_f,
#       species, gtf_dt_f, gsea_dir,
#       zoom_res, exc_regex, zoom_min_n_cl, min_cells,
#       not_ok_re, min_cpm_go, max_zero_p, gsea_cut, n_cores)
#   }
# }
# 
# impute_clusters_for_small_samples <- function(hmny_dt, zoom_res, zoom_sel_cls, 
#   hmny_sub_f, sce_all_f, hvg_pcs_f, boost_f, pairwise_f, imputed_f, exc_regex, 
#   n_hvgs, n_dims, min_n_cl, n_train, n_cores = 4, overwrite = TRUE) {
#   if ( all(file.exists(hvg_pcs_f, boost_f, pairwise_f, imputed_f)) & !overwrite ) {
#     message("already done!")
#     return(NULL)
#   }
#   # get all cells and selected cells
#   selected_dt = hmny_dt %>% .[ cluster %in% zoom_sel_cls ]
#   clusts_dt   = fread(hmny_sub_f) %>% .[ !is.na(UMAP1) ] %>% 
#     .[, cluster := get(paste0("RNA_snn_res.", zoom_res)) %>% factor ] %>% 
#     .[, .(sample_id, cell_id, cluster) ]
#   clust_ns    = clusts_dt[, .N, by = cluster ]
#   clusts_dt   = clusts_dt[ cluster %in% clust_ns[ N >= min_n_cl ]$cluster ] %>% 
#     .[, cluster   := cluster %>% fct_drop ]
#   all_ids     = selected_dt$cell_id
#   sel_ids     = clusts_dt$cell_id
# 
#   # get HVGs
#   hvg_pcs_dt  = get_hvgs_dt(hvg_pcs_f, sce_all_f, all_ids, sel_ids, what = "pca",
#     exc_regex = exc_regex, n_hvgs = n_hvgs, n_dims = n_dims, 
#     n_cores = n_cores, overwrite = overwrite)
# 
#   # extract data
#   data_ls     = load_train_test_data(clusts_dt, hvg_pcs_dt, min_n_cl, n_train)
#   train_dt    = data_ls$train
#   valid_dt    = data_ls$valid
#   valid_rest  = data_ls$valid_rest
#   test_dt     = data_ls$test
# 
#   # run xgboost
#   boost_obj   = run_xgboost(train_dt, valid_dt, n_cores)
#   saveRDS(boost_obj, file = boost_f, compress = FALSE)  
# 
#   # predict on validation data
#   valid_all   = rbind(valid_dt, valid_rest)
#   pred_valid  = get_pred_valid(boost_obj, valid_all)
#   conf_boost  = calc_confuse_xgboost_dt(pred_valid)
# 
#   pairwise_dt = calc_pairwise_cluster_predictions(pred_valid, conf_boost)
#   fwrite(pairwise_dt, file = pairwise_f)
# 
#   # make predictions
#   cl_ls       = levels(train_dt$cluster)
#   preds_dt    = predict_w_boost(test_dt, boost_obj, cl_ls)
#   preds_dt[, .(n_na = .N, n_conf = sum(max_p > 0.5)), by = sample_id ] %>% 
#     .[ order(n_conf) ]
# 
#   # join with already-clustered data
#   selected_tmp  = selected_dt[, .(sample_id, cell_id)] %>% 
#     merge(clusts_dt, by = c("sample_id", "cell_id"), all.x = TRUE)
#   imputed_dt    = add_predictions(selected_tmp, preds_dt, cl_ls, p_cutoff = 0.5) %>% 
#     .[, .(sample_id, cell_id, cluster = cl_new)]
#   fwrite(imputed_dt, file = imputed_f)
# }
# 
# get_hvgs_dt <- function(save_f, sce_f, all_ids, sel_ids, what = c("pca", "hvgs"), 
#   exc_regex = NULL, n_hvgs = 2000, n_dims = 50, n_cores = 4, overwrite = FALSE) {
#   if (file.exists(save_f) & overwrite == FALSE) {
#     message('already done!')
#     hvg_pcs_dt  = fread(save_f)
# 
#     return(hvg_pcs_dt)
#   }
#   # check inputs
#   what        = match.arg(what)
# 
#   message('subsetting to HVGs')
#   message("  extracting ", what, " values")
#   # load sce, restrict to ok cells
#   message('  setting up cluster')
#   plan("multicore", workers = n_cores)
#   options( future.globals.maxSize = 2^36 )
# 
#   message('  loading sce')
#   sce_sel     = readRDS(sce_f) %>% .[, sel_ids ]
# 
#   # exclude genes if requested
#   if (!is.null(exc_regex)) {
#     exc_idx     = rownames(sce_sel) %>% str_detect(exc_regex)
#     exc_gs_str  = rowData(sce_sel)$symbol[ exc_idx ] %>% paste0(collapse = " ")
#     sprintf("    excluding %d genes: %s", sum(exc_idx), exc_gs_str) %>% message
#     sce_sel     = sce_sel[ !exc_idx, ]
#   }
# 
#   # turn into seurat object
#   message('  converting to Seurat object')
#   seu_sel     = Seurat::CreateSeuratObject(
#     counts      = counts(sce_sel),
#     meta.data   = data.frame(colData(sce_sel)),
#     project     = "MS2"
#     )
#   rm(sce_sel); gc()
#   
#   # run Seurat pipeline, plus clustering    
#   message('  finding HVGs')
#   seu_sel     = NormalizeData(seu_sel, verbose = FALSE ) %>% 
#     FindVariableFeatures( nfeatures = n_hvgs, verbose = FALSE )
#   var_feats   = VariableFeatures(seu_sel)
# 
#   # now use these genes for all cells
#   message('  creating seurat object for all data')
#   sce_all     = readRDS(sce_f) %>% .[, all_ids ]
#   seu_all     = Seurat::CreateSeuratObject(
#     counts      = counts(sce_all),
#     meta.data   = data.frame(colData(sce_all)),
#     project     = "MS2"
#     )
#   rm(sce_all); gc()
# 
#   # process this data
#   message('  normalize all data')
#   seu_all     = NormalizeData(seu_all, verbose = FALSE ) %>% 
#     ScaleData( verbose = FALSE ) %>% 
#     RunPCA( features = var_feats, verbose = FALSE, npcs = n_dims ) %>% 
#     identity()
# 
#   # switch cluster off
#   plan("sequential")
# 
#   if (what == "pca") {
#     # now use these genes for all cells
#     message('  extract PCs, save')
#     save_dt   = Embeddings(seu_all, reduction = "pca") %>% 
#       as.data.table(keep.rownames = "cell_id")    
#   } else if (what == "hvgs") {
#     # now use these genes for all cells
#     message('  extract HVGs, save')
#     save_dt   = GetAssayData(seu_all, slot = "data", assay = "RNA") %>% 
#       .[ var_feats, ] %>% t %>% 
#       as.data.table(keep.rownames = "cell_id")    
#   }
#   fwrite(save_dt, file = save_f)
# 
#   message('done!')
# 
#   return(save_dt)
# }
# 
# load_train_test_data <- function(clusts_dt, hvg_pcs_dt, min_cells, n_train) {
#   # check inputs
#   assert_that( all(clusts_dt$cell_id %in% hvg_pcs_dt$cell_id) )
#   clusts_na   = merge(clusts_dt, hvg_pcs_dt, by = "cell_id", all = TRUE) %>% 
#     .[, sample_id := NULL ]
# 
#   # which clusters are too small to bother with?
#   ns_dt       = clusts_na[ !is.na(cluster) ] %>% .[, .N, by = cluster]
#   keep_cl     = ns_dt[ N >= min_cells ]$cluster %>% as.character
# 
#   # get train data, balance samples
#   train_dt    = clusts_na[ cluster %in% keep_cl ] %>% 
#     .[, .SD[ sample(ceiling(min(.N, n_train) / 2)) ], by = cluster ] %>% 
#     .[, cluster := cluster %>% fct_drop ]
#   
#   # get validation data
#   valid_dt    = clusts_na[ cluster %in% keep_cl ] %>%   
#     .[ !(cell_id %in% train_dt$cell_id) ] %>% 
#     .[, .SD[ sample(ceiling(min(.N, n_train) / 2)) ], by = cluster ] %>% 
#     .[, cluster := cluster %>% fct_drop ]
#   valid_rest  = clusts_na[ cluster %in% keep_cl ] %>% 
#     .[ !(cell_id %in% c(train_dt$cell_id, valid_dt$cell_id)) ] %>% 
#     .[, cluster := cluster %>% fct_drop ]
# 
#   # get test data
#   test_dt     = clusts_na[ is.na(cluster) ]
# 
#   # some checks
#   assert_that( all( levels(valid_dt$cluster) == levels(train_dt$cluster) ) )
#   chk_dt_1    = clusts_na[ is.na(cluster) | (cluster %in% keep_cl) ]
#   chk_dt_2    = rbind(train_dt, valid_dt, valid_rest, test_dt)
#   assert_that( nrow(chk_dt_1) == nrow(chk_dt_2) )
#   assert_that( all( sort(chk_dt_1$cell_id) == sort(chk_dt_2$cell_id) ) )
# 
#   return(list(
#     train       = train_dt, 
#     valid       = valid_dt, 
#     valid_rest  = valid_rest, 
#     test        = test_dt
#   ))
# }
# 
# run_xgboost <- function(train_dt, valid_dt, n_cores) {
#   # convert training data to expected format
#   train_vars  = colnames(train_dt) %>% setdiff(c("cluster", "cell_id"))
#   assert_that( length(train_vars) > 0 )
#   train_mat   = train_dt[, c("cell_id", train_vars), with = FALSE] %>% 
#     as.matrix( rownames = "cell_id" )
#   train_cl    = train_dt$cluster
#   train_y     = train_cl %>% as.integer %>% `-`(1)
#   # weights_v   = 1 / table(train_y) * 1000
#   # weights_y   = weights_v[ train_y + 1 ]
#   assert_that(
#     min(train_y) == 0,
#     max(train_y) + 1 == length(levels(train_dt$cluster))
#   )
# 
#   # convert validation data to expected format
#   valid_mat   = valid_dt[, c("cell_id", train_vars), with = FALSE] %>% 
#     as.matrix( rownames = "cell_id" )
#   valid_cl    = valid_dt$cluster
#   valid_y     = valid_cl %>% as.integer %>% `-`(1)
# 
#   # set up parameters for xgboost
#   dtrain      = xgb.DMatrix( data = train_mat, label = train_y )
#   dvalid      = xgb.DMatrix( data = valid_mat, label = valid_y )
#   evals       = list( train = dtrain, eval = dvalid )
#   param       = xgb.params(
#     nthread     = n_cores, 
#     objective   = "multi:softprob", 
#     num_class   = max(train_y) + 1
#   )
# 
#   # run boost
#   boost_obj   = xgb.train(params = param, data = dtrain, nrounds = 100, 
#     evals = evals, early_stopping_rounds = 5, verbose = 2)
# 
#   # # run standard boost
#   # boost_obj   = xgboost(
#   #   data = train_mat, label = train_y, 
#   #   weight = weights_y,
#   #   objective = "multi:softprob", num_class = max(train_y) + 1,
#   #   nthread = n_cores, nrounds = 10, verbose = 2)
# 
#   return(boost_obj)
# }
# 
# get_pred_valid <- function(boost_obj, valid_dt) {
#   # get validation data
#   train_vars  = colnames(valid_dt) %>% setdiff(c("cluster", "cell_id"))
#   assert_that( length(train_vars) > 0 )
#   valid_mat   = valid_dt[, c("cell_id", train_vars), with = FALSE] %>% 
#     as.matrix( rownames = "cell_id" )
# 
#   # get probabilities for each predicted cluster
#   probs_mat   = predict(boost_obj, valid_mat, reshape = TRUE)
#   assert_that(
#     length(levels(valid_dt$cluster)) == ncol(probs_mat),
#     length(rownames(valid_mat)) == nrow(probs_mat)
#   )
#   probs_mat   = probs_mat %>% 
#     set_colnames(levels(valid_dt$cluster)) %>% 
#     set_rownames(rownames(valid_mat))
# 
#   # prediction for each cell
#   pred_valid  = data.table(
#     cell_id     = rownames(probs_mat),
#     cl_true     = valid_dt$cluster,
#     cl_pred     = colnames(probs_mat)[ apply(probs_mat, 1, which.max) ],
#     p_pred      = apply(probs_mat, 1, max)
#     ) %>% cbind(probs_mat)
# 
#   return(pred_valid)
# }
# 
# calc_confuse_xgboost_dt <- function(pred_valid) {
#   n_cats    = length(unique(pred_valid$cl_true))
#   confuse_dt  = pred_valid[, .N, by=.(cl_true, cl_pred)] %>%
#     .[, N_true  := sum(N), by = cl_true] %>%
#     .[, prop    := N / N_true ] %>%
#     .[, logit   := qlogis((N+1) / (N_true + n_cats))] %>%
#     .[, H       := -sum(prop*log2(prop)), by = cl_true ]
# 
#   return(confuse_dt)
# }
# 
# calc_pairwise_cluster_predictions <- function(pred_valid, confuse_dt, n_cores = 8) {
#   # list all types, make pairs
#   message('  making pairs to compare')
#   cl_list     = levels(pred_valid$cl_true)
#   pairs_list  = t(combn(cl_list, 2)) %>% asplit(1)
# 
#   # calculate for every pair
#   message('  calculating measures for ', length(pairs_list), ' pairs', sep='')
#   bpparam     = MulticoreParam(workers = n_cores, progressbar = TRUE, tasks = length(pairs_list))
#   on.exit(bpstop(bpparam))
#   pairwise_dt = bplapply(seq_along(pairs_list), .x_entropy_fn, pairs_list, pred_valid,
#     BPPARAM = bpparam) %>% rbindlist
# 
#   # get sccaf value, add to dt
#   message('  calculating SCCAF scores')
#   sccaf_dt    = .calc_sccaf_scores(confuse_dt, pairs_list)
#   assert_that(
#     all(pairwise_dt$cl1 == sccaf_dt$cl1),
#     all(pairwise_dt$cl2 == sccaf_dt$cl2)
#     )
#   pairwise_dt[, sccaf_val := sccaf_dt$sccaf_val ]
# 
#   # duplicate to have both directions
#   plot_dt   = rbind(
#     pairwise_dt[, .(cl1, cl2, error)],
#     pairwise_dt[, .(cl1 = cl2, cl2 = cl1, error)]
#     )
# 
#   # make into matrix
#   plot_wide = dcast(plot_dt, cl1 ~ cl2, value.var = 'error', fill = 0)
#   h_mat     = plot_wide[, -'cl1', with=FALSE] %>%
#     as.matrix %>% set_rownames(plot_wide$cl1)
#   assert_that( all(rownames(h_mat) == colnames(h_mat)) )
# 
#   # calculate nice ordering, add as factor level
#   seriate_obj = seriate(h_mat)
#   seriate_ord = get_order(seriate_obj, 1)
#   new_levels  = rownames(h_mat)[seriate_ord]
#   pairwise_dt[, cl1 := factor(cl1, levels = new_levels) ]
#   pairwise_dt[, cl2 := factor(cl2, levels = new_levels) ]
# 
#   return(pairwise_dt)
# }
# 
# .x_entropy_fn <- function(i, pairs_list, pred_valid) {
#   # restrict to these
#   pair    = pairs_list[[i]]
#   tmp_dt  = pred_valid[ (cl_true %in% pair) & (cl_pred %in% pair) ] %>%
#     melt(id = c("cell_id", "cl_true", "cl_pred", "p_pred"),
#       variable.name = "cl", value.name = "p") %>%
#     .[ cl %in% pair ]
#   if (nrow(tmp_dt) == 0) {
#     return(data.table())
#   }
#   tmp_dt  = tmp_dt %>%
#     .[, which_p := ifelse(cl_true == cl, 'prob_T_raw',
#       'prob_F_raw') %>% factor(levels = c('prob_T_raw', 'prob_F_raw')) ] %>%
#     .[, correct := cl_true == cl_pred ] %>%
#     dcast( cell_id + cl_true + cl_pred + correct ~ which_p,
#       value.var = 'p')
#   # add cols of zeros if necessary
#   if ( !('prob_F_raw' %in% names(tmp_dt)) )
#     tmp_dt[, prob_F_raw := 0 ]
#   if ( !('prob_T_raw' %in% names(tmp_dt)) )
#     tmp_dt[, prob_T_raw := 0 ]
# 
#   # record prob_correct; prob_other. normalize too
#   tmp_dt  = tmp_dt %>%
#     .[, prob_T_norm := prob_T_raw / (prob_T_raw + prob_F_raw) ] %>%
#     .[, prob_F_norm := prob_F_raw / (prob_T_raw + prob_F_raw) ] %>%
#     .[, H_F_raw   := ifelse(prob_F_raw == 0, 0, -log2(prob_F_raw)) ] %>%
#     .[, H_T_raw   := ifelse(prob_T_raw == 0, 0, -log2(prob_T_raw)) ] %>%
#     .[, H_F_norm  := ifelse(prob_F_norm == 0, 0, -log2(prob_F_norm)) ] %>%
#     .[, H_T_norm  := ifelse(prob_T_norm == 0, 0, -log2(prob_T_norm)) ] %>%
#     .[, H_raw     := correct*H_T_raw + (1-correct)*H_F_raw ] %>%
#     .[, H_norm    := correct*H_T_norm + (1-correct)*H_F_norm ]
# 
#   # summarize: mean x-entropy (raw); mean x-entropy (norm); mean accuracy; N_{combns}
#   summary_dt  = tmp_dt[,
#     .(
#       cl1     = pair[[1]], 
#       cl2     = pair[[2]],
#       H_raw   = mean(H_raw),
#       H_norm  = mean(H_norm),
#       error   = mean(!correct),
#       N_1     = sum( cl_true == pair[[1]] ),
#       N_2     = sum( cl_true == pair[[2]] ),
#       N       = .N
#       )
#     ]
# 
#   return(summary_dt)
# }
# 
# .calc_sccaf_scores <- function(confuse_dt, pairs_list) {
#   sccaf_wide  = confuse_dt %>% 
#     dcast.data.table(cl_true ~ cl_pred, value.var='prop', fill=0)
#   sccaf_mat   = sccaf_wide %>% as.matrix(rownames = "cl_true")
#   sccaf_vals_dt   = lapply(pairs_list,
#     function(p) {
#       v1  = sccaf_mat[p[[1]], p[[2]]] / sccaf_mat[p[[1]], p[[1]]]
#       v2  = sccaf_mat[p[[2]], p[[1]]] / sccaf_mat[p[[2]], p[[2]]]
#       sccaf_val = max(v1, v2)
#       if (is.infinite(sccaf_val))
#         sccaf_val   = NA
#       return( data.table(cl1=p[[1]], cl2=p[[2]], sccaf_val=sccaf_val) )
#     }) %>% rbindlist
# }
# 
# calc_merged_dt <- function(pairwise_dt, labels_dt, save_dir, v_cut, 
#   merge_var = c("error", "H_norm", "sccaf_val")) {
#   # check inputs
#   merge_var   = match.arg(merge_var)
# 
#   # # which to keep?
#   # dont_merge  = mrg_custom %>% unlist
#   # assert_that( all(dont_merge %in% pairwise_dt$conos1) )
# 
#   # turn into matrix
#   merge_wide  = rbind(
#     pairwise_dt[ cl1 != cl2, .(cl1, cl2, var = get(merge_var))],
#     pairwise_dt[ cl1 != cl2, .(cl2 = cl1, cl1 = cl2, var = get(merge_var))]
#     ) %>%
#     # .[ !(cl1 %in% dont_merge) & !(cl2 %in% dont_merge) ] %>%
#     dcast.data.table( cl1 ~ cl2, value.var = "var", fill = 0)
#   merge_mat   = merge_wide %>% as.matrix(rownames = "cl1")
#   assert_that( all(rownames(merge_mat) == colnames(merge_mat)) )
#   diag(merge_mat) = 1
# 
#   # try out different clustering approaches
#   method_list = names(v_cut)
#   merged_dt   = method_list %>% 
#     lapply(.calc_one_merge, merge_mat, merge_var, v_cut, save_dir) %>%
#     rbindlist
# 
#   # # add back in the preserved ones
#   # custom_dt   = .calc_merge_custom(method_list, merge_var, mrg_custom)
#   # merged_dt   = rbind(merged_dt, custom_dt) %>%
#   #   setorder('method', 'cl_merge', 'cl')
# 
#   # add some labels
#   merged_dt   = merge(merged_dt, labels_dt, by='cl') %>%
#     setcolorder(c("metric", "hclust_link", "type_broad", "cl_merge", "N_merge", "cl"))
# 
#   return(merged_dt)
# }
# 
# .calc_one_merge <- function(m, merge_mat, merge_var, v_cut, save_dir) {
#   # do hierarchical clustering
#   merge_clust = hclust(as.dist(1 - merge_mat), method = m)
#   merge_dend  = as.dendrogram(merge_clust)
#   merge_cut   = cutree(merge_clust, h = 1 - v_cut[[m]])
# 
#   # assemble into data.table
#   dt = data.table(
#     metric      = merge_var,
#     hclust_link = m,
#     cl          = names(merge_cut),
#     merge_clust = merge_cut
#     ) %>%
#     .[, N_merge := .N, by = merge_clust] %>%
#     .[, cl_merge := sort(.SD$cl)[[1]], by = merge_clust]
# 
#   # save dendrogram
#   png(file.path(save_dir, sprintf('merging_dend_%s.png', m)),
#     h=5, w=10, units='in', res=300)
#   plot(merge_dend)
#   abline(h = 1 - v_cut[[m]], col = 'blue')
#   dev.off()
# 
#   return(dt)
# }
# 
# predict_w_boost <- function(query_dt, boost_obj, cl_ls) {
#   if (nrow(query_dt) == 0) {
#     return(data.table(
#       sample_id = vector("character", 0), 
#       cell_id   = vector("character", 0), 
#       max_p     = vector("numeric", 0), 
#       cl_pred   = vector("character", 0)
#     ))
#   }
#   # make predictions
#   train_vars  = colnames(query_dt) %>% setdiff(c("cluster", "cell_id"))
#   test_mat    = query_dt[, c("cell_id", train_vars), with = FALSE] %>% 
#     as.matrix( rownames = "cell_id" )
#   prob_test   = predict(boost_obj, test_mat) %>% 
#     matrix(nrow = nrow(test_mat), byrow = TRUE)
# 
#   # assemble results
#   preds_dt    = data.table(
#     sample_id   = rownames(test_mat) %>% str_extract("^.+(?=:)"),
#     cell_id     = rownames(test_mat),
#     max_p       = apply(prob_test, 1, max),
#     cl_pred     = apply(prob_test, 1, which.max) %>% cl_ls[ . ]
#     )
# 
#   return(preds_dt)
# }
# 
# add_predictions <- function(clusts_dt, preds_dt, cl_ls, p_cutoff = 0.5) {
#   if (nrow(preds_dt) == 0) {
#     adjusted_dt = copy(clusts_dt) %>% 
#       .[, .(sample_id, cell_id, cl_orig = cluster, cl_pred = NA, cl_new = cluster)]
#     return(adjusted_dt)
#   }
#   # what to add?
#   to_merge    = preds_dt[ max_p > p_cutoff, .(cell_id, cl_pred)]
# 
#   # add these to combined
#   adjusted_dt = clusts_dt[, .(sample_id, cell_id, cl_orig = cluster)] %>% 
#     merge(to_merge, by = "cell_id", all.x = TRUE)
#   assert_that( all( is.na(adjusted_dt$cl_orig) + is.na(adjusted_dt$cl_pred) >= 1 ) )
# 
#   # make column with new clusters
#   adjusted_dt = adjusted_dt[, cl_new := is.na(cl_orig) %>% 
#       ifelse(cl_pred, as.character(cl_orig)) %>% factor(levels = cl_ls) ]
# 
#   return(adjusted_dt)
# }
# 
# load_cell_expression <- function(sce_f, sel_dt, hmny_dt) {
#   sce       = sce_sub_f %>% readRDS %>% logNormCounts
# 
#   assert_that( all( sel_dt$gene_id %in% rownames(sce) ) )
#   sel_idx   = rownames(sce) %in% sel_dt$gene_id
#   assert_that( sum( sel_idx ) == nrow(sel_dt) )
# 
#   cell_exp_dt = logcounts(sce)[ sel_idx, ] %>% 
#     as.data.table( keep.rownames = "gene_id" ) %>% 
#     melt.data.table( id = "gene_id", var = "cell_id", val = "logcount" ) %>% 
#     merge( sel_dt, by = "gene_id" ) %>% 
#     merge( hmny_dt[, .(cell_id, UMAP1, UMAP2)], by = "cell_id" ) %>% 
#     .[, max_val   := quantile(logcount, 0.99), by = gene_id] %>% 
#     .[ is.na(max_val) | max_val == 0, max_val := max(logcount), by = gene_id ] %>% 
#     .[, norm_val  := logcount %>% `/`(max_val) %>% pmin(1)] %>% 
#     .[, symbol := symbol %>% factor( levels = sel_dt$symbol ) ]
# 
#   return(cell_exp_dt)
# }
# 
# plot_selected_genes_umap <- function(sel_dt) {
#   g = ggplot(sel_dt) +
#     aes( x = UMAP1, y = UMAP2, colour = norm_val ) +
#     geom_point( size = 0.1 ) +
#     scale_colour_viridis( breaks = pretty_breaks(),
#       guide = guide_legend( override.aes = list(size = 3) )) +
#     facet_wrap( ~ symbol ) +
#     theme_bw() +
#     theme( axis.text = element_blank(), panel.grid = element_blank(),
#       strip.background = element_rect( fill = "white" ) ) +
#     labs( colour = "log expression\n(max val. = 1)" )
# }
>>>>>>> dev
