# source('scripts/label_celltypes.R')
suppressPackageStartupMessages({
  library('RColorBrewer')
  library("BiocParallel")
  library('circlize')
  library('magrittr')
  library('data.table')
  library('stringr')
  library('assertthat')
  library('viridis')
  library('scales')
  library('ggplot2')
  library('patchwork')
  library('forcats')
  library('readxl')

  library('future')
  library('SingleCellExperiment')

  library('scater')
  library('Seurat')

  library('ComplexHeatmap')
  library('seriation')
  library('purrr')
  library('xgboost')
  library('ggrepel')
})


train_celltype_labeller <- function(sce_f, hvgs_xgb_f, xgb_f, allowed_f,
  clusters_dt, meta_dt, clust_var = "cluster", use_all_samples = FALSE, meta_vars = NULL,
  min_n_cl = 200, n_train = 1000, n_dims = 50, sel_gs = NULL, n_hvgs = 2000,
  seed = 123, n_cores = 4) {
  
  # randomly sample evenly across lesion types, until we have >= 500 of each type
  if (use_all_samples) {
    message('  using all samples')
    samples_dt  = copy(clusters_dt)
  } else {
    message('  getting representative subset of cells for training')
    set.seed(seed)
    samples_dt  = .get_representative_subset(clusters_dt, meta_dt,
      clust_var = clust_var, meta_vars = meta_vars,
      n_per_type = min_n_cl, min_n = 10)
    message(sprintf('    %d cells chosen, split like so:', nrow(samples_dt)))
    samples_dt[, .N, by = 'cluster'] %>% .[ order(-N) ] %>% print
  }
  # now subset
  sel_ids     = samples_dt$cell_id
  sel_cls     = clusters_dt[ cell_id %in% sel_ids ] %>%
    .[, .(sample_id, cell_id, cluster = get(clust_var))]

  # get HVGs for these
  message('  calculate HVGs on these')
  hvgs_mat    = .get_hvgs_mat(hvgs_xgb_f, sce_f, sel_ids, what = "hvgs",
    sel_gs = sel_gs, n_hvgs = n_hvgs, n_dims = n_dims, n_cores = n_cores)
  hvg_gs      = hvgs_mat %>% colnames %>% setdiff(c("cluster", "cell_id"))

  # make data for xgboost
  message('  split into train/test')
  set.seed(seed)
  data_broad  = .load_train_test_data(sel_cls, hvgs_mat, min_n_cl, n_train)
  train_dt    = data_broad$train
  valid_dt    = data_broad$valid

  # run xgboost
  message('  run XGBoost')
  set.seed(seed)
  xgb_obj     = .run_boost_watchlist(train_dt, valid_dt, n_cores)
  xgb_obj$cl_lvls = levels(train_dt$cluster)
  message('    saving')
  saveRDS(xgb_obj, file = xgb_f, compress = FALSE)
  allowed_dt  = data.table( cluster = xgb_obj$cl_lvls )
  fwrite(allowed_dt, file = allowed_f)

  # predict on validation data
  message('  print some outputs to show performance on validation data')
  valid_all   = rbind(data_broad$valid, data_broad$valid_rest)
  pred_valid  = .get_pred_valid(xgb_obj, valid_all)
  conf_dt     = .calc_confuse_xgboost_dt(pred_valid)
  conf_tmp    = .calc_confuse_xgboost_dt(pred_valid[ p_pred > 0.5 ])
  conf_dt[ (cl_pred != cl_true) & (prop > 0.1) ] %>%
    .[, .(cl_true, cl_pred, N_true, pct_pred = round(prop * 100, 1))] %>%
    .[ order(-pct_pred) ] %>%
    print
  conf_tmp[ (cl_pred != cl_true) & (prop > 0.01) ] %>%
    .[, .(cl_true, cl_pred, N_true, pct_pred = round(prop * 100, 1))] %>%
    .[ order(-pct_pred) ] %>%
    print
  conf_tmp[ (cl_pred == cl_true) ] %>% .[ order(prop) ] %>%
    .[, .(cl_true, N_true, pct_true = round(prop * 100, 1))] %>% print

  message('done.')
}

label_celltypes_with_xgboost <- function(xgb_f, sces_yaml_f, integration_f, qc_sample_stats_f, 
  hvg_mat_f, guesses_f, custom_labels_f, exclude_mito, sel_res,  gene_var = c("gene_id", "ensembl_id"), 
  min_pred = 0.5, min_cl_prop = 0.5, min_cl_size = 100, n_cores = 4) {
  
  exclude_mito = as.logical(exclude_mito)
  
  # check inputs
  if(custom_labels_f != ""){
    cust_lbls = fread(custom_labels_f)
    # check if cell_ids match cell_ids in sce_f
    sce = readRDS(sce_f)

    sel_res_str = paste0("RNA_snn_res.", sel_res)

    sce_cells_df = as.data.table(colData(sce)) %>%
    .[, .(sample_id, cell_id, UMAP1, UMAP2, cl_hmny = get(sel_res_str))]

    cell_ids_olap = intersect(sce_cells_df$cell_id, cust_lbls$cell_id)

    assert_that(
      length(cell_ids_olap) != 0, 
      msg = "values in cell_id column of the 'custom_labels' file don't match cell_id values in the sce object"
      )
    
    cell_ids_df = sce_cells_df %>%
    merge(cust_lbls, by = 'cell_id', all.x = TRUE, all.y = FALSE)
    
    # write file with celltype annotations
    fwrite(cell_ids_df, file = guesses_f)
    # write empty hvg_mat_f
    file.create(hvg_mat_f)
  
  }else{
  # check inputs
  assert_that( file.exists(xgb_f) )
  gene_var    = match.arg(gene_var)

  # load XGBoost object
  message('  loading XGBoost classifier')
  xgb_obj     = readRDS(xgb_f)
  hvgs        = xgb_obj$feature_names

  # get values for these genes in new datasets
  message('  getting lognorm counts of HVGs')
  
  hvg_mat     = .calc_logcounts(
    hvg_mat_f, sces_yaml_f, qc_sample_stats_f, gene_var,
    hvgs, exclude_mito, n_cores = n_cores
    )
  
  assert_that( is(hvg_mat, "sparseMatrix") )
  #assert_that( all( colnames(hvg_mat) == xgb_obj$feature_names ) )

  # predict for new data
  message('  predicting celltypes for all cells')
  preds_dt    = .predict_on_new_data(xgb_obj, hvg_mat, min_pred)

  # label harmony clusters
  message('  predicting majority celltype for each cluster')
  int_dt     = .load_clusters(integration_f)
  guesses_dt  = .apply_labels_by_cluster(int_dt, preds_dt, min_cl_prop, min_cl_size)

  # save
  message('  saving results')
  fwrite(guesses_dt, file = guesses_f)
  message('done.')
  }
}

save_subset_sces <- function(sce_f, guesses_f, sel_res_cl, subset_df_f, 
  sce_ls_concat, subset_names_concat, allowed_cls_f, custom_labels_f, n_cores = 4) {
  message('saving sce subsets')
  # unpack inputs
  sce_ls        = sce_ls_concat %>% str_split(" ") %>% unlist
  subset_names  = subset_names_concat %>% str_split(" ") %>% unlist
  assert_that( length(sce_ls) == length(subset_names) )
  assert_that( all( str_detect(sce_ls, subset_names) ) )
  names(sce_ls) = subset_names

  # get specifications
  subsets_dt  = fread(subset_df_f)
  message('  subset specifications:')
  print(subsets_dt)
  assert_that( length(setdiff(subset_names, subsets_dt$subset_name)) == 0 )

  if(custom_labels_f != ""){
    custom_labels = fread(custom_labels_f)
    allowed_cls = custom_labels$label %>% unique
    guess_col_name = 'label'
  }else{
    allowed_cls = allowed_cls_f %>% fread(sep = ",") %>% .$cluster
    guess_col_name = paste0("cl_pred_", sel_res_cl)
  }
  
  assert_that( length(setdiff(subsets_dt$guess, allowed_cls)) == 0 )

  
  # load guesses
  guesses_all = fread(guesses_f)
  assert_that( guess_col_name %in% names(guesses_all) )
  guesses_dt  = guesses_all[, .(cell_id, guess = get(guess_col_name)) ]

  # load sce
  sce         = readRDS(sce_f)
  assert_that( all(guesses_dt$cell_id == colnames(sce)) )

  # get subsets
  for (nn in subset_names) {
    message('  getting subset for ', nn)
    # where are they?
    sel_types   = subsets_dt[ subset_name == nn ]$guess
    sel_ids     = guesses_dt[ guess %in% sel_types ]$cell_id

    # take subset, save
    message('    saving')
    sce_subset  = sce[, sel_ids]
    saveRDS(sce_subset, file = sce_ls[[ nn ]])
  }
}

.get_representative_subset <- function(cl_dt, meta_dt, clust_var = "cluster", 
  meta_vars = NULL, n_per_type = 100, min_n = 10) {
  # initialize
  cl_tmp      = copy(cl_dt) %>%
    .[, .(sample_id, cell_id, cluster = get(clust_var))]
  sample_list = NULL
  cl_sub      = data.table(cell_id = character(0), sample_id = character(0),
    cluster = character(0))

  # define metadata combinations we want to balance
  if (is.null(meta_vars)) {
    meta_track  = meta_dt[, .(sample_id, combn_var = 'dummy')]
  } else {
    meta_track  = meta_dt[, c('sample_id', meta_vars), with = FALSE] %>%
      .[, combn_var := do.call(paste, c(.SD, sep = "_")), 
        .SDcols = meta_vars, by = meta_vars ]    
  }
  props_all  = meta_track[, .N, by = combn_var] %>%
    .[, p_all := N / sum(N) ]

  # add samples one by one until we have at least that many per conos cluster
  totals_dt   = cl_tmp[, .N, by = .(sample_id, cluster)] %>%
    .[ N > min_n ]
  types_list  = cl_tmp[, .N, .(cluster)] %>% 
    .[ order(N) ] %>% 
    use_series('cluster') %>% as.character

  # loop
  for (this_type in types_list) {
    n_type    = cl_sub[ cluster == this_type ] %>% nrow
    n_total   = totals_dt[ cluster == this_type ]$N %>% sum
    while (n_type < min(n_per_type, n_total)) {
      # which samples would help?
      sample_opts = cl_tmp[ cluster == this_type, .N, by = sample_id] %>% 
        setorder(-N) %>% .[ N > min_n ] %>% use_series("sample_id")

      # pick one which improves metadata representation
      sel_sample  = .pick_next_sample(meta_track, props_all, sample_list, sample_opts)

      # add to list
      sample_list = c(sample_list, sel_sample)
      cl_sub      = rbind(cl_sub, cl_tmp[sample_id == sel_sample])
      cl_tmp      = cl_tmp[!(sample_id %in% sample_list)]

      # update count
      n_type      = cl_sub[ cluster == this_type ] %>% nrow
    }
  }

  # check worked
  ns_all      = totals_dt[, .(n_all = sum(.N)), by = cluster]
  ns_sub      = cl_sub[, .(n_sub = .N), by = cluster]
  check_dt    = merge(ns_all, ns_sub, by = 'cluster', all = TRUE) %>% 
    .[, n_per_type  := n_per_type ] %>% 
    .[, is_ok       := n_sub >= pmin(n_per_type, n_all) ]
  assert_that( all(check_dt$is_ok) )

  return(cl_sub)
}

.pick_next_sample <- function(meta_track, props_all, sample_list, sample_opts) {
  # to start pick one at random
  if (is.null(sample_list))
    return(sample(sample_opts, 1))

  # otherwise calc current props
  props_now   = meta_track[ sample_id %in% sample_list ] %>%
    .[, .N, by = combn_var ] %>%
    .[, p_now   := N / sum(N) ]

  # combine with target props
  props_now   = merge(props_now, props_all, by = 'combn_var', all.y = TRUE) %>%
    .[ is.na(p_now), p_now := 0 ] %>%
    .[, p_delta := p_all - p_now ]

  # which is most out of whack, in the samples where we see this celltype?
  vars_ok     = meta_track[ sample_id %in% sample_opts ]$combn_var %>% unique
  sel_val     = props_now[ order(-p_delta) ] %>% 
    .[ combn_var %in% vars_ok ] %>%
    use_series('combn_var') %>% .[[1]]

  # pick a sample at random from these
  sel_sample  = meta_track[ (combn_var == sel_val) & (sample_id %in% sample_opts) ]$sample_id %>%
    sample(1)

  return(sel_sample)
}

.get_hvgs_mat <- function(hvgs_xgb_f, sce_f, sel_ids, what = c("pca", "hvgs"), 
  sel_gs = NULL, n_hvgs = 2000, n_dims = 50, n_cores = 4, overwrite = FALSE) {
  if (file.exists(hvgs_xgb_f) & overwrite == FALSE) {
    message('  already done')
    hvgs_mat    = readRDS(hvgs_xgb_f)

    return(hvgs_mat)
  }
  # check inputs
  what        = match.arg(what)

  message('  subsetting to HVGs')
  # load sce, restrict to ok cells
  message('    setting up cluster')
  plan("multicore", workers = n_cores)
  options( future.globals.maxSize = 2^35 )

  message('    loading sce')
  sce_sel     = readRDS(sce_f) %>% .[, sel_ids ]

  # restrict to specified genes
  if (!is.null(sel_gs)) {
    ref_gs      = rowData(sce_sel)[[ names(sel_gs) ]]
    sel_gs_v    = unlist(sel_gs)
    assert_that( all(sel_gs_v %in% ref_gs) )
    sel_idx     = ref_gs %in% sel_gs_v
    sce_sel     = sce_sel[ sel_idx, ]
  }
  # get counts
  counts_mat  = counts(sce_sel)
  if (!is.null(sel_gs)) {
    rownames(counts_mat) = ref_gs[ sel_idx ]
  }  

  # turn into seurat object
  message('    converting to Seurat object')
  seu_obj     = Seurat::CreateSeuratObject(
    counts      = counts_mat,
    meta.data   = data.frame(colData(sce_sel)),
    project     = "MS2"
    )
  rm(sce_sel); gc()
  
  # run Seurat pipeline, plus clustering    
  message('    finding HVGs')
  seu_obj     = NormalizeData(seu_obj, verbose = FALSE ) %>% 
    FindVariableFeatures( nfeatures = n_hvgs, verbose = FALSE )
  var_feats   = VariableFeatures(seu_obj)
  seu_obj     = seu_obj %>% 
    ScaleData( verbose = FALSE ) %>% 
    RunPCA( features = var_feats, verbose = FALSE, npcs = n_dims ) %>% 
    identity()

  # switch cluster off
  plan("sequential")

  if (what == "pca") {
    # now use these genes for all cells
    message('    extract PCs, save')
    save_mat    = Embeddings(seu_obj, reduction = "pca")
  } else if (what == "hvgs") {
    # now use these genes for all cells
    message('    extract HVGs, save')
    save_mat    = GetAssayData(seu_obj, slot = "data", assay = "RNA") %>% 
      .[ var_feats, ] %>% t
  }
  saveRDS(save_mat, file = hvgs_xgb_f)

  message('  done')

  return(save_mat)
}

.calc_logcounts <- function(hvg_mat_f, sces_yaml_f, qc_sample_stats_f, gene_var, hvgs, 
  exclude_mito, n_cores = 4, overwrite = FALSE) {
  if (file.exists(hvg_mat_f) & overwrite == FALSE) {
    message('    already done!')
    hvg_mat   = readRDS(hvg_mat_f)

    return(hvg_mat)
  }
  # load sce, restrict to ok cells
  message('    setting up cluster')
  plan("multicore", workers = n_cores)
  options( future.globals.maxSize = 2^35 )

  qc_dt      = fread(qc_sample_stats_f)
  ok_samples = qc_dt[bad_sample == FALSE, sample_id]
  
  sce_fs = yaml::read_yaml(sces_yaml_f)
  
  bpparam     = MulticoreParam(workers = n_cores, tasks = length(ok_samples))
  on.exit(bpstop(bpparam))
  
  message('    getting HVG matrix for each sample')
  hvg_mats = ok_samples %>% bplapply(
    FUN = .get_one_hvg_mat, BPPARAM = bpparam,
    sce_fs = sce_fs, hvgs = hvgs,
    gene_var = gene_var, exclude_mito = exclude_mito
    )

  message('    making one large HVG matrix')
  
  hvg_mat = do.call('cbind', hvg_mats) %>% t()
  # switch cluster off
  message('    saving')
  saveRDS(hvg_mat, file = hvg_mat_f)

  plan("sequential")
  message('done!')

  return(hvg_mat)
}


.get_one_hvg_mat <- function(sel_s, sce_fs, hvgs, gene_var, exclude_mito){
  
  message('    loading sce for sample ', sel_s)
  sce_f = sce_fs[[sel_s]]
  sce   = readRDS(sce_f)
  
  # get selected genes
  ref_gs      = rowData(sce)[[ gene_var ]]
  if (gene_var == "gene_id")
    ref_gs = ref_gs %>% str_replace("_ENSG", "-ENSG")
  
  # add zero counts matrix for missing hvgs if there are any
  if(!all(hvgs %in% ref_gs)){
    missing_hvgs = setdiff(hvgs, ref_gs)
    n_cols = ncol(sce)
    missing_mat = matrix(0, length(missing_hvgs), n_cols)
    rownames(missing_mat) = missing_hvgs
    colnames(missing_mat) = colnames(sce)
    # convert to sparse
    missing_mat = as(missing_mat, 'dgTMatrix')
    
    # get counts
    counts_mat  = counts(sce)
    rownames(counts_mat) = ref_gs
    
    counts_mat = rbind(counts_mat, missing_mat)
  }else{
    
    counts_mat  = counts(sce)
    rownames(counts_mat) = ref_gs
    
  }
  # select highly variable genes
  message('    creating normalized HVG matrix for sample ', sel_s)
  
  # normalize
  coldata      = colData(sce) %>% as.data.table()
  hvg_mat      = counts_mat[hvgs, ]
  norm_hvg_mat = normalize_hvg_mat(
    hvg_mat, coldata, exclude_mito
  )
  
  return(norm_hvg_mat)
  
}

.load_train_test_data <- function(clusts_dt, hvgs_mat, min_cells, n_train) {
  # check inputs
  assert_that( all(clusts_dt$cell_id %in% rownames(hvgs_mat)) )
  clusts_na   = hvgs_mat %>% as.data.table( keep.rownames = 'cell_id' ) %>% 
    merge(clusts_dt, by = "cell_id", all = TRUE) %>% 
    .[, sample_id := NULL ]

  # which clusters are too small to bother with?
  ns_dt       = clusts_na[ !is.na(cluster) ] %>% .[, .N, by = cluster]
  keep_cl     = ns_dt[ N >= min_cells ]$cluster %>% as.character

  # get train data, balance samples
  train_dt    = clusts_na[ cluster %in% keep_cl ] %>% 
    .[, .SD[ sample(ceiling(min(.N, n_train) / 2)) ], by = cluster ] %>% 
    .[, cluster := cluster %>% fct_drop ]
  
  # get validation data
  valid_dt    = clusts_na[ cluster %in% keep_cl ] %>%   
    .[ !(cell_id %in% train_dt$cell_id) ] %>% 
    .[, .SD[ sample(ceiling(min(.N, n_train) / 2)) ], by = cluster ] %>% 
    .[, cluster := cluster %>% fct_drop ]
  valid_rest  = clusts_na[ cluster %in% keep_cl ] %>% 
    .[ !(cell_id %in% c(train_dt$cell_id, valid_dt$cell_id)) ] %>% 
    .[, cluster := cluster %>% fct_drop ]

  # get test data
  test_dt     = clusts_na[ is.na(cluster) ]

  # some checks
  assert_that( all( levels(valid_dt$cluster) == levels(train_dt$cluster) ) )
  chk_dt_1    = clusts_na[ is.na(cluster) | (cluster %in% keep_cl) ]
  chk_dt_2    = rbind(train_dt, valid_dt, valid_rest, test_dt)
  assert_that( nrow(chk_dt_1) == nrow(chk_dt_2) )
  assert_that( all( sort(chk_dt_1$cell_id) == sort(chk_dt_2$cell_id) ) )

  return(list(
    train       = train_dt, 
    valid       = valid_dt, 
    valid_rest  = valid_rest, 
    test        = test_dt
  ))
}

.run_boost_watchlist <- function(train_dt, valid_dt, n_cores) {
  # convert training data to expected format
  assert_that( !is.null(train_dt$cluster) )
  train_vars  = colnames(train_dt) %>% setdiff(c("cluster", "cell_id"))
  assert_that( length(train_vars) > 0 )
  train_mat   = train_dt[, c("cell_id", train_vars), with = FALSE] %>% 
    as.matrix( rownames = "cell_id" )
  train_cl    = train_dt$cluster
  train_y     = train_cl %>% as.integer %>% `-`(1)
  # weights_v   = 1 / table(train_y) * 1000
  # weights_y   = weights_v[ train_y + 1 ]
  assert_that(
    min(train_y) == 0,
    max(train_y) + 1 == length(levels(train_dt$cluster))
  )

  # convert validation data to expected format
  valid_mat   = valid_dt[, c("cell_id", train_vars), with = FALSE] %>% 
    as.matrix( rownames = "cell_id" )
  valid_cl    = valid_dt$cluster
  valid_y     = valid_cl %>% as.integer %>% `-`(1)

  # blah
  dtrain      = xgb.DMatrix( data = train_mat, label = train_y )
  dvalid      = xgb.DMatrix( data = valid_mat, label = valid_y )
  watchlist   = list( train = dtrain, test = dvalid )

  # run boost
  boost_obj   = xgb.train( data = dtrain, watchlist = watchlist, 
    objective = "multi:softprob", num_class = max(train_y) + 1,
    nrounds = 100, early_stopping_rounds = 5, 
    nthread = n_cores, verbose = 2)

  # # run standard boost
  # boost_obj   = xgboost(
  #   data = train_mat, label = train_y, 
  #   weight = weights_y,
  #   objective = "multi:softprob", num_class = max(train_y) + 1,
  #   nthread = n_cores, nrounds = 10, verbose = 2)

  return(boost_obj)
}

.get_pred_valid <- function(boost_obj, valid_dt) {
  # get validation data
  train_vars  = colnames(valid_dt) %>% setdiff(c("cluster", "cell_id"))
  assert_that( length(train_vars) > 0 )
  valid_mat   = valid_dt[, c("cell_id", train_vars), with = FALSE] %>% 
    as.matrix( rownames = "cell_id" )

  # get probabilities for each predicted cluster
  probs_mat   = predict(boost_obj, valid_mat, reshape = TRUE)
  assert_that(
    length(levels(valid_dt$cluster)) == ncol(probs_mat),
    length(rownames(valid_mat)) == nrow(probs_mat)
  )
  probs_mat   = probs_mat %>% 
    set_colnames(levels(valid_dt$cluster)) %>% 
    set_rownames(rownames(valid_mat))

  # prediction for each cell
  pred_valid  = data.table(
    cell_id     = rownames(probs_mat),
    cl_true     = valid_dt$cluster,
    cl_pred     = colnames(probs_mat)[ apply(probs_mat, 1, which.max) ],
    p_pred      = apply(probs_mat, 1, max)
    ) %>% cbind(probs_mat)

  return(pred_valid)
}

.calc_confuse_xgboost_dt <- function(pred_valid) {
  n_cats    = length(unique(pred_valid$cl_true))
  confuse_dt  = pred_valid[, .N, by=.(cl_true, cl_pred)] %>%
    .[, N_true  := sum(N), by = cl_true] %>%
    .[, prop    := N / N_true ] %>%
    .[, logit   := qlogis((N+1) / (N_true + n_cats))] %>%
    .[, H       := -sum(prop*log2(prop)), by = cl_true ]

  return(confuse_dt)
}

.predict_on_new_data <- function(xgb_obj, hvg_mat, min_pred) {
  # unpack
  xgb_names   = xgb_obj$cl_lvls

  # get probabilities for each predicted cluster
  probs_mat   = predict(xgb_obj, as.matrix(hvg_mat), reshape = TRUE)
  assert_that(
    length(rownames(hvg_mat)) == nrow(probs_mat)
  )
  probs_mat   = probs_mat %>%
    set_colnames( xgb_names ) %>%
    set_rownames( rownames(hvg_mat) )

  # prediction for each cell
  preds_dt    = data.table(
    cell_id     = rownames(probs_mat),
    cl_pred_raw = colnames(probs_mat)[ apply(probs_mat, 1, which.max) ],
    p_pred      = apply(probs_mat, 1, max)
    ) %>% cbind(probs_mat) %>% 
    .[, cl_pred_naive := (p_pred > min_pred) %>% ifelse(cl_pred_raw, "unknown") %>% 
      factor(levels = c(xgb_names, "unknown"))  ]

  return(preds_dt)
}

.load_clusters <- function(cls_f) {
  cls_dt      = cls_f %>% fread(na.strings = "") %>% .[ !is.na(UMAP1) ]
  cl_cols     = colnames(cls_dt) %>% str_subset("RNA_snn_res")
  cls_dt      = cls_dt %>% 
    .[, c("sample_id", "cell_id", "UMAP1", "UMAP2", cl_cols), with = FALSE]

  return(cls_dt)
}

.apply_labels_by_cluster <- function(int_dt, preds_dt, min_cl_prop, min_cl_size) {
  # melt clusters
  non_cl_vars = c("sample_id", "cell_id", "UMAP1", "UMAP2")
  int_cls    = int_dt %>% 
    melt.data.table( id = non_cl_vars, var = "res_int", val = "cl_int")

  # exclude tiny clusters
  int_ns     = int_cls[, .(N_cl = .N), by = .(res_int, cl_int) ]
  keep_cls    = int_ns[ N_cl >= min_cl_size ]
  if ( nrow(keep_cls) > 0 ) {
    message("  excluding some clusters bc they are tiny:")
    int_ns[ N_cl < min_cl_size ] %>% .[ order(res_int, cl_int) ] %>% print
    int_cls  = int_cls %>% merge(int_ns, by = c("res_hmny", "cl_hmny")) %>% 
      .[, N_cl := NULL ]
  }

  # match these up to predictions, calculate proportions for each cluster
  match_dt    = preds_dt[, .(cell_id, cl_pred = cl_pred_naive)] %>% 
    merge(hmny_cls, by = "cell_id") %>% 
    .[, .N,                 by = .(res_hmny, cl_hmny, cl_pred)] %>% 
    .[, prop := N / sum(N), by = .(res_hmny, cl_hmny) ] %>% 
    setorder(res_hmny, cl_hmny, -prop)

  # take top prediction for each cluster
  match_lu    = match_dt[, .SD[1], by = .(res_hmny, cl_hmny)] %>% 
    .[ (cl_pred != "unknown") & (prop > min_cl_prop) ]

  # add these results to original cluster labels
  guesses_dt  = match_lu[, .(res_hmny, cl_hmny, cl_pred, prop_pred = prop)] %>% 
    merge(hmny_cls, by = c("res_hmny", "cl_hmny"), all.y = TRUE) %>% 
    merge(preds_dt[, .(cell_id, cl_pred_raw, cl_pred_naive, p_pred)], by = "cell_id") %>%
    setcolorder( c(non_cl_vars, "cl_pred_raw", "cl_pred_naive", "p_pred") ) %>% 
    dcast.data.table( sample_id + cell_id + UMAP1 + UMAP2 + 
      cl_pred_raw + cl_pred_naive + p_pred ~ res_hmny, 
      value.var = c("cl_hmny", "cl_pred", "prop_pred") )

  # broad_short     = c(`OPCs + COPs`='opc_cop', `Oligodendrocytes`='oligo', 
  #   `Astrocytes`='astro', `Microglia`='micro', `Excitatory neurons`='excit',
  #   `Inhibitory neurons`='inhib', `Endo + Peri`='endo_peri', `T cells`='t_cells', 
  #   `B cells`='b_cells', unknown = "?")
  # guesses_dt[, .(
  #     cl_hmny   = cl_hmny_RNA_snn_res.2, 
  #     pred_cl   = broad_short[ cl_pred_RNA_snn_res.2 ] %>% factor(levels = broad_short), 
  #     pred_raw  = broad_short[ cl_pred_raw ] %>% factor(levels = broad_short) %>% fct_drop
  #   ) ] %>% 
  #   .[, .N, by = .(cl_hmny, pred_cl, pred_raw)] %>%
  #   .[ !(pred_raw %in% c("t_cells", "b_cells")) ] %>% 
  #   dcast.data.table( cl_hmny + pred_cl ~ pred_raw, fill = 0 ) %>% 
  #   .[ order(-`?`) ]

  return(guesses_dt)
}

# code for Rmd
get_guesses_melt <- function(guesses_dt, res_ls, cl_lvls, min_cl_size) {
  # define some useful variables
  id_vars       = c("sample_id", "cell_id", "UMAP1", "UMAP2", 
    "cl_pred_raw", "p_pred", "cl_pred_naive")
  measure_vars  = c("cl_hmny", "cl_pred", "prop_pred")

  # split out by resolution
  guesses_melt  = guesses_dt[, -setdiff(id_vars, "cell_id"), with = FALSE ] %>%
    melt.data.table( id = "cell_id", 
      measure = patterns("^cl_hmny_", "^cl_pred_", "^prop_pred_") ) %>% 
    set_colnames(c("cell_id", "res_idx", measure_vars)) %>% 
    .[, res     := sort(res_ls)[ res_idx ] %>% factor(levels = res_ls) ] %>% 
    .[, cl_pred := cl_pred %>% factor(levels = cl_lvls) ]

  # exclude tiny clusters
  cl_ns         = guesses_melt[, .(N_cl = .N), by = .(res, cl_hmny)]
  keep_cls      = cl_ns[ N_cl >= min_cl_size ]
  guesses_melt  = guesses_melt %>% merge(keep_cls[, -"N_cl"], by = c("res", "cl_hmny"))

  # add useful labels back in
  guesses_melt  = guesses_dt[, id_vars, with = FALSE] %>% 
    merge(guesses_melt, id = "cell_id") %>% 
    .[, cl_pred_raw := cl_pred_raw %>% factor(levels = cl_lvls) ]

  return(guesses_melt)
}

calc_confuse_dt <- function(cl1_dt, cl2_dt, cl1, cl2) {
  assert_that( cl1 %in% names(cl1_dt) )
  assert_that( cl2 %in% names(cl2_dt) )
  assert_that( !is.null(levels(cl1_dt[[ cl1 ]])) )
  assert_that( !is.null(levels(cl2_dt[[ cl2 ]])) )
  assert_that( cl2 %in% names(cl2_dt) )
  assert_that(
    "cell_id" %in% names(cl1_dt),
    "cell_id" %in% names(cl2_dt)
  )

  confuse_dt = merge(
    cl1_dt[, .(cell_id, cl1 = get(cl1)) ], 
    cl2_dt[, .(cell_id, cl2 = get(cl2)) ], 
    by = "cell_id") %>% 
    .[, .N, by = .(cl1, cl2) ]
  # sort factor levels
  lvls_cl1    = confuse_dt$cl1 %>% levels
  if (is.null(lvls_cl1)) {
    lvls_cl1    = confuse_dt$cl1 %>% unique %>% sort
    confuse_dt  = confuse_dt[, cl1 := factor(cl1, levels = lvls_cl1)]
  }
  lvls_cl2    = confuse_dt$cl2 %>% levels
  if (is.null(lvls_cl2)) {
    lvls_cl2    = confuse_dt$cl2 %>% unique %>% sort
    confuse_dt  = confuse_dt[, cl2 := factor(cl2, levels = lvls_cl2)]
  }

  match_dt    = expand.grid(
    cl1  = unique(confuse_dt$cl1), 
    cl2  = unique(confuse_dt$cl2)
    )
  confuse_dt  = confuse_dt %>% 
    merge( match_dt, by = c("cl1", "cl2"), all = TRUE ) %>% 
    .[ is.na(N), N := 0 ] %>%
    .[, N0        := N + 1 ] %>%
    .[, log_N     := log(N0) ] %>%
    .[, p_cl1     := N0 / sum(N0), by = cl1 ] %>%
    .[, log_p_cl1 := log(p_cl1) ] %>% 
    .[, p_cl2     := N0 / sum(N0), by = cl2 ] %>%
    .[, log_p_cl2 := log(p_cl2) ]

  return(confuse_dt)
}

plot_cluster_comparison_heatmap <- function(confuse_dt, cl1, cl2, 
  plot_var = c("log_p_cl1", "p_cl1", "log_p_cl2", "p_cl2", "N", "log_N"),
  do_sort = c("no", "hclust", "seriate"), order_cl1 = NULL, order_cl2 = NULL, 
  lbl_threshold = 0.05) {
  # check inputs
  plot_var    = match.arg(plot_var)
  do_sort     = match.arg(do_sort)
  if (!is.null(order_cl1))
    assert_that( all(unique(confuse_dt$cl1) %in% order_cl1) )
  if (!is.null(order_cl2))
    assert_that( all(unique(confuse_dt$cl2) %in% order_cl2) )

  # don't make any changes!
  copy_dt     = copy(confuse_dt)
  if (!is.null(order_cl1))
    copy_dt     = copy_dt[, cl1 := factor(cl1, levels = order_cl1)]
  if (!is.null(order_cl2))
    copy_dt     = copy_dt[, cl2 := factor(cl2, levels = order_cl2)]

  # decide what to plot
  if (plot_var == "N") {
    data_wide   = dcast(copy_dt, cl1 ~ cl2, value.var = "N", fill = 0)
    value_name  = "no. cells"

    # define colours
    max_val     = copy_dt$N %>% max
    chunk_opts  = c(5e2, 1e3, 5e3, 1e4)
    best_chunk  = ((max_val / 10) / chunk_opts) %>% `<`(1) %>% which %>% min %>% chunk_opts[ . ]
    n_brks      = seq(0, ceiling(max_val / best_chunk) * best_chunk, by = best_chunk)
    n_labs      = n_brks
    mat_cols  = cols_fn(seq(min(n_brks), max(n_brks), best_chunk), 
      res = best_chunk, pal = "viridis", pal_dir = 1, range = "natural")
    lgd       = list(title = "no. cells\nin sample", 
      at = n_brks, labels = n_labs)

  } else if (plot_var == "log_N") {
    data_wide   = dcast(copy_dt, cl1 ~ cl2, value.var = "log_N")

    # define log breaks
    log_brks  = c(0, 3, 10, 30, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5) %>% 
      `+`(1) %>% log
    log_labs  = c("0", "3", "10", "30", "100", "300", 
      "1k", "3k", "10k", "30k", "100k", "300k")
    assert_that( length(log_brks) == length(log_labs) )

    # truncate to observed range
    max_val   = copy_dt$log_N %>% max
    assert_that( max_val <= max(log_brks) )
    max_brk   = (log_brks <= max_val) %>% which %>% max
    log_brks  = log_brks[seq.int(max_brk)]
    log_labs  = log_labs[seq.int(max_brk)]

    # get colours
    res       = 0.1
    mat_cols  = cols_fn(seq(min(log_brks), max_val, res), 
      res = res, pal = "viridis", pal_dir = 1, range = "natural")
    lgd       = list(title = "no. cells", 
      at = log_brks, labels = log_labs)

  } else if (plot_var == "p_cl1") {
    data_wide   = dcast(copy_dt, cl1 ~ cl2, value.var = "p_cl1")

    # define percentage breaks
    pct_brks  = seq(0, 1, 0.2)
    pct_labs  = sprintf("%.0f%%", 100 * pct_brks)
    res       = 0.1
    mat_cols  = cols_fn(seq(min(pct_brks), max(pct_brks), res), 
      res = res, pal = "viridis", pal_dir = 1, range = "natural")
    lgd       = list(title = "cluster\nproportion\n(rows sum to 1)", 
      at = pct_brks, labels = pct_labs)

  } else if (plot_var == "log_p_cl1") {
    data_wide   = dcast(copy_dt, cl1 ~ cl2, value.var = "log_p_cl1")

    # define colours
    pct_brks  = c(0.001, 0.002, 0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1) %>% log
    pct_labs  = c("0.1%", "0.2%", "0.4%", "1%", "2%", "4%", "10%", "20%", "40%", "100%")
    res       = 0.1
    mat_cols  = cols_fn(seq(min(pct_brks), max(pct_brks), res), 
      res = res, pal = "viridis", pal_dir = 1, range = "natural")
    lgd       = list(title = "cluster\nproportion\n(rows sum to 1)", 
      at = pct_brks, labels = pct_labs)

  } else if (plot_var == "p_cl2") {
    data_wide   = dcast(copy_dt, cl1 ~ cl2, value.var = "p_cl2")

    # define percentage breaks
    pct_brks  = seq(0, 1, 0.2)
    pct_labs  = sprintf("%.0f%%", 100 * pct_brks)
    res       = 0.1
    mat_cols  = cols_fn(seq(min(pct_brks), max(pct_brks), res), 
      res = res, pal = "viridis", pal_dir = 1, range = "natural")
    lgd       = list(title = "original\nproportion\n(cols sum to 1)", 
      at = pct_brks, labels = pct_labs)

  } else if (plot_var == "log_p_cl2") {
    data_wide   = dcast( copy_dt, cl1 ~ cl2, value.var = "log_p_cl2" )

    # define colours
    pct_brks  = c(0.001, 0.002, 0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1) %>% log
    pct_labs  = c("0.1%", "0.2%", "0.4%", "1%", "2%", "4%", "10%", "20%", "40%", "100%")
    res       = 0.1
    mat_cols  = cols_fn(seq(min(pct_brks), max(pct_brks), res), 
      res = res, pal = "viridis", pal_dir = 1, range = "natural")
    lgd       = list(title = "original\nproportion\n(cols sum to 1)", 
      at = pct_brks, labels = pct_labs)

  }

  # add text annotations
  if ( plot_var %in% c("p_cl1", "log_p_cl1", "p_cl2", "log_p_cl2")) {
    sel_cl    = str_extract(plot_var, "cl[0-9]")
    # define annotations
    txt_mat   = dcast( copy_dt, cl1 ~ cl2, value.var = paste0("p_", sel_cl) ) %>% 
      as.matrix(rownames = "cl1")
    .lbl_fn <- function(j, i, x, y, width, height, fill) {
      val     = txt_mat[i, j]
      if (val < lbl_threshold) {
        s       = ""
      } else {
        # n_dec   = ifelse( log10(val) > 1, 0, 1)
        # s       = paste0("%.", n_dec, "f%%") %>% sprintf(val)
        s       = sprintf("%.0f%%", 100 * val)
      }
      return(grid.text(s, x, y, gp = gpar(fontsize = 6)))
    }

  } else if ( plot_var %in% c("N", "log_N")) {
    # define annotations
    txt_mat   = dcast( copy_dt, cl1 ~ cl2, value.var = "N" ) %>% 
      as.matrix(rownames = "cl1")
    .lbl_fn <- function(j, i, x, y, width, height, fill)
      grid.text(sprintf("%s", signif(txt_mat[i, j], 2)), 
        x, y, gp = gpar(fontsize = 8))
  }
  
  # turn into matrix
  data_mat    = data_wide %>% as.matrix( rownames = "cl1" )

  # make data for annotations
  rows_dt     = copy_dt[, .(total_cl1 = sum(N)), by = cl1 ] %>%
    .[, log_total_cl1 := log(total_cl1) ] %>%
    setkey("cl1") %>% .[ rownames(data_mat) ]
  cols_dt     = copy_dt[, .(total_cl2 = sum(N)), by = cl2 ] %>%
    .[, log_total_cl2 := log(total_cl2) ] %>%
    setkey("cl2") %>% .[ colnames(data_mat) ]
  assert_that(
    nrow(data_mat) == nrow(rows_dt),
    ncol(data_mat) == nrow(cols_dt)
    )

  # do annotations
  # if (plot_var == "log_N") {
    log_brks  = c(0, 1, 10, 1e2, 1e3, 1e4, 1e5) %>% 
      `+`(1) %>% log
    log_labs  = c("0", "1", "10", "100", "1k", "10k", "100k")
    res       = 0.1
    log_cols  = cols_fn(seq(min(log_brks), max(log_brks), res), 
      res = res, pal = "magma", pal_dir = 1, range = "natural")

    # label with broad types
    row_annots  = rowAnnotation(
      `cl1 total`  = rows_dt$log_total_cl1,
      col = list(`cl1 total` = log_cols),
      annotation_name_side = "top",
      annotation_legend_param = list(
        `cl1 total` = list(at = log_brks, labels = log_labs)
        )
      )
    col_annots  = HeatmapAnnotation(
      `cl2 total`  = cols_dt$log_total_cl2,
      col = list(`cl2 total` = log_cols),
      annotation_name_side = "right",
      annotation_legend_param = list(
        `cl2 total` = list(at = log_brks, labels = log_labs)
        )
      )

  # } else {
  #   row_annots  = rowAnnotation(
  #     `orig. only`  = rows_dt$N_orig_only,
  #     col = list(
  #       `orig. only` = cols_fn(rows_dt$N_orig_only,
  #         res = max(rows_dt$N_orig_only) / 8, 
  #         pal = "Greys", pal_dir = 1, range = "has_zero")
  #        ),
  #     annotation_name_side = "top"
  #     )
  #   col_annots  = HeatmapAnnotation(
  #     `new only`  = cols_dt$N_cl_only,
  #     col = list(
  #       `new only`  = cols_fn(cols_dt$N_cl_only,
  #         res = max(cols_dt$N_cl_only) / 8, 
  #         pal = "Greys", pal_dir = 1, range = "has_zero")
  #       ),
  #     annotation_name_side = "right"
  #     )
  # }

  # # split by cell type
  # if (which_type == "type_fine") {
  #   orig_lvls   = fine_ord
  #   row_split   = fine_split[ rownames(data_mat) ] %>% factor(levels = broad_ord)
  # } else if (which_type == "type_broad") {
  #   orig_lvls   = broad_ord
  #   row_split   = NULL
  # }

  # # put in nice order, always order by log_N
  # set.seed(20220602)
  
  # # put matrix in nice order
  # order_dt    = copy_dt %>% 
  #   dcast.data.table(cl1 ~ cl2, value.var = "log_N", fill = 0) %>% 
  #   melt( id = "cl1", var = "cluster" ) %>% 
  #   .[, cl1 := factor(cl1, levels = orig_lvls) ] %>% 
  #   .[ order(cluster, cl1) ] %>%
  #   .[, smth_score := ksmooth(as.numeric(cl1), value, 
  #     kernel = "normal", x.points = as.numeric(cl1))$y, by = cluster ] %>%
  #   .[, .SD[ which.max(smth_score) ], by = cluster ] %>%
  #   .[ order(cl1) ]
  # assert_that( all( sort(order_dt$cluster) == colnames(data_mat) ) )
  # put in nice order

  if (do_sort == "no") {
    rows_flag   = FALSE
    cols_flag   = FALSE
    row_order   = NULL
    col_order   = NULL

  } else if (do_sort == "hclust") {
    # define vars
    rows_flag   = TRUE
    cols_flag   = TRUE
    row_order   = NULL
    col_order   = NULL

  } else if (do_sort == "seriate") {
    # do seriate
    data_min    = data_mat %>% as.vector %>% min(na.rm = TRUE)
    data_mat[is.na(data_mat)] = data_mat
    seriate_obj = seriate(data_mat - data_min, method = "BEA")

    # define vars
    rows_flag   = FALSE
    cols_flag   = FALSE
    row_order   = get_order(seriate_obj, 1)
    col_order   = get_order(seriate_obj, 2)
  }

  # heatmap
  hm_obj      = Heatmap(
    matrix = data_mat, col = mat_cols, 
    cell_fun = .lbl_fn,
    cluster_rows = rows_flag, cluster_columns = cols_flag,
    row_order = row_order, column_order = col_order,
    # row_names_gp = gpar(fontsize  =  8), column_names_gp = gpar(fontsize  =  8),
    row_title = cl1, column_title = cl2,
    left_annotation = row_annots, top_annotation = col_annots,
    heatmap_legend_param = lgd,
    row_names_side = "left", column_names_side = "top",
    na_col = "grey"
    )

  return(hm_obj)
}

plot_qc_by_cluster <- function(clusts_dt, qc_melt, x_lab) {
  plot_dt     = merge(qc_melt, clusts_dt, by = "cell_id")

  cl_lvls     = levels(plot_dt$cluster) %>% setdiff("unknown")
  cl_cols     = seq_along( cl_lvls ) %>% 
    rep(nice_cols, times = 10)[ . ] %>% 
    setNames( cl_lvls ) %>% c(unknown = "grey")

  # define breaks
  log_brks    = c(1e1, 2e1, 5e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4) %>%
    log10
  log_labs    = c("10", "20", "50", "100", "200", "500",
    "1k", "2k", "5k", "10k", "20k", "50k")
  logit_brks  = c(1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.10, 0.30,
    0.50, 0.70, 0.90, 0.97, 0.99) %>% qlogis
  logit_labs  = c("0.01%", "0.03%", "0.1%", "0.3%", "1%", "3%", "10%", "30%",
    "50%", "70%", "90%", "97%", "99%")
  splice_brks = seq(0.1, 0.9, 0.1) %>% qlogis
  splice_labs = (100*seq(0.1, 0.9, 0.1)) %>% as.integer %>% sprintf("%d%%", .)

  # plot
  g = ggplot(plot_dt) + aes( x = cluster, y = qc_val, fill = cluster ) +
    geom_violin() +
    scale_fill_manual( values = cl_cols, guide = "none" ) +
    facet_grid( qc_full ~ ., scales = 'free_y' ) +
    facetted_pos_scales(
      y = list(
        qc_full == "library size"    ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full == "no. of features" ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full == "mito pct"        ~
          scale_y_continuous(breaks = logit_brks, labels = logit_labs),
        qc_full == "spliced pct"     ~
          scale_y_continuous(breaks = splice_brks, labels = splice_labs)
        )
      ) +
    theme_bw() +
    theme(
      axis.text.x       = element_text( angle = 90, hjust = 1, vjust = 0.5 ),
      panel.grid        = element_blank(),
      strip.background  = element_rect( fill = 'white')
      ) +
    labs( x = x_lab, y = "QC metric" )

  return(g)
}

plot_cluster_entropies <- function(input_dt, what = c("norm", "raw")) {
  # check inputs
  what      = match.arg(what)

  # calculate mixing
  ns_dt     = input_dt %>% 
    .[, .N, by = .(sample_id, cluster) ] %>% 
    .[, p_sample  := N / sum(N), by = sample_id ] %>% 
    .[, p_cluster := N / sum(N), by = cluster ] %>% 
    .[, p_cl_norm := p_sample / sum(p_sample), by = cluster ]

  # calculate entropy
  entropy_dt  = ns_dt[, .(
      h_cl_raw      = -sum(p_cluster * log2(p_cluster), na.rm = TRUE), 
      h_cl_norm     = -sum(p_cl_norm * log2(p_cl_norm), na.rm = TRUE), 
      max_pct_raw   = 100 * max(p_cluster), 
      max_pct_norm  = 100 * max(p_cl_norm), 
      N             = sum(N)
    ), by = cluster ]
  labels_dt   = entropy_dt[ order(cluster) ]

  # get nice colours
  cl_ls     = entropy_dt$cluster %>% unique %>% sort
  cl_cols   = nice_cols[ seq_along(cl_ls) ] %>% setNames(cl_ls)

  # plot
  g = ggplot(entropy_dt) + 
    aes_string( x = paste0('h_cl_', what), y = paste0('max_pct_', what) ) +
    geom_smooth( method = "lm", formula = y ~ x, se = FALSE, colour = "grey" ) +
    geom_text_repel( data = labels_dt, aes(label = cluster), size = 3, 
      min.segment.length = 0, max.overlaps = Inf, box.padding = 0.5 ) +
    geom_point( shape = 21, aes( size = sqrt(N), fill = cluster ) ) +
    scale_x_continuous(breaks = pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = pretty_breaks(n = 3)) +
    scale_fill_manual( values = cl_cols, guide = "none" ) +
    expand_limits( y = 0 ) +
    scale_size(
      range   = c(1, 8),
      breaks  = c(2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4) %>% sqrt, 
      labels  = c('200', '500', '1k', '2k', '5k', '10k', '20k', '50k')
      ) +
    theme_bw() + 
    theme( panel.grid = element_blank() ) +
    labs(
      x     = 'entropy (high when clusters even across samples)',
      y     = 'max. pct. of one sample (high when concentrated in few samples)',
      size  = 'total # cells'
    )

  return(g)
}
