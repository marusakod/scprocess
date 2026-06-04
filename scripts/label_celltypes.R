
suppressPackageStartupMessages({
  library('RColorBrewer')
  library("BiocParallel")
  library('BiocStyle')
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
  library('anndataR')

  library('future')
  library('SingleCellExperiment')

  library('scater')
  library('Seurat')

  library('ComplexHeatmap')
  library('seriation')
  library('purrr')
  library('xgboost')
  library('ggrepel')

  library('Matrix')
  library('yaml')
  library('kableExtra')
})


label_with_xgboost_one_batch <- function(sel_batch, batch_var, model_name, xgb_f, xgb_cls_f,
  adata_f, pred_f) {
  # check inputs
  assert_that( file.exists(xgb_f) )

  # load XGBoost object
  message('  loading XGBoost classifier')
  xgb_obj     = readRDS(xgb_f)
  xgb_cls_dt  = fread(xgb_cls_f)
  hvgs        = variable.names(xgb_obj)
  
  # get values for these genes in new datasets
  message('  getting counts for HVGs')
  counts_mat  = .get_counts_mat(adata_f)
  hvg_mat     = .normalize_hvg_mat(counts_mat, hvgs)

  # predict for new data
  message('  predicting celltypes for all cells')
  preds_dt    = .predict_on_new_data(xgb_obj, xgb_cls_dt, hvg_mat)

  # add labels
  preds_dt    = preds_dt %>%
    .[, (batch_var) := sel_batch ] %>% 
    .[, labeller    := "scprocess"] %>% 
    .[, model       := model_name ]

  # save
  message('  saving results')
  fwrite(preds_dt, file = pred_f)
  message('done.')
}

.get_counts_mat <- function(adata_f) {
  sce    = read_h5ad(adata_f, as = 'SingleCellExperiment')
  counts = assay(sce, 'X')
  rownames(counts) = str_replace(rownames(counts), "_ENSG", "-ENSG")

  return(counts)  
}

.normalize_hvg_mat = function(counts_mat, hvgs, scale_f = 10000) {
  
  if (!all(hvgs %in% rownames(counts_mat))) {
    warning("not all HVGs present")
    missing_gs  = setdiff(hvgs, rownames(counts_mat))
    n_cols      = ncol(counts_mat)
    missing_mat = matrix(0, length(missing_gs), n_cols)
    rownames(missing_mat) = missing_gs
    colnames(missing_mat) = colnames(counts_mat)

    # convert to sparse, add to counts
    missing_mat = as(missing_mat, 'TsparseMatrix')
    counts_mat  = rbind(counts_mat, missing_mat)
  }

  hvg_mat   = counts_mat[hvgs, ]
  lib_sizes = colSums(counts_mat)

  norm_mat  = sweep(hvg_mat, 2, lib_sizes, FUN = "/")
  norm_mat  = norm_mat * scale_f
  log_mat   = log1p(norm_mat) %>% t()

  return(log_mat)
}

.predict_on_new_data <- function(xgb_obj, allow_dt, hvg_mat, chunk_size = 10000) {
  num_chunks = ceiling(nrow(hvg_mat) / chunk_size)
  idx_vec = rep(1:num_chunks, each = chunk_size, length.out = nrow(hvg_mat))
  cell_chunks= split(rownames(hvg_mat), idx_vec)

  probs_mat_ls = cell_chunks %>% lapply(function(cells_sub){
    sub_hvg_mat = hvg_mat[cells_sub, ] %>% as.matrix()
    # get probabilities for each cluster
    sub_probs_mat = predict(xgb_obj, sub_hvg_mat, reshape = TRUE)
    return(sub_probs_mat)
  })
  
  # merge all predictions
  probs_mat = do.call('rbind', probs_mat_ls)
    
  assert_that(
    nrow(hvg_mat) == nrow(probs_mat)
  )
  
  probs_mat   = probs_mat %>%
    set_colnames( allow_dt$cluster ) %>%
    set_rownames( rownames(hvg_mat) )

  # make data.table with predictions
  pred_dt     = data.table(
    cell_id         = rownames(probs_mat),
    predicted_label = colnames(probs_mat)[ apply(probs_mat, 1, which.max) ],
    probability     = apply(probs_mat, 1, max)
  )

  return(pred_dt)
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
    int_cls  = int_cls %>% merge(int_ns, by = c("res_int", "cl_int")) %>%
      .[, N_cl := NULL ]
  }

  # match these up to predictions, calculate proportions for each cluster
  match_dt    = preds_dt[, .(cell_id, cl_pred = cl_pred_naive)] %>%
    merge(int_cls, by = "cell_id") %>%
    .[, .N,                 by = .(res_int, cl_int, cl_pred)] %>%
    .[, prop := N / sum(N), by = .(res_int, cl_int) ] %>%
    setorder(res_int, cl_int, -prop)

  # take top prediction for each cluster
  match_lu    = match_dt[, .SD[1], by = .(res_int, cl_int)] %>%
    .[ (cl_pred != "ambiguous") & (prop > min_cl_prop) ]

  # add these results to original cluster labels
  guesses_dt  = match_lu[, .(res_int, cl_int, cl_pred, prop_pred = prop)] %>%
    merge(int_cls, by = c("res_int", "cl_int"), all.y = TRUE) %>%
    merge(preds_dt[, .(cell_id, cl_pred_raw, cl_pred_naive, p_pred)], by = "cell_id") %>%
    setcolorder( c(non_cl_vars, "cl_pred_raw", "cl_pred_naive", "p_pred") ) %>%
    dcast.data.table( sample_id + cell_id + UMAP1 + UMAP2 +
      cl_pred_raw + cl_pred_naive + p_pred ~ res_int,
      value.var = c("cl_int", "cl_pred", "prop_pred") )

  return(guesses_dt)
}

plot_umap_labels <- function(umap_dt, clust_dt, name) {
  # join umap and clusters
  assert_that(
    all(c('UMAP1', 'UMAP2') %in% names(umap_dt)),
    'cell_id' %in% names(umap_dt),
    'cell_id' %in% names(clust_dt),
    'cluster' %in% names(clust_dt)
  )

  # define cluster name
  plot_dt     = merge(umap_dt, clust_dt, by = 'cell_id', all.x = TRUE) %>%
    .[, .(
      UMAP1   = rescale(UMAP1, to = c(0.05, 0.95)),
      UMAP2   = rescale(UMAP2, to = c(0.05, 0.95)),
      cluster
      )] %>% .[, cluster := factor(cluster) %>% fct_infreq ]

  # tweak if "ambiguous"
  if ("ambiguous" %in% levels(plot_dt$cluster) )
    plot_dt     = plot_dt[, cluster := cluster %>% fct_relevel("ambiguous", after = Inf) ]

  # define colours
  cl_lvls     = levels(plot_dt$cluster)
  cl_cols     = seq_along( cl_lvls ) %>% rep(nice_cols, times = 10)[ . ] %>% setNames( cl_lvls )
  if ("ambiguous" %in% cl_lvls)
    cl_cols[[ "ambiguous" ]] = "grey"

  # make plot
  plot_dt     = plot_dt[ sample(.N, .N) ]
  g = ggplot(plot_dt) +
    aes( x = UMAP1, y = UMAP2, colour = cluster ) +
    geom_point(size = 0.1) +
    scale_colour_manual( values = cl_cols, guide = guide_legend(override.aes = list(size = 3)) ) +
    scale_x_continuous( breaks = pretty_breaks(), limits = c(0, 1) ) +
    scale_y_continuous( breaks = pretty_breaks(), limits = c(0, 1) ) +
    theme_bw() +
    theme( panel.grid = element_blank(), aspect.ratio = 1, 
      axis.ticks = element_blank(), axis.text = element_blank() ) +
    labs( colour = name )

  return(g)
}

calc_labels_table <- function(guesses_dt) {
  ns_dt   = guesses_dt %>%
    .[, .(cell_id, unagg = predicted_label_naive, agg = predicted_label_agg)] %>%
    melt(id = "cell_id", variable.name = "label_method", value.name = "label") %>%
    .[, .N, by = .(label_method, label)] %>%
    dcast( label ~ label_method, value.var = "N", fill = 0) %>%
    .[ order(-agg, -unagg) ] %>% setcolorder(c("label", "agg", "unagg"))
  # put in nice order
  if ("ambiguous" %in% ns_dt$label)
    ns_dt     = rbind( ns_dt[ label != "ambiguous" ], ns_dt[ label == "ambiguous" ] )
  # change names
  ns_dt     = ns_dt %>% set_colnames(c("predicted\nlabel", "no. cells,\naggregated", "no. cells, not\naggregated"))

  return(ns_dt)
}
