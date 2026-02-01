
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
  library('zellkonverter')

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
  preds_dt    = .predict_on_new_data(xgb_obj, xgb_cls_dt, hvg_mat, min_pred)

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
  sce    = readH5AD(adata_f)
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

.predict_on_new_data <- function(xgb_obj, allow_dt, hvg_mat, min_pred, chunk_size = 10000) {
  browser()
  # predict on chunks of cells for efficiency
  num_chunks = ceiling(nrow(hvg_mat/chunk_size))
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

# code for Rmd
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
        x, y, gp = gpar(fontsize = 10))
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
      annotation_name_side = "bottom",
      annotation_legend_param = list(
        `cl1 total` = list(at = log_brks, labels = log_labs)
        )
      )
    col_annots  = HeatmapAnnotation(
      `cl2 total`  = cols_dt$log_total_cl2,
      col = list(`cl2 total` = log_cols),
      annotation_name_side = "left",
      annotation_legend_param = list(
        `cl2 total` = list(at = log_brks, labels = log_labs)
        )
      )
    
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

plot_umap_cluster <- function(umap_dt, clust_dt, name) {
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
