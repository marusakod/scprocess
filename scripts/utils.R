# utils.R
suppressPackageStartupMessages({
  library('magrittr')
  library('data.table')
  library('assertthat')
  library('forcats')
  library('stringr')
  library('rhdf5')

  library('RColorBrewer')
  library('circlize')
  library('viridis')
  library('scales')

  library('ggplot2')
  library("ggbeeswarm")
  library('ggrepel')
  library('ggh4x')
  library("patchwork")

  library('ComplexHeatmap')
  library("seriation")
  library("strex")

  library('readxl')
})

# nice_cols   = CATALYST:::.cluster_cols
nice_cols   = c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6",
  "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A",
  "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4",
  "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3",
  "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
  )

# function for making nice colours
cols_fn <- function(mat, res, pal, pal_dir, range=c('has_zero', 'symmetric', 'natural')) {
  # check inputs
  stopifnot(pal_dir %in% c(-1,1))
  range     = match.arg(range)
  mat       = mat[!is.na(mat)]
  stopifnot(length(mat) > 0)

  # check values
  n_vals    = length(unique(mat))
  if (n_vals == 1) {
    pal_cols  = .get_pal_cols(pal, 3)
    cols    = pal_cols[[1]]
    names(cols) = unique(mat)
    return(cols)
  }

  # check inputs
  sgn     = 1
  if (range=='has_zero') {
    assert_that( all(mat >= 0) | all(mat <= 0) )
    if (all(mat <= 0) ) {
      sgn   = -1
      mat   = mat * -1
    }
    max_val   = mat %>% as.vector %>% max %>% `/`(res) %>% ceiling %>% `*`(res)
    max_val   = round(max_val, 3)
    min_val   = 0
  } else if (range=='symmetric') {
    max_val   = mat %>% as.vector %>% abs %>% max %>% `/`(res) %>% ceiling %>% `*`(res)
    max_val   = round(max_val, 3)
    min_val   = -max_val
  } else if (range=='natural') {
    max_val   = mat %>% as.vector %>% max %>% `/`(res) %>% ceiling %>% `*`(res)
    min_val   = mat %>% as.vector %>% min %>% `/`(res) %>% floor %>% `*`(res)
    max_val   = round(max_val, 3)
    min_val   = round(min_val, 3)
  }

  # define colours
  seq_vals  = seq(min_val, max_val, res)
  n_cols    = length(seq_vals)
  pal_cols  = .get_pal_cols(pal, n_cols)
  if ( length(pal_cols) < n_cols ) {
    n_cols    = length(pal_cols)
    seq_vals  = seq(min_val, max_val, length.out=n_cols)
  }

  # define colours
  if (pal_dir == -1)
    pal_cols  = rev(pal_cols)

  # make colour function
  if (sgn == 1)
    cols  = colorRamp2(seq_vals, pal_cols)
  else
    cols  = colorRamp2(-seq_vals, rev(pal_cols))

  return(cols)
}

.get_pal_cols <- function(pal_str = c('viridis', 'magma', 'Blues', 'BuGn',
  'BuPu', 'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd', 'PuBu', 'PuBuGn',
  'PuRd', 'Purples', 'RdPu', 'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd',
  'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn',
  'Spectral', 'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2',
  'Set1', 'Set2', 'Set3'), n_cols) {
  pal_str   = match.arg(pal_str)
  if ( pal_str == 'viridis' ) {
    pal_cols  = str_remove(viridis(n_cols), 'FF$')
  } else if ( pal_str == 'magma' ) {
    pal_cols  = str_remove(magma(n_cols), 'FF$')
  } else {
    suppressWarnings({pal_cols = brewer.pal(n_cols, pal_str)})
  }

  return(pal_cols)
}

load_gene_biotypes <- function(gtf_dt_f) {
    # .[, .(ensembl_id = gene_id, symbol = gene_name, gene_type)] %>%
    # unique %>%
    # .[, .(gene_id = paste0(symbol, '_', ensembl_id), symbol, gene_type)]
  biotypes_dt = fread(gtf_dt_f) %>%
    .[, .(gene_id, symbol, gene_type)]
  assert_that( all(table(biotypes_dt$gene_id) == 1) )

  return(biotypes_dt)
}

print_top_markers <- function(mkrs_dt, min_cpm = 50, top_n = 10, max_fdr = 0.05, 
  order_by    = c("pval", "lfc")) {
  order_by    = match.arg(order_by)
  show_dt     = mkrs_dt[ !str_detect(gene_type, "(lincRNA|lncRNA|pseudogene)") ] %>%
    .[ (FDR < max_fdr) & (logFC > 0) & (logcpm.sel > log(min_cpm + 1)) ]
  if (order_by == "pval") {
    show_dt     = show_dt %>% 
      .[ order(cluster, PValue, -logFC) ]
  } else if (order_by == 'lfc') {
    show_dt     = show_dt %>% 
      .[ order(cluster, -logFC) ]    
  }
  show_dt %>% 
    .[, .SD[ 1:min(top_n, .N) ], by = cluster ] %>% 
    .[, .(cluster, symbol, 
      CPM = logcpm.sel %>% exp %>% `-`(1) %>% signif(2) %>% round,
      log2fc = logFC %>% round(1), FDR = FDR %>% signif(2)
    )]
}

print_sel_gene <- function(mkrs_dt, sel_ls) {
  mkrs_dt %>% 
    .[ symbol %in% sel_ls ] %>% .[, symbol := symbol %>% factor(sel_ls) ] %>% 
    .[ order(symbol, log(PValue)*sign(logFC)) ] %>% 
    .[, .(symbol, gene_type, cluster, 
      CPM = logcpm.sel %>% exp %>% `-`(1) %>% signif(2) %>% round,
      log2fc = logFC %>% round(1), FDR = FDR %>% signif(2)
    )]
}

.get_h5_mx <- function(af_mat_f, sel_s) {
  # get this file
  h5_filt   = H5Fopen(af_mat_f, flags = "H5F_ACC_RDONLY" )

  # get indices of barcodes
  mat       = sparseMatrix(
    i = as.vector(h5_filt$matrix$indices +1),
    p = as.vector(h5_filt$matrix$indptr),
    x = as.vector(h5_filt$matrix$data),
    repr = "C",
    dims = h5_filt$matrix$shape
  ) %>% as("TsparseMatrix")

  # add names
  bcs           = h5_filt$matrix$barcodes
  colnames(mat) = paste0(sel_s, bcs)
  rownames(mat) = as.character(h5_filt$matrix$features$name)

  return(mat)
}

# sum spliced, unspliced and ambiguous counts for same gene
.sum_SUA <- function(sua_mat) {
  types = c('_S$', '_U$', '_A')
  
  mats  = lapply(types, function(t) sua_mat[grepl(t, rownames(sua_mat)), ])
  # check if symbols are all in the same order
  genes = lapply(mats, function(m) rownames(m) %>% str_before_first(pattern = '_'))

  assert_that(
    sapply(genes, function(gs) identical(gs, genes[[1]])) %>% all(), 
    msg = "gene names in split matrices don't match"
  )

  # remove suffixes from rownames
  mats = lapply(mats, function(m) {
    rownames(m) = str_before_first(rownames(m), pattern = '_')
    m
  })
  mats_sum = mats[[1]] + mats[[2]] + mats[[3]]
  return(mats_sum)
}


calc_confuse_dt <- function(cl1_dt, cl2_dt, cl1, cl2, min_cl2_p = NULL) {
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

  # aggregate if requested
  if (!is.null(min_cl2_p)) {
    cl1_max_ps  = copy(confuse_dt) %>% .[, p_cl2 := N / sum(N), by = cl2 ] %>%
      .[, .(max_p = max(p_cl2)), by = cl1]
    cl1_to_exc  = cl1_max_ps[ max_p < min_cl2_p ]$cl1 %>% as.character
    confuse_dt  = confuse_dt %>%
      .[ !(cl1 %in% cl1_to_exc) ] %>% .[, cl1 := cl1 %>% fct_drop ]
  }

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
    .[, p_cl1     := N / sum(N), by = cl1 ] %>%
    .[, log_p_cl1 := log(N0 / sum(N0)), by = cl1 ] %>%
    .[, p_cl2     := N / sum(N), by = cl2 ] %>%
    .[, log_p_cl2 := log(N0 / sum(N0)), by = cl2 ]

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
    data_mat[is.na(data_mat)] = data_min
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
    row_title = cl1, column_title = cl2,
    left_annotation = row_annots, top_annotation = col_annots,
    heatmap_legend_param = lgd,
    row_names_side = "left", column_names_side = "top",
    na_col = "grey"
    )

  return(hm_obj)
}
