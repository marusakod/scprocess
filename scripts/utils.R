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

.get_alevin_mx <- function(af_mat_f, sel_s) {
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
