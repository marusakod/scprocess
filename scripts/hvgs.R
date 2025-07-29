# hvgs.R
suppressPackageStartupMessages({
  library("DESeq2")
})

calc_vst_obj <- function(pb, edger_dt) {
  # get counts
  empty_mat = assay(pb)
  suppressMessages({
    dds       = DESeqDataSetFromMatrix(countData = empty_mat,
      colData = data.frame(dummy = rep(1, ncol(empty_mat))), design = ~ 1)
  })

  # calculate vst
  vst_obj     = tryCatch({
    vst(dds, blind = TRUE)
    }, error = function(cond) {
    varianceStabilizingTransformation(dds, blind = TRUE)
    })

  # add edger stats
  edger_tmp   = copy(edger_dt) %>% 
    .[, .(gene_id, log2fc.empty = logFC, pval.empty = PValue, padj.empty = FDR)]
  row_dt      = rowData(vst_obj) %>% as.data.frame %>% as.data.table(keep.rownames = "gene_id")
  row_dt      = row_dt %>% merge(edger_tmp, by = "gene_id", all.x = TRUE) %>% 
    setkey("gene_id") %>% .[ rownames(vst_obj) ]
  rowData(vst_obj) = as(row_dt, "DataFrame")

  return(vst_obj)
}

plot_hvg_stats_vs_empty_log2fc <- function(hvgs_dt, edger_dt, n_top = 10) {
  col_vals  = c(
    hvg_clean = "#1965B0",
    hvg_dirty = "#DC050C",
    not_clean = "grey",
    not_dirty = "#FB8072"
  )
  status_labs = c(
    hvg_clean = "highly variable, \"clean\"",
    hvg_dirty = "highly variable, \"ambient\"",
    not_clean = "not variable, \"clean\"",
    not_dirty = "not variable, \"ambient\""
  )
  plot_dt   = merge(hvgs_dt, edger_dt, by = "gene_id") %>% 
    # .[, is_dirty  := (FDR < 0.01) & (logFC > 0) ] %>% 
    # .[ is.na(FDR), is_dirty := FALSE ] %>% 
    .[, .(gene_id, log2fc = logFC, padj = FDR, 
      hv_n = highly_variable_nbatches, mean_var = variances_norm,
      is_hvg = highly_variable, is_ambient, 
      status = ifelse(highly_variable, 
        ifelse(is_ambient, "hvg_dirty", "hvg_clean"),
        ifelse(is_ambient, "not_dirty", "not_clean") ) %>% factor(levels = names(status_labs))
      )]
  labels_dt = plot_dt[ order(-hv_n, -mean_var) ] %>% head(n_top) %>% 
    .[, symbol := gene_id %>% str_extract(".+(?=_ENS)") ]

  g = ggplot(plot_dt[ order(-status) ]) + 
    aes( x = log2fc, y = mean_var, colour = status ) +
    geom_hline( yintercept = 0, linewidth = 0.1, colour = 'grey20', alpha = 0.5 ) +
    geom_vline( xintercept = 0, linewidth = 0.1, colour = 'grey20', alpha = 0.5 ) +
    geom_point( size = 0.2, alpha = 0.5, show.legend = TRUE ) +
    geom_label_repel( data = labels_dt, aes( label = symbol ), 
      size = 3, max.overlaps = Inf, show.legend = FALSE ) +
    scale_x_continuous( breaks = pretty_breaks() ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    scale_colour_manual( values = col_vals, breaks = names(status_labs), labels = status_labs, 
      drop = FALSE, guide = guide_legend(override.aes = list(alpha = 1, size = 2)) ) +
    theme_classic( base_size = 14 ) +
    theme( plot.caption = element_text(vjust = -0.5) ) +
    labs(
      x       = "log2fc of \"empty\" drops vs all cells",
      y       = "mean standardized var. across samples",
      colour  = "HVG classification\nof gene",
      caption = sprintf("Labelled genes are top %d most variable.", n_top)
    )

  return(g)
}

plot_ambient_gene_calculations <- function(edger_dt, max_padj = 0.01, n_top = 10, 
  min_cpm_empty = 0) {
  # what is pval for max_padj?
  pval_1      = edger_dt[ FDR <= max_padj ] %>% .$PValue %>% max
  pval_2      = edger_dt[ FDR > max_padj ] %>% .$PValue %>% min
  max_pval    = mean(pval_1, pval_2)

  # take top 10 both x and y
  if (min_cpm_empty > 0) {
    tmp_dt      = edger_dt %>% .[ mean_logcpm.empties > log(min_cpm_empty + 1)]
  } else {
    tmp_dt      = copy(edger_dt)
  }
  
  # choose what to label
  labels_dt   = rbind(
    tmp_dt %>% .[ order(-abs(logFC)) ] %>% .[ 1:n_top ],
    tmp_dt %>% .[ order(PValue) ] %>% .[ 1:n_top ]
    ) %>% unique %>% 
    .[, symbol := gene_id %>% str_extract(".+(?=_ENS)") ]

  g = ggplot(tmp_dt) +
    aes( x = logFC, y = -log10(PValue) ) +
    geom_hline( yintercept = -log10(max_pval), colour = 'grey', linetype = "dashed" ) +
    geom_point( size = 0.1 ) +
    geom_label_repel( data = labels_dt, aes( label = symbol ), size = 3, max.overlaps = Inf ) +
    scale_x_continuous( breaks = pretty_breaks() ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    theme_classic( base_size = 14 ) +
    theme( plot.caption = element_text(vjust = -0.5) ) +
    labs(
      x       = "log2fc of \"empty\" drops vs all cells",
      y       = "-log10( nominal p-value of \"empty\" drops vs all cells )",
      caption = sprintf("Dotted horizontal line shows cutoff for adjusted p-value < %.1e.", 
        max_padj)
    )

  return(g)
}

plot_heatmap_of_ambient_profiles <- function(vst_mat, top_var = c("var", "mean", 
  "log2fc.empty", "pval.empty"), n_top = 30) {
  top_var     = match.arg(top_var)

  # calc
  vst_mat     = assay(vst_obj)
  if (top_var == "var") {
    top_vals    = rowVars(vst_mat)
  } else if (top_var == "mean") {
    top_vals    = rowMeans(vst_mat)
  } else if (top_var == "log2fc.empty") {
    top_vals    = rowData(vst_obj)$log2fc.empty %>% setNames(rownames(vst_obj))
  } else if (top_var == "pval.empty") {
    rows_dt     = rowData(vst_obj) %>% as.data.frame %>% as.data.table
    rows_dt     = rows_dt[ !is.na(log2fc.empty) ] %>% .[ log2fc.empty > 0 ]
    top_vals    = rows_dt$pval.empty %>% setNames(rows_dt$gene_id) %>% `*`(-1)
  }
  top_gs      = top_vals %>% sort %>% tail(n_top) %>% names

  # get these genes
  plot_mat    = vst_mat[ top_gs, ]
  row_labs    = top_gs %>% str_extract(".+(?=_ENS)")
  res         = 0.5
  mat_cols    = cols_fn(seq(max(5, min(plot_mat)), ceiling(max(plot_mat)), res),
    res = res, pal = "viridis", pal_dir = 1, range = "natural")
  lgd         = list(title = "log2cpm-like values\nfrom DESeq2::vst")

  # heatmap
  hm_obj      = Heatmap(
    matrix = plot_mat, col = mat_cols,
    row_labels = row_labs, row_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = lgd,
    row_names_side = "left", column_names_side = "top",
    na_col = "grey"
    )

  return(hm_obj)
}
