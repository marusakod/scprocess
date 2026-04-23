# CLUSTER-LEVEL PLOT FUNCTIONS -------------------------------------------------

# UMAP coloured by cluster identity, with centroid labels.
plot_cluster_umap <- function(umap_dt, col_pal) {
  centroids = copy(umap_dt) %>%
    .[, `:=`(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = cluster] %>%
    unique()

  ggplot(umap_dt, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    geom_scattermore(pointsize = 2.5, pixels = c(1024, 1024)) +
    theme_classic() +
    labs(x = 'UMAP 1', y = 'UMAP 2', color = "cluster") +
    theme(
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      axis.title   = element_text(size = FONT_AXIS),
      legend.text  = element_text(size = FONT_SMALL),
      legend.title = element_text(size = FONT_TEXT),
      plot.title   = element_text(size = FONT_TITLE)
    ) +
    guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) +
    geom_text(data = centroids, aes(label = cluster),
              size = 4, vjust = 2, color = "black", fontface = 'bold') +
    scale_color_manual(values = col_pal, guide = 'none')
}

# PAGA graph: nodes sized by cell count, edges weighted by connectivity.
# col_pal: named colour vector for clusters (names = cluster levels).
plot_cluster_paga <- function(pos_dt, mat_dt, weight_cut = 0.2, col_pal) {
  cl_ls = levels(mat_dt$cluster)
  assert_that(all( mat_dt$cluster == colnames(mat_dt)[!colnames(mat_dt) %in% c('cluster', 'n_cells')] ))
  assert_that(all( mat_dt$cluster == pos_dt$cluster ))

  edges_dt = copy(mat_dt) %>%
    .[, c('cluster', cl_ls), with = FALSE] %>%
    setnames('cluster', 'cluster_i') %>%
    melt(id = 'cluster_i', variable.name = 'cluster_j', value.name = 'weight') %>%
    merge(pos_dt[, .(paga1_i = paga1, paga2_i = paga2, cluster_i = cluster)],
          by = 'cluster_i', all.y = TRUE, all.x = FALSE) %>%
    merge(pos_dt[, .(paga1_j = paga1, paga2_j = paga2, cluster_j = cluster)],
          by = 'cluster_j', all.y = TRUE, all.x = FALSE)

  pos_dt = merge(pos_dt, mat_dt[, .(cluster, n_cells)], by = 'cluster',
                 all.x = TRUE, all.y = FALSE)

  ggplot() +
    geom_segment(
      data = edges_dt[weight > weight_cut, ],
      aes(x = paga1_i, xend = paga1_j, y = paga2_i, yend = paga2_j,
          alpha = weight, linewidth = weight),
      lineend = 'round'
    ) +
    geom_point(
      data = pos_dt,
      aes(x = paga1, y = paga2, fill = cluster, size = sqrt(n_cells)),
      colour = 'black', shape = 21
    ) +
    scale_fill_manual(values = col_pal, guide = "none") +
    scale_linewidth_continuous(breaks = scales::pretty_breaks(), range = c(weight_cut, 4), limits = c(weight_cut, 1)) +
    scale_alpha_continuous(breaks = scales::pretty_breaks(), range = c(weight_cut, 1),  limits = c(weight_cut, 1)) +
    scale_size_continuous(
      range  = c(1, 12),
      breaks = c(1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 5e4) %>% sqrt,
      labels = c('100', '300', '1k', '3k', '10k', '30k', '50k')
    ) +
    ggrepel::geom_label_repel(
      data = pos_dt,
      aes(x = paga1, y = paga2, label = cluster),
      size = 4, max.overlaps = Inf, force_pull = 0.5
    ) +
    theme_bw() +
    theme(
      axis.text    = element_blank(),
      panel.grid   = element_blank(),
      axis.ticks   = element_blank(),
      axis.title   = element_text(size = FONT_AXIS),
      legend.title = element_text(size = FONT_TEXT),
      legend.text  = element_text(size = FONT_SMALL)
    ) +
    labs(
      x           = 'PAGA 1', y = 'PAGA 2',
      fill        = 'cluster',
      size        = '# cells\nin cluster',
      alpha       = 'PAGA\nconnectivity',
      linewidth   = 'PAGA\nconnectivity'
    )
}

# Cell function for heatmap significance annotations (.lbl_fn).
# Placed at module level so it can be reused by both log2fc and logcpm branches.
.heatmap_lbl_fn <- function(fdr_mat) {
  function(j, i, x, y, width, height, fill)
    grid.text(sprintf("%s", fdr_mat[i, j]), x, y, gp = gpar(fontsize = 9))
}

# Heatmap of marker genes across clusters.
# show: 'log2fc' (default) or 'logcpm' — controls which matrix and legend are drawn.
# Row and column dendrograms are computed on the log2fc matrix; when show='logcpm'
# the same ordering is applied to the logcpm matrix so gene order is consistent.
plot_heatmap <- function(mkrs_dt, hmap_gs, max_fc = 4,
                         cluster_rows = FALSE, cluster_columns = FALSE,
                         show = 'log2fc') {
  show = match.arg(show, c('log2fc', 'logcpm'))

  mkrs_sel = mkrs_dt[symbol %in% hmap_gs]

  # build matrices once
  log2fc_mat = mkrs_sel[, .(symbol, cluster, log2fc)] %>%
    data.table::dcast(symbol ~ cluster, value.var = 'log2fc') %>%
    tibble::column_to_rownames('symbol') %>%
    as.matrix() %>%
    .[hmap_gs, ]

  fdr_mat = mkrs_sel %>%
    .[, signif := fcase(
      FDR <= 0.001 & log2fc > 0, '***',
      FDR <= 0.01  & log2fc > 0, '**',
      FDR <= 0.05  & log2fc > 0, '*',
      default = '')] %>%
    .[, .(symbol, cluster, signif)] %>%
    data.table::dcast(symbol ~ cluster, value.var = 'signif') %>%
    tibble::column_to_rownames('symbol') %>%
    as.matrix() %>%
    .[hmap_gs, ]

  logcpm_mat = mkrs_sel %>%
    .[, logcpm := log(CPM + 10)] %>%
    .[, .(symbol, cluster, logcpm)] %>%
    data.table::dcast(symbol ~ cluster, value.var = 'logcpm') %>%
    tibble::column_to_rownames('symbol') %>%
    as.matrix() %>%
    .[hmap_gs, ]

  # common heatmap params
  hm_params = list(
    cluster_rows      = cluster_rows,
    cluster_columns   = cluster_columns,
    column_names_gp   = gpar(fontsize = FONT_SMALL),
    row_names_gp      = gpar(fontsize = FONT_SMALL),
    show_heatmap_legend = FALSE,
    row_names_side    = "left",
    column_names_side = "top",
    row_dend_side     = "right",
    column_dend_side  = "top",
    na_col            = "grey"
  )

  # choose matrix/colours/legend based on show
  if (show == 'log2fc') {
    fc_cols  = cols_fn(seq(-max_fc, max_fc, 1), res = 0.1, pal = "RdBu", pal_dir = -1, range = "natural")
    lgd      = Legend(
      title      = "log2fc in\ncluster",
      col_fun    = fc_cols,
      at         = c(-max_fc, 0, max_fc),
      title_gp   = gpar(fontsize = FONT_SMALL, fontface = 'bold'),
      labels_gp  = gpar(fontsize = 10),
      grid_width = unit(0.5, "cm"),
      legend_width  = unit(6, "cm"),
      legend_height = unit(4, "cm")
    )
    hm_obj = do.call(Heatmap, c(list(matrix = log2fc_mat, col = fc_cols,
                                     cell_fun = .heatmap_lbl_fn(fdr_mat)), hm_params))
  } else {
    # compute row/col order from log2fc heatmap so gene order is preserved
    ref_hm  = do.call(Heatmap, c(list(matrix = log2fc_mat,
                                      col = cols_fn(seq(-max_fc, max_fc, 1), res = 0.1,
                                                    pal = "RdBu", pal_dir = -1, range = "natural"),
                                      cell_fun = .heatmap_lbl_fn(fdr_mat)), hm_params))
    row_ord = if (cluster_rows)    row_dend(ref_hm)    else FALSE
    col_ord = if (cluster_columns) column_dend(ref_hm) else FALSE

    cpm_cols = cols_fn(seq(log(0 + 10), log(10000 + 10), 0.1), res = 0.1,
                       pal = "viridis", pal_dir = 1, range = "natural")
    lgd      = Legend(
      title      = "avg. CPM\nin cluster",
      col_fun    = cpm_cols,
      at         = log(c(0, 10, 100, 1000, 10000) + 10),
      labels     = c('0', '10', '100', '1k', '10k'),
      title_gp   = gpar(fontsize = FONT_SMALL, fontface = 'bold'),
      labels_gp  = gpar(fontsize = 10),
      grid_width = unit(0.5, "cm"),
      legend_width  = unit(6, "cm"),
      legend_height = unit(4, "cm")
    )
    hm_obj = do.call(Heatmap, c(list(matrix = logcpm_mat, col = cpm_cols,
                                     cluster_rows = row_ord,
                                     cluster_columns = col_ord), hm_params[!names(hm_params) %in%
                                       c('cluster_rows', 'cluster_columns', 'cell_fun')]))
  }

  draw(hm_obj, padding = unit(c(5, 5, 5, 10), "mm"), annotation_legend_list = lgd)
}
