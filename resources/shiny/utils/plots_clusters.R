# CLUSTER-LEVEL PLOT FUNCTIONS -------------------------------------------------

# UMAP coloured by cluster identity, with filled-circle centroid markers.
# centroids:  data.table (cluster, UMAP_1, UMAP_2) — pre-computed from full data.
# repel_pos:  data.table (cluster, repel_1, repel_2) — low-density label positions.
# label_mode: "centroid" (default) places numbered circles on centroids;
#             "repel" places them in low-density space with a connecting segment.
# Returns a patchwork composite: UMAP on the left, numbered-circle legend on the right.
plot_cluster_umap <- function(umap_dt, col_pal, centroids,
                              repel_pos  = NULL,
                              label_mode = c("centroid", "repel")) {
  label_mode = match.arg(label_mode)

  centroids = copy(centroids) %>% .[order(cluster)]
  centroids[, label_num := seq_len(.N)]

  # Pick black or white number/stroke based on fill luminance (BT.709 coefficients)
  fill_hex = col_pal[as.character(centroids$cluster)]
  rgb_m    = col2rgb(fill_hex) / 255
  lum      = 0.2126 * rgb_m[1, ] + 0.7152 * rgb_m[2, ] + 0.0722 * rgb_m[3, ]
  centroids[, text_color := ifelse(lum > 0.5, "black", "white")]

  ct_dark  = centroids[text_color == "black"]
  ct_light = centroids[text_color == "white"]

  # -- UMAP panel (no legend) --------------------------------------------------
  umap_p = ggplot(umap_dt, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    geom_scattermore(pointsize = 2.5, pixels = c(1024, 1024)) +
    theme_classic() +
    labs(x = 'UMAP 1', y = 'UMAP 2') +
    theme(
      aspect.ratio = 1,
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      axis.title   = element_text(size = FONT_AXIS),
      plot.title   = element_text(size = FONT_TITLE)
    ) +
    scale_color_manual(values = col_pal, guide = 'none') +
    scale_fill_manual(values = col_pal,  guide = 'none')

  if (label_mode == "centroid" || is.null(repel_pos)) {
    # Numbered circles placed directly on centroids
    umap_p = umap_p +
      geom_point(data = ct_dark,  aes(fill = cluster),
                 shape = 21, size = 9, stroke = 0.8, color = "black") +
      geom_point(data = ct_light, aes(fill = cluster),
                 shape = 21, size = 9, stroke = 0.8, color = "white") +
      geom_text(data = ct_dark,  aes(label = label_num), size = 3.5, color = "black", fontface = "bold") +
      geom_text(data = ct_light, aes(label = label_num), size = 3.5, color = "white", fontface = "bold")
  } else {
    # Small centroid dot + segment to low-density label position
    rc = merge(centroids, repel_pos, by = "cluster")
    rc_dark  = rc[text_color == "black"]
    rc_light = rc[text_color == "white"]

    umap_p = umap_p +
      geom_point(data = centroids, aes(x = UMAP_1, y = UMAP_2),
                 shape = 16, size = 2, color = "grey30") +
      geom_segment(data = rc,
                   aes(x = UMAP_1, y = UMAP_2, xend = repel_1, yend = repel_2),
                   color = "grey30", linewidth = 0.35, lineend = "round") +
      geom_point(data = rc_dark,  aes(x = repel_1, y = repel_2, fill = cluster),
                 shape = 21, size = 9, stroke = 0.8, color = "black") +
      geom_point(data = rc_light, aes(x = repel_1, y = repel_2, fill = cluster),
                 shape = 21, size = 9, stroke = 0.8, color = "black") +
      geom_text(data = rc_dark,  aes(x = repel_1, y = repel_2, label = label_num),
                size = 3.5, color = "black", fontface = "bold") +
      geom_text(data = rc_light, aes(x = repel_1, y = repel_2, label = label_num),
                size = 3.5, color = "white", fontface = "bold")
  }

  # Size the legend on the active graphics device, independently of the UMAP.
  device_size = grDevices::dev.size("in")
  legend = cluster_umap_legend(centroids, col_pal, device_size)
  patchwork::wrap_plots(
    umap_p, patchwork::wrap_elements(full = legend$grob),
    widths = grid::unit(c(1, legend$width), c("null", "in"))
  )
}

# Keep specification order down columns. Wrapped labels contribute their full
# line count when deciding whether one column fits at the current device size.
cluster_umap_legend <- function(centroids, col_pal, device_size) {
  n_labels = nrow(centroids)
  font_size = max(9, min(12, round(12 - 0.06 * max(0, n_labels - 15))))
  legend_width = max(1.5, device_size[1] * 0.45)
  available_height = max(1, device_size[2] - 0.65)
  labels = as.character(centroids$cluster)

  max_columns = max(1, min(n_labels, floor(legend_width / 1.4)))
  fits = FALSE
  for (font_size in seq(font_size, 9)) {
    for (n_columns in seq_len(max_columns)) {
      column_width = legend_width / n_columns
      wrap_chars = max(8, floor((column_width - 0.35) * 72 / (font_size * 0.52)))
      wrapped = stringr::str_wrap(labels, width = wrap_chars)
      row_heights = pmax(0.20, (stringr::str_count(wrapped, "\\n") + 1) *
                          font_size / 72 * 1.25) + 0.045
      column = ceiling(seq_len(n_labels) / ceiling(n_labels / n_columns))
      column_heights = tapply(row_heights, column, sum)
      fits = max(column_heights) <= available_height
      if (fits) break
    }
    if (fits) break
  }

  # Vertically centre the legend; all positions and sizes are in physical units.
  top = (device_size[2] + max(column_heights)) / 2 - 0.15
  y = numeric(n_labels)
  for (i in unique(column)) {
    rows = which(column == i)
    y[rows] = top - cumsum(row_heights[rows]) + row_heights[rows] / 2
  }
  x = (column - 1) * column_width + 0.12
  colours = unname(col_pal[labels])
  children = list(
    grid::textGrob("Cluster", x = grid::unit(0.02, "in"),
                   y = grid::unit(top + 0.24, "in"), just = "left",
                   gp = grid::gpar(fontsize = font_size + 1, fontface = "bold")),
    grid::pointsGrob(x = grid::unit(x, "in"), y = grid::unit(y, "in"),
                     pch = 21, size = grid::unit(0.18, "in"),
                     gp = grid::gpar(fill = colours, col = colours)),
    grid::textGrob(centroids$label_num, x = grid::unit(x, "in"),
                   y = grid::unit(y, "in"),
                   gp = grid::gpar(fontsize = min(8, font_size - 1),
                                   col = centroids$text_color)),
    grid::textGrob(wrapped, x = grid::unit(x + 0.15, "in"),
                   y = grid::unit(y, "in"), just = "left",
                   gp = grid::gpar(fontsize = font_size, lineheight = 1.25))
  )
  list(grob = grid::gTree(children = do.call(grid::gList, children)),
       width = legend_width, n_columns = n_columns, font_size = font_size)
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
