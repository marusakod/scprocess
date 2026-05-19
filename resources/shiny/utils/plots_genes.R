# GENE EXPRESSION PLOT FUNCTIONS -----------------------------------------------

# UMAP coloured by expression of a single gene.
plot_gene_umap <- function(gene_umap_dt, gene, pal_name) {
  palette = palette_map[[pal_name]]

  if (max(gene_umap_dt$gene_exp) == 0) {
    orig_brks_r = 0
    trans       = log1p(0)
  } else {
    trans_brks  = pretty(gene_umap_dt$gene_exp, n = 5)
    orig_brks_r = unique(c(c(1, 2, 5), round(expm1(trans_brks), digits = -1)))
    trans       = log1p(orig_brks_r)
  }

  p = ggplot(gene_umap_dt, aes(x = UMAP_1, y = UMAP_2, color = gene_exp)) +
    geom_point(size = 1.5, shape = 16) +
    theme_classic() +
    labs(
      x     = 'UMAP 1', y = 'UMAP 2',
      color = "Reads per 10k",
      title = paste0("Selected gene: ", gene)
    ) +
    theme(
      axis.text          = element_blank(),
      axis.ticks         = element_blank(),
      axis.title         = element_text(size = FONT_AXIS),
      legend.text        = element_text(size = FONT_SMALL),
      legend.title       = element_text(size = FONT_TEXT),
      plot.title         = element_text(size = FONT_TITLE),
      legend.key.width   = unit(0.8, "cm"),
      legend.key.height  = unit(1.5, "cm")
    )

  if (!pal_name %in% c('greytored', 'greytopurple')) {
    p = p + scale_color_viridis_c(
      option    = palette$option,
      direction = palette$direction,
      breaks    = trans,
      labels    = orig_brks_r
    )
  } else {
    p = p + scale_color_gradientn(
      colours = palette$cols,
      breaks  = trans,
      labels  = orig_brks_r
    )
  }

  return(p)
}

# Pseudobulk dotplot for a single gene, with optional subsetting, faceting,
# x-axis reordering, and FDR significance annotations.
plot_gene_dotplot <- function(dotplot_dt, subset, x_axis, fill, fill_name, fill_pal,
                              facet_by, gene, anno_dt = NULL, note = NULL, order = FALSE) {
  # apply subsetting
  for (nn in names(subset)) {
    if (!is.null(subset[[ nn ]])) {
      keep_idx   = dotplot_dt[[ nn ]] %in% subset[[ nn ]]
      dotplot_dt = dotplot_dt[ keep_idx ]
    }
  }
  if ('cluster' %in% names(subset))
    dotplot_dt = dotplot_dt[, cluster := fct_drop(cluster)]

  show_ps = (!is.null(anno_dt)) & (x_axis == 'cluster') & (facet_by %in% c("none", ""))

  # compute axis order
  means = dotplot_dt %>%
    .[, .(mean_gene_exp = mean(gene_exp)), by = x_axis] %>%
    tibble::deframe() %>%
    sort(decreasing = TRUE)
  ord = if (order) names(means) else levels(dotplot_dt[[x_axis]])

  p = ggplot(dotplot_dt,
             aes(x    = forcats::fct_relevel(get(x_axis), ord),
                 y    = log(gene_exp + 10),
                 fill = get(fill),
                 size = n_cells)) +
    geom_quasirandom(colour = 'black', shape = 21) +
    scale_y_continuous(breaks = LOG_BRKS, labels = LOG_LABS) +
    expand_limits(y = log(10 + 100)) +
    scale_size_continuous(trans = 'identity', range = c(1, 8), breaks = pretty_breaks()) +
    expand_limits(size = 0)

  if (show_ps) {
    anno_dt  = anno_dt[ cluster %in% levels(dotplot_dt$cluster) ]
    range_y  = range(log(dotplot_dt$gene_exp + 10), na.rm = TRUE)
    buff     = (range_y[2] - range_y[1]) * 0.05
    anno_ypos           = tapply(log(dotplot_dt$gene_exp + 10), dotplot_dt$cluster, max) + buff
    assert_that( all(names(anno_ypos) == anno_dt$cluster) )
    anno_dt$ypos        = anno_ypos

    p = p + geom_text(
      data        = anno_dt,
      aes(x = forcats::fct_relevel(cluster, ord), y = ypos, label = label),
      inherit.aes = FALSE, size = 5
    )
  }

  if (!is.null(facet_by) & !(facet_by %in% c("none", ""))) {
    p = p + facet_grid(as.formula(sprintf("~ %s", facet_by)), scales = 'free_x', space = 'free_x')
  }

  p = p +
    theme_scjoin_dots() +
    labs(
      x        = NULL,
      y        = "counts per million",
      size     = "# cells",
      fill     = fill_name,
      title    = paste0("Selected gene: ", gene),
      subtitle = note
    ) +
    scale_fill_manual(values = fill_pal) +
    guides(
      fill = guide_legend(order = 1, override.aes = list(size = 4)),
      size = guide_legend(order = 2)
    )

  if (show_ps) {
    p = p +
      labs(caption = "Annotations show adjusted p-values of expression of gene in cluster relative to average across all other clusters.") +
      theme(plot.caption = element_text(size = 10, colour = "grey60"))
  }

  return(p)
}
