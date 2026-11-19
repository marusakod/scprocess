
suppressPackageStartupMessages({
  library('magrittr')
  library('data.table')
  library('assertthat')
  library('scales')
  library('ggplot2')
  library('forcats')
})


plot_umap_labels <- function(umap_dt, clust_dt, name, label_order = NULL, label_colours = NULL,
  label_clusters = FALSE) {
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

  # Keep the canonical relative order among labels actually present, rather
  # than unioning in the full (possibly much larger) label_order -- otherwise
  # colours get spread across every canonical label, not just the ones shown.
  if (!is.null(label_order)) {
    present_order = intersect(label_order, levels(plot_dt$cluster))
    plot_dt[, cluster := factor(as.character(cluster), levels = unique(c(present_order, levels(cluster))))]
  }

  # tweak if "ambiguous"
  if ("ambiguous" %in% levels(plot_dt$cluster) )
    plot_dt     = plot_dt[, cluster := cluster %>% fct_relevel("ambiguous", after = Inf) ]

  # define colours
  cl_lvls     = levels(plot_dt$cluster)
  cl_cols     = resolve_label_colours(cl_lvls, label_colours)
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

  # Print each cluster's name at its median position, so that clusters sharing
  # a repeated colour can still be told apart
  if (label_clusters) {
    centres_dt = plot_dt[, .(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2)), by = cluster]
    g = g + geom_text(data = centres_dt, aes(label = cluster), colour = 'black', size = 2.5)
  }

  return(g)
}

calc_labels_table <- function(guesses_dt, label_order = NULL) {
  ns_dt   = guesses_dt %>%
    .[, .(cell_id, unagg = predicted_label_naive, agg = predicted_label_agg)] %>%
    melt(id = "cell_id", variable.name = "label_method", value.name = "label") %>%
    .[, .N, by = .(label_method, label)] %>%
    dcast( label ~ label_method, value.var = "N", fill = 0) %>%
    .[ order(-agg, -unagg) ] %>% setcolorder(c("label", "agg", "unagg"))
  if (!is.null(label_order)) {
    ns_dt[, label := factor(label, levels = unique(c(label_order, as.character(label))))]
    setorder(ns_dt, label)
  }

  # Keep ambiguous calls after the named labels.
  if ("ambiguous" %in% ns_dt$label)
    ns_dt     = rbind( ns_dt[ label != "ambiguous" ], ns_dt[ label == "ambiguous" ] )

  # Use descriptive headings in the report table.
  ns_dt     = ns_dt %>% set_colnames(c("predicted\nlabel", "no. cells,\naggregated", "no. cells, not\naggregated"))

  return(ns_dt)
}


get_label_views <- function(guesses_dt, style) {
  views = list(Labels = list(
    labels_dt = guesses_dt[, .(cell_id, predicted_label_naive, predicted_label_agg)],
    order = style$order, colours = style$colours))
  if (all(c('coarse_predicted_label_naive', 'coarse_predicted_label_agg') %in% names(guesses_dt))) {
    views[['Grouped labels']] = list(
      labels_dt = guesses_dt[, .(cell_id, predicted_label_naive = coarse_predicted_label_naive,
        predicted_label_agg = coarse_predicted_label_agg)],
      order = style$coarse_order, colours = style$coarse_colours)
  }
  return(views)
}

cat_label_views_heading <- function(title, views, level) {
  cat(strrep('#', level), ' ', title, if (length(views) > 1) ' {.tabset}', '\n\n', sep = '')
}

cat_label_view_tab <- function(view_name, views, level) {
  if (length(views) > 1) cat(strrep('#', level + 1), ' ', view_name, '\n\n', sep = '')
}
