
suppressPackageStartupMessages({
  library('magrittr')
  library('data.table')
  library('assertthat')
  library('scales')
  library('ggplot2')
  library('forcats')
})


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
