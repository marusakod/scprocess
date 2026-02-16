# integration.R

suppressPackageStartupMessages({
  library('magrittr')
  library('data.table')
  library('forcats')
  library('stringr')
  library('assertthat')
  library('BiocParallel')
  library('Seurat')
  library('harmony')
  library('scran')
  library('uwot')
  library('future')
  library('ggplot.multistats')
  library('zellkonverter')
})



plot_umap_density <- function(input_dt) {
  # eps         = 0.001
  umap_dt     = copy(input_dt) %>%
    .[, .(
      UMAP1 = rescale(UMAP1, to = c(0.05, 0.95)),
      UMAP2 = rescale(UMAP2, to = c(0.05, 0.95))
      )]

  g = ggplot(umap_dt) + aes( x = UMAP1, y = UMAP2 ) +
    geom_bin2d(bins = 50) +
    scale_fill_distiller( palette = "RdBu", trans = "log10" ) +
    scale_x_continuous( breaks = pretty_breaks(), limits = c(0, 1) ) +
    scale_y_continuous( breaks = pretty_breaks(), limits = c(0, 1) ) +
    theme_bw() +
    theme( panel.grid = element_blank(),
      legend.title.position = "left", legend.position = "bottom",
      axis.ticks = element_blank(), axis.text = element_blank(),
      aspect.ratio = 1 )

  return(g)
}

plot_umap_doublets <- function(input_dt) {
  dbl_dt      = copy(input_dt) %>%
    .[, .(
      UMAP1 = rescale(UMAP1, to = c(0.05, 0.95)),
      UMAP2 = rescale(UMAP2, to = c(0.05, 0.95)),
      is_dbl
      )]

  g = ggplot(dbl_dt) +
    aes(x = UMAP1, y = UMAP2, z = is_dbl*1) +
    stat_summary_hex(fun = 'mean', bins = 30) +
    scale_fill_distiller(palette = 'PiYG', limits = c(0, 1),
      breaks = pretty_breaks()) +
    labs(fill = 'mean doublets\nby bin') +
    scale_x_continuous( breaks = pretty_breaks() ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    theme_bw() +
    theme( panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), aspect.ratio = 1 )

  return(g)
}

plot_doublet_clusters <- function(hmny_dbl, dbl_cl_prop) {
  # calculate doublet cluster sizes and proportion doublets
  dbl_clusts  = hmny_dbl %>%
    .[, .(
      dbl_cl_size    = .N,
      dbl_cl_prop    = sum(is_dbl) / .N * 100
      ), by = dbl_cluster ]

  log_brks    = c(1e1, 2e1, 5e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4) %>%
    log10
  log_labs    = c("10", "20", "50", "100", "200", "500",
    "1k", "2k", "5k", "10k", "20k", "50k")

  # plot doublet cluster log size and doublet proportions against each other
  g = ggplot(dbl_clusts) +
    aes( x = log10(dbl_cl_size), y = dbl_cl_prop ) +
    geom_hline( yintercept = dbl_cl_prop * 100, linetype = "dashed", colour = 'black' ) +
    geom_point() +
    scale_x_continuous( breaks = log_brks, labels = log_labs ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    expand_limits( y = c(0, 100) ) +
    theme_classic() +
    labs( x = "# cells in cluster", y = "doublet pct. of cluster")

  return(g)
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
      )] %>% .[, cluster := factor(cluster) ]
  # plot_dt     = rbind(plot_dt[ is.na(cluster) ], plot_dt[ !is.na(cluster) ])
  plot_dt     = plot_dt[ sample(.N, .N) ]

  # add labels to clusters
  cl_labels   = plot_dt[, .( N = .N ), by = cluster ] %>%
    .[ order(cluster) ] %>%
    .[, cl_label  := sprintf("%s (%d)", cluster, signif(N, 2)) ]
  label_lu    = cl_labels$cl_label %>% setNames(cl_labels$cluster)

  # define colours
  cl_cols     = seq_along( label_lu ) %>% rep(nice_cols, times = 10)[ . ] %>%
    setNames( label_lu )
  n_col       = 4
  n_rows_lgd  = ceiling(length(label_lu) / n_col)

  # make plot
  g = ggplot(plot_dt) +
    aes( x = UMAP1, y = UMAP2, colour = label_lu[cluster] ) +
    geom_point(size = 0.1) +
    scale_colour_manual( values = cl_cols, 
      guide = guide_legend(override.aes = list(size = 3), nrow = n_rows_lgd) ) +
    scale_x_continuous( breaks = pretty_breaks(), limits = c(0, 1) ) +
    scale_y_continuous( breaks = pretty_breaks(), limits = c(0, 1) ) +
    theme_bw() +
    theme( panel.grid = element_blank(), aspect.ratio = 1,
      legend.title.position = "left", legend.position = "bottom",
      axis.ticks = element_blank(), axis.text = element_blank() ) +
    labs( colour = name )

  return(g)
}

plot_cluster_entropies <- function(input_dt, batch_var, what = c("norm", "raw")) {
  # check inputs
  what      = match.arg(what)

  # calculate mixing
  ns_dt     = input_dt %>% copy %>% 
    .[, .N, by = .(batch_var, cluster) ] %>%
    .[, p_sample  := N / sum(N), by = batch_var ] %>%
    .[, p_cluster := N / sum(N), by = cluster ] %>%
    .[, p_cl_norm := p_sample / sum(p_sample), by = cluster ]

  # calculate entropy
  entropy_dt  = ns_dt[, .(
      h_cl_raw      = -sum(p_cluster * log2(p_cluster), na.rm = TRUE),
      h_cl_norm     = -sum(p_cl_norm * log2(p_cl_norm), na.rm = TRUE),
      max_pct_raw   = 100 * max(p_cluster),
      max_pct_norm  = 100 * max(p_cl_norm),
      N             = sum(N)
    ), by = cluster ]
  labels_dt   = entropy_dt[ order(cluster) ]

  # get nice colours
  cl_ls     = entropy_dt$cluster %>% unique %>% sort
  cl_cols   = nice_cols[ seq_along(cl_ls) ] %>% setNames(cl_ls)

  # plot
  g = ggplot(entropy_dt) +
    aes_string( x = paste0('h_cl_', what), y = paste0('max_pct_', what) ) +
    geom_smooth( method = "lm", formula = y ~ x, se = FALSE, colour = "grey" ) +
    geom_text_repel( data = labels_dt, aes(label = cluster), size = 3,
      min.segment.length = 0, max.overlaps = Inf, box.padding = 0.5 ) +
    geom_point( shape = 21, aes( size = sqrt(N), fill = cluster ) ) +
    scale_x_continuous(breaks = pretty_breaks(n = 3)) +
    scale_y_continuous(breaks = pretty_breaks(n = 3)) +
    scale_fill_manual( values = cl_cols, guide = "none" ) +
    expand_limits( x = 0, y = 0 ) +
    scale_size(
      range   = c(1, 8),
      breaks  = c(2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4) %>% sqrt,
      labels  = c('200', '500', '1k', '2k', '5k', '10k', '20k', '50k')
      ) +
    theme_bw() +
    theme( panel.grid = element_blank() ) +
    labs(
      x     = 'entropy',
      y     = sprintf('max. pct. of one %s', batch_var),
      size  = 'total # cells'
    )

  return(g)
}

plot_cluster_qc_distns <- function(qc_melt, clust_dt, name, min_cl_size = 1e2) {
  # edit cluster names
  if (is.numeric(name)) {
    # exclude tiny clusters
    cl_n_dt     = clust_dt[, .N, by = cluster]
    cls_tiny    = cl_n_dt$N < min_cl_size
    if (any(cls_tiny)) {
      clust_dt    = clust_dt[ cluster %in% cl_n_dt[ !cls_tiny ]$cluster ]
    }

    # change clusters to characters
    clust_dt    = clust_dt %>% .[, cluster := fct_infreq(cluster) ]

    # update name
    name      = sprintf('res = %s', name)
  }

  # join
  plot_dt   = merge(qc_melt, clust_dt, by = 'cell_id')

  # define colours
  cl_cols     = seq_along( levels(plot_dt$cluster) ) %>%
    rep(nice_cols, times = 10)[ . ] %>% setNames( levels(plot_dt$cluster) )

  # define breaks
  log_brks    = c(1e1, 2e1, 5e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5) %>%
    log10
  log_labs    = c("10", "20", "50", "100", "200", "500",
    "1k", "2k", "5k", "10k", "20k", "50k", "100k", "200k", "500k")
  logit_brks  = c(1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.10, 0.30,
    0.50, 0.70, 0.90, 0.97, 0.99) %>% qlogis
  logit_labs  = c("0.01%", "0.03%", "0.1%", "0.3%", "1%", "3%", "10%", "30%",
    "50%", "70%", "90%", "97%", "99%")
  splice_brks = seq(0.1, 0.9, 0.1) %>% qlogis
  splice_labs = (100*seq(0.1, 0.9, 0.1)) %>% as.integer %>% sprintf("%d%%", .)

  # plot
  g = ggplot(plot_dt) + aes( x = cluster, y = qc_val, fill = cluster ) +
    geom_violin(kernel = 'rectangular', adjust = 0.5, scale = 'width', width = 0.8, colour = NA) +
    scale_fill_manual( values = cl_cols, guide = "none" ) +
    facet_grid( qc_full ~ ., scales = 'free_y' ) +
    facetted_pos_scales(
      y = list(
        qc_full == "no. of UMIs"    ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full == "no. of genes" ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full == "mito pct."        ~
          scale_y_continuous(breaks = logit_brks, labels = logit_labs),
        qc_full == "spliced pct."     ~
          scale_y_continuous(breaks = logit_brks, labels = logit_labs)
        )
      ) +
    theme_bw() +
    theme(
      axis.text.x       = element_text( angle = 90, hjust = 1, vjust = 0.5 ),
      panel.grid        = element_blank(),
      strip.background  = element_rect( fill = 'white')
      ) +
    labs( x = name, y = "QC metric" )

  return(g)
}

plot_marker_dotplot <- function(exp_dt, clust_dt, name, markers_dt, min_cl_size = 1e2) {
  # new or old?
  is_new    = is.numeric(name)

  # fiddle with clusters
  if (is_new) {
    # exclude tiny clusters
    cl_n_dt     = clust_dt[, .N, by = cluster]
    cls_tiny    = cl_n_dt$N < min_cl_size
    if (any(cls_tiny)) {
      clust_dt    = clust_dt[ cluster %in% cl_n_dt[ !cls_tiny ]$cluster ]
    }

    # change clusters to characters
    n_digits    = clust_dt$cluster %>% as.integer %>% max(na.rm = TRUE) %>%
      log10 %>% ceiling
    clust_dt    = clust_dt %>%
      .[, cluster := sprintf("cl%%0%dd", n_digits) %>%
      sprintf(as.integer(cluster)) %>% factor ]

    # update name
    name      = sprintf('res = %s', name)
  }

  # calc expression per cluster
  p_count   = 1
  dots_dt   = merge(clust_dt, exp_dt, by = 'cell_id') %>%
    .[, .(
      count         = sum(count),
      prop_detected = mean(count > 0),
      lib_size      = sum(lib_size)),
      by = c('cluster', 'gene_id')] %>%
    .[, CPM       := ifelse(lib_size == 0, 0, count / lib_size * 1e6) ] %>%
    .[, log10cpm  := log10( CPM + p_count ) ] %>%
    .[, symbol    := gene_id %>% str_extract('^[^_]+') ] %>%
    merge(markers_dt, by = 'symbol') %>%
    .[, symbol    := symbol %>% factor(levels = markers_dt$symbol) ]

  if (is_new) {
    # put matrix in nice order
    order_dt    = dots_dt %>%
      dcast(symbol ~ cluster, value.var = 'log10cpm', fill = 0) %>%
      melt( id = 'symbol', var = 'cluster' ) %>%
      .[, symbol := factor(symbol, levels = markers_dt$symbol) ] %>%
      .[ order(cluster, symbol) ] %>%
      .[, smth_score := ksmooth(as.numeric(symbol), value,
        kernel = 'normal', x.points = as.numeric(symbol))$y, by = cluster ] %>%
      .[, .SD[ which.max(smth_score) ], by = cluster ] %>%
      .[ order(symbol) ]
    assert_that( all( sort(order_dt$cluster) == levels(clust_dt$cluster) ) )

    # add this order to dots
    dots_dt     = dots_dt[, cluster := factor(cluster, levels = order_dt$cluster) ]
  }

  # dotplot
  brk_vals  = log10(p_count + c(0, 10^seq(0, 4)))
  brk_labs  = c('0', '1', '10', '100', '1k', '10k')
  plot_dt   = dots_dt[ prop_detected > 0.01 ] %>% setorder( CPM )
  g = ggplot(plot_dt) +
    aes(x = cluster, y = fct_rev(symbol), fill = log10cpm, size = prop_detected) +
    geom_point(shape = 21) +
    scale_fill_viridis(breaks = brk_vals, labels = brk_labs, option = 'magma') +
    scale_size(range = c(0, 4), breaks = pretty_breaks()) +
    expand_limits(size = 0, fill = log10(p_count)) +
    facet_grid( marker_cat ~ ., scales = 'free_y', space = 'free_y' ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5, hjust = 1),
      panel.grid  = element_blank(),
      panel.spacing     = unit(0.2, "lines"),
      strip.text.y      = element_text(size = 8),
      strip.background  = element_rect( fill = 'white')
    ) +
    labs(
      x = name, y = 'marker gene',
      fill = 'CPM', size = 'propn. cells\nexpressing'
      )

  return(g)
}

make_clean_sce_from_h5ad <- function(sel_batch, adata_f, sce_f){
  sce = readH5AD(adata_f)
  # rename X assay to counts
  assayNames(sce)[assayNames(sce) == "X"] = "counts"
  
  message('saving sce object for ', sel_batch)
  saveRDS(sce, file=sce_f, compress = FALSE)
}