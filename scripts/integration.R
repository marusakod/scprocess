# integration.R


suppressPackageStartupMessages({
  library('magrittr')
  library('data.table')
  library('forcats')
  library('stringr')
  library('assertthat')

  library('BiocParallel')
  library('Seurat')
  # library('glmGamPoi')
  library('harmony')
  library('scran')
  # library('reticulate')
  # use_condaenv("scanpy")
  # sc        = reticulate::import("scanorama")
  library('uwot')
  library('future')
  library('ggplot.multistats')

  # library('ComplexHeatmap')
  # library('seriation')
  # library('ggh4x')
})

run_harmony <- function(sce_all_f, keep_f, dbl_f,
  exc_regex, n_hvgs, n_dims, dbl_res, dbl_cl_prop, theta, res_ls_concat,
  harmony_f, hvgs_f, sce_clean_f, n_cores = 4) {
  # unpack inputs
  res_ls      = res_ls_concat %>% str_split(" ") %>% unlist %>% as.numeric

  message('running Harmony')
  # load sce, restrict to ok cells
  message('  setting up cluster')
  plan("multicore", workers = n_cores)
  options( future.globals.maxSize = 2^35 )

  message('  loading relevant cell ids')
  keep_ids    = fread(keep_f) %>% .$cell_id
  dbl_ids     = fread(dbl_f) %>% .[ class == "doublet" ] %>% .$cell_id
  load_ids    = c(keep_ids, dbl_ids)

  message('  loading sce')
  assert_that( length(keep_ids) > 0 )
  sce_dbl     = readRDS(sce_all_f) %>% .[, load_ids ]
  assert_that("sample_id" %in% colnames(colData(sce_dbl)))

  # exclude genes if requested
  if (!is.null(exc_regex)) {
    exc_idx     = rownames(sce_dbl) %>% str_detect(exc_regex)
    exc_gs      = rowData(sce_dbl)$symbol[ exc_idx ]
    sprintf("    excluding %d genes: %s", sum(exc_idx), paste0(exc_gs, collapse = " ")) %>%
      message
    sce_dbl     = sce_dbl[ !exc_idx, ]
  }

  # turn into seurat object
  message('  prepping seurat object')
  suppressWarnings({
    seu_dbl     = Seurat::CreateSeuratObject(
      counts      = counts(sce_dbl),
      meta.data   = data.frame(colData(sce_dbl)),
      project     = "dummy"
      )
  })

  # run harmony including doublets
  message('  running Harmony to find more doublets')
  hmny_dbl    = .run_one_harmony(seu_dbl, n_hvgs, n_dims, theta = 0, dbl_res)
  dbl_data    = .calc_dbl_data(hmny_dbl, dbl_ids, dbl_res, dbl_cl_prop)

  # make new seurat object without these doublets
  message('  running on clean data')
  ok_ids      = setdiff(keep_ids, dbl_data[ in_dbl_cl == TRUE ]$cell_id)
  # seu_ok      = seu_dbl[ , ok_ids ]
  suppressWarnings({
    seu_ok      = Seurat::CreateSeuratObject(
      counts      = seu_dbl[['RNA']]$counts[ , ok_ids ],
      meta.data   = data.frame(colData(sce_dbl)[ ok_ids, ]),
      project     = "dummy"
      )
  })
  rm(seu_dbl); gc()
  hmny_ok     = .run_one_harmony(seu_ok, n_hvgs, n_dims, theta = theta, res_ls)

  # join these together
  hmny_dt     = merge(hmny_ok, dbl_data, by = c("sample_id", "cell_id"), all = TRUE)

  # save outputs
  fwrite(hmny_dt, file = harmony_f)

  # do HVGs
  message('  saving HVGs')
  seu_ok      = seu_ok %>%
    NormalizeData( verbose = FALSE ) %>%
    FindVariableFeatures( nfeatures = n_hvgs, verbose = FALSE ) %>%
    ScaleData( verbose = FALSE )
  hvgs_dt   = seu_ok@assays$RNA@meta.data %>%
    as.data.table(keep.rownames = "gene_id") %>%
    .[, is_hvg := vf_vst_counts_variable ]
  fwrite(hvgs_dt, file = hvgs_f)

  # make sce file with only ok cells, and nice clusters
  message('  making nice clean sce')
  sce       = .annotate_sce_w_harmony(sce_dbl, hmny_ok)
  message('    saving')
  saveRDS(sce, file = sce_clean_f, compress = FALSE)
  message('done!')
}

.run_one_harmony <- function(seu_obj, n_hvgs, n_dims, theta = 0, res_ls) {
  message('    normalizing and finding HVGs')
  seu_obj   = seu_obj %>%
    NormalizeData( verbose = FALSE ) %>%
    FindVariableFeatures( nfeatures = n_hvgs, verbose = FALSE ) %>%
    ScaleData( verbose = FALSE )

  # check whether we have one or more values of batch
  n_samples       = seu_obj$sample_id %>% unique %>% length
  is_one_batch    = n_samples == 1
  this_reduction  = ifelse(is_one_batch, "pca", "harmony")
  if (is_one_batch) {
    stop("only one sample; doesn't make sense to do this")
    # run Seurat pipeline, plus clustering
    warning('    only one sample; not doing Harmony')
    message('    running UMAP')
    seu_obj   = seu_obj %>%
      RunPCA( npcs = n_dims, verbose = FALSE ) %>%
      RunUMAP( reduction = this_reduction, dims = 1:n_dims, verbose = FALSE  )
  } else {
    # run Seurat pipeline, plus clustering
    message('    running Harmony and UMAP')
    seu_obj   = seu_obj %>%
      RunPCA( npcs = n_dims, verbose = FALSE ) %>%
      RunHarmony( c("sample_id"), theta = theta, verbose = TRUE ) %>%
      RunUMAP( reduction = this_reduction, dims = 1:n_dims, verbose = FALSE  )
  }

  message('    finding clusters')
  seu_obj   = seu_obj %>%
    FindNeighbors( reduction = this_reduction, dims = 1:n_dims, verbose = FALSE  ) %>%
    FindClusters( resolution = res_ls, verbose = FALSE ) %>%
    identity()

  # switch cluster off
  plan("sequential")

  message('    recording clusters')
  clusts_dt   = seu_obj[[]] %>% as.data.table %>% .[, integration := "harmony" ]
  cl_vs       = names(clusts_dt) %>% str_subset('RNA_snn_res\\.')
  clusts_dt   = clusts_dt[, c('integration', 'cell_id', 'sample_id', cl_vs), with = FALSE ]
  for (cl_v in cl_vs) {
    clusts_dt[[ cl_v ]] = clusts_dt[[ cl_v ]] %>% fct_infreq %>% as.integer %>%
      sprintf("cl%02d", .)
  }

  # get embeddings
  message('    extracting other outputs')
  hmny_pca_dt = Embeddings(seu_obj, reduction = this_reduction) %>%
    as.data.table( keep.rownames = "cell_id" ) %>%
    setnames(names(.), names(.) %>% str_replace("_(?=[0-9]$)", "_0") %>%
      str_replace("harmony", "hmny_pca"))
  umap_dt     = Embeddings(seu_obj, reduction = "umap") %>%
    as.data.table( keep.rownames = "cell_id" ) %>%
    setnames(names(.), names(.) %>% str_replace("umap_", "UMAP"))

  # join together
  hmny_dt     = clusts_dt %>% merge(umap_dt, by = "cell_id") %>%
    merge(hmny_pca_dt, by = "cell_id")

  return(hmny_dt)
}

.calc_dbl_data <- function(hmny_dbl, dbl_ids, dbl_res, dbl_cl_prop) {
  # assemble useful doublet data
  dbl_data = hmny_dbl %>%
    .[, .(cell_id, sample_id, dbl_UMAP1 = UMAP1, dbl_UMAP2 = UMAP2,
      dbl_cluster = get(paste0("RNA_snn_res.", dbl_res)))] %>%
    .[, is_dbl := cell_id %in% dbl_ids]

  # calculate proportion doublets per dbl_cluster
  dbl_data = dbl_data %>%
    .[, dbl_prop  := sum(is_dbl) / .N, by = dbl_cluster ] %>%
    .[, in_dbl_cl := dbl_prop > dbl_cl_prop ]

  return(dbl_data)
}

.annotate_sce_w_harmony <- function(sce_dbl, hmny_ok) {
  # restrict to just ok cells
  sce       = sce_dbl[ , hmny_ok$cell_id ]

  # get useful harmony variables
  hmny_vs   = c('UMAP1', 'UMAP2', str_subset(names(hmny_ok), "RNA_snn_res"))

  # add these to sce object
  for (v in hmny_vs) {
    if (str_detect(v, "RNA_snn_res")) {
      colData(sce)[[ v ]] = hmny_ok[[ v ]] %>% factor
    } else {
      colData(sce)[[ v ]] = hmny_ok[[ v ]]
    }
  }

  return(sce)
}

plot_umap_density <- function(input_dt) {
  # eps         = 0.001
  umap_dt     = copy(input_dt) %>%
    .[, .(
      UMAP1 = rescale(UMAP1, to = c(0.05, 0.95)),
      UMAP2 = rescale(UMAP2, to = c(0.05, 0.95))
      )]
    # .[, UMAP1_trunc := UMAP1 %>% pmin(quantile(UMAP1, 1 - eps)) %>%
    #   pmax(quantile(UMAP1, eps)) %>% rescale(c(0.05, 0.95)) ] %>%
    # .[, UMAP2_trunc := UMAP2 %>% pmin(quantile(UMAP2, 1 - eps)) %>%
    #   pmax(quantile(UMAP2, eps)) %>% rescale(c(0.05, 0.95)) ]

  g = ggplot(umap_dt) + aes( x = UMAP1, y = UMAP2 ) +
    geom_bin2d(bins = 50) +
    scale_fill_distiller( palette = "RdBu", trans = "log10" ) +
    scale_x_continuous( breaks = pretty_breaks(), limits = c(0, 1) ) +
    scale_y_continuous( breaks = pretty_breaks(), limits = c(0, 1) ) +
    theme_bw() +
    theme( panel.grid = element_blank(), axis.text = element_blank(), aspect.ratio = 1 )

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
    theme( panel.grid = element_blank(), axis.text = element_blank(), aspect.ratio = 1 )

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

  # define colours
  if (name == 'type_broad') {
    cl_cols     = broad_cols
  } else {
    cl_cols     = seq_along( levels(plot_dt$cluster) ) %>%
      rep(nice_cols, times = 10)[ . ] %>%
      setNames( levels(umap_dt$cluster) )
  }

  # make plot
  g = ggplot(plot_dt) +
    aes( x = UMAP1, y = UMAP2, colour = cluster ) +
    geom_point(size = 0.1) +
    scale_colour_manual( values = cl_cols, guide = guide_legend(override.aes = list(size = 3)) ) +
    scale_x_continuous( breaks = pretty_breaks(), limits = c(0, 1) ) +
    scale_y_continuous( breaks = pretty_breaks(), limits = c(0, 1) ) +
    theme_bw() +
    theme( panel.grid = element_blank(), axis.text = element_blank(), aspect.ratio = 1 ) +
    labs( colour = name )

  return(g)
}

plot_cluster_entropies <- function(input_dt, what = c("norm", "raw")) {
  # check inputs
  what      = match.arg(what)

  # calculate mixing
  ns_dt     = input_dt %>%
    .[, .N, by = .(sample_id, cluster) ] %>%
    .[, p_sample  := N / sum(N), by = sample_id ] %>%
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
    expand_limits( y = 0 ) +
    scale_size(
      range   = c(1, 8),
      breaks  = c(2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4) %>% sqrt,
      labels  = c('200', '500', '1k', '2k', '5k', '10k', '20k', '50k')
      ) +
    theme_bw() +
    theme( panel.grid = element_blank() ) +
    labs(
      x     = 'entropy (high when clusters even across samples)',
      y     = 'max. pct. of one sample (high when concentrated in few samples)',
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
  if (name == 'type_broad') {
    cl_cols     = broad_cols
  } else {
    cl_cols     = seq_along( levels(plot_dt$cluster) ) %>%
      rep(nice_cols, times = 10)[ . ] %>%
      setNames( levels(plot_dt$cluster) )
  }

  # define breaks
  log_brks    = c(1e1, 2e1, 5e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4) %>%
    log10
  log_labs    = c("10", "20", "50", "100", "200", "500",
    "1k", "2k", "5k", "10k", "20k", "50k")
  logit_brks  = c(1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.10, 0.30,
    0.50, 0.70, 0.90, 0.97, 0.99) %>% qlogis
  logit_labs  = c("0.01%", "0.03%", "0.1%", "0.3%", "1%", "3%", "10%", "30%",
    "50%", "70%", "90%", "97%", "99%")
  splice_brks = seq(0.1, 0.9, 0.1) %>% qlogis
  splice_labs = (100*seq(0.1, 0.9, 0.1)) %>% as.integer %>% sprintf("%d%%", .)

  # plot
  g = ggplot(plot_dt) + aes( x = cluster, y = qc_val, fill = cluster ) +
    geom_violin() +
    scale_fill_manual( values = cl_cols, guide = "none" ) +
    facet_grid( qc_full ~ ., scales = 'free_y' ) +
    facetted_pos_scales(
      y = list(
        qc_full == "library size"    ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full == "no. of features" ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full == "mito pct."        ~
          scale_y_continuous(breaks = logit_brks, labels = logit_labs),
        qc_full == "spliced pct."     ~
          scale_y_continuous(breaks = splice_brks, labels = splice_labs)
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
