#libraries
suppressMessages({
  library('data.table')
  library('stringr')
  library('magrittr')
  library('assertthat')

  library('SingleCellExperiment')
  library('scater')
  library('Matrix')
  library('tidyverse')

  library('scDblFinder')
  library('BiocParallel')
  RhpcBLASctl::omp_set_num_threads(1L)
  library('DelayedArray')
  library('HDF5Array')
  library('rhdf5')

  library('ggplot2')
  library('ggplot.multistats')
  library('viridis')
  library('scales')
  library('patchwork')
})


main_doublet_id <- function(sel_sample, sce_f, sample_stats_f = NULL, sample_var = 'sample_id', ambient_method, dbl_f, dimred_f, min_feats = 100, min_cells = 100){
  message('running scDblFinder')
  
  # if ambient method is cellbender exclude bad samples
  smpl_status = FALSE

  if(ambient_method == 'cellbender'){
  # loading file with bad bender samples
  message(' loading cellbender sample stats file')
  sample_stats_df = fread(sample_stats_f)
  smpl_status = unique(sample_stats_df[get(sample_var) == sel_sample, bad_sample])

  if(smpl_status){
    message('  sample ', sel_sample, ' has been excluded. Saving empty results file')
    file.create(dbl_f)
    file.create(dimred_f)
    message('done!')

    return(NULL)
  }

  }else{

  # load sce file
  message('  loading sce object')
  sce         = sce_f %>% readRDS
  sample_idx  = sce[[sample_var]] == sel_sample
  sce         = sce[, sample_idx ]

  # exclude tiny cells
  message('  filtering out cells with low counts')
  assert_that( 'detected' %in% colnames(colData(sce)) )
  keep_idx    = sce$detected >= min_feats
  if (sum(keep_idx) < min_cells)

  assert_that( sum(keep_idx) >= min_cells,
    msg = "insufficient cells to run scDblFinder :(")
  message(sprintf('    keeping %d / %d cells (%.0f%%)',
    sum(keep_idx), length(keep_idx), 100 * sum(keep_idx) / length(keep_idx)))
  sce         = sce[, keep_idx]

  message('  running scDblFinder')
  dbl_dt      = scDblFinder(sce, returnType = 'table',
    multiSampleMode = 'singleModel', verbose = FALSE ) %>%
    as.data.table(keep.rownames = 'cell_id') %>%
    .[, (sample_var) := sel_sample ] %>% setcolorder(sample_var) %>%
    # keep just real cells
    .[type == 'real']

  # check if class is available from demultiplexing
  if('demux_class' %in% colnames(colData(sce))){
    #extract demux_class
    demux_dt = colData(sce) %>% as.data.table(keep.rownames = 'cell_id') %>%
      .[, c("cell_id", sample_var, "demux_class"), with = FALSE]

    # define dt to combine scdbl and demux class
    dbl_dict = data.table(
      class        = c('doublet', 'singlet', 'singlet', 'doublet', 'singlet', 'doublet'),
      demux_class  = c('negative', 'negative', 'doublet', 'singlet', 'singlet', 'doublet'),
      dbl_class    = c('doublet', 'negative', 'doublet', 'doublet', 'singlet', 'doublet')
    )

    dbl_dt = dbl_dt %>%
      merge(demux_dt, by = c('cell_id', sample_var), all.x = TRUE, all.y = FALSE) %>%
      merge(dbl_dict, by = c('class', 'demux_class')) %>%
      setnames("class", "scdbl_class")
  }else{
    dbl_dt = dbl_dt %>% setnames("class", "dbl_class")
  }

  setkeyv(dbl_dt, "cell_id")
  dbl_dt      = dbl_dt[ colnames(sce) ]

  message('  running PCA')
  dimred_dt   = .calc_one_dimred(sce, sel_sample)

  # check they match
  assert_that( all(sort(dbl_dt$cell_id) == sort(colnames(sce))) )
  assert_that( all(sort(dimred_dt$cell_id) == sort(dbl_dt$cell_id)) )

  # save
  message('  saving results')
  fwrite(dbl_dt, file = dbl_f)
  fwrite(dimred_dt, file = dimred_f)
  message('done!')

  return(NULL)
  }
}

.calc_one_dimred <- function(sce, sel_sample) {
  # run PCA on this sce
  sce       = sce %>% logNormCounts %>% runPCA
  pca_dt    = reducedDim(sce, "PCA") %>%
    as.data.table %>%
    .[,1:2] %>% set_colnames(c('pc1', 'pc2'))
  dimred_dt = data.table(
    cell_id   = colnames(sce),
    sample_id = sel_sample
    ) %>% cbind(pca_dt)

  return(dimred_dt)
}


combine_scDblFinder_outputs <- function(dbl_fs_f, sample_var, combn_dbl_f, combn_dimred_f, demux_type, n_cores) {
  setDTthreads(n_cores)

  # unpack some inputs
  dbl_fs_dt   = fread(dbl_fs_f)
  samples     = dbl_fs_dt[, get(sample_var)]
  dbl_fs      = dbl_fs_dt$dbl_f
  dimred_fs   = dbl_fs_dt$dimred_f

  # check some inputs
  assert_that(
    length(samples) == length(dbl_fs),
    length(samples) == length(dimred_fs)
  )
  assert_that(
    all(str_detect(dbl_fs, samples)),
    all(str_detect(dimred_fs, samples))
  )

  # get relevant files, exclude nulls
  dbl_ls      = lapply(dbl_fs, fread)
  dbl_ls      = dbl_ls[ sapply(dbl_ls, nrow) > 0 ]

  # get common columns, that we want
  first_cols  = c(sample_var, 'cell_id', 'dbl_class')
  exc_cols    = c('type', 'src')
  col_counts  = lapply(dbl_ls, colnames) %>% unlist %>% table
  keep_cols   = names(col_counts)[ col_counts == length(dbl_ls) ]
  keep_cols   = first_cols %>% c(setdiff(keep_cols, first_cols)) %>% setdiff(exc_cols)
  dbl_dt      = dbl_ls %>% lapply(function(d) d[, keep_cols, with = FALSE ] ) %>%
    rbindlist %>% setkey('cell_id')

  # save
  fwrite(dbl_dt, file = combn_dbl_f)

  # get dimred files, join and save
  dimred_dt   = lapply(dimred_fs, fread) %>% rbindlist
  assert_that( nrow(dimred_dt) == nrow(dbl_dt) )
  fwrite(dimred_dt, file = combn_dimred_f)
}

# summary plots
scdblfinder_diagnostic_plot <- function(s, sc_dbl_dt, dimred_dt) {
  # restrict to sample
  plot_dt   = merge(
    dimred_dt[ sample_id == s ],
    sc_dbl_dt,
    by = 'cell_id') %>%
    .[dbl_class != 'negative', is_doublet := dbl_class == 'doublet' ]

  # calc prop doublet
  prop_dbl  = mean(plot_dt$is_doublet)

  # calc cutoff
  if ( 'doublet' %in% plot_dt$dbl_class ) {
    cut_val   = mean(c(
      plot_dt[ dbl_class == 'singlet' ]$score %>% max,
      plot_dt[ dbl_class == 'doublet' ]$score %>% min
      ))
  } else {
    cut_val   = 1
  }

  # density of cells
  g_pca_dens = ggplot(plot_dt) +
    aes(x = pc1, y = pc2) +
    geom_hex(bins = 30) +
    scale_fill_distiller(palette = 'RdBu', trans = 'log10',
      limits = c(1,250)) +
    scale_x_continuous( breaks = pretty_breaks() ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      aspect.ratio = 1,
      legend.position = 'right'
      ) +
    labs(
      title   = paste0(s, ' (', round(100*prop_dbl), '% doublets, ',
        nrow(plot_dt),' cells)'),
      fill  = 'count'
      )

  # hex: doublet proportion
  g_pca_dbl = ggplot(plot_dt) +
    aes(x = pc1, y = pc2, z = is_doublet*1) +
    stat_summary_hex(fun = 'mean', bins = 30) +
    scale_fill_distiller(palette = 'PiYG', limits = c(0,1),
      breaks = pretty_breaks()) +
    labs(fill = 'mean doublets\nby bin') +
    scale_x_continuous( breaks = pretty_breaks() ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      title = element_text(''),
      axis.text = element_blank(),
      aspect.ratio = 1,
      legend.position = 'right'
      )

  # hex: doublet scores
  g_pca_score = ggplot(plot_dt) +
    aes(x = pc1, y = pc2, z = score) +
    stat_summary_hex(fun = 'mean', bins = 30) +
    scale_fill_viridis(breaks = pretty_breaks(), limits = c(0,1)) +
    labs(fill = 'mean score \nby bin') +
    scale_x_continuous( breaks = pretty_breaks() ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      title = element_text(''),
      axis.text = element_blank(),
      aspect.ratio = 1,
      legend.position = 'right'
      )

  # density of cells
  g_umap_dens = ggplot(plot_dt) +
    aes(x = umap1, y = umap2) +
    geom_hex(bins = 30) +
    scale_fill_distiller(palette = 'RdBu', trans = 'log10',
      limits = c(1,250)) +
    labs(fill = 'count') +
    scale_x_continuous( breaks = pretty_breaks() ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      aspect.ratio = 1,
      legend.position = 'right')

  # hex: doublet proportion
  g_umap_dbl = ggplot(plot_dt) +
    aes(x = umap1, y = umap2, z = is_doublet*1) +
    stat_summary_hex(fun = 'mean', bins = 30) +
    scale_fill_distiller(palette = 'PiYG', limits = c(0,1),
      breaks = pretty_breaks()) +
    labs(fill = 'mean doublets\nby bin') +
    scale_x_continuous( breaks = pretty_breaks() ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      aspect.ratio = 1,
      legend.position = 'right'
      )

  # hex: doublet scores
  g_umap_score = ggplot(plot_dt) +
    aes(x = umap1, y = umap2, z = score) +
    stat_summary_hex(fun = 'mean', bins = 30) +
    scale_fill_viridis(breaks = pretty_breaks(), limits = c(0,1)) +
    labs(fill = 'mean score \nby bin') +
    scale_x_continuous( breaks = pretty_breaks() ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      aspect.ratio = 1,
      legend.position = 'right'
      )

  # hist: doublet scores
  g_hist = ggplot(plot_dt) +
    aes( x = score ) +
    geom_histogram( boundary = 0, binwidth = 0.02 ) +
    geom_vline( xintercept = cut_val ) +
    scale_x_continuous( breaks = pretty_breaks(), limits = c(0,1) ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      aspect.ratio = 1,
      legend.position = 'right'
    )

  # aggregate plots
  g = list(g_pca_dens, g_pca_dbl, g_pca_score, g_hist,
    g_umap_dens, g_umap_dbl, g_umap_score) %>%
    wrap_plots(ncol = 4)
}
