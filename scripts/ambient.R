# QC of cellbender input parameters

suppressPackageStartupMessages({
  library('data.table')
  library('magrittr')
  library('Matrix')
  library('assertthat')
  library('DropletUtils')

  library('gridExtra')
  library('tidyverse')
  library('Rfast')
  #library('Seurat')
  library('ggbeeswarm')
  library('BiocParallel')
  library('robustbase')

  library("decontX")
  library("strex")
  library("testit")
  library('yaml')
})

# get matrix with (decontx cleaned) cells and a list of barcodes called as cells
get_cell_mat_and_barcodes <- function(out_mat_f, out_bcs_f, out_dcx_f = NULL, sel_s, af_mat_f, knee_f, 
  cell_calls_method = c('barcodeRanks', 'emptyDrops'), ncores = 4, niters = 1000, hvg_n = 2000,
  ambient_method = c('none', 'decontx')) {
  
  # check inputs
  call_m  = match.arg(cell_calls_method)
  amb_m   = match.arg(ambient_method)

  # load alevin matrix
  af_mat  = .get_alevin_mx(af_mat_f = af_mat_f, sel_s = '')
  
  # load knee
  knee_dt = fread(knee_f)
 
  message('Looking for cells and empty droplets in sample ', sel_s)
  # call cells and empties
  cell_empty_bcs = call_cells_and_empties(
    af_mat  = af_mat,
    knee_dt = knee_dt,
    ncores  = ncores,
    n_iters = niters,
    call_method = call_m
    )
  
  # get matrix with cells
  cell_mat = af_mat[, cell_empty_bcs[["cells"]]]

  # do ambient correction
  if (ambient_method == 'none') {

    message("Ambient RNA removal method is 'none'. Saving uncorrected count matrix for ", sel_s )
    # save list of retained barcodes
    cell_bcs_dt = data.table(barcode = colnames(cell_mat))
    fwrite(cell_bcs_dt, out_bcs_f, col.names = FALSE, quote = FALSE)
    # save matrix
    write10xCounts(out_mat_f, cell_mat, version = "3", overwrite = TRUE)

  } else {
    message("Running 'decontx' for ", sel_s)
    # get empty matrix
    empty_mat = af_mat[, cell_empty_bcs[["empty"]]]
    dcx_res   = decontX(x = cell_mat, background = empty_mat, varGenes = hvg_n)

    clean_mat = dcx_res$decontXcounts %>% round()
    # remove barcodes with no remaining counts
    keep_idx = which(Matrix::colSums(clean_mat) != 0)
    message("Keeping ", length(keep_idx), " barcodes with non-zero counts in decontaminated matrix")
    clean_mat = clean_mat[, keep_idx]
    
    message("Saving decontaminated count matrix for ", sel_s)
    # save list of retained barcodes
    cell_bcs_dt = data.table(barcode = colnames(clean_mat))
    fwrite(cell_bcs_dt, out_bcs_f, col.names = FALSE, quote = FALSE)
    # save matrix
    write10xCounts(out_mat_f, clean_mat, version = "3", overwrite = TRUE)

    # save contamination proportions and decontX cluster assignment
    dcx_clust  = dcx_res$z %>% sprintf("cl%02d", .) %>% factor
    dcx_contam = dcx_res$contamination
    
    dcx_params = data.table(
      barcode = colnames(clean_mat),
      pct.contamination = dcx_contam[keep_idx],
      dcx_cluster = dcx_clust[keep_idx]
      )

    fwrite(dcx_params, out_dcx_f)

    message("'decontx' completed for ", sel_s)

  }

  return(NULL)
}

call_cells_and_empties <- function(af_mat, knee_dt, ncores = 4, n_iters = 1000, fdr_thr = 0.001,
  call_method = c('barcodeRanks', 'emptyDrops')) {
  
  # get empty droplets and their indices
  empty_bcs   = knee_dt %>%
    .[in_empty_plateau == TRUE, barcode]  
  empty_idx = which(colnames(af_mat) %in% empty_bcs)

  if (call_method == 'barcodeRanks') {
    exp_cells = knee_dt$expected_cells %>% unique
    cell_bcs =  knee_dt$barcode[1: exp_cells]
  } else {
    # get retain (threshold for the total umi count above which all barcodes are assumed to contain cells)
    knee_1 = knee_dt$knee1 %>% unique
    # get sum of s+u+a counts (instead of this maybe try removing genes with v low expression)
    # doing this because otherwise takes too long
    af_mat_sum =.sum_SUA(af_mat)
  
    bpparam = MulticoreParam(workers = ncores, progressbar = TRUE)
    
    # call cells
    edrops_res = emptyDrops(
      m           = af_mat_sum,
      niters      = n_iters,
      BPPARAM     = bpparam,
      known.empty = empty_idx,
      retain      = knee_1  
      )

    # get cell barcodes
    cell_bcs = edrops_res %>% as.data.table(keep.rownames = 'barcode') %>%
      .[FDR <= fdr_thr, barcode]
  }

  # return a list with cell and empty barcodes 
  empty_cell_bcs_ls = list(
    empty = empty_bcs,
    cells = cell_bcs
    )
  
  return(empty_cell_bcs_ls)
}


save_barcode_qc_metrics <- function(af_h5_f, amb_out_yaml, out_qc_f, ambient_method) {
  # read in the yaml file
  amb_yaml = yaml::read_yaml(amb_out_yaml)

  # get counts matrix
  if (ambient_method == 'cellbender') {
    amb_mat_f   = amb_yaml$raw_counts_f
  } else {
    amb_mat_f   = amb_yaml$filt_counts_f
  }

  # get alevin counts
  af_mat      = .get_alevin_mx(af_h5_f, '')
  af_dt       = .get_usa_dt(af_mat, prefix = 'pre')

  # get pre and post ambient removal stats
  if (ambient_method == 'cellbender') {
    # get ranks
    ranks       = barcodeRanks(af_mat) %>% as.data.frame %>%  
      as.data.table( keep.rownames = "barcode") %>% .[ order(rank) ]

    # get cellbender matrix
    cb_mat      = .get_alevin_mx(amb_mat_f, '')

    # sum s/u/a counts for each barcode
    cb_dt       = .get_usa_dt(cb_mat, prefix = 'post')

    # merge together
    qc_dt       = merge(af_dt, cb_dt, by = 'barcode')

  } else if (ambient_method == 'decontx') {

    # get decontx matrix
    dcx_mat     = .get_alevin_mx(amb_mat_f, '')

    # sum s/u/a counts for each barcode
    dcx_dt      = .get_usa_dt(dcx_mat, prefix = 'post')

    # merge together
    qc_dt       = merge(af_dt, dcx_dt, by = 'barcode', all.x = TRUE)

  } else if (ambient_method == "none") {

    # just use full matrix
    qc_dt       = af_dt

  }

  fwrite(qc_dt, out_qc_f)

  return(NULL)
}

.get_usa_dt <- function(mat, prefix) {
  # get counts
  usa_ls      = c("S", "U", "A")
  sums_ls     = usa_ls %>% sapply( function(x)
    colSums(mat[ str_detect(rownames(mat), paste0("_", x, "$")), ]) )
  usa_dt      = sums_ls %>% set_colnames(paste0(prefix, "_", usa_ls)) %>%
    as.data.table(keep.rownames = "barcode")
  return(usa_dt)
}

get_bender_log <- function(f, sample, sample_var) {
  ll =  read_lines(f, n_max = 25)
  .get_match <- function(ll, pat) {
    ll %>% str_subset(pat) %>% str_match(pat) %>% .[, 3] %>% as.integer()
  }
  bender_dt = data.table(
    cb_prior_cells              = .get_match(ll, '(Prior on counts for cells is )([0-9]+)'),
    cb_prior_empty              = .get_match(ll, '(Prior on counts [infor]{2,3} empty droplets is )([0-9]+)'),
    excluding_bc_w_counts_below = .get_match(ll, '(Excluding barcodes with counts below )([0-9]+)'),
    used_probable_cell_barcodes = .get_match(ll, '(Using )([0-9]+)( probable cell barcodes)'),
    used_additional_barcodes    = .get_match(ll, '(plus an additional )([0-9]+)( barcodes)'),
    used_empty_droplets         = .get_match(ll, '(and )([0-9]+)( empty droplets)')
  ) %>%
    .[, (sample_var) := sample]
  assert_that( nrow(bender_dt) == 1 )

  return(bender_dt)
}

plot_reads_removed_as_ambient <- function(usa_dt_ls, ok_bcs_ls) {
  logit_brks  = c(1e-4, 1e-3, 1e-2, 0.10, 0.50, 0.90, 0.99, 0.999) %>% qlogis
  logit_labs  = c("0.01%", "0.1%", "1%", "10%", "50%", "90%", "99%", "99.9%")

  # prepare results
  plot_dt   = lapply(names(usa_dt_ls), function(nn) {
    usa_tmp   = usa_dt_ls[[ nn ]]
    bcs_tmp   = ok_bcs_ls[[ nn ]]
    plot_tmp  = usa_tmp[ barcode %in% bcs_tmp ] %>% 
      melt( id = "barcode", measure.vars = measure(value.name, transcript, sep = "_") ) %>% 
      .[, .(total_pre = sum(pre), total_post = sum(post)), by = barcode ] %>% 
      .[, logit_removed := qlogis(1 - total_post / total_pre) ] %>% 
      .[, sample_id := nn ]

    return( plot_tmp )
  }) %>% rbindlist

  # do plot
  g = ggplot(plot_dt) +
    aes( x = sample_id, 
      y = logit_removed %>% pmin(tail(logit_brks, 1)) %>% pmax(head(logit_brks, 1)) ) +
    geom_violin(colour = NA, fill = 'grey60',
      kernel = 'rectangular', adjust = 0.1, scale = 'width', width = 0.8) +
    scale_y_continuous( breaks = logit_brks, labels = logit_labs ) +
    theme_classic() +
    theme(
     axis.text.x = element_text( angle = -45, hjust = 0, vjust = 0.5), 
     plot.margin = unit(c(1, 4, 1, 1), "lines")) +
    labs( y = "pct. reads removed as ambient",x = NULL)

  return(g)
}




plot_spliced_vs_umis <- function(ss, usa_dt, ok_bcs, total_inc) {
  pscount     = 10
  PCT_BRKS    = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.5, 0.9, 0.97, 0.99, 0.997, 0.999) %>% qlogis
  PCT_LABS    = c('0.1%', '0.3%', '1%', '3%', '10%', '50%', '90%', '97%', '99%',
    '99.7%', '99.9%')
  LOG_BRKS    = c(pscount + c(0, 10, 30, 100, 300, 1000, 3000, 1e4, 3e4, 1e5, 3e5)) %>% log10
  LOG_LABS    = c('0', '10', '30', '100', '300', '1k', '3k', '10k', '30k', '100k', '300k')
  COL_VALS    = c(
    cell                        = 'grey80', 
    "no, used for\nambient"     = "#FB8072", 
    "no, excluded\ncompletely"  = "grey30")

  # arrange data
  plot_dt     = usa_dt[, .(
    barcode,
    pre_total     = pre_S + pre_U + pre_A,
    post_total    = post_S + post_U + post_A,
    pre_logit_S   = qlogis((pre_S + 1)/(pre_S + pre_U + 2)),
    post_logit_S  = qlogis((post_S + 1)/(post_S + post_U + 2))
    )] %>% 
    melt( measure.vars = measure(stage, value.name, pattern = "(pre|post)_(.+)")) %>% 
    .[, stage := factor(stage, levels = c("pre", "post")) ] %>% 
    .[ total > 0 ] %>% 
    .[ order(total) ]

  # to keep
  inc_dt      = plot_dt[ stage == "pre" ] %>% .[ order(-total) ] %>% .[ 1:total_inc ]
  plot_dt     = plot_dt %>% 
    .[, status   := ifelse( barcode %in% ok_bcs, names(COL_VALS)[[1]], 
      ifelse(barcode %in% inc_dt$barcode, names(COL_VALS)[[2]], names(COL_VALS)[[3]])) %>% 
      factor( levels = names(COL_VALS) ) ]

  # make plot
  g_dens  = ggplot(plot_dt) +
    aes( y = logit_S, x = log10(total + pscount) ) +
    geom_bin2d( bins = c(80, 50) ) +
    scale_x_continuous( breaks = LOG_BRKS, labels = LOG_LABS ) + expand_limits( x = log10(pscount) ) +
    scale_y_continuous( breaks = PCT_BRKS, labels = PCT_LABS ) +
    scale_fill_distiller( palette = "RdBu", trans = "log10" ) +
    facet_grid( . ~ stage ) +
    theme_classic() + 
    labs( y = "spliced pct.", x = 'total UMIs', colour = "called as cell?", title = ss )

  # make plot
  g_dots  = ggplot(plot_dt) +
    aes( y = logit_S, x = log10(total + pscount), colour = status ) +
    geom_point( size = 0.1 ) +
    scale_x_continuous( breaks = LOG_BRKS, labels = LOG_LABS ) + expand_limits( x = log10(pscount) ) +
    scale_y_continuous( breaks = PCT_BRKS, labels = PCT_LABS ) +
    scale_color_manual( values = COL_VALS,
      guide = guide_legend(override.aes = list(size = 3, alpha = 1)) ) +
    facet_grid( . ~ stage ) +
    theme_classic() + 
    labs( y = "spliced pct.", x = 'total UMIs', colour = "called as cell?" )

  # join together
  g = g_dens / g_dots + plot_layout( axes = "collect" )

  return(g)
}

calc_ambient_exclusions <- function(stats_dt, run_var) {
  exc_dt  = stats_dt[, .(get(run_var), total_droplets, kept_droplets, 
    pct_kept = round(prop_kept_by_cb * 100, 1), bad_run)] %>% 
    setnames("V1", run_var)

  return(exc_dt)
}

get_usa_dt <- function(usa_f, min_umi = 10) {
  usa_dt      = usa_f %>% fread
  row_sums    = usa_dt %>% .[, 2:4, with = FALSE ] %>% as.matrix %>% rowSums
  keep_rows   = row_sums >= min_umi
  return(usa_dt[ keep_rows ])
}
