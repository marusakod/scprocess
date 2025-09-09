# QC of cellbender input parameters

suppressPackageStartupMessages({
  library('data.table')
  library('magrittr')
  library('rhdf5')
  library('Matrix')
  library('assertthat')
  library('DropletUtils')

  library('gridExtra')
  library('tidyverse')
  library('Rfast')
  library('Seurat')
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

# sum spliced, unspliced and ambiguous counts for same gene
.sum_SUA <- function(sua_mat) {
  types = c('_S$', '_U$', '_A')
  
  mats  = lapply(types, function(t) sua_mat[grepl(t, rownames(sua_mat)), ])
  # check if symbols are all in the same order
  genes = lapply(mats, function(m) rownames(m) %>% str_before_first(pattern = '_'))

  assert_that(
    sapply(genes, function(gs) identical(gs, genes[[1]])) %>% all(), 
    msg = "gene names in split matrices don't match"
  )

  # remove suffixes from rownames
  mats = lapply(mats, function(m) {
    rownames(m) = str_before_first(rownames(m), pattern = '_')
    m
  })
  mats_sum = mats[[1]] + mats[[2]] + mats[[3]]
  return(mats_sum)
}

.get_alevin_mx <- function(af_mat_f, sel_s) {
  # get this file
  h5_filt   = H5Fopen(af_mat_f, flags = "H5F_ACC_RDONLY" )

  # get indices of barcodes
  mat       = sparseMatrix(
    i = as.vector(h5_filt$matrix$indices +1),
    p = as.vector(h5_filt$matrix$indptr),
    x = as.vector(h5_filt$matrix$data),
    repr = "C",
    dims = h5_filt$matrix$shape
  ) %>% as("TsparseMatrix")

  # add names
  bcs           = h5_filt$matrix$barcodes
  colnames(mat) = paste0(sel_s, bcs)
  rownames(mat) = as.character(h5_filt$matrix$features$name)

  return(mat)
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

# find slope at first inflection and total droplets included & expected_cells/total ratio
get_knee_params <- function(knee_f, sample_var) {
  ranks_df  = fread(knee_f)
  total_thr = unique(ranks_df$total_droplets_included) %>% log10()
  
  inf1 = unique(ranks_df$inf1)
  
  # get x coordinate of inf1
  infl1_idx = which.min( abs(ranks_df$total - inf1) )[1]
  
  # get x coordinate of inf1
  inf_1_x = ranks_df[ infl1_idx, rank ] %>%
    log10()
  
  # fit curve to all points
  ranks_df = ranks_df %>% 
    .[total > 5] %>%
    .[, `:=`(
      ranks_log = log10(rank),
      total_log = log10(total)
    )
    ] %>%
    unique
  
  fit = smooth.spline(x = ranks_df$ranks_log, y = ranks_df$total_log)
  fitted.vals = 10^fitted(fit)
  
  # get value of the first derivative at total included and inflection1
  d1       = predict(fit, deriv=1)
  d1_inf   = d1$y[ which.min(abs(d1$x - inf_1_x))[1] ]
  d1_total = d1$y[ which.min(abs(d1$x - total_thr))[1] ]
  
  keep_cols = c(sample_var, 'knee1', 'inf1', 'knee2', 'inf2', 'total_droplets_included', 'expected_cells')
  
  final = ranks_df %>%
    .[, ..keep_cols] %>%
    unique() %>%
    .[, `:=`(
      slope_inf1 = d1_inf,
      slope_total_included = d1_total
    )]%>% 
    .[, `:=`(
      slope_ratio = abs(slope_total_included) / abs(slope_inf1),
      expected_total_ratio = expected_cells / total_droplets_included
    )]
  
  return(final)
}


plot_barcode_ranks_w_params <- function(knee_fs, ambient_knees_df, sample_var, bender_priors_df = NULL, show_lines = TRUE) {
  
  s_ord = names(knee_fs)
  
  # Add knee and inflection to params
  knee_data = lapply(s_ord, function(s) {
    knee_f = knee_fs[[s]]
    x = fread(knee_f) %>% as.data.table
    x %>%
      .[, .(n_bc = .N), by = .(lib_size = total)] %>%
      .[order(-lib_size)] %>%
      .[, bc_rank := cumsum(n_bc)] %>%
      .[, (sample_var) := s]
  }) %>% rbindlist()
  
  knee_vars = c(sample_var, 'knee1', 'inf1', 'knee2', 'inf2',
                'total_droplets_included', 'expected_cells')
  
  lines_knees = ambient_knees_df %>% as.data.table %>% 
    .[ get(sample_var) %in% s_ord, ..knee_vars] %>%
    setnames( c("inf1", "inf2"), c("shin1", "shin2") ) %>% 
    melt(id.vars = sample_var) %>%
    .[, `:=`(
      axis = fifelse(variable %in% c('knee1', 'knee2', 'shin1', 'shin2'), 'y', 'x'),
      type = fifelse(variable %in% c('knee1', 'knee2', 'shin1', 'shin2'),
                     'cellbender intermediate\nparameter', 
                     'cellbender input\nparameter')
    )]
  
  if ( is.null(bender_priors_df) ) {
    lines_priors = NULL
  } else {
    prior_vars = c(sample_var, 'cb_prior_cells', 'cb_prior_empty')
    lines_priors = bender_priors_df %>% as.data.table() %>%
      .[get(sample_var) %in% s_ord, ..prior_vars] %>%
      melt(id.vars= sample_var) %>%
      .[, `:=`(
        axis = 'y',
        type = 'cellbender prior\nparameter'
      )]
  }
  
  lines = list(lines_knees, lines_priors) %>% rbindlist()
  
  hlines = lines %>% filter(axis == 'y')
  vlines = lines %>% filter(axis == 'x')
  
  # plot everything above low count threshold
  p_labels =  c("1", "10", "100", "1k", "10k", "100k", "1M")
  p_breaks =  c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6)
  
  # set factor levels for sample_var so the samples will appear in the right order in the plot
  knee_data[[sample_var]] = factor(knee_data[[sample_var]], levels = s_ord)
  if (show_lines) {
    hlines[[sample_var]] = factor(hlines[[sample_var]], levels = s_ord)
    vlines[[sample_var]] = factor(vlines[[sample_var]], levels = s_ord)
  }
  
  p = ggplot(knee_data) +
    aes(x = bc_rank, y = lib_size) +
    geom_line(linewidth = 0.3, color = '#283747') +
    facet_wrap( ~ get(sample_var), ncol = 4 ) +
    scale_x_log10(labels = p_labels, breaks = p_breaks) +
    scale_y_log10(labels = p_labels, breaks = p_breaks) +
    scale_color_manual(
      values = c("#7c4b73", "#88a0dc", "#ab3329"),
      breaks = c('cellbender input\nparameter', 'cellbender intermediate\nparameter',
                 'cellbender prior\nparameter')) +
    theme_classic(base_size = 9) +
    theme(legend.position = 'none') +
    labs(x = 'barcode rank', y = 'library size', color = NULL)
  
  # add lines only if show_lines is TRUE
  if (show_lines) {
    p = p +
      geom_hline(data = hlines,
                 mapping = aes(yintercept = value, color = type)) +
      geom_vline(data = vlines,
                 mapping = aes(xintercept = value, color = type)) +
      geom_text_repel(data = hlines, mapping = aes(y = value, x = 10, label = variable),
                      size = 2.5) +
      geom_text_repel(data = vlines, mapping = aes(x = value, y = 100, label = variable),
                      size = 2.5, angle = 90)
  }
  
  return(p)
}



find_outlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

# boxplots of log ratios of slopes and barcode percents
# log: get outliers on the log10 scale
plot_amb_params_dotplot <- function(params_qc, sample_var, scales = 'fixed') {
  all_scales_opts = c('fixed', 'free')
  scale = match.arg(scales, all_scales_opts)
  
  # get outliers
  outliers_te = params_qc$expected_total_ratio %>%
    set_names(params_qc[[sample_var]]) %>%
    log10() %>%
    find_outlier(.) %>% .[.] %>%
    names()
  
  outliers_slopes = params_qc$slope_ratio %>%
    set_names(params_qc[[sample_var]]) %>%
    log10() %>%
    find_outlier(.) %>% .[.] %>%
    names()
  
  keep_cols = c(sample_var, 'slope_ratio', 'expected_total_ratio')
  plot_df = params_qc %>% 
    .[, ..keep_cols] %>%
    melt(id.vars = sample_var) %>%
    .[, `:=`(
      is.outlier_slope = fifelse(get(sample_var) %in% outliers_slopes, TRUE, FALSE), 
      is.outlier_te    = fifelse(get(sample_var) %in% outliers_te, TRUE, FALSE)
    )
     ]
  
  outlier_df_slope = copy(plot_df) %>%
    .[is.outlier_slope == TRUE & variable == 'slope_ratio']
    
  outlier_df_te = copy(plot_df) %>% 
    .[is.outlier_te == TRUE & variable == 'expected_total_ratio'] 
  
  pl =  ggplot(plot_df, aes(x = variable, y = value) ) +
    geom_quasirandom( fill = 'grey', shape = 21, size = 3 ) +
    labs(x = NULL, y = 'ratio') +
    ggrepel::geom_text_repel(data = outlier_df_slope, mapping = aes(label = get(sample_var))) +
    ggrepel::geom_text_repel(data = outlier_df_te, mapping = aes(label = get(sample_var))) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12)) +
    scale_y_log10()
  
  if (scales == 'free') {
    pl =  pl + facet_wrap(~variable, scales ='free') +
      theme(axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            strip.text = element_text(size = 12))
  }
  
  return(pl)
}

get_amb_sample_level_qc <- function(qc, sel_s, amb_method = c('cellbender', 'decontx')) {

  amb = match.arg(amb_method)

  sum_qc = qc %>%
    as_tibble() %>%
    dplyr::select(-barcode) %>%
    rowwise()

  if (amb == 'cellbender') {
  sum_qc = sum_qc %>%
    mutate(af_all = sum(af_S, af_U, af_A),
           cb_all = sum(cb_S, cb_U, cb_A)) %>%
    colSums()

  smpl_qc = c((sum_qc['af_S']/sum_qc['af_all']) *100,
                 (sum_qc['cb_S']/sum_qc['cb_all'])*100,
                 sum_qc['af_all'],
                 sum_qc['af_all'] - sum_qc['cb_all'])

  } else {
    sum_qc = sum_qc %>%
      mutate(af_all = sum(af_S, af_U, af_A),
             dcx_all = sum(dcx_S, dcx_U, dcx_A)) %>%
      colSums()

    smpl_qc <- c((sum_qc['af_S']/sum_qc['af_all']) *100,
                   (sum_qc['dcx_S']/sum_qc['dcx_all'])*100,
                   sum_qc['af_all'],
                   sum_qc['af_all'] - sum_qc['dcx_all'])


  }

  names(smpl_qc) = c('af_spliced_pct', 'amb_spliced_pct', 'af_all', 'removed')
  smpl_qc$sample_id = sel_s

  return(smpl_qc)
}

get_amb_sample_qc_outliers <- function(qc_df, var1, var2) {
  bivar =  qc_df %>% dplyr::select(all_of(c(var1, var2)))

  mcd = robustbase::covMcd(bivar)
  chi_thr = chi_threshold <- qchisq(0.95, df = 2)
  outliers_df = qc_df[which(mcd$mah > chi_threshold), ]

  return(outliers_df)
}

make_amb_sample_qc_oulier_plots <- function(qc_df, var1, var2, outliers_df,
  x_title, y_title, y_thr = NULL, x_thr = NULL) {

  p = ggplot(qc_df, aes(x = get(var1), y = get(var2))) +
    geom_point(shape = 21, fill = 'grey', color = 'black') +
    labs(x = x_title,
         y = y_title) +
    theme_classic()

  if (nrow(outliers_df) != 0) {
    p = p + geom_text_repel(data = outliers_df, mapping = aes(label = sample_id), size = 3)
  }

  if (!is.null(y_thr)) {
    p = p + geom_hline(yintercept = y_thr, linewidth = 0.2, linetype = 'dashed')
  }

  if (!is.null(x_thr)) {
    p = p + geom_vline(xintercept = x_thr, linewidth = 0.2, linetype = 'dashed')
  }

  return(p)
}

get_usa_dt <- function(usa_f, min_umi = 10) {
  usa_dt      = usa_f %>% fread
  row_sums    = usa_dt %>% .[, 2:4, with = FALSE ] %>% as.matrix %>% rowSums
  keep_rows   = row_sums >= min_umi
  return(usa_dt[ keep_rows ])
}

plot_qc_metrics_split_by_cells_empties <- function(rna_knee_dfs, 
  metric = c("umis", "splice_pct"), sample_var = "sample_id", min_umis = 10) {
  metric    = match.arg(metric)

  # get cells and empties
  plot_dt   = knee_fs %>% lapply(function(f) {
    tmp_dt = fread(f) %>% 
      .[rank <= expected_cells | in_empty_plateau == TRUE ] %>% 
      .[, `:=`(
        umis       = log10(total), 
        splice_pct = qlogis((spliced + 1) / (spliced + unspliced + 2)), 
        what       = fifelse(rank <= expected_cells, "cell", "empty")
      )] %>%
      .[, c(sample_var, 'barcode', 'rank', 'umis', 'splice_pct', 'what'), with = FALSE]
    }) %>% rbindlist

  # plot these
  if (metric == "umis") {
    y_brks    = c(1e0, 1e1, 3e1, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6) %>% log10
    y_labs    = c("1", "10", "30", "100", "300", "1k", "3k", "10k", "30k", "100k", 
      "300k", "1M")
    y_title   = "UMIs"
  } else if (metric == "splice_pct") {
    y_brks    = c(0.01, 0.03, 0.1, 0.3, 0.5, 0.7, 0.9, 0.97, 0.99) %>% qlogis
    y_labs    = c("1%", "3%", "10%", "30%", "50%", "70%", "90%", "97%", "99%")
    y_title   = "spliced pct."
  }
  g = ggplot(plot_dt) + 
    aes( fill = what, x = get(sample_var), y = get(metric) ) +
    geom_violin( colour = NA,
      kernel = 'rectangular', adjust = 0.1, scale = 'width') +
    scale_y_continuous( breaks = y_brks, labels = y_labs ) +
    scale_fill_manual( values = c(cell = "#1965B0", empty = "grey") ) +
    theme_classic() + 
    theme( axis.text.x = element_text( angle = -45, hjust = 0, vjust = 0.5 ) ) +
    labs( y = y_title, fill = "what does\nthe barcode\ncontain?" )

  return(g)
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
    theme( axis.text.x = element_text( angle = -45, hjust = 0, vjust = 0.5) ) +
    labs( y = "pct. reads removed as ambient")

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

calc_ambient_exclusions <- function(stats_dt, sample_var) {
  exc_dt  = stats_dt[, .(get(sample_var), total_droplets, kept_droplets, 
    pct_kept = round(prop_kept_by_cb * 100, 1), bad_sample)] %>% 
    setnames("V1", sample_var)

  return(exc_dt)
}
