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
get_cell_mat_and_barcodes <- function(out_mat_f, out_bcs_f, out_dcx_f = NULL, sel_s, af_mat_f,
                                      knee_1 = NULL, knee_2 = NULL, inf_1 = NULL,
                                      total_included = NULL , exp_cells = NULL,
                                      cell_calls_method = c('barcodeRanks', 'emptyDrops'),
                                      ncores = 4, niters = 1000, hvg_n = 2000,
                                      ambient_method = c('none', 'decontx')){
  # check inputs
  call_m = match.arg(cell_calls_method)
  amb_m = match.arg(ambient_method)

  # load alevin matrix
  af_mat = .get_alevin_mx(af_mat_f = af_mat_f, sel_s = '')

  # call cells and empties
  cell_empty_bcs = call_cells_and_empties(af_mat = af_mat,
                                          knee_1 = knee_1,
                                          knee_2 = knee_2,
                                          inf_1 = inf_1,
                                          total_included = total_included,
                                          exp_cells = exp_cells,
                                          ncores = ncores,
                                          n_iters = niters,
                                          call_method = call_m)

  message('Cells and empty droplets found for ', sel_s)

  # get matrix with cells
  cell_bc_df = data.frame(barcode = cell_empty_bcs[["cells"]])
  cell_mat = af_mat[, cell_bc_df$barcode]

  # do ambient correction
  if(ambient_method == 'none'){

    message("Ambient RNA removal method is 'none'. Saving uncorrected count matrix for ", sel_s )
    # return uncorrected matrix with just cells
    write.table(cell_bc_df, out_bcs_f, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = ',')
    write10xCounts(out_mat_f, cell_mat, version = "3", overwrite = TRUE)


  }else{

    message("Running 'decontx' for ", sel_s)
    # get empty matrix
    empty_mat = af_mat[, cell_empty_bcs[["empty"]]]
    dcx_res = decontX(x = cell_mat, background = empty_mat, varGenes = hvg_n)

    # save decontaminated matrix
    clean_mat = dcx_res$decontXcounts %>% round()
    write10xCounts(out_mat_f, clean_mat, version = "3", overwrite = TRUE)

    # save barcodes
    clean_bcs_df = data.frame(barcode = colnames(clean_mat))
    write.table(clean_bcs_df, out_bcs_f, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = ',')

    # add save contamination proportions and decontx cluster assignment
    dcx_params = data.frame(barcode = colnames(clean_mat),
                            pct.contamination = dcx_res$contamination,
                            dcx_cluster = dcx_res$z %>% sprintf("cl%02d", .) %>% factor)

    fwrite(dcx_params, out_dcx_f)

    message("'decontx' completed for ", sel_s)

  }

  return(NULL)
}


call_cells_and_empties <- function(af_mat,
                       knee_1 = NULL, knee_2 = NULL, inf_1 = NULL,
                       total_included = NULL , exp_cells = NULL,
                       ncores = 4, n_iters = 1000, fdr_thr = 0.001,
                       call_method = c('barcodeRanks', 'emptyDrops')){

  # get barcode ranks
  ranks = barcodeRanks(af_mat) %>% as.data.frame() %>%
    rownames_to_column('barcode') %>%
    arrange(rank)

  # get empty droplets
    empty_locs = .get_empty_plateau(ranks_df = ranks,
                                    inf1 = inf_1,
                                    total_inc = total_included,
                                    knee2 = knee_2)

 print(empty_locs)

  # get indices of empty barcodes
  empty_bcs   = ranks %>%
    filter(between(rank, empty_locs[1], empty_locs[2])) %>%
    pull(barcode)
  empty_idx = which(colnames(af_mat) %in% empty_bcs)

  if(call_method == 'barcodeRanks'){
    assert_that(!is.null(exp_cells))
    cell_bcs = ranks$barcode[1: exp_cells]


  }else{

    # get sum of s+u+a counts (instead of this maybe try removing genes with v low expression)
    # doing this because othewise takes too long
    af_mat_sum =.sum_SUA(af_mat)

    bpparam = MulticoreParam(workers = ncores, progressbar = TRUE)
    # call cells

    edrops_res = emptyDrops(m = af_mat_sum,
                            niters = n_iters,
                            BPPARAM = bpparam,
                            known.empty = empty_idx,
                            retain =  knee_1)

    # get cell barcodes
    cell_bcs = edrops_res %>% as.data.frame() %>%
      filter(FDR <= fdr_thr) %>%
      rownames()

  }

  # return a list with both cell barcodes and empty barcodes
  empty_cell_bcs_ls = list(empty = empty_bcs,
                           cells = cell_bcs)
  return(empty_cell_bcs_ls)

}



# sum spliced, unspliced and ambiguous counts for same gene
.sum_SUA <- function(sua_mat){
  types = c('_S$', '_U$', '_A')
  mats = lapply(types, function(t) sua_mat[grepl(t, rownames(sua_mat)), ])
  # check if symbols are all in the same order
  genes = lapply(mats, function(m) rownames(m) %>% str_before_first(pattern = '_'))

  assert("gene names in split matrices don't match",
         sapply(genes, function(gs) identical(gs, genes[[1]])) %>% all())

  # remove suffixes from rownames
  mats = lapply(mats, function(m){
    rownames(m) = str_before_first(rownames(m), pattern = '_')
    m
  })

  mats_sum = mats[[1]] + mats[[2]] + mats[[3]]
  return(mats_sum)

}



.get_alevin_mx <- function(af_mat_f, sel_s){
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


.get_empty_plateau <- function(ranks_df, inf1, total_inc, knee2){

  # calculate empty plateau

   infl1_idx = which.min( abs(ranks_df$total - inf1) )[1]
   infl1_x   = ranks_df[ infl1_idx, "rank" ]

  empty_start = ranks_df %>% mutate(n = 1:nrow(.)) %>%
    filter( between(rank, infl1_x, total_inc)) %>%
    pull(n) %>%
    log10 %>%
    mean() %>%
    10^.

  empty_end = unique(ranks_df[ranks_df$total == knee2, "rank"])
  empty_plateau = c(empty_start, empty_end) %>% as.numeric

  return(empty_plateau)
}


#
# save_cellbender_qc_metrics <- function(knee_data_f, af_h5_f, cb_h5_f, cb_qc_f, do_cellbender) {
#   # get expected cells from cellbender yaml
#   knee_dt     = knee_data_f %>% fread
#   n_exp_cells = knee_dt$expected_cells %>% unique
#   exp_cells   = knee_dt$barcode[ 1:n_exp_cells ]
#
#   # get alevin counts
#   af_dt       = .get_alevin_USA_sums(af_h5_f, exp_cells)
#   if (as.logical(do_cellbender)) {
#     cb_dt       = .get_cellbender_USA_sums(cb_h5_f, exp_cells)
#     assert_that( all( sort(af_dt$barcode) == sort(cb_dt$barcode) ))
#
#     # join together
#     cb_qc_dt    = merge(af_dt, cb_dt, by = "barcode") %>%
#       setkey("barcode") %>% .[ exp_cells ]
#   } else {
#     cb_qc_dt    = af_dt %>%
#       setkey("barcode") %>% .[ exp_cells ]
#   }
#   fwrite(cb_qc_dt, file = cb_qc_f)
# }
#
#
# .get_alevin_USA_sums <- function(af_h5_f, exp_cells) {
#   sce_af      = read10xCounts(af_h5_f)
#   af_idx      = sce_af$Barcode %in% exp_cells
#   mat_af      = counts(sce_af)[, af_idx ] %>%
#     as("sparseMatrix") %>%
#     set_colnames(sce_af$Barcode[ af_idx ])
#
#   af_dt       = .get_usa_dt(mat_af, "af")
#
#   return(af_dt)
# }



save_barcode_qc_metrics <- function(af_h5_f, amb_out_yaml, out_qc_f, expected_cells, ambient_method) {

  # read in the yaml file
  amb_yaml = yaml::read_yaml(amb_out_yaml)

  # get the right ambient matrix file path from yaml
  if(ambient_method == 'cellbender'){
    # get unfiltered cellbender matrix
    amb_mat_f = amb_yaml$cb_full_f
  }else if(ambient_method == 'decontx'){
    amb_mat_f = amb_yaml$dcx_filt_f
  }else{
    amb_mat_f = amb_yaml$cell_filt_f
  }


  # get alevin counts
  af_mat = .get_alevin_mx(af_h5_f, '')

  if(ambient_method == 'cellbender'){
    # get ranks
    ranks = barcodeRanks(af_mat) %>% as.data.frame() %>%
      rownames_to_column('barcode') %>%
      arrange(rank)

    # get expected barcodes
    exp_cells = ranks$barcode[1:expected_cells]

    # get cellbender matrix
    cb_mat = .get_alevin_mx(amb_mat_f, '')

    # restrict to barcodes that we care about
    af_filt_mat = af_mat[, exp_cells]
    cb_filt_mat = cb_mat[, exp_cells]

    # sum s/u/a counts for each barcode
    af_dt = .get_usa_dt(af_filt_mat, prefix = 'af')
    cb_dt = .get_usa_dt(cb_filt_mat, prefix = 'cb')

    # merge together
    qc_dt = merge(af_dt, cb_dt, by = 'barcode') %>%
      setkey("barcode") %>% .[ exp_cells ]

  }else if(ambient_method == 'decontx'){

    # get decontx matrix
    dcx_mat = .get_alevin_mx(amb_mat_f, '')
    # get same barcodes from raw alevin mat
    af_filt_mat = af_mat[, colnames(dcx_mat)]

    # sum s/u/a counts for each barcode
    af_dt = .get_usa_dt(af_filt_mat, prefix = 'af')
    dcx_dt = .get_usa_dt(dcx_mat, prefix = 'dcx')

    # merge together
    qc_dt = merge(af_dt, dcx_dt, by = 'barcode')

  }else{

    cell_mat = .get_alevin_mx(amb_mat_f, '')
    qc_dt = .get_usa_dt(cell_mat, prefix = 'af')

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


# get_bender_log_file <- function(sample_name, bender_dir) {
#   b_dir_full <- file.path(bender_dir, paste0('bender_', sample_name))
#   log_f <- list.files(b_dir, pattern = '.log', full.names = TRUE)
#   return(log_f)
# }

# test comment

get_bender_log <- function(f, sample) {
  ll <- read_lines(f, n_max = 25)
  .get_match <- function(ll, pat) {
    ll %>% str_subset(pat) %>% str_match(pat) %>% .[, 3] %>% as.integer()
  }
  bender_df <- tibble(
    sample_id                   = sample,
    cb_prior_cells              = .get_match(ll, '(Prior on counts for cells is )([0-9]+)'),
    cb_prior_empty              = .get_match(ll, '(Prior on counts [infor]{2,3} empty droplets is )([0-9]+)'),
    excluding_bc_w_counts_below = .get_match(ll, '(Excluding barcodes with counts below )([0-9]+)'),
    used_probable_cell_barcodes = .get_match(ll, '(Using )([0-9]+)( probable cell barcodes)'),
    used_additional_barcodes    = .get_match(ll, '(plus an additional )([0-9]+)( barcodes)'),
    used_empty_droplets         = .get_match(ll, '(and )([0-9]+)( empty droplets)')
  )
  assert_that( nrow(bender_df) == 1 )

  return(bender_df)
}

# find slope at first inflection and total droplets included & expected_cells/total ratio
get_knee_params <- function(ranks_df) {
  total_thr = unique(ranks_df$total_droplets_included) %>%
    log10()
  inf1 = unique(ranks_df$inf1)

  # get x coordinate of inf1
  infl1_idx = which.min( abs(ranks_df$total - inf1) )[1]

  # get x coordinate of inf1
  inf_1_x = ranks_df[ infl1_idx, "rank" ] %>%
    log10()

  # fit curve to all points
  ranks_df = ranks_df %>% filter(total > 5) %>%
    mutate(ranks_log  = log10(rank),
           total_log = log10(total)) %>%
    distinct()

  fit = smooth.spline(x = ranks_df$ranks_log, y = ranks_df$total_log)
  fitted.vals = 10^fitted(fit)

  # get value of the first derivative at total included and inflection1
  d1 = predict(fit, deriv=1)
  d1_inf = d1$y[ which.min(abs(d1$x - inf_1_x))[1] ]
  d1_total = d1$y[ which.min(abs(d1$x - total_thr))[1] ]

  # return(c(d1_inf, d1_total) %>% setNames(c('slope_inf1', 'slope_total_included')))
  final = ranks_df %>% dplyr::select(sample_id, knee1, inf1, knee2, inf2, total_droplets_included, expected_cells) %>%
    distinct()

  final = final %>%
    mutate(slope_inf1 = d1_inf,
           slope_total_included = d1_total) %>%
    mutate(slope_ratio = abs(slope_total_included)/abs(slope_inf1),
           expected_total_ratio = expected_cells/total_droplets_included) %>%
    as_tibble()

  return(final)
}

make_br_plot1 <- function(ranks) { # ranks == output of calc_cellbender_params()
  # add knee and inflection to params
  knee_data <- ranks %>%
    as_tibble() %>%
    group_by(lib_size = total) %>%
    summarize(n_bc = n()) %>%
    arrange(desc(lib_size)) %>%
    mutate(bc_rank = cumsum(n_bc))

  lines <- ranks[, 4:11] %>% distinct() %>%
    reshape2::melt(.) %>%
    mutate(axis = c(rep('y', 4), rep('x', 2)),
           type = c(rep('cellbender intermediate\nparameter', 4),
                    rep('cellbender input\nparameter', 2)))

  hlines <- lines %>% filter(axis == 'y')
  vlines <- lines %>% filter(axis == 'x')
  # plot everything above low count threshold

  p_labels =  c("10", "100", "1k", "10k", "100k", "1M")
  p_breaks =  c(10, 100, 1000, 10000, 100000, 1000000)

  p <- ggplot(knee_data,
              aes(x = bc_rank, y = lib_size)) +
    geom_line(linewidth = 0.5, color = '#283747') +
    scale_x_log10(labels = p_labels,
                  breaks = p_breaks) +
    scale_y_log10(labels = p_labels,
                  breaks = p_breaks) +
    theme_light() +
    geom_hline(data = hlines,
               mapping = aes(yintercept = value, color = type)) +
    geom_vline(data = vlines,
               mapping = aes(xintercept = value, color = type)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = 'none') +
    labs(x = 'barcode rank',
         y = 'library size',
         color = NULL,
         title = unique(ranks$sample_id),
         subtitle = unique(ranks$study_id)) +
    scale_color_manual(values = c('#E74C3C', '#5DADE2'),
                       breaks = c('cellbender input\nparameter',
                                  'cellbender intermediate\nparameter')) +
    geom_text_repel(data = hlines, mapping = aes(y = value, x = 10,  label = variable)) +
    geom_text_repel(data = vlines, mapping = aes(x = value, y = 100, label = variable), angle = 90)

  return(p)
}



plot_barcode_ranks_w_params <- function(knees, ambient_knees_df, bender_priors_df = NULL) {
  # get sample order
  s_ord <- names(knees)

  # add knee and inflection to params
  knee_data <- lapply(s_ord, function(s) {
    x <- knees[[ s ]]
    x %>%
      as_tibble() %>%
      group_by(lib_size = total) %>%
      summarize(n_bc = n()) %>%
      arrange(desc(lib_size)) %>%
      mutate(bc_rank = cumsum(n_bc)) %>%
      mutate(sample_id = s)
    }) %>% bind_rows()

  knee_vars <- c('sample_id', 'knee1', 'inf1', 'knee2', 'inf2',
    'total_droplets_included', 'expected_cells')
  lines_knees <- ambient_knees_df[ambient_knees_df$sample_id %in% s_ord, knee_vars] %>%
    reshape2::melt(., id = "sample_id") %>%
    mutate(
      axis = case_when(
        variable %in% c('knee1', 'knee2', 'inf1', 'inf2') ~ 'y',
        TRUE ~ 'x'),
      type = case_when(
        variable %in% c('knee1', 'knee2', 'inf1', 'inf2') ~ 'cellbender intermediate\nparameter',
        TRUE ~ 'cellbender input\nparameter')
    )
  if ( is.null(bender_priors_df) ) {
    lines_priors <- NULL
  } else {
    prior_vars <- c('sample_id', 'cb_prior_cells', 'cb_prior_empty')
    lines_priors <- bender_priors_df[ bender_priors_df$sample_id %in% s_ord, prior_vars] %>%
      reshape2::melt(., id = "sample_id") %>%
      mutate( axis = 'y', type = 'cellbender prior\nparameter')
  }

  lines <- list(lines_knees, lines_priors) %>% bind_rows

  hlines <- lines %>% filter(axis == 'y')
  vlines <- lines %>% filter(axis == 'x')

  # plot everything above low count threshold
  p_labels =  c("1", "10", "100", "1k", "10k", "100k", "1M")
  p_breaks =  c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6)

  # set factor levels for sample_id so the samples will appear in the right order in the plot
  knee_data$sample_id <- factor( knee_data$sample_id, levels = s_ord )
  hlines$sample_id <- factor(hlines$sample_id, levels = s_ord )
  vlines$sample_id <- factor(vlines$sample_id, levels = s_ord )

  p <- ggplot(knee_data) +
    aes( x = bc_rank, y = lib_size ) +
    geom_line(linewidth = 0.3, color = '#283747') +
    geom_hline(data = hlines,
               mapping = aes(yintercept = value, color = type)) +
    geom_vline(data = vlines,
               mapping = aes(xintercept = value, color = type)) +
    geom_text_repel(data = hlines, mapping = aes(y = value, x = 10,  label = variable),
      size = 2.5) +
    geom_text_repel(data = vlines, mapping = aes(x = value, y = 100, label = variable),
      size = 2.5, angle = 90) +
    facet_wrap( ~ sample_id, ncol = 4 ) +
    scale_x_log10(labels = p_labels,
                  breaks = p_breaks) +
    scale_y_log10(labels = p_labels,
                  breaks = p_breaks) +
    scale_color_manual(
      # values = c("#FAD510", "#CB2314", "#273046"),
      values = c("#7c4b73", "#88a0dc", "#ab3329"),
      breaks = c('cellbender input\nparameter', 'cellbender intermediate\nparameter',
        'cellbender prior\nparameter')) +
    theme_classic(base_size = 9) +
    theme( legend.position = 'none' ) +
    labs(x = 'barcode rank', y = 'library size', color = NULL)

  return(p)
}

find_outlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

# boxplots of log ratios of slopes and barcode percents
# log: get outliers on the log10 scale
make_amb_params_dotplot <- function(params_qc, log = FALSE, scales = 'fixed') {
  all_scales_opts = c('fixed', 'free')
  scale = match.arg(scales, all_scales_opts)

  # get outliers
  outliers_te = params_qc$expected_total_ratio %>%
    set_names(params_qc$sample_id)

  outliers_slopes = params_qc$slope_ratio %>%
    set_names(params_qc$sample_id)

  if (log == TRUE) {
    outliers_te = outliers_te %>% log10()
    outliers_slopes = outliers_slopes %>% log10()
  }

  outliers_te = outliers_te %>%
    find_outlier(.) %>%
    .[.] %>%
    names()

  outliers_slopes = outliers_slopes %>%
    find_outlier(.) %>%
    .[.] %>%
    names()

  plot_df = params_qc %>% dplyr::select(sample_id, slope_ratio, expected_total_ratio) %>%
    reshape2::melt(., id = "sample_id") %>%
    mutate(is.outlier_slope = case_when(sample_id %in% outliers_slopes ~ TRUE,
                                  TRUE ~ FALSE)) %>%
    mutate(is.outlier_te = case_when(sample_id %in% outliers_te ~ TRUE,
         TRUE ~ FALSE))

  outlier_df_slope = plot_df %>% filter(is.outlier_slope == TRUE) %>%
    filter(variable == 'slope_ratio')
  outlier_df_te = plot_df %>% filter(is.outlier_te == TRUE) %>%
    filter(variable == 'expected_total_ratio')

 pl =  ggplot(plot_df, aes(x = variable, y = value) ) +
    geom_quasirandom( fill = 'grey', shape = 21, size = 3 ) +
    labs(x = NULL, y = 'ratio') +
    ggrepel::geom_text_repel(data = outlier_df_slope, mapping = aes(label = sample_id)) +
    ggrepel::geom_text_repel(data = outlier_df_te, mapping = aes(label = sample_id)) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12))

  if (log == TRUE) {
    pl  = pl + scale_y_log10()
  }

  if (scales == 'free') {
   pl =  pl + facet_wrap(~variable, scales ='free') +
      theme(axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            strip.text = element_text(size = 12))
  }

 return(pl)
}



get_amb_sample_level_qc <- function(qc, sel_s, amb_method = c('cellbender', 'decontx')){

  amb = match.arg(amb_method)

  sum_qc = qc %>%
    as_tibble() %>%
    dplyr::select(-barcode) %>%
    rowwise()

  if(amb == 'cellbender'){
  sum_qc = sum_qc %>%
    mutate(af_all = sum(af_S, af_U, af_A),
           cb_all = sum(cb_S, cb_U, cb_A)) %>%
    colSums()

  smpl_qc = c((sum_qc['af_S']/sum_qc['af_all']) *100,
                 (sum_qc['cb_S']/sum_qc['cb_all'])*100,
                 sum_qc['af_all'],
                 sum_qc['af_all'] - sum_qc['cb_all'])

  }else{
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




get_amb_sample_qc_outliers <- function(qc_df, var1, var2){
  bivar =  qc_df %>% dplyr::select(all_of(c(var1, var2)))

  mcd = robustbase::covMcd(bivar)
  chi_thr = chi_threshold <- qchisq(0.95, df = 2)
  outliers_df = qc_df[which(mcd$mah > chi_threshold), ]

  return(outliers_df)
}



make_amb_sample_qc_oulier_plots <- function(qc_df, var1, var2, outliers_df,
                                           x_title, y_title, y_thr = NULL, x_thr = NULL){

  p = ggplot(qc_df, aes(x = get(var1), y = get(var2))) +
    geom_point(shape = 21, fill = 'grey', color = 'black') +
    labs(x = x_title,
         y = y_title) +
    theme_classic()

  if(nrow(outliers_df) != 0){
    p = p + geom_text_repel(data = outliers_df, mapping = aes(label = sample_id), size = 3)
  }

  if(!is.null(y_thr)){
    p = p + geom_hline(yintercept = y_thr, linewidth = 0.2, linetype = 'dashed')
  }

  if(!is.null(x_thr)){
    p = p + geom_vline(xintercept = x_thr, linewidth = 0.2, linetype = 'dashed')
  }

  return(p)

}


