# alevin_fry.R

suppressPackageStartupMessages({
  library("magrittr")
  library("fishpond")
  library("SingleCellExperiment")
  library("DropletUtils")
  library("tidyverse")
  library("data.table")
  library("parallel")
  library("testit")
  library('strex')
  library('BiocParallel')
})

# load counts data into sce object
save_alevin_h5_ambient_params <- function(run, fry_dir, h5_f, cb_yaml_f, knee_data_f, 
  run_var, knee1, shin1, knee2, shin2, exp_cells, total_included, low_count_thr) {
  # load the data, save to h5
  bender_ps = save_alevin_h5_knee_params_df(run, fry_dir, h5_f, knee_data_f, 
    hto_mat = 0, run_var, knee1, shin1, knee2, shin2, 
    exp_cells, total_included, low_count_thr)

  # write these parameters to yaml file
  con_obj     = file(cb_yaml_f)
  writeLines(c(
    sprintf("run: %s",                          run),
    sprintf("cb_total_droplets_included: %.0f", unique(bender_ps$total_droplets_included) ),
    sprintf("cb_expected_cells: %.0f",          unique(bender_ps$expected_cells) ),
    sprintf("cb_low_count_threshold: %.0f",     unique(bender_ps$low_count_threshold) ),
    sprintf("knee1: %.0f",                      unique(bender_ps$knee1)),
    sprintf("shin1: %.0f",                      unique(bender_ps$shin1)),
    sprintf("knee2: %.0f",                      unique(bender_ps$knee2)),
    sprintf("shin2: %.0f",                      unique(bender_ps$shin2))
  ), con = con_obj)
  close(con_obj)
}

save_alevin_h5_knee_params_df <- function(run, fry_dir, h5_f, knee_data_f, 
  hto_mat = 0, run_var, knee1 = '', shin1 = '', knee2 = '', shin2 ='',
  exp_cells ='', total_included ='', low_count_thr ='') {
  # load the data
  if (hto_mat) {
    sce = loadFry(fry_dir)
    mat = counts(sce)
  } else {
    # get sce object
    sce = loadFry(
      fry_dir,
      outputFormat = list(S = c("S"), U = c("U"), A = c("A"))
      )

    # convert to matrix
    mat = assayNames(sce) %>% lapply(function(n) {
      mat       = assay(sce, n)
      rownames(mat) = paste0(rownames(mat), "_", n)
      return(mat)
      }) %>% do.call(rbind, .)
  }

  # remove zero cols
  mat             = mat[, colSums(mat) > 0]
  message("number of barcodes kept: ", ncol(mat))

  # save to h5 file
  write10xCounts(h5_f, mat, version = "3", overwrite = TRUE)

  # convert custom knees, shins and cellbender params to integers
  knee1           = as.integer(knee1)
  shin1           = as.integer(shin1)
  knee2           = as.integer(knee2)
  shin2           = as.integer(shin2)
  exp_cells       = as.integer(exp_cells)
  total_included  = as.integer(total_included)
  low_count_thr   = as.integer(low_count_thr)

  # check if low count threshold is defined
  if (is.na(low_count_thr)) {
    low_count_thr = 'shin2'
  }

  # estimate ambient(cellbender) parameters, write to csv
  bender_ps   = calc_ambient_params(split_mat = mat, run = run,
    knee1 = knee1, shin1 = shin1, knee2 = knee2, shin2 = shin2,
    run_var = run_var, low_count_threshold = low_count_thr, 
    expected_cells = exp_cells, total_included = total_included )

  # add spliced stats if not hto
  if (hto_mat == 0) {
    # get spliced / unspliced values
    splice_dt = data.table(
      barcode   = colnames(sce), 
      spliced   = colSums(assay(sce, "S")), 
      unspliced = colSums(assay(sce, "U"))
    )
    # add to bender_ps
    bender_ps = merge(bender_ps, splice_dt, by = "barcode") %>% .[ order(rank) ]
  }

  fwrite(bender_ps, file = knee_data_f)

  return(bender_ps)
}


# split_mat: matrix with split spliced, unspliced, ambiguous counts (could also be a normal matrix)
# this_run: sample_id
# min_umis_empty: library size of droplets that are definitely below the second knee and inflection
# min_umis_cells: minimum library size expected for cell containing droplets
# rank_empty_plateau: rank of any barcode within the empty droplet plateau
# low_count_threshold: 'shin2', 'knee2' or a specific library_size;
# low count threshold can be equal to second knee or second inflection or can be set manually

calc_ambient_params <- function(split_mat, run, min_umis_empty = 5, 
  min_umis_cells = NULL, rank_empty_plateau = NULL, low_count_threshold = 'shin2', 
  expected_cells = NA, total_included = NA, run_var= "sample_id",
  knee1 = NA, shin1 = NA, knee2 = NA, shin2 = NA) {
  # some checks on inputs
  if ( class(low_count_threshold) == 'character') {
    # check is the value of low_count_threshold is valid
    assert('low_count_threshold needs to be either "knee2", "shin2" or an integer',
      any(low_count_threshold %in% c('knee2', 'shin2')))
  }

  # get first knee
  knee1_ls  = .get_knee_and_shin_1(split_mat, min_umis_cells, knee1, shin1, knee2)

  # get second knee
  knee2_ls  = .get_knee_and_shin_2(split_mat, knee1_ls$ranks_dt, rank_empty_plateau, 
    min_umis_empty, knee1_ls$shin1_x, knee2, shin2)

  # get parameters
  params_ls = .get_params_ls(knee1_ls, knee2_ls, low_count_threshold = low_count_threshold,
    expected_cells = expected_cells, total_included = total_included)

  # return a dataframe with ranks and all parameters
  bender_ps = knee1_ls$ranks_dt %>%
    .[, (run_var) := run] %>%
    .[, `:=`(
      knee1                   = knee1_ls$sel_knee[ 'knee' ],
      shin1                   = knee1_ls$sel_knee[ 'shin' ],
      knee2                   = knee2_ls$sel_knee[ 'knee' ],
      shin2                   = knee2_ls$sel_knee[ 'shin' ],
      total_droplets_included = params_ls$total_included,
      low_count_threshold     = params_ls$lc,
      expected_cells          = params_ls$expected_cells
    )]

  # label cells in empty plateau
  bender_ps = .get_empty_plateau(
    knee_df        = bender_ps, 
    shin1          = knee1_ls$sel_knee[ 'shin' ], 
    total_included = params_ls$total_included, 
    knee2          = knee2_ls$sel_knee[ 'knee' ]
  )

  return(bender_ps)
}


.get_empty_plateau <- function(knee_df, shin1, total_included, knee2) {
  shin1_idx = which.min(abs(knee_df$total - shin1))[1]
  shin1_x   = knee_df[shin1_idx, rank] 

  empty_start = copy(knee_df)[, n := .I] %>%
    .[rank %between% c(shin1_x, total_included), n] %>%  
    log10() %>%
    mean() %>%
    (function(x) 10^x)() 

  empty_end = copy(knee_df)[total == knee2, unique(rank)]  

  knee_df[, in_empty_plateau := fifelse(rank %between% c(empty_start, empty_end), TRUE, FALSE)]

  return(knee_df)
}


.get_knee_and_shin_1 <- function(split_mat, min_umis_cells, knee1 = NA, shin1 = NA, knee2 = NA) {
  # check if custom knees and shins are defined
  if (all(sapply(c(knee1, shin1, knee2), function(p) !is.na(p)))) {
    # get knee
    low       = median(c(knee2, shin1))
    ranks_obj = barcodeRanks( split_mat, lower = low )
    sel_knee  = c(
      shin      = shin1,
      knee      = knee1
    )

  } else if (!is.null(min_umis_cells)) {
    # if min_umis_cells is specified use it as 'lower' parameter in barcodeRanks()
    # to find the first knee and inflection
    # min_umis_cells should be below expected first knee and inflection and ideally
    # above the second knee (on y axis)
    ranks_obj   = barcodeRanks( split_mat, lower = min_umis_cells )
    sel_knee    = c(
      shin    = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$inflection)))),
      knee    = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$knee))))
    )

  } else {
    # if min_umis_cells is not specified different 'lower' parameters are tested
    # and knee and inflection point with most votes are selected
    ranks_ls    = lapply(seq(1000, 100, by = -100), function(x) barcodeRanks(split_mat, lower = x) )

    # which parameters do these point to?
    ranks_ls_knees_and_shins = lapply(ranks_ls,
      function(x) paste0( as.integer(round(as.numeric(as.character(metadata(x)$inflection)))), 
        '_', as.integer(round(as.numeric(as.character(metadata(x)$knee)))) )
    ) %>% unlist()

    # pick the cutpoint with the largest number of "votes"
    votes_tbl   = ranks_ls_knees_and_shins %>% table()
    sel_cut     = names(votes_tbl)[ which.max(votes_tbl) ]

    # extract knee parameters
    sel_knee    = sel_cut %>%
      strsplit(., split ='_') %>% unlist() %>%
      as.integer() %>%
      setNames(c('shin', 'knee'))

    # get ranks object with selected knee and inflection
    ranks_obj = ranks_ls[[ which(ranks_ls_knees_and_shins == sel_cut)[1] ]]
  }

  # convert rankings object to data.frame
  ranks_dt = ranks_obj %>% as.data.frame() %>%
    as.data.table(keep.rownames = TRUE) %>%
    setnames("rn", "barcode") %>%
    .[order(rank)]

  # get x coordinates of selected inflection point (as.character is used because
  # as.integer() only returns the wrong number)
  shin1_idx   = which.min( abs(ranks_dt$total - sel_knee[1]) )[1]
  shin1_x     = ranks_dt[ shin1_idx, rank ]

  # put list of outputs together
  return(list(ranks_dt = ranks_dt, sel_knee = sel_knee, shin1_x = shin1_x))
}


.get_knee_and_shin_2 <- function(split_mat, ranks_dt, rank_empty_plateau,
  min_umis_empty, shin1_x, knee2, shin2) {
  # if rank_empty_plateau is specified use it to select barcodes for second call to barcodeRanks()
  # rank_empty_plateau should ideally be above expected second knee and inflection (on y axis)
  if (all(sapply(c(knee2, shin2), function(p) !is.na(p)))) {
    # find umi values in ranks_dt closest to predefined knee and shin

    shin2_idx = which.min( abs(ranks_dt$total - shin2) )[1]
    shin2_corr = ranks_dt[ shin2_idx, total]

    knee2_idx = which.min( abs(ranks_dt$total - knee2) )[1]
    knee2_corr = ranks_dt[ knee2_idx, total]

    sel_knee    = c(
      shin        = shin2_corr,
      knee        = knee2_corr
    )

  } else if (!is.null(rank_empty_plateau)) {
    # restrict to barcodes below (with higher ranks) specified threshold
    ranks_smol  = ranks_dt[rank > rank_empty_plateau, barcode]

    # use barcodeRanks to find knee
    ranks_obj   = barcodeRanks(split_mat[, ranks_smol], lower = min_umis_empty)
    sel_knee    = c(
      shin        = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$inflection)))),
      knee        = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$knee))))
    )
  } else {
    # if rank_empty_plateaus is not specified, filter barcodes based on multiple
    # different thresholds and select knee+inflection with most votes
    cuts        = .calc_small_knee_cuts_ls(ranks_dt, min_umis_empty, shin1_x)

    # rerun barcode ranks excluding everything above cuts
    ranks_ls    = lapply(cuts, function(this_cut) {
      # restrict to this cut
      ranks_smol  = ranks_dt[rank > this_cut, barcode]

      # run barcodeRanks, extract knee values
      ranks_obj   = barcodeRanks( split_mat[, ranks_smol], lower = min_umis_empty )
      shin2       = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$inflection))))
      knee2       = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$knee))))

      return( paste0(shin2, "_", knee2) )
    }) %>% unlist()

    # find the most consistent second knee and inflection
    shin_tbl    = str_before_first(ranks_ls, pattern = '_') %>% table()
    sel_i       = names(shin_tbl)[ which.max(shin_tbl) ]

    # different cuts often give identical inflection points but varying knees
    # from the knees corresponding to identical inflection points pick the one
    # closest to the median
    match_idx   = grepl(paste0('^', sel_i, '_'), ranks_ls)
    match_ks    = ranks_ls[ match_idx ] %>%
      str_after_last(., pattern = '_') %>%
      as.numeric()
    med_val     = median(match_ks)
    sel_k       = match_ks[ which.min(abs(match_ks - med_val))[[1]] ]

    # we have a knee!
    sel_knee    = c(
      shin        = as.integer(sel_i),
      knee        = as.integer(sel_k)
    )
  }

  # get rank corresponding to the second knee
  knee2_x = ranks_dt[total == sel_knee['knee'], rank] %>%
    unique()

  return(list(sel_knee = sel_knee, knee2_x = knee2_x))
}


.calc_small_knee_cuts_ls <- function(ranks_dt, min_umis_empty, shin_x) {
  # pick a 'total' value below which there are still enough data points to run
  # barcodeRanks(), barcodeRanks() need as least 3 unique 'total' (library size) values
  last        = tail(unique(ranks_dt[total > min_umis_empty, total]), n = 3)[1]
  last_x      = ranks_dt[total == last, rank] %>% .[1]

  # get what is in the middle of the last barcode and first inflection point on the log scale
  # we want to land approximatelly at the end of the empty_droplet_plateau
  middle = copy(ranks_dt) %>%
  .[, n:= 1:.N] %>%
  .[rank %between% c(shin_x, last_x), n] %>%
  log10() %>%
  mean() %>%
  10^.

  # pick 10 values (including infection one and middle value) to be used to filter
  # barcodes for second knee and inflection detection

  cuts      = 10^seq(log10(shin_x), log10(middle), length.out = 10)

  return( cuts )
}


.get_params_ls <- function(knee1_ls, knee2_ls, low_count_threshold, expected_cells = NA, 
  total_included = NA) {
  # unpack some things
  ranks_dt  = knee1_ls$ranks_dt
  shin1_x   = knee1_ls$shin1_x
  knee2_x   = knee2_ls$knee2_x

  if (is.na(expected_cells)) {
    # expected cells at first inflection point
    expected_cells  = shin1_x
  }

  if (is.na(total_included)) {
    # get total_droplets_included: halfway between 1st inflection point and
    # second knee on the log10 scale
    total_included = copy(ranks_dt) %>%
      .[, n:= 1:.N] %>%
      .[ rank %between% c(shin1_x, knee2_x), n ] %>%
      log10() %>% mean() %>% 10^. %>% round()
  }

  # get low count threshold
  if ( "character" %in% class(low_count_threshold) ) {
    if (low_count_threshold == 'knee2') {
      lc    = knee2_ls$sel_knee[ "knee" ]
    } else if (low_count_threshold == "shin2") {
      lc    = knee2_ls$sel_knee[ "shin" ]
    }
  } else {
    # if low_count_threshold is an integer, check that it's bigger than
    # total_droplets_included and expected_cells
    expected_x = ranks_dt[which.min(abs(rank - expected_cells)), total]
    total_x    = ranks_dt[which.min(abs(rank - total_included)), total]
    assert('low count threshold exceeds expected_cells and/or total_droplets_included',
           (low_count_threshold < expected_x) & (low_count_threshold < total_x))

    # it's ok, so we use it
    lc          = low_count_threshold
  }

  return(list(
    expected_cells  = expected_cells,
    total_included  = total_included,
    lc              = lc
  ))
}

# find slope at first inflection and total droplets included & expected_cells/total ratio
get_knee_params <- function(knee_f, sample_var) {
  ranks_df  = fread(knee_f)
  total_thr = unique(ranks_df$total_droplets_included) %>% log10()
  
  # get x coordinate of shin1
  shin1     = unique(ranks_df$shin1)
  shin1_row = which.min( abs(ranks_df$total - shin1) )[1]
  
  # get x coordinate of shin1
  shin1_x   = ranks_df[ shin1_row, rank ] %>% log10
  
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
  d1_shin  = d1$y[ which.min(abs(d1$x - shin1_x))[1] ]
  d1_total = d1$y[ which.min(abs(d1$x - total_thr))[1] ]
  
  keep_cols = c(sample_var, 'knee1', 'shin1', 'knee2', 'shin2', 'total_droplets_included', 'expected_cells')
  
  final = ranks_df %>%
    .[, ..keep_cols] %>%
    unique() %>%
    .[, `:=`(
      slope_shin1 = d1_shin,
      slope_total_included = d1_total
    )]%>% 
    .[, `:=`(
      slope_ratio = abs(slope_total_included) / abs(slope_shin1),
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
  
  knee_vars = c(sample_var, 'knee1', 'shin1', 'knee2', 'shin2',
    'total_droplets_included', 'expected_cells')
  
  lines_knees = ambient_knees_df %>% as.data.table %>% 
    .[ get(sample_var) %in% s_ord, ..knee_vars] %>%
    setnames( "total_droplets_included", "empty_plateau_middle" ) %>% 
    setnames( "expected_cells", "expected_cells" ) %>% 
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

plot_qc_metrics_split_by_cells_empties <- function(rna_knee_fs, 
  sample_var = "sample_id", min_umis = 10, n_cores) {
  
  # setup cluster
  bpparam = MulticoreParam(workers = n_cores, progressbar = FALSE)
  
  # get cells and empties
  plot_dt   = rna_knee_fs %>% bplapply(function(f) {
    tmp_dt = fread(f) %>% 
      .[rank <= expected_cells | in_empty_plateau == TRUE ] %>% 
      .[, `:=`(
        umis       = log10(total), 
        splice_pct = qlogis((spliced + 1) / (spliced + unspliced + 2)), 
        what       = fifelse(rank <= expected_cells, "cell", "empty")
      )] %>%
      .[, c(sample_var, 'barcode', 'umis', 'splice_pct', 'what'), with = FALSE]
    }, BPPARAM = bpparam) %>% rbindlist %>%
    melt(id.vars = c(sample_var, "barcode", "what"), measure.vars = c("umis", "splice_pct"),
      variable.name = "qc_metric", value.name = "qc_val"
    ) %>%
    .[, qc_metric := fcase(
      qc_metric == 'umis', 'no. of UMIs', 
      qc_metric == 'splice_pct', "spliced pct.", 
      default = NA_character_
    )]

  # breaks and labs
    umis_brks    = c(1e0, 1e1, 3e1, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6) %>% log10
    umis_labs    = c("1", "10", "30", "100", "300", "1k", "3k", "10k", "30k", "100k", "300k", "1M")
  
    splice_brks  = c(0.01, 0.03, 0.1, 0.3, 0.5, 0.7, 0.9, 0.97, 0.99) %>% qlogis
    splice_labs = c("1%", "3%", "10%", "30%", "50%", "70%", "90%", "97%", "99%")
    
  g_violin = ggplot() +
    geom_violin( data = plot_dt[ !is.na(qc_val) ],
                 aes( x = get(sample_var), y = qc_val, fill = what), colour = NA, 
                 kernel = 'rectangular', adjust = 0.1, scale = 'width', width = 0.8) +
    facet_grid( . ~ qc_metric, scales = 'free', space = 'free_y' ) +
    facetted_pos_scales(
      y = list(
        qc_metric == "no. of UMIs"     ~
          scale_y_continuous(breaks = umis_brks, labels = umis_labs),
        qc_metric == "spliced pct."    ~
          scale_y_continuous(breaks = splice_brks, labels = splice_labs)
      )
    ) +
    scale_fill_manual( values = c(cell = "#1965B0", empty = "grey") ) +
    coord_flip() +
    theme_classic() +
    labs( x = NULL, y = NULL, fill = "what does\nthe barcode\nrepresent?" ) +
    theme(panel.spacing = unit(1, "lines"))
  
  return(g_violin)
}
