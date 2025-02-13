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
})


# load counts data into sce object
save_alevin_h5_ambient_params <- function(sample, fry_dir, h5_f, cb_yaml_f, knee_data_f, demux_type,
                                          knee1, inf1, knee2, inf2, exp_cells, total_included, low_count_thr){
  # load the data, save to h5
  bender_ps = save_alevin_h5_knee_params_df(sample, fry_dir, h5_f, knee_data_f, hto_mat = 0, demux_type, 
                                            knee1, inf1, knee2, inf2, exp_cells, total_included, low_count_thr)

  # write these parameters to yaml file
  con_obj     = file(cb_yaml_f)
  writeLines(c(
    sprintf("sample: %s",                     sample),
    sprintf("total_droplets_included: %.0f",  unique(bender_ps$total_droplets_included) ),
    sprintf("expected_cells: %.0f",           unique(bender_ps$expected_cells) ),
    sprintf("low_count_threshold: %.0f",      unique(bender_ps$low_count_threshold) ),
    sprintf("knee_1: %.0f",                   unique(bender_ps$knee1)),
    sprintf("inflection_1: %.0f",             unique(bender_ps$inf1)),
    sprintf("knee_2: %.0f",                   unique(bender_ps$knee2)),
    sprintf("inflection_2: %.0f",             unique(bender_ps$inf2))
  ), con = con_obj)
  close(con_obj)
}


save_alevin_h5_knee_params_df <- function(sample, fry_dir, h5_f,  knee_data_f, hto_mat = 0, demux_type,
                                          knee1 = '', inf1 = '', knee2 = '', inf2 ='',
                                          exp_cells ='', total_included ='', low_count_thr =''){
 # load the data
  if(hto_mat){
    sce = loadFry(fry_dir)
    mat = counts(sce)
  }else{
     sce = loadFry(
      fry_dir,
      outputFormat = list(S = c("S"), U = c("U"), A = c("A"))
      )
      
  mat = assayNames(sce) %>% lapply(function(n) {
    mat       = assay(sce, n)
    rownames(mat) = paste0(rownames(mat), "_", n)
    return(mat)
  }) %>% do.call(rbind, .)
  }

  # remove zero cols
  mat = mat[, colSums(mat) > 0]
  message("number of barcodes kept: ", ncol(mat))

  # save to h5 file
  write10xCounts(h5_f, mat, version = "3", overwrite = TRUE)

  # convert custom knees, shins and cellbender params to integers
  knee1           = as.integer(knee1)
  inf1            = as.integer(inf1)
  knee2           = as.integer(knee2)
  inf2            = as.integer(inf2)
  exp_cells       = as.integer(exp_cells)
  total_included  = as.integer(total_included)
  low_count_thr   = as.integer(low_count_thr)

  # check if low count threshold is defined
  if(is.na(low_count_thr)){
    low_count_thr = 'inf2'
  }

  # estimate ambient(cellbender) parameters, write to csv
  bender_ps   = calc_ambient_params(
    split_mat = mat,
    sel_s = sample,
    knee1 = knee1,
    inf1 = inf1,
    knee2 = knee2,
    inf2 = inf2,
    multiplexing = demux_type,
    low_count_threshold = low_count_thr,
    expected_cells = exp_cells,
    total_included = total_included
  )

  fwrite(bender_ps, file = knee_data_f)

  return(bender_ps)

}


# split_mat: matrix with split spliced, unspliced, ambiguous counts (could also be a normal matrix)
# sample: sample_id
# min_umis_empty: library size of droplets that are definitely below the second knee and inflection
# min_umis_cells: minimum library size expected for cell containing droplets
# rank_empty_plateau: rank of any barcode within the empty droplet plateau
# low_count_threshold: 'inf2', 'knee2' or a specific library_size;
# low count threshold can be equal to second knee or second inflection or can be set manually

calc_ambient_params <- function(split_mat, sel_s, min_umis_empty = 5, min_umis_cells = NULL,
                                   rank_empty_plateau = NULL, low_count_threshold = 'inf2',
                                   expected_cells = NA, total_included =NA, multiplexing = "",
                                   knee1 = NA, inf1 = NA, knee2 = NA, inf2 = NA) {
  # some checks on inputs
  if ( class(low_count_threshold) == 'character') {
    # check is the value of low_count_threshold is valid
    assert('low_count_threshold needs to be either "knee2", "inf2" or an integer',
           any(low_count_threshold %in% c('knee2', 'inf2')))
  }

  # get first knee
  knee1_ls  = .get_knee_and_inf_1(
    split_mat,
    min_umis_cells,
    knee1,
    inf1,
    knee2
    )


  # get second knee
  knee2_ls  = .get_knee_and_inf_2(
    split_mat,
    knee1_ls$ranks_dt,
    rank_empty_plateau,
    min_umis_empty,
    knee1_ls$infl1_x,
    knee2,
    inf2
    )

  # get parameters
  params_ls = .get_params_ls(
    knee1_ls,
    knee2_ls,
    low_count_threshold =low_count_threshold,
    expected_cells = expected_cells,
    total_included = total_included
    )

  # return a dataframe with ranks and all parameters
  if(multiplexing != ""){
    run = 'pool_id'
  }else{
    run = 'sample_id'
  }

  bender_ps = knee1_ls$ranks_dt %>%
    .[, (run):= sel_s] %>%
    .[, `:=`(
      knee1                   = knee1_ls$sel_knee[ 'knee' ],
      inf1                    = knee1_ls$sel_knee[ 'inflection' ],
      knee2                   = knee2_ls$sel_knee[ 'knee' ],
      inf2                    = knee2_ls$sel_knee[ 'inflection' ],
      total_droplets_included = params_ls$total_included,
      low_count_threshold     = params_ls$lc,
      expected_cells          = params_ls$expected_cells
    )]


  return(bender_ps)
}

.get_knee_and_inf_1 <- function(split_mat, min_umis_cells, knee1 = NA, inf1 = NA, knee2 = NA) {
  # check if custom knees and shins are defined
  if(all(sapply(c(knee1, inf1, knee2), function(p) !is.na(p)))){
      low = median(c(knee2, inf1))

      ranks_obj = barcodeRanks( split_mat, lower = low )
      sel_knee    = c(
      inflection  = inf1,
      knee        = knee1
    )

  }else if(!is.null(min_umis_cells)) {
    # if min_umis_cells is specified use it as 'lower' parameter in barcodeRanks()
    # to find the first knee and inflection
    # min_umis_cells should be below expected first knee and inflection and ideally
    # above the second knee (on y axis)
    ranks_obj   = barcodeRanks( split_mat, lower = min_umis_cells )
    sel_knee    = c(
      inflection  = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$inflection)))),
      knee        = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$knee))))
    )

  } else {
    # if min_umis_cells is not specified different 'lower' parameters are tested
    # and knee and inflection point with most votes are selected
    ranks_ls    = lapply(seq(1000, 100, by = -100),
                         function(x) barcodeRanks(split_mat, lower = x) )

    # which parameters do these point to?
    ranks_ls_knees_and_inf   = lapply(ranks_ls,
                                      function(x) paste0(as.integer(round(as.numeric(as.character(metadata(x)$inflection)))), '_',
                                                         as.integer(round(as.numeric(as.character(metadata(x)$knee)))))
    ) %>% unlist()

    # pick the cutpoint with the largest number of "votes"
    votes_tbl   = ranks_ls_knees_and_inf %>% table()
    sel_cut     = names(votes_tbl)[ which.max(votes_tbl) ]

    # extract knee parameters
    sel_knee    = sel_cut %>%
      strsplit(., split ='_') %>% unlist() %>%
      as.integer() %>%
      setNames(c('inflection', 'knee'))

    # get ranks object with selected knee and inflection
    ranks_obj = ranks_ls[[ which(ranks_ls_knees_and_inf == sel_cut)[1] ]]
  }

  # convert rankings object to data.frame
  ranks_dt = ranks_obj %>% as.data.frame() %>%
    as.data.table(keep.rownames = TRUE) %>%
    setnames("rn", "barcode") %>%
    .[order(rank)]

  # get x coordinates of selected inflection point (as.character is used because
  # as.integer() only returns the wrong number)
  infl1_idx = which.min( abs(ranks_dt$total - sel_knee[1]) )[1]
  infl1_x   = ranks_dt[ infl1_idx, rank ]

  # put list of outputs together
  return(list(ranks_dt = ranks_dt, sel_knee = sel_knee, infl1_x = infl1_x))
}

.get_knee_and_inf_2 <- function(split_mat, ranks_dt, rank_empty_plateau,
                                min_umis_empty, infl1_x,
                                knee2, inf2) {

  # if rank_empty_plateau is specified use it to select barcodes for second call to barcodeRanks()
  # rank_empty_plateau should ideally be above expected second knee and inflection (on y axis)
  if(all(sapply(c(knee2, inf2), function(p) !is.na(p)))){
    # find umi values in ranks_dt closest to predefined knee and shin

    inf2_idx = which.min( abs(ranks_dt$total - inf2) )[1]
    inf2_corr = ranks_dt[ inf2_idx, total]

    knee2_idx = which.min( abs(ranks_dt$total - knee2) )[1]
    knee2_corr = ranks_dt[ knee2_idx, total]

      sel_knee    = c(
      inflection  = inf2_corr,
      knee        = knee2_corr
    )

  }else if (!is.null(rank_empty_plateau)) {
    # restrict to barcodes below (with higher ranks) specified threshold
    ranks_smol  = ranks_dt[rank > rank_empty_plateau, barcode]


    # use barcodeRanks to find knee
    ranks_obj   = barcodeRanks(split_mat[, ranks_smol], lower = min_umis_empty)
    sel_knee    = c(
      inflection  = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$inflection)))),
      knee        = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$knee))))
    )
  } else {
    # if rank_empty_plateaus is not specified, filter barcodes based on multiple
    # different thresholds and select knee+inflection with most votes
    cuts        = .calc_small_knee_cuts_ls(ranks_dt, min_umis_empty, infl1_x)

    # rerun barcode ranks excluding everything above cuts
    ranks_ls    = lapply(cuts, function(this_cut) {
      # restrict to this cut
      ranks_smol  = ranks_dt[rank > this_cut, barcode]

      # run barcodeRanks, extract knee values
      ranks_obj   = barcodeRanks( split_mat[, ranks_smol], lower = min_umis_empty )
      inf_2       = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$inflection))))
      knee_2      = as.integer(round(as.numeric(as.character(metadata(ranks_obj)$knee))))

      return( paste0(inf_2, "_", knee_2) )
    }) %>% unlist()

    # find the most consistent second knee and inflection
    infl_tbl    = str_before_first(ranks_ls, pattern = '_') %>% table()
    sel_i       = names(infl_tbl)[ which.max(infl_tbl) ]

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
      inflection  = as.integer(sel_i),
      knee        = as.integer(sel_k)
    )
  }

  # get rank corresponding to the second knee
  knee2_x = ranks_dt[total == sel_knee['knee'], rank] %>%
    unique()

  return(list(sel_knee = sel_knee, knee2_x = knee2_x))
}



.calc_small_knee_cuts_ls <- function(ranks_dt, min_umis_empty, infl1_x) {
  # pick a 'total' value below which there are still enough data points to run
  # barcodeRanks(), barcodeRanks() need as least 3 unique 'total' (library size) values
  last        = tail(unique(ranks_dt[total > min_umis_empty, total]), n = 3)[1]
  last_x      = ranks_dt[total == last, rank] %>% .[1]

  # get what is in the middle of the last barcode and first inflection point on the log scale
  # we want to land approximatelly at the end of the empty_droplet_plateau
  middle = copy(ranks_dt) %>%
  .[, n:= 1:.N] %>%
  .[rank %between% c(infl1_x, last_x), n] %>%
  log10() %>%
  mean() %>%
  10^.

  # pick 10 values (including infection one and middle value) to be used to filter
  # barcodes for second knee and inflection detection

  cuts      = 10^seq(log10(infl1_x), log10(middle), length.out = 10)

  return( cuts )
}

.get_params_ls <- function(knee1_ls, knee2_ls,
                          low_count_threshold, expected_cells = NA, total_included = NA) {
  # unpack some things
  ranks_dt  = knee1_ls$ranks_dt
  infl1_x   = knee1_ls$infl1_x
  knee2_x   = knee2_ls$knee2_x

  if(is.na(expected_cells)){
  # expected cells at first inflection point
  expected_cells  = infl1_x
  }

  if(is.na(total_included)){
  # get total_droplets_included: halfway between 1st inflection point and
  # second knee on the log10 scale
  total_included = copy(ranks_dt) %>%
  .[, n:= 1:.N] %>%
  .[rank %between% c(infl1_x, knee2_x), n] %>%
  log10() %>%
  mean() %>%
  10^. %>%
  round()

  }

  # get low count threshold
  if ( "character" %in% class(low_count_threshold) ) {
    if (low_count_threshold == 'knee2') {
      lc    = knee2_ls$sel_knee[ "knee" ]
    } else if (low_count_threshold == "inf2") {
      lc    = knee2_ls$sel_knee[ "inflection" ]
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



