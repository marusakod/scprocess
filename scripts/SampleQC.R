suppressPackageStartupMessages({
  library('magrittr')
  library('data.table')
  library('forcats')
  library('assertthat')
  library('stringr')

  library('ggplot2')
  library('scales')
  library('patchwork')
  library('ComplexHeatmap')
  library('ggbeeswarm')
  library('ggh4x')

  library('Matrix')
  library('SingleCellExperiment')
  # library('SampleQC')

  library('BiocParallel')
})

# define some breaks
log_brks    = c(1e1, 2e1, 5e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4) %>%
  log10
log_labs    = c("10", "20", "50", "100", "200", "500",
  "1k", "2k", "5k", "10k", "20k", "50k")
logit_brks  = c(1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.10, 0.30,
  0.50, 0.70, 0.90, 0.97, 0.99, 0.999) %>% qlogis
logit_labs  = c("0.01%", "0.03%", "0.1%", "0.3%", "1%", "3%", "10%", "30%",
  "50%", "70%", "90%", "97%", "99%", "99.9%")
splice_brks = c(1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.10, 0.50,
  0.90, 0.97, 0.99, 0.999) %>% qlogis
splice_labs = c("0.01%", "0.03%", "0.1%", "0.3%", "1%", "3%", "10%", "50%",
  "90%", "97%", "99%", "99.9%")
prob_brks   = c(0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999, 0.999999) %>% qlogis
prob_labs   = c("50%", "90%", "99%", "99.9%", "99.99%", "99.999%", "99.9999%")


main_qc <- function(sce_f, dbl_f, dbl_smpl_var, 
  hard_min_counts, hard_min_feats, hard_max_mito,
  min_counts, min_feats, min_mito, max_mito, min_splice, max_splice, min_cells, filter_bender, amb_method,
  qc_f, keep_f) {
  # check inputs
  filter_bender = as.logical(filter_bender)
  # load doublet info
  dbl_dt        = dbl_f %>% fread %>%
    .[, c('cell_id', dbl_smpl_var, 'dbl_class'), with = FALSE]
  sngl_ids      = dbl_dt[ dbl_class == 'singlet' ]$cell_id

  if(amb_method == 'cellbender'){
    annot_vars = c('bender_n_ok', 'bender_n_used',
      'bender_prop_ok', 'bender_logit_ok')
  }else{
    annot_vars = NULL
  }
  # load sce file, get QC data
  sce         = readRDS(sce_f) %>% .[, dbl_dt$cell_id ]
  qc_all      = make_qc_dt(colData(sce),
    sample_var  = 'sample_id',
    qc_names    = c('log_counts', 'log_feats', 'logit_mito', 'logit_spliced'),
    annot_vars  = annot_vars
    )

  # some checks
  assert_that( all(sort(qc_all$cell_id) == sort(dbl_dt$cell_id)))
  assert_that( all( !is.na(qc_all$logit_splice) ) )

  # restrict to singlets
  qc_all        = qc_all[ cell_id %in% sngl_ids ]


  if(amb_method == 'cellbender'){
  # processing / calculations
  bender_chks   = calc_bender_chks(qc_all)
  qc_all        = qc_all %>%
    merge(bender_chks[, .(sample_id, bender_ok)], by = "sample_id")
  }

  # restrict to samples with enough cells
  qc_dt         = qc_all %>%
    .[ log_counts >= log10(hard_min_counts) ] %>%
    .[ log_feats  >= log10(hard_min_feats)  ] %>%
    .[ logit_mito <  qlogis(hard_max_mito)  ]

  # filter on splicing
  keep_dt       = qc_dt %>%
    .[ log_counts    >= log10(min_counts)  ] %>%
    .[ log_feats     >= log10(min_feats)   ] %>%
    .[ logit_mito    >  qlogis(min_mito)   ] %>%
    .[ logit_mito    <  qlogis(max_mito)   ] %>%
    .[ logit_spliced >  qlogis(min_splice) ] %>%
    .[ logit_spliced <  qlogis(max_splice) ]

  if(amb_method == 'cellbender'){
  if (filter_bender)
    keep_dt       = keep_dt[ bender_ok == TRUE ]
  }


  # exclude samples with v few cells
  n_dt          = keep_dt[, .N, by = sample_id]
  keep_s        = n_dt[ N >= min_cells ]$sample_id
  keep_dt       = keep_dt[ sample_id %in% keep_s ]
  keep_ids      = keep_dt$cell_id
  assert_that( length(keep_ids) > 0 )

  # check that that worked
  assert_that( all(keep_dt[, .N, by = sample_id]$N >= min_cells) )
  assert_that( all(dbl_dt[ cell_id %in% keep_ids ]$dbl_class == "singlet") )

  # record which kept
  qc_all    = qc_all %>%
    .[, keep_hard := cell_id %in% qc_dt$cell_id ] %>%
    .[, keep      := cell_id %in% keep_ids ]

  # save outputs
  fwrite(qc_all, file = qc_f)
  fwrite(keep_dt[, .(sample_id, cell_id)], file = keep_f)
}

calc_bender_chks <- function(qc_all, mad_cutoff = 3) {
  logits_dt   = qc_all[, .(sample_id, bender_prop_ok, bender_logit_ok)] %>% unique
  bender_chks = logits_dt %>%
    .[, med_logit := median(bender_logit_ok) ] %>%
    .[, mad_logit := mad(bender_logit_ok) ] %>%
    .[, bender_ok := abs(bender_logit_ok - med_logit) < (mad_cutoff * mad_logit) ]

  return(bender_chks)
}

plot_qc_ranges_marginals <- function(qc_input, s_lvls, qc_names, qc_lu, thrshlds_dt, amb_method) {
  # melt, add names
  qc_melt   = copy(qc_input) %>%
    melt(measure = qc_names, val = 'qc_val', var = 'qc_var') %>%
    .[, qc_full   := qc_lu[ qc_var ] ] %>%
    .[, qc_var    := factor(qc_var, levels = qc_names) ] %>%
    .[, qc_full   := fct_reorder(qc_full, as.integer(qc_var)) ]
  hlines_dt = thrshlds_dt %>% copy %>%
    .[, qc_full   := qc_lu[ qc_var ] ] %>%
    .[, qc_var    := factor(qc_var, levels = qc_names) ] %>%
    .[, qc_full   := fct_reorder(qc_full, as.integer(qc_var)) ]

  # calculate medians etc
if(amb_method == 'cellbender'){
  qc_meds = qc_melt %>%
  .[, .(
    log10_N = log10(.N),
    bender_logit_ok =  unique(bender_logit_ok),
    q50 = median(qc_val, na.rm = TRUE),
    q10 = quantile(qc_val, 0.1, na.rm = TRUE),
    q90 = quantile(qc_val, 0.9, na.rm = TRUE),
    q025 = quantile(qc_val, 0.025, na.rm = TRUE),
    q975 = quantile(qc_val, 0.975, na.rm = TRUE)
  ),
  by = c('sample_id', 'qc_var', 'qc_full')]
  }else{
    qc_meds = qc_melt %>%
  .[, .(
    log10_N = log10(.N),
    q50 = median(qc_val, na.rm = TRUE),
    q10 = quantile(qc_val, 0.1, na.rm = TRUE),
    q90 = quantile(qc_val, 0.9, na.rm = TRUE),
    q025 = quantile(qc_val, 0.025, na.rm = TRUE),
    q975 = quantile(qc_val, 0.975, na.rm = TRUE)
  ),
  by = c('sample_id', 'qc_var', 'qc_full')]
  }


  # bar width
  bar_w     = 0.4
  n_dt      = qc_meds[, .(sample_id, `no. cells` = log10_N) ] %>% unique %>%
    melt.data.table( id = "sample_id", var = "var", val = "value")

  # put in nice order
  n_dt      = n_dt %>%
    .[, sample_id := factor(sample_id, levels = rev(s_lvls)) ]
  qc_melt   = qc_melt %>%
    .[, sample_id := factor(sample_id, levels = rev(s_lvls)) ]
  qc_meds   = qc_meds %>%
    .[, sample_id := factor(sample_id, levels = rev(s_lvls)) ]

  # make plot
  g_n = ggplot( n_dt ) +
    geom_point( aes( y = value, x = as.integer(sample_id) ),
      size = 4, shape = 21, fill = 'grey60') +
    scale_x_continuous( breaks = seq.int(length(s_lvls)), labels = levels(n_dt$sample_id) ) +
    facet_grid( . ~ var, scales = 'free', space = 'free_y' ) +
    scale_y_continuous(breaks = log_brks, labels = log_labs) +
    expand_limits( y = log10(c(1e3, 1e4)) ) +
    coord_flip( xlim = c(0.5, length(s_lvls) + 0.5), expand = FALSE ) +
    theme_classic() +
    theme(
      axis.text.x       = element_text( angle = 90, hjust = 1, vjust = 0.5 ),
      strip.text.y      = element_blank()
      ) +
    labs( x = NULL, y = 'sample' )

  g_violin = ggplot() +
    geom_violin( data = qc_melt[ !is.na(qc_val) ],
      aes( x = sample_id, y = qc_val ), colour = NA, fill = 'grey60',
      kernel = 'rectangular', adjust = 0.1, scale = 'width') +
    geom_hline( data = hlines_dt, aes( yintercept = cut_point ),
      colour = 'black', linetype = 'dashed', size = 0.5, alpha = 0.5 ) +
    facet_grid( . ~ qc_full, scales = 'free', space = 'free_y' ) +
    scale_x_discrete( breaks = levels(qc_meds$sample_id), drop = FALSE ) +
    facetted_pos_scales(
      y = list(
        qc_full == "library size"     ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full == "no. of features"  ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full == "mito pct."         ~
          scale_y_continuous(breaks = logit_brks, labels = logit_labs),
        qc_full == "spliced pct."      ~
          scale_y_continuous(breaks = splice_brks, labels = splice_labs)
        )
      ) +
    coord_flip() +
    theme_classic() +
    theme(
      axis.ticks.y      = element_blank(),
      axis.text.y       = element_blank(),
      axis.text.x       = element_text( angle = 90, hjust = 1, vjust = 0.5 )
      ) +
    labs( x = NULL, y = NULL )

  g = g_n + g_violin + plot_layout(widths = c(1, 5))

  return(g)
}

plot_qc_ranges_pairwise <- function(qc_input, qc_names, qc_lu, thrshlds_dt, amb_method) {
  # add logit ok to names
  if(amb_method == 'cellbender'){
  qc_names  = c(qc_names, 'bender_logit_ok')
  qc_lu     = c(qc_lu, bender_logit_ok = "bender cell pct.")
  }

  # calc medians etc
  qc_meds   = qc_input %>%
    melt(measure = qc_names, val = 'qc_val', var = 'qc_var') %>%
    .[, .(
      log10_N   = log10(.N),
      q50       = median(qc_val, na.rm = TRUE),
      q10       = quantile(qc_val, 0.1, na.rm = TRUE),
      q90       = quantile(qc_val, 0.9, na.rm = TRUE),
      q025      = quantile(qc_val, 0.025, na.rm = TRUE),
      q975      = quantile(qc_val, 0.975, na.rm = TRUE)
      ),
      by = c('sample_id', 'qc_var')] %>%
    .[, qc_full := qc_lu[ qc_var ] ]

    if(amb_method == 'cellbender'){
    qc_meds = qc_meds %>%
    .[ qc_var == "bender_logit_ok", q50 := pmin(q50, qlogis(0.999)) ]
    }

  # make pairwise plot
  pairs_dt  = merge(qc_meds, qc_meds,
    by = c('sample_id', 'log10_N'), allow.cartesian = TRUE) %>%
    .[, qc_var.x  := factor(qc_var.x, levels = qc_names) ] %>%
    .[, qc_full.x := fct_reorder(qc_full.x, as.integer(qc_var.x)) ] %>%
    .[, qc_var.y  := factor(qc_var.y, levels = qc_names) ] %>%
    .[, qc_full.y := fct_reorder(qc_full.y, as.integer(qc_var.y)) ] %>%
    .[ as.integer(qc_var.x) > as.integer(qc_var.y) ]

  scales_x_ls = list(
        qc_full.x == "library size"    ~
          scale_x_continuous(breaks = log_brks, labels = log_labs),
        qc_full.x == "no. of features" ~
          scale_x_continuous(breaks = log_brks, labels = log_labs),
        qc_full.x == "mito pct."        ~
          scale_x_continuous(breaks = logit_brks, labels = logit_labs),
        qc_full.x == "spliced pct."     ~
          scale_x_continuous(breaks = splice_brks, labels = splice_labs)
  )

  scales_y_ls = list(
        qc_full.y == "library size"    ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full.y == "no. of features" ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full.y == "mito pct."        ~
          scale_y_continuous(breaks = logit_brks, labels = logit_labs),
        qc_full.y == "spliced pct."     ~
          scale_y_continuous(breaks = splice_brks, labels = splice_labs)
  )

  if(amb_method == 'cellbender'){
    # add bender cell pct. to scales
    scales_x_ls = c(scales_x_ls,
    qc_full.x == "bender cell pct." ~
          scale_x_continuous(breaks = logit_brks, labels = logit_labs)
          )

    scales_y_ls = c(scales_y_ls,
    qc_full.y == "bender cell pct." ~
          scale_y_continuous(breaks = logit_brks, labels = logit_labs)
    )

  }


  # make plot
  g = ggplot(pairs_dt) +

    geom_point(
      aes( x = q50.x, y = q50.y, size = log10_N ),
      colour = 'black', shape = 21, alpha = 0.8
      ) +
    scale_x_continuous( breaks = pretty_breaks() ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    scale_size( range = c(1, 4), breaks = log_brks, labels = log_labs ) +
    guides(fill = guide_legend(override.aes = list(size = 3, shape = 21) ) ) +
    facet_grid( qc_full.y ~ qc_full.x, scales = 'free' ) +
    facetted_pos_scales(
      x = scales_x_ls,
      y = scales_y_ls
      ) +
    theme_bw() +
    theme(
      axis.text.x       = element_text( angle = 90, hjust = 1, vjust = 0.5 ),
      panel.grid        = element_blank(),
      strip.background  = element_rect( fill = 'white')
      ) +
    labs( x = NULL, y = NULL, size = 'no. cells\nin sample' )

  return(g)
}

make_qc_dt_file <- function(sce_f, qc_f, overwrite=FALSE) {
  # check if already done
  if (file.exists(qc_f) & overwrite == FALSE)
    return(fread(qc_f))

  # load sce file
  sce         = sce_f %>% readRDS

  # calculate reads
  all_counts  = sce %>% counts %>% Matrix::colSums(.)

  # calculate mito reads
  mt_idx      = rowData(sce)$symbol %>%
    str_detect('^MT-')
  assert_that( sum(mt_idx) == 13 )
  mt_counts   = sce[mt_idx, ] %>% counts %>% Matrix::colSums(.)

  # calculate no. features
  all_feats   = sce %>% counts %>% `>`(0) %>% Matrix::colSums(.)

  # assemble
  qc_dt = data.table(
    cell_id     = colnames(sce),
    all_counts  = all_counts,
    all_feats   = all_feats,
    mito_counts = mt_counts
    )
  fwrite(qc_dt, file=qc_f)

  return(qc_dt)
}

get_cols_dt <- function(sce_f, cols_f, overwrite=FALSE) {
  # check if already done
  if (file.exists(cols_f) & overwrite == FALSE)
    return(fread(cols_f))

  # load sce file
  sce     = sce_f %>% readRDS

  # assemble
  cols_dt   = colData(sce) %>% as.data.frame %>%
    as.data.table(keep.rownames = 'cell_id')
  fwrite(cols_dt, file=cols_f)

  return(cols_dt)
}

plot_qc_metric_scatter <- function(dt, qc_names, qc_lu, thrshlds_dt) {
  melt_dt = dt %>%
    melt(measure=qc_names, value.name='qc_val', variable.name='qc_var') %>%
    .[, qc_full := qc_lu[ qc_var ] ] %>%
    .[, qc_var  := factor(qc_var, levels = qc_names) ] %>%
    .[, qc_full := fct_reorder(qc_full, as.integer(qc_var)) ]

  # make some lines
  hlines_dt = thrshlds_dt[, .(qc_var.y = qc_var, cut.y = cut_point)] %>%
    .[, qc_y      := qc_lu[ qc_var.y ] ] %>%
    .[, qc_var.y  := factor(qc_var.y, levels = qc_names) ] %>%
    .[, qc_y      := fct_reorder(qc_y, as.integer(qc_var.y)) ] %>%
    .[, dummy     := "dummy" ]
  vlines_dt = thrshlds_dt[, .(qc_var.x = qc_var, cut.x = cut_point)] %>%
    .[, qc_x      := qc_lu[ qc_var.x ] ] %>%
    .[, qc_var.x  := factor(qc_var.x, levels = qc_names) ] %>%
    .[, qc_x      := fct_reorder(qc_x, as.integer(qc_var.x)) ] %>%
    .[, dummy     := "dummy" ]

  # remove duplications
  lines_dt  = merge(hlines_dt, vlines_dt, by = "dummy", allow.cartesian = TRUE) %>%
    .[ as.integer(qc_x) > as.integer(qc_y) ]

  # get what to plot
  plot_dt = merge(
    melt_dt[, .(cell_id, qc_x = qc_full, val_x = qc_val)],
    melt_dt[, .(cell_id, qc_y = qc_full, val_y = qc_val)],
    by = 'cell_id', allow.cartesian = TRUE
    ) %>% .[ as.integer(qc_x) > as.integer(qc_y) ]

  # get sample name
  sel_s   = unique(dt$sample_id)
  g = ggplot(plot_dt) + aes( x=val_x, y=val_y ) +
    geom_bin2d() + scale_fill_distiller( palette='RdBu', trans='log10' ) +
    geom_hline( data = lines_dt, aes( yintercept = cut.y ),
      color = 'black', linetype = "dashed", size = 0.5, alpha = 0.5 ) +
    geom_vline( data = lines_dt, aes( xintercept = cut.x ),
      color = 'black', linetype = "dashed", size = 0.5, alpha = 0.5 ) +
    facet_grid( qc_y ~ qc_x, scales='free' ) +
    theme_classic() +
    theme(
      panel.grid        = element_blank(),
      strip.background  = element_rect(fill = 'white')
    ) +
    facetted_pos_scales(
      x = list(
        qc_x == "library size"    ~
          scale_x_continuous(breaks = log_brks, labels = log_labs),
        qc_x == "no. of features" ~
          scale_x_continuous(breaks = log_brks, labels = log_labs),
        qc_x == "mito pct."        ~
          scale_x_continuous(breaks = logit_brks, labels = logit_labs),
        qc_x == "spliced pct."     ~
          scale_x_continuous(breaks = splice_brks, labels = splice_labs)
        ),
      y = list(
        qc_y == "library size"    ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_y == "no. of features" ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_y == "mito pct."        ~
          scale_y_continuous(breaks = logit_brks, labels = logit_labs),
        qc_y == "spliced pct."     ~
          scale_y_continuous(breaks = splice_brks, labels = splice_labs)
        )
      ) +
    labs(
      x     = 'QC metric 1',
      y     = 'QC metric 2',
      fill  = 'no. cells',
      title   = sel_s
      )

  return(g)
}

plot_totals_split_by_meta <- function(pre_dt, post_dt, meta_dt) {
  # define age splits
  age_breaks  = c(0, 40, 50, 60, 70, 80, 100)
  age_labels  = paste0('<=', age_breaks[-1])
  yrs_breaks  = c(0, 10, 20, 30, 40, 60)
  yrs_labels  = paste0(yrs_breaks[-length(yrs_breaks)], ' to ', yrs_breaks[-1])

  # load metadata
  meta_dt     = copy(meta_dt) %>%
    .[, lesion_type := lesion_type %>%
      fct_recode(`WM (ctrl)` = "WM", `GM (ctrl)` = "GM")] %>%
    .[, age_cat     := cut(age_at_death, breaks = age_breaks,
      labels = age_labels) %>% factor(levels = age_labels) ] %>%
    .[, yrs_w_ms    := cut(years_w_ms, breaks = yrs_breaks,
      labels = yrs_labels) %>% factor(levels = yrs_labels) ] %>%
    .[, .(sample_id, subject_id, matter, lesion_type, diagnosis,
      sex, age_cat, yrs_w_ms, pmi_cat, brain_bank = sample_source, seq_pool)]

  # join to keep totals
  pre_n_dt    = pre_dt[, .(pre_n = .N), by = sample_id]
  post_n_dt   = post_dt[, .(post_n = .N), by = sample_id]
  keep_n_dt   = merge(pre_n_dt, post_n_dt, by = 'sample_id', all.x = TRUE) %>%
    .[ is.na(post_n), post_n := 0 ] %>%
    merge(meta_dt, by = 'sample_id')

  # pick which to show, get their values
  meta_vars   = c('matter', 'lesion_type', 'diagnosis', 'sex', 'age_cat',
    'yrs_w_ms', 'pmi_cat', 'brain_bank', 'seq_pool')
  level_ord   = meta_vars %>%
    lapply(function(v) levels(factor(keep_n_dt[[v]]))) %>%
    do.call(c, .)

  # melt by variable
  tmp_dt      = keep_n_dt %>%
    melt(id = c('sample_id', 'subject_id', 'pre_n', 'post_n'),
      measure = meta_vars, variable.name = 'meta_var', value.name = 'meta_val') %>%
    .[, meta_val := factor(meta_val, levels=level_ord)] %>%
    melt(measure = c('pre_n', 'post_n'),
      value.name = 'n_cells', variable.name = 'qc_status') %>%
    .[n_cells > 0, .(
      n_cells   = round( sum(n_cells) / 1e3 ) %>% as.integer,
      n_samples = .N,
      n_donors  = length(unique(subject_id))
      ), by = .(meta_var, meta_val, qc_status)] %>%
    melt(measure = c('n_cells', 'n_samples', 'n_donors'),
      value.name = 'count_val', variable.name = 'count_var')

  # calculate lost cells etc
  lost_dt     = merge(
    tmp_dt[qc_status == 'pre_n', .(meta_var, meta_val, count_var, pre_val=count_val)],
    tmp_dt[qc_status == 'post_n', .(meta_var, meta_val, count_var, post_val=count_val)],
    by = c('meta_var', 'meta_val', 'count_var')
    ) %>%
    .[, lost_val := pre_val - post_val]
  all_dt      = rbind(
    tmp_dt[qc_status == 'post_n', .(meta_var, meta_val, qc_status = 'kept', count_var, count_val)],
    lost_dt[, .(meta_var, meta_val, qc_status = 'excluded', count_var, count_val=lost_val)]
    ) %>% .[, qc_status := qc_status %>% factor(levels=c('excluded', 'kept')) ] %>%
    .[, count_var := count_var %>% fct_recode(`#k cells` = 'n_cells',
      `# samples` = 'n_samples', `# donors` = 'n_donors') ]

  # plot
  g = ggplot(all_dt[count_val > 0]) +
    aes( y = fct_rev(meta_val), x = count_val,
      label = count_val, fill = qc_status ) +
    geom_col(colour = NA) +
    geom_text(colour = 'black', position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = c(kept = 'grey50', excluded = 'grey80'),
      breaks = c('kept', 'excluded')) +
    scale_x_continuous(breaks = pretty_breaks()) +
    facet_grid( meta_var ~ count_var, scales='free', space='free_y') +
    theme_classic() +
    labs(y = NULL, x = NULL, fill = 'QC status',
      title = 'Summary of cells, samples and donors retained')

  return(g)
}



plot_qc_summary_heatmap <- function(qc_stats, meta_input) {
  # make matrix of z-scores
  stats_tmp = copy(qc_stats) %>%
    .[ sample_id %in% meta_input$sample_id ] %>%
    .[, z := scale(med_val), by = qc_var ] %>%
    .[ qc_var %in% c('logit_mito', 'logit_splice'), z := z * -1 ]

  # make matrix of z-scores
  z_wide    = stats_tmp %>%
    dcast.data.table( sample_id ~ qc_var, value.var = 'z')
  z_mat     = z_wide[, -'sample_id', with = FALSE] %>%
    as.matrix %>% set_rownames(z_wide$sample_id)

  # make matrix of text labels
  lab_wide  = stats_tmp %>%
    .[ qc_var == "log10_N",
      lab := sprintf("%0.1fk", 10^med_val / 1e3) ] %>%
    .[ qc_var == "log_counts",
      lab := sprintf("%0.1fk", 10^med_val / 1e3) ] %>%
    .[ qc_var == "log_feats",
      lab := sprintf("%0.1fk", 10^med_val / 1e3) ] %>%
    .[ qc_var == "logit_mito",
      lab := sprintf("%0.1f%%", plogis(med_val) * 1e2) ] %>%
    .[ qc_var == "logit_splice",
      lab := sprintf("%0.0f%%", 1/(1 + 2^med_val) * 1e2) ] %>%
    dcast.data.table( sample_id ~ qc_var, value.var = 'lab')
  lab_mat   = lab_wide[, -'sample_id', with = FALSE] %>%
    as.matrix %>% set_rownames(lab_wide$sample_id)

  # lots of colours
  z_cols    = cols_fn(seq(-2, 2, 0.2), 0.2, 'RdBu',
    pal_dir = -1, range = 'symmetric')

  # define function for labelling
  labelling_fn <- function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%s", lab_mat[i, j]), x, y,
      gp = gpar(fontsize = 6))
  }

  # make sample annotations
  rows_dt   = copy(meta_input) %>% setkey('sample_id') %>%
    .[rownames(z_mat)]

  # do colours for patients
  pats_dt     = meta_input[, .(subject_n = .N), by = subject_id] %>%
    .[ subject_n > 1, .(subject_id, subject_n) ] %>% .[ order(-subject_n)]
  pat_list    = pats_dt$subject_id
  pat_cols    = nice_cols[seq.int(length(pat_list))] %>% setNames(pat_list)
  donor_vals  = rows_dt$subject_id
  donor_vals[ !(donor_vals %in% names(pat_cols)) ] = NA
  les_vals    = rows_dt$lesion_type %>% fct_drop %>% levels
  les_cols    = brewer.pal(length(les_vals), 'BrBG') %>%
    setNames(les_vals)

  # make col annotations
  row_annots  = rowAnnotation(
    lesion    = rows_dt$lesion_type,
    diagnosis = rows_dt$diagnosis,
    donor     = donor_vals,
    sex       = rows_dt$sex,
    age       = rows_dt$age_at_death,
    pmi_cat   = rows_dt$pmi_cat,
    col     = list(
      lesion    = les_cols,
      diagnosis = disease_cols,
      donor     = pat_cols,
      sex       = c(M = 'black', F = 'white'),
      age       = cols_fn(rows_dt$age_at_death, 10, 'Greys',
        pal_dir = 1, range='natural'),
      pmi_cat   = c(up_to_6H = '#f7fcf5', `6H_to_12H` = '#41ab5d',
        over_12H = '#005a32')
    ),
    annotation_name_side='top', show_legend = c(donor = FALSE)
  )

  # do column titles
  var_lookup  = c(
    log10_N       = "no. cells",
    log_counts    = "library size",
    log_feats     = "no. features",
    logit_mito    = "mito pct.",
    logit_splice  = "pct. spliced"
  )

  # do heatmap for these genes
  name_str  = "median\nQC value\n(z-scored,\nsign flipped:\n+ve = good)"
  hm_obj    = Heatmap(
    matrix = z_mat, col = z_cols, cell_fun = labelling_fn,
    column_labels = var_lookup[ colnames(z_mat) ],
    cluster_rows = TRUE, cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 7),
    name = name_str, heatmap_legend_param = list(title = name_str),
    left_annotation = row_annots,
    show_row_names = FALSE, show_column_names = FALSE,
    # row_split = row_split, cluster_row_slices = FALSE,
    row_names_side = "left", column_names_side = "top")

  return(hm_obj)
}

calc_qc_summary <- function(qc_dt, kept_dt) {
  qc_summary  = merge(
    qc_dt[, .(n_pre_QC = .N), by = .(sample_id)],
    kept_dt[, .(n_post_QC = .N), by = .(sample_id)],
    by = 'sample_id', all = TRUE) %>%
    .[ is.na(n_post_QC), n_post_QC := 0 ] %>%
    .[, n_excluded    := n_pre_QC - n_post_QC ] %>%
    .[, pct_excluded := round(100*(1 - n_post_QC / n_pre_QC), 1) ] %>%
    .[ order(-pct_excluded, -n_post_QC) ]
}




############## FUNCTIONS FORM SampleQC package

# from make_qc_dt.R
make_qc_dt <- function(qc_df, sample_var = 'sample_id',
                       qc_names = c('log_counts', 'log_feats', 'logit_mito'), annot_vars = NULL) {

  # some checks
  if ( 'DFrame' %in% class(qc_df) )
    qc_df      = as.data.frame(qc_df)
  assert_that( is.data.frame(qc_df), msg = "qc_df must be a data.frame" )

  assert_that( sample_var %in% colnames(qc_df),
               msg = sprintf("%s is listed as variable for samples but is not in data.frame",
                             sample_var))

  reserved_ns  = c('sample_id', 'group_id', 'cell_id')
  assert_that( length(intersect(annot_vars, reserved_ns)) == 0,
               msg = paste0("The following variable names are reserved and cannot be used ",
                            "as annot_vars:\n", paste(reserved_ns, collapse = ", ")))

  assert_that( all(annot_vars %in% names(qc_df)),
               msg = sprintf("the following variables are listed in annot_vars but not in qc_df:\n%s",
                             paste(setdiff(annot_vars, names(qc_df)), collapse = ", ")))

  # set up qc_dt
  qc_dt   = .init_qc_dt(qc_df, sample_var)

  # add known metrics
  if ('log_counts' %in% qc_names) {
    qc_dt   = .add_log_counts(qc_dt, qc_df)
  }
  if ('log_feats' %in% qc_names) {
    qc_dt   = .add_log_feats(qc_dt, qc_df)
  }
  if ('logit_mito' %in% qc_names) {
    qc_dt   = .add_logit_mito(qc_dt, qc_df)
  }
  if ('splice_ratio' %in% qc_names) {
    qc_dt   = .add_splice_ratio(qc_dt, qc_df)
  }

  # add unknown metrics
  qc_dt   = .add_unknown_metrics(qc_dt, qc_df, qc_names)

  # add some useful annotations
  qc_dt   = .add_qc_annots(qc_dt)

  # add specified annotation variables
  qc_dt   = .add_annot_vars(qc_dt, qc_df, annot_vars)

  # put in nice order
  setcolorder(qc_dt, c('cell_id', 'sample_id', qc_names))

  # double-check everything is ok
  .check_qc_dt(qc_dt, qc_names, annot_vars)

  return(qc_dt)
}


.init_qc_dt <- function(qc_df, sample_var) {
  # add cell identifiers
  if ('cell_id' %in% colnames(qc_df)) {
    qc_dt   = data.table(cell_id = qc_df$cell_id)
  } else if ( !is.null(rownames(qc_df)) ) {
    qc_dt   = data.table(cell_id = rownames(qc_df))
  } else {
    stop("input data.frame must have either rownames or 'cell_id' as a column")
  }
  assert_that( length(unique(qc_dt$cell_id)) == nrow(qc_dt),
               msg = "cell identifiers are not unique")

  # add sample identifiers
  qc_dt[, sample_id := qc_df[[sample_var]] ]

  # check no missing values or NAs
  assert_that( all(!is.na(qc_dt$cell_id)), msg = "missing values in cell_id")
  assert_that( all(!is.na(qc_dt$sample_id)), msg = "missing values in sample_id")

  return(qc_dt)
}


.add_log_counts <- function(qc_dt, qc_df) {
  # what names do we have, and want?
  df_names  = colnames(qc_df)
  valid_ns  = c('log_counts', 'total', 'sum', 'nCount_RNA')

  # check which are present
  here_ns   = vapply(valid_ns, function(v) v %in% df_names, logical(1))
  assert_that( sum(here_ns) >= 1,
               msg = paste0(
                 "no valid column present for log_counts\n",
                 paste0("valid columns are: ", paste(valid_ns, collapse = ", "))
               ))
  to_use    = valid_ns[here_ns][[1]]

  # add values
  if (to_use %in% 'log_counts') {
    qc_dt[, log_counts := qc_df[[ to_use ]] ]

  } else if (to_use %in% c('total', 'sum', 'nCount_RNA')) {
    assert_that( all(qc_df[[ to_use ]] > 0) )
    qc_dt[, log_counts := log10(qc_df[[ to_use ]]) ]

  } else {
    stop("log_counts requested but required variables not present")

  }

  # do some checks
  assert_that( "log_counts" %in% names(qc_dt) )
  assert_that( !any(is.na(qc_dt$log_counts)),
               msg = "some log_counts values are NA")
  assert_that( !any(is.infinite(qc_dt$log_counts)),
               msg = "some log_counts values are infinite")
  assert_that( all(qc_dt$log_counts >= 0),
               msg = "some log_counts values are <= 0")

  return(qc_dt)
}


.add_log_feats <- function(qc_dt, qc_df) {
  # what names do we have, and want?
  df_names  = colnames(qc_df)
  valid_ns  = c('log_feats', 'detected', 'nFeature_RNA')

  # check which are present
  here_ns   = vapply(valid_ns, function(v) v %in% df_names, logical(1))
  assert_that( sum(here_ns) >= 1,
               msg = paste0(
                 "no valid column present for log_feats\n",
                 paste0("valid columns are: ", paste(valid_ns, collapse = ", "))
               ))
  to_use    = valid_ns[here_ns][[1]]

  # add values
  if (to_use %in% 'log_feats') {
    qc_dt[, log_feats := qc_df[[ to_use ]] ]

  } else if (to_use %in% c('detected', 'nFeature_RNA')) {
    assert_that( all(qc_df[[ to_use ]] > 0) )
    qc_dt[, log_feats := log10(qc_df[[ to_use ]]) ]

  } else {
    stop("log_feats requested but required variables not present")

  }

  # do some checks
  assert_that( "log_feats" %in% names(qc_dt) )
  assert_that( !any(is.na(qc_dt$log_feats)),
               msg = "some log_feats values are NA")
  assert_that( !any(is.infinite(qc_dt$log_feats)),
               msg = "some log_feats values are infinite")
  assert_that( all(qc_dt$log_feats >= 0),
               msg = "some log_feats values are <= 0")

  return(qc_dt)
}


.add_logit_mito <- function(qc_dt, qc_df) {
  # what names do we have, and want?
  df_names  = colnames(qc_df)

  # add logit-transformed mitochondrial proportion to qc_dt
  if ('logit_mito' %in% df_names) {
    qc_dt[, logit_mito  := qc_df$logit_mito ]

  } else if ( ('subsets_mito_sum' %in% df_names) & ('total' %in% df_names) ) {
    qc_dt[, logit_mito  := qlogis( (qc_df$subsets_mito_sum + 1) / (qc_df$total + 2) ) ]

  } else if ( ('subsets_mt_sum' %in% df_names) & ('total' %in% df_names) ) {
    qc_dt[, logit_mito  := qlogis( (qc_df$subsets_mt_sum + 1) / (qc_df$total + 2) ) ]

  } else if ( ('percent.mt' %in% df_names) & ('nCount_RNA' %in% df_names) ) {
    total_counts  = qc_df$nCount_RNA
    mt_counts     = qc_df$nCount_RNA * qc_df$percent.mt / 100
    assert_that( all(abs(mt_counts - round(mt_counts, 0)) < 1e-10) )
    qc_dt[, logit_mito  := qlogis( (mt_counts + 1) / (total_counts + 2) ) ]

  } else if ( ('mito_prop' %in% df_names) & ('log_counts' %in% df_names) ) {
    total_counts  = 10^qc_df$log_counts
    mt_counts     = qc_df$mito_prop * total_counts
    assert_that( all(abs(mt_counts - round(mt_counts, 0)) < 1e-8) )
    qc_dt[, logit_mito  := qlogis( (mt_counts + 1) / (total_counts + 2) ) ]

  } else {
    stop("logit_mito requested but required variables not present")
  }

  # do some checks
  assert_that( "logit_mito" %in% names(qc_dt) )
  assert_that( !any(is.na(qc_dt$logit_mito)),
               msg = "some logit_mito values are NA")
  assert_that( !any(is.infinite(qc_dt$logit_mito)),
               msg = "some logit_mito values are infinite")

  return(qc_dt)
}


.add_splice_ratio <- function(qc_dt, qc_df) {
  # what names do we have, and want?
  df_names  = colnames(qc_df)

  # add logit-transformed mitochondrial proportion to qc_dt
  if ('splice_ratio' %in% df_names) {
    qc_dt[, splice_ratio  := qc_df$splice_ratio ]

  } else if ( ('total_spliced' %in% df_names) & ('total_unspliced' %in% df_names) ) {
    qc_dt[, splice_ratio  := log2( (qc_df$total_spliced + 1) / (qc_df$total_unspliced + 1) ) ]

  } else {
    stop("logit_mito requested but required variables not present")

  }

  # do some checks
  assert_that( "splice_ratio" %in% names(qc_dt) )
  assert_that( !any(is.na(qc_dt$logit_mito)),
               msg = "some logit_mito values are NA")
  assert_that( !any(is.infinite(qc_dt$logit_mito)),
               msg = "some logit_mito values are infinite")

  return(qc_dt)
}


list_known_metrics <- function() {
  return(c('log_counts', 'log_feats', 'logit_mito', 'splice_ratio'))
}


.add_unknown_metrics <- function(qc_dt, qc_df, qc_names)  {
  # anything to add?
  to_add  = setdiff(qc_names, list_known_metrics())
  if ( length(to_add) == 0 )
    return(qc_dt)

  # add them
  message("adding the following metrics that are not known to `SampleQC`:")
  message(paste(to_add, collapse = ", "))
  for (v in to_add) {
    assert_that( v %in% names(qc_df), msg = paste0(v, " missing from qc_df"))
    set(qc_dt, i = NULL, v, qc_df[[v]])
    assert_that( !any(is.na(qc_dt$v)), msg = paste0("NA values for ", v))
    assert_that( !any(is.infinite(qc_dt$v)), msg = paste0("infinite values for ", v))
  }

  return(qc_dt)
}


.add_qc_annots <- function(qc_dt) {
  # add annotations for sample size
  qc_dt[, log_N  := log10(.N), by='sample_id']

  # and factor version
  N_cuts    = c(1,100,200,400,1000,2000,4000,10000,20000,40000,Inf)
  N_labs    = paste0('<=', N_cuts[-1])
  qc_dt[, N_cat  := factor(
    cut(10^log_N, breaks = N_cuts, labels = N_labs),
    levels = N_labs), by = 'sample_id']

  # add annotations relating to library sizes
  if ('log_counts' %in% names(qc_dt) ) {
    # add median log counts per sample
    qc_dt[, med_counts  := median(log_counts), by='sample_id']

    # put mito level into categories
    counts_cuts = c(1,100,300,1000,3000,10000,30000, Inf)
    counts_labs = paste0('<=', counts_cuts[-1])
    qc_dt[, counts_cat  := factor(
      cut(10^med_counts, breaks = counts_cuts, labels = counts_labs),
      levels = counts_labs), by = 'sample_id']
  }

  # add annotations relating to features
  if ('log_feats' %in% names(qc_dt) ) {
    # add median log feats per sample
    qc_dt[, med_feats   := median(log_feats), by='sample_id']

    # put mito level into categories
    feats_cuts = c(1,100,300,1000,3000,10000,30000, Inf)
    feats_labs = paste0('<=', feats_cuts[-1])
    qc_dt[, feats_cat  := factor(
      cut(10^med_feats, breaks = feats_cuts, labels = feats_labs),
      levels = feats_labs), by = 'sample_id']
  }

  # add annotations relating to mitochondrial proportions
  if ('logit_mito' %in% names(qc_dt) ) {
    # add median mito proportion
    qc_dt[, med_mito  := median(plogis(logit_mito)), by='sample_id']

    # put mito level into categories
    mito_cuts   = c(0,0.01,0.05,0.1,0.2,0.5,1)
    mito_labs   = paste0('<=', mito_cuts[-1])
    qc_dt[, mito_cat  := factor(
      cut(med_mito, breaks = mito_cuts, labels = mito_labs),
      levels = mito_labs), by = 'sample_id']
  }

  # add annotations relating to mitochondrial proportions
  if ('splice_ratio' %in% names(qc_dt) ) {
    # add median mito proportion
    qc_dt[, med_splice  := median(plogis(splice_ratio)), by='sample_id']

    # put mito level into categories
    splice_cuts = c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 1)
    splice_labs = paste0('<=', splice_cuts[-1])
    qc_dt[, splice_cat  := factor(
      cut(med_splice, breaks = splice_cuts, labels = splice_labs),
      levels = splice_labs), by = 'sample_id']
  }

  return(qc_dt)
}


.add_annot_vars <- function(qc_dt, qc_df, annot_vars) {
  # check all present
  assert_that( all(annot_vars %in% names(qc_df)) )
  # add them
  for (v in annot_vars)
    qc_dt[[v]] = qc_df[[v]]

  # check that they all sample level
  for (v in annot_vars) {
    check_dt  = qc_dt[, c('sample_id', v), with = FALSE] %>%
      .[, .N, by = c('sample_id', v) ]
    assert_that( nrow(check_dt) == length(unique(check_dt$sample_id)),
                 msg = paste0("annotation variable ", v, " has more than one value per\n",
                              "sample (should be sample-level only)"))
  }

  return(qc_dt)
}


.check_qc_dt <- function(qc_dt, qc_names, annot_vars) {
  # unpack
  col_names   = colnames(qc_dt)

  # check specific names
  if ('log_counts' %in% col_names)
    assert_that( all(qc_dt$log_counts >= 0) )
  if ('log_feats' %in% col_names)
    assert_that( all(qc_dt$log_feats >= 0) )
  if ('logit_mito' %in% col_names)
    assert_that( all(is.finite(qc_dt$logit_mito)) )
  if ('splice_ratio' %in% col_names)
    assert_that( all(is.finite(qc_dt$splice_ratio)) )

  # check qc metrics and annotations for NAs
  for (n in qc_names) {
    assert_that( all(!is.na(qc_dt[[n]])) )
  }
  annots_auto   = c(
    "med_counts", "counts_cat",
    "med_feats", "feats_cat",
    "med_mito", "mito_cat",
    "med_splice", "splice_cat",
    "log_N", "N_cat")
  for (n in c(annots_auto, annot_vars)) {
    if ( n %in% names(qc_dt) )
      assert_that( all(!is.na(qc_dt[[n]])),
                   msg = paste0('NA present in an annotation variable, ', n) )
  }
}


#
