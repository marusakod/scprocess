# METADATA / PREVALENCE PLOT FUNCTIONS -----------------------------------------

# Stacked bar chart: proportion of each cluster within each level of meta_var.
# cl_ord: if a valid cluster name, rows are sorted by that cluster's proportion;
#         otherwise rows are ordered by the variable's factor levels.
plot_metadata_barplot <- function(meta_dt, meta_var, col_pal, cl_ord = 'no') {
  bar_dt = copy(meta_dt) %>%
    .[, .(n = sum(n_cells)), by = c('cluster', meta_var)]

  cl_ls = levels(bar_dt$cluster)

  if (cl_ord %in% cl_ls) {
    bar_y_ord = bar_dt %>%
      .[, tot  := sum(n), by = meta_var] %>%
      .[, prop := n / tot] %>%
      .[cluster == cl_ord] %>%
      setorder(-prop) %>%
      .[[meta_var]] %>%
      as.character()
    missing   = setdiff(bar_dt[[meta_var]], bar_y_ord)
    bar_y_ord = c(bar_y_ord, missing)
    bar_dt[[meta_var]] = factor(bar_dt[[meta_var]], levels = bar_y_ord)
  }

  y_counts = bar_dt %>%
    .[, .(n = sum(n)), by = meta_var] %>%
    .[, n := (n / 1e3) %>% signif(2) %>% round]
  recode_vector = sprintf("%s (~%sk)", y_counts[[meta_var]], y_counts$n) %>%
    setNames(y_counts[[meta_var]])
  levels(bar_dt[[meta_var]]) = recode_vector[ levels(bar_dt[[meta_var]]) ]

  ggplot(bar_dt) +
    aes(
      y    = forcats::fct_rev(get(meta_var)),
      x    = n * 100,
      fill = forcats::fct_rev(cluster)
    ) +
    geom_bar(position = 'fill', stat = 'identity') +
    theme_classic() +
    labs(x = 'pct. of condition', y = NULL, fill = NULL) +
    theme(
      axis.text.x  = element_text(size = FONT_SMALL),
      axis.text.y  = element_text(size = FONT_TEXT),
      axis.title   = element_text(size = FONT_AXIS),
      legend.text  = element_text(size = FONT_SMALL)
    ) +
    scale_fill_manual(values = col_pal) +
    scale_x_continuous(breaks = pretty_breaks()) +
    guides(fill = guide_legend(reverse = TRUE, ncol = 1))
}

# Build a CLR-transformed composition data.table for use in plot_metadata_dotplot.
make_composition_dt <- function(meta_dt, sample_col = 'sample_id', zero_const = 1,
                                prior_size = 2e1) {
  meta_dt   = meta_dt %>% setnames(., sample_col, 'sample')
  meta_vars = names(meta_dt) %>% setdiff(c("cell_id", "cluster", "n_cells"))
  meta_tmp  = meta_dt %>%
    .[, .(log10_N_sample = log10(sum(n_cells))), by = meta_vars]

  prior_dt = meta_dt %>%
    .[, .(N_all = sum(n_cells)), by = cluster] %>%
    .[, prop0 := N_all / sum(N_all)] %>%
    .[, N0    := prop0 * prior_size]

  match_dt = expand.grid(
    sample  = unique(meta_dt$sample),
    cluster = unique(meta_dt$cluster)
  ) %>% as.data.table

  props_dt = meta_dt %>%
    .[, .(N = sum(n_cells)), by = .(cluster, sample)] %>%
    merge(match_dt, by = c('cluster', 'sample'), all.y = TRUE) %>%
    merge(prior_dt, by = 'cluster') %>%
    .[, log10_N         := log10(N + 1)] %>%
    .[is.na(N), log10_N := 0] %>%
    .[, N               := N * 1.0000001] %>%
    .[is.na(N), N       := (runif(.N) * 0.9 * zero_const) + (0.1 * zero_const)] %>%
    .[, prop_sample      := N / sum(N),                    by = sample] %>%
    .[, prop0            := (N + N0) / (sum(N) + prior_size), by = .(sample)] %>%
    .[, log10_sample     := log10(prop0)] %>%
    .[, clr              := log10_sample - mean(log10_sample, na.rm = TRUE), by = sample] %>%
    .[, log10_sample_z   := (log10_sample - median(log10_sample, na.rm = TRUE)) /
                             mad(log10_sample, na.rm = TRUE), by = cluster] %>%
    .[, prop_cl          := N / sum(N), by = cluster] %>%
    merge(meta_tmp, by = "sample")

  return(props_dt)
}

# Composition dotplot: one dot per sample, with median + 5th–95th percentile
# segment per cluster, faceted by facet_by.
plot_metadata_dotplot <- function(props_dt, leg_title, x_var, facet_by, col_var, col_pal,
                                  order = FALSE, n_cols = 2) {
  stats_dt = props_dt %>%
    .[, .(
      med_log10_sample = median(log10_sample, na.rm = TRUE),
      q05_log10_sample = quantile(log10_sample, 0.05, na.rm = TRUE),
      q95_log10_sample = quantile(log10_sample, 0.95, na.rm = TRUE)
    ),
    by = c(x_var, facet_by, 'cluster')]

  dot_size  = 6
  base_size = 20

  if (order) {
    means = stats_dt %>%
      .[, .(mean = mean(med_log10_sample)), by = x_var] %>%
      tibble::deframe() %>%
      sort(decreasing = TRUE)
    ord = names(means)
  } else {
    ord = levels(stats_dt[[x_var]])
  }

  min_y = min(PCT_BRKS)

  ggplot() +
    geom_quasirandom(
      data = props_dt[ sample(.N, .N) ],
      aes(x     = forcats::fct_relevel(get(x_var), ord),
          group = forcats::fct_relevel(get(x_var), ord),
          y     = log10_sample %>% pmax(min_y),
          fill  = get(col_var)),
      shape = 21, size = dot_size
    ) +
    geom_segment(
      data = stats_dt,
      aes(x    = forcats::fct_relevel(get(x_var), ord),
          xend = forcats::fct_relevel(get(x_var), ord),
          y    = q05_log10_sample %>% pmax(min_y),
          yend = q95_log10_sample %>% pmax(min_y)),
      linewidth = 1, colour = "black"
    ) +
    geom_point(
      data = stats_dt,
      aes(x = forcats::fct_relevel(get(x_var), ord),
          y = med_log10_sample %>% pmax(min_y)),
      colour = "black", size = 10
    ) +
    scale_x_discrete(breaks = ord) +
    scale_y_continuous(breaks = PCT_BRKS, labels = PCT_LABS) +
    scale_fill_manual(values = col_pal,
      guide = guide_legend(override.aes = list(size = 5, alpha = 1), nrow = 5)) +
    facet_wrap(~ get(facet_by), scales = "free_y", ncol = n_cols) +
    theme_classic(base_size = base_size) +
    theme(
      axis.text.x    = element_text(angle = -45, hjust = 0, vjust = 0.5),
      strip.text     = element_text(size = FONT_AXIS),
      axis.title.y   = element_text(size = FONT_AXIS),
      axis.text.y    = element_text(size = FONT_TEXT),
      legend.text    = element_text(size = FONT_TEXT),
      legend.title   = element_text(size = FONT_AXIS),
      plot.margin    = margin(0, 5, 0, 0.5, "cm"),
      legend.position  = 'bottom',
      legend.direction = 'horizontal'
    ) +
    labs(y = "pct. of sample", x = NULL, fill = leg_title)
}

# Faceted density (bin2d) UMAP plots, one panel per level of meta_var.
plot_metadata_density <- function(meta_dt, meta_var, x_range, y_range, n_cols = NULL) {
  plot_ratio = 1.5
  n_panels   = meta_dt[[ meta_var ]] %>% unique %>% length

  if (is.na(n_cols) | is.null(n_cols)) {
    n_cols = sqrt(n_panels * plot_ratio) %>% ceiling
  } else {
    n_cols = min(n_cols, n_panels)
  }
  n_rows = (n_panels / n_cols) %>% ceiling

  ggplot(meta_dt, aes(x = UMAP_1, y = UMAP_2)) +
    geom_bin2d(bins = 40,
      aes(fill = after_stat(density) %>% pmin(0.01) %>% pmax(0.0001))) +
    scale_fill_distiller(palette = "RdBu", trans = "log10", limits = c(0.0001, 0.01),
      breaks = c(0.0001, 0.001, 0.01), labels = c("0.01%", "0.1%", "1%")) +
    facet_wrap(~ get(meta_var), ncol = n_cols, nrow = n_rows) +
    coord_cartesian(xlim = x_range, ylim = y_range) +
    theme_classic() +
    theme(
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      panel.grid   = element_blank(),
      axis.title   = element_text(size = FONT_TEXT),
      strip.text   = element_text(size = FONT_AXIS),
      legend.text  = element_text(size = FONT_SMALL),
      legend.title = element_text(size = FONT_TEXT),
      aspect.ratio = 1
    ) +
    labs(fill = "pct. of\npanel", x = 'UMAP 1', y = 'UMAP 2')
}
