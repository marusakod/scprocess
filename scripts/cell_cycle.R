suppressPackageStartupMessages({
  library(data.table)
})


load_cell_cycle_report_data <- function(scores_f, origin_f, sensitivity_f,
  diagnostics_f, ellipse_f, coldata_f, marker_expression_f) {

  # Load saved outputs, retain QC-passed cells, and prepare shared plot inputs.
  scores <- fread(scores_f)
  origin <- fread(origin_f)
  origin_sensitivity <- fread(sensitivity_f)
  diagnostics <- fread(diagnostics_f)
  display_sample <- cell_cycle_display_sample(
    diagnostics[['inclusion_probability']], seed=origin[['seed']][1])
  blob <- read_cell_cycle_blob(ellipse_f)
  coldata <- fread(coldata_f, select=c('cell_id', 'keep'))
  scores_qc <- merge(scores, coldata[keep == TRUE], by='cell_id')
  scores_qc <- merge(
    scores_qc,
    diagnostics[, .(cell_id, kde_density, inclusion_probability)],
    by='cell_id')
  scores_qc[, approximate_phase := assign_approximate_cell_cycle_phase(tricycle_theta)]
  marker_expression <- read_cell_cycle_marker_expression(
    marker_expression_f, scores_qc)
  marker_annotation <- cell_cycle_marker_annotation()
  if (nrow(marker_expression)) {
    marker_expression <- merge(marker_expression, marker_annotation, by='gene')
    marker_expression[, gene := factor(gene, levels=marker_annotation[['gene']])]
  }
  expression_smooth <- smooth_periodic_expression(marker_expression)

  pc_x_limits <- range(c(diagnostics[['tricycle_pc1']], origin[['tricycle_pc1_origin']]),
                       finite=TRUE)
  pc_y_limits <- range(c(diagnostics[['tricycle_pc2']], origin[['tricycle_pc2_origin']]),
                       finite=TRUE)
  pc_x_limits <- pc_x_limits + c(-1, 1) * max(diff(pc_x_limits) * 0.025, 1e-6)
  pc_y_limits <- pc_y_limits + c(-1, 1) * max(diff(pc_y_limits) * 0.025, 1e-6)
  origin_x <- origin[['tricycle_pc1_origin']][1]
  origin_y <- origin[['tricycle_pc2_origin']][1]
  pc_corners <- CJ(x=pc_x_limits, y=pc_y_limits)
  phase_ray_radius <- max(sqrt(
    (pc_corners[['x']] - origin_x)^2 + (pc_corners[['y']] - origin_y)^2))
  phase_boundaries <- data.table(theta=c(0.25, 0.5, 1, 1.75) * pi)
  phase_boundaries[, `:=`(
    x=origin_x,
    y=origin_y,
    xend=origin_x + phase_ray_radius * cos(theta),
    yend=origin_y + phase_ray_radius * sin(theta))]
  theta_breaks <- seq(0, 2, 0.5)
  theta_labels <- c('0', expression(pi/2), expression(pi),
                    expression(3*pi/2), expression(2*pi))
  cyclic_colours <- hsv(seq(0, 1, length.out=13), s=0.7, v=0.8)
  phase_colours <- c(
    'G1/G0'='#6A3D9A', 'G1/S'='#1F78B4', 'S'='#33A02C', 'G2/M'='#E31A1C')
  gene_colours <- c(
    'CDKN1A'='#636363',
    'PCNA'='#A1D99B', 'MCM4'='#74C476', 'MCM5'='#31A354', 'TYMS'='#006D2C',
    'TOP2A'='#FEE6CE', 'CCNB1'='#FDD0A2', 'CDK1'='#FDAE6B',
    'TPX2'='#FD8D3C', 'SMC2'='#F16913', 'MKI67'='#D94801', 'UBE2C'='#8C2D04')
  gene_legend_labels <- setNames(
    sprintf('%s · %s', marker_annotation[['marker_phase']], marker_annotation[['gene']]),
    marker_annotation[['gene']])
  density_percentile_levels <- as.character(sort(unique(origin_sensitivity[['density_percentile']])))
  origin_is_manual <- identical(origin[['origin_method']][1], 'manual')
  origin_sensitivity_points <- if (origin_is_manual) {
    origin_sensitivity
  } else {
    origin_sensitivity[selected == FALSE]
  }
  percentile_caption <- sprintf(
    'Density percentile: %s; full weight: %s of cells',
    format_cell_cycle_density_percentile(origin[['density_percentile']][1]),
    format_cell_cycle_density_percentile(origin[['full_weight_percent']][1]))
  origin_caption <- if (origin_is_manual) {
    sprintf('Manual origin: (%g, %g); %s',
      origin[['tricycle_pc1_origin']][1], origin[['tricycle_pc2_origin']][1],
      percentile_caption)
  } else {
    percentile_caption
  }

  # Per-cell tables, a one-row origin table, marker curves and shared plot styles.
  list(
    scores_qc = scores_qc,
    origin = origin,
    origin_sensitivity = origin_sensitivity,
    diagnostics = diagnostics,
    display_sample = display_sample,
    blob = blob,
    marker_expression = marker_expression,
    marker_annotation = marker_annotation,
    expression_smooth = expression_smooth,
    pc_x_limits = pc_x_limits,
    pc_y_limits = pc_y_limits,
    phase_boundaries = phase_boundaries,
    theta_breaks = theta_breaks,
    theta_labels = theta_labels,
    cyclic_colours = cyclic_colours,
    phase_colours = phase_colours,
    gene_colours = gene_colours,
    gene_legend_labels = gene_legend_labels,
    density_percentile_levels = density_percentile_levels,
    origin_sensitivity_points = origin_sensitivity_points,
    origin_caption = origin_caption)
}


plot_cell_cycle_expression <- function(report_data) {
  if (nrow(report_data$expression_smooth[['curves']])) {
    curve_dt <- merge(report_data$expression_smooth[['curves']], report_data$marker_annotation, by='gene')
    curve_dt[, gene := factor(gene, levels=report_data$marker_annotation[['gene']])]
    curve_dt[, theta_over_pi := theta / pi]
    theta_projection_dt <- .evenly_spaced_cell_cycle_rows(report_data$scores_qc, 50000L)
    g_theta_projection <- ggplot(
        theta_projection_dt,
        aes(tricycle_pc1, tricycle_pc2, colour=tricycle_theta)) +
      geom_point(size=0.25, alpha=0.6) +
      .cell_cycle_phase_layer(report_data$phase_boundaries) +
      scale_colour_gradientn(
        colours=report_data$cyclic_colours, limits=c(0, 2*pi),
        breaks=seq(0, 2*pi, pi/2), labels=report_data$theta_labels, name='theta') +
      .cell_cycle_origin_layer(report_data$origin) +
      coord_equal(xlim=report_data$pc_x_limits, ylim=report_data$pc_y_limits, expand=FALSE) +
      theme_bw() +
      labs(x='tricycle PC1', y='tricycle PC2')

    g_expression_curves <- ggplot() +
      geom_line(
        data=curve_dt,
        aes(theta_over_pi, mean_scaled, colour=gene), linewidth=0.9) +
      geom_vline(xintercept=c(0.25, 0.5, 1, 1.75),
                 linewidth=0.25, linetype='dotted', colour='grey50') +
      scale_colour_manual(values=report_data$gene_colours, labels=report_data$gene_legend_labels) +
      scale_x_continuous(breaks=report_data$theta_breaks, labels=report_data$theta_labels,
                         limits=c(0, 2), expand=c(0, 0)) +
      scale_y_continuous(breaks=pretty_breaks()) +
      coord_cartesian(ylim=c(0, 1.05)) +
      theme_bw() +
      labs(x=expression('tricycle theta / '*pi),
           y='pooled mean expression (curve maximum = 1)',
           colour='marker')
    g_theta_projection + g_expression_curves +
      plot_layout(widths=c(1, 1.4))
  } else {
    cat('No marker-expression curves could be estimated.')
  }
}


plot_cell_cycle_pooled_projection <- function(report_data) {
  g_pooled <- ggplot(report_data$diagnostics, aes(tricycle_pc1, tricycle_pc2)) +
    geom_hex(
      aes(colour=after_stat(count), fill=after_scale(colour)),
      bins=180, linewidth=0) +
    geom_path(
      data=report_data$origin_sensitivity[order(density_percentile)],
      aes(tricycle_pc1_origin, tricycle_pc2_origin), inherit.aes=FALSE,
      linewidth=1.3, colour='white') +
    geom_path(
      data=report_data$origin_sensitivity[order(density_percentile)],
      aes(tricycle_pc1_origin, tricycle_pc2_origin), inherit.aes=FALSE,
      linewidth=0.55, colour='grey50') +
    geom_point(
      data=report_data$origin_sensitivity_points,
      aes(
        tricycle_pc1_origin, tricycle_pc2_origin,
        fill=factor(density_percentile, levels=report_data$density_percentile_levels)),
      inherit.aes=FALSE, shape=21, size=2.5, stroke=0.5, colour='grey20') +
    geom_path(data=report_data$blob$ellipse, colour='white', linewidth=1.1) +
    geom_path(data=report_data$blob$ellipse, colour='#D55E00', linewidth=0.6) +
    .cell_cycle_phase_layer(report_data$phase_boundaries) +
    .cell_cycle_origin_layer(report_data$origin) +
    scale_colour_viridis_c(
      trans='log10', name='cell count',
      guide=guide_colourbar(order=1, direction='horizontal')) +
    scale_fill_discrete(
      palette=function(n) grey.colors(n, start=0.88, end=0.28),
      name='density percentile', limits=report_data$density_percentile_levels, drop=FALSE,
      labels=function(x) format_cell_cycle_density_percentile(as.numeric(x)),
      guide=guide_legend(order=2, direction='horizontal', nrow=1)) +
    coord_equal(xlim=report_data$pc_x_limits, ylim=report_data$pc_y_limits, expand=FALSE) +
    theme_bw() +
    labs(
      title='All QC-passed cells', x='tricycle PC1', y='tricycle PC2',
      subtitle=report_data$blob$label, caption=report_data$origin_caption)

  g_retained <- ggplot(
      report_data$diagnostics[report_data$display_sample$retained],
      aes(tricycle_pc1, tricycle_pc2)) +
    geom_point(size=0.5, alpha=0.7, colour='#1F4E79') +
    .cell_cycle_phase_layer(report_data$phase_boundaries) +
    .cell_cycle_origin_layer(report_data$origin) +
    coord_equal(xlim=report_data$pc_x_limits, ylim=report_data$pc_y_limits, expand=FALSE) +
    theme_bw() +
    labs(
      title=sprintf(
        'Display sample (%s of %s cells)',
        scales::comma(sum(report_data$display_sample$retained)),
        scales::comma(nrow(report_data$diagnostics))),
      x='tricycle PC1', y='tricycle PC2')

  g_pooled + g_retained + plot_layout(guides='collect') &
    theme(legend.position='bottom')
}


plot_cell_cycle_spatial_expression <- function(report_data) {
  report_data$marker_expression = copy(report_data$marker_expression)

  if (nrow(report_data$marker_expression)) {
    report_data$marker_expression[, expression_q98 := {
      positive <- log_normalized_expression[log_normalized_expression > 0]
      cap <- if (length(positive)) as.numeric(quantile(positive, 0.98)) else 1
      pmin(log_normalized_expression / max(cap, 1e-12), 1)
    }, by=gene]
    gene_summaries <- report_data$marker_expression[, .(
      label=sprintf('%s (%s > 0)', gene[1], scales::comma(sum(log_normalized_expression > 0)))
    ), by=gene]
    gene_labels <- setNames(gene_summaries[['label']], as.character(gene_summaries[['gene']]))
    set.seed(20230308)
    report_data$marker_expression <- report_data$marker_expression[sample.int(.N)]

    ggplot(report_data$marker_expression,
           aes(tricycle_pc1, tricycle_pc2, colour=expression_q98)) +
      geom_point(size=0.18, alpha=0.8) +
      .cell_cycle_phase_layer(report_data$phase_boundaries) +
      .cell_cycle_origin_layer(report_data$origin) +
      facet_wrap(~gene, ncol=4, labeller=as_labeller(gene_labels)) +
      scale_colour_viridis_c(option='magma', limits=c(0, 1),
                             name='expression /\ngene q98') +
      coord_equal(xlim=report_data$pc_x_limits, ylim=report_data$pc_y_limits, expand=FALSE) +
      theme_bw() +
      theme(
        panel.grid=element_blank(),
        strip.background=element_rect(fill='white', colour='grey50')) +
      labs(x='tricycle PC1', y='tricycle PC2')
  } else {
    cat('No canonical marker-expression values were available for plotting.')
  }
}


render_cell_cycle_sample_diagnostics <- function(report_data) {
  for (sample in sort(unique(report_data$scores_qc$sample_id))) {
    cat('### ', sample, '\n\n')
    sample_dt <- report_data$scores_qc[sample_id == sample]
    sample_plot_dt <- .evenly_spaced_cell_cycle_rows(sample_dt, 20000L)
    g_projection <- ggplot(sample_plot_dt,
      aes(
        tricycle_pc1, tricycle_pc2, colour=tricycle_theta,
        size=inclusion_probability)) +
      geom_point(alpha=0.55) +
      .cell_cycle_phase_layer(report_data$phase_boundaries) +
      scale_colour_gradientn(
        colours=report_data$cyclic_colours, limits=c(0, 2*pi),
        breaks=seq(0, 2*pi, pi/2), labels=report_data$theta_labels, name='theta') +
      scale_size_continuous(range=c(0.2, 2), limits=c(0, 1), guide='none') +
      .cell_cycle_origin_layer(report_data$origin) +
      coord_equal(xlim=report_data$pc_x_limits, ylim=report_data$pc_y_limits, expand=FALSE) +
      theme_bw() + labs(x='tricycle PC1', y='tricycle PC2')
    g_theta <- ggplot(sample_dt, aes(tricycle_theta / pi)) +
      geom_density(fill='steelblue', alpha=0.35) +
      scale_x_continuous(breaks=report_data$theta_breaks, labels=report_data$theta_labels, limits=c(0, 2)) +
      theme_bw() + labs(x=expression('tricycle theta / '*pi), y='density')
    print(g_projection + g_theta + plot_layout(widths=c(1.5, 1)))
    cat('\n\n')
  }
}


make_cell_cycle_phase_summary <- function(report_data) {
  phase_dt <- report_data$scores_qc[, .(n_cells=.N), by=.(sample_id, approximate_phase)]
  phase_dt[, percent_cells := 100 * n_cells / sum(n_cells), by=sample_id]
  cycling_phase_dt <- phase_dt[approximate_phase != 'G1/G0']
  cycling_phase_dt[, approximate_phase := droplevels(approximate_phase)]
  phase_plot <- ggplot(cycling_phase_dt, aes(percent_cells, sample_id, fill=approximate_phase)) +
    geom_quasirandom(
      shape=21, size=3, width=0.25, orientation='y', colour='grey20') +
    scale_y_discrete(limits=rev) +
    scale_x_log10(
      breaks=c(0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100),
      labels=label_number(suffix='%')) +
    scale_fill_manual(values=report_data$phase_colours) + theme_bw() +
    labs(x='percent of QC-passed cells (log scale)', y='sample',
         fill='approximate phase')

  phase_table <- dcast(cycling_phase_dt, sample_id ~ approximate_phase,
    value.var = 'percent_cells', fill = 0)
  list(plot = phase_plot, table = phase_table)
}


.cell_cycle_origin_layer <- function(origin) {
  geom_point(
    data = origin, aes(tricycle_pc1_origin, tricycle_pc2_origin),
    inherit.aes = FALSE, shape = 4, size = 4.5, stroke = 1.5, colour = 'black')
}


.cell_cycle_phase_layer <- function(phase_boundaries) {
  geom_segment(
    data = phase_boundaries, aes(x, y, xend = xend, yend = yend),
    inherit.aes = FALSE, linewidth = 0.5, linetype = 'dotted', colour = 'grey30')
}


.evenly_spaced_cell_cycle_rows <- function(cells_dt, maximum) {
  if (nrow(cells_dt) <= maximum) return(cells_dt)
  cells_dt[unique(round(seq(1, nrow(cells_dt), length.out = maximum)))]
}


format_cell_cycle_density_percentile <- function(x) {
  vapply(x, function(value) {
    paste0(format(value, trim=TRUE, scientific=FALSE, digits=15), '%')
  }, character(1))
}


cell_cycle_display_sample <- function(probabilities, seed, n_draws=100L) {
  stopifnot(length(probabilities) > 0L, all(is.finite(probabilities)),
            all(probabilities > 0 & probabilities <= 1),
            length(n_draws) == 1L, is.finite(n_draws),
            n_draws >= 1, n_draws == as.integer(n_draws))

  # Enlarge the display without altering relative weights or the fitted origin.
  display_probabilities <- probabilities / max(probabilities)
  expected <- sum(display_probabilities)
  retained <- withr::with_seed(seed, {
    best <- NULL
    best_error <- Inf
    for (i in seq_len(n_draws)) {
      draw <- runif(length(probabilities)) < display_probabilities
      error <- abs(sum(draw) - expected)
      if (error < best_error) {
        best <- draw
        best_error <- error
      }
    }
    best
  })
  list(retained=retained, probabilities=display_probabilities, expected_cells=expected)
}


cell_cycle_origin_caption <- function(origin) {
  if (origin$origin_method == 'manual') {
    sprintf('Manual origin: (%s, %s)', signif(origin$tricycle_pc1_origin, 3),
            signif(origin$tricycle_pc2_origin, 3))
  } else {
    sprintf('%s line fraction: %s; origin: (%s, %s)',
            if (origin$origin_method == 'ellipse_auto') 'Automatic' else 'Selected',
            signif(origin$origin_line_fraction, 3), signif(origin$tricycle_pc1_origin, 3),
            signif(origin$tricycle_pc2_origin, 3))
  }
}


estimate_cell_cycle_blob <- function(diagnostics_f, out_f) {
  # Fit once for origin selection and the report's reference ellipse.
  points <- as.matrix(fread(diagnostics_f, select=c('tricycle_pc1', 'tricycle_pc2')))
  summary <- data.table(status='unavailable', reason='', n_candidate_cells=nrow(points),
    coverage=0.99, alpha=0.5, centre_pc1=NA_real_, centre_pc2=NA_real_,
    cov_11=NA_real_, cov_12=NA_real_, cov_22=NA_real_,
    n_outside=NA_integer_, prop_outside=NA_real_,
    outside_mean_pc1=NA_real_, outside_mean_pc2=NA_real_,
    robustbase_version=as.character(utils::packageVersion('robustbase')))

  # A degenerate projection can still be reported without an ellipse.
  tryCatch({
    if (nrow(points) < 5L || any(!is.finite(points))) {
      stop('At least five finite QC-passed projections are required')
    }
    fit <- robustbase::covMcd(points, alpha=0.5, nsamp='deterministic',
                             raw.only=FALSE, use.correction=TRUE)
    if (any(!is.finite(fit$center)) || any(!is.finite(fit$cov)) ||
        any(eigen(fit$cov, symmetric=TRUE)$values <= 0)) {
      stop('The robust covariance is singular or non-finite')
    }
    distance2 <- mahalanobis(points, fit$center, fit$cov)
    outside <- distance2 > qchisq(0.99, df=2)
    summary[, `:=`(status='ok', centre_pc1=fit$center[1], centre_pc2=fit$center[2],
      cov_11=fit$cov[1, 1], cov_12=fit$cov[1, 2], cov_22=fit$cov[2, 2],
      n_outside=sum(outside), prop_outside=mean(outside))]
    if (any(outside)) {
      outside_mean <- colMeans(points[outside, , drop=FALSE])
      summary[, `:=`(outside_mean_pc1=outside_mean[1], outside_mean_pc2=outside_mean[2])]
    }
  }, error=function(e) {
    summary[, reason := conditionMessage(e)]
    warning('Cell-cycle blob reference unavailable: ', conditionMessage(e), call.=FALSE)
  })
  fwrite(summary, out_f)
  invisible(summary)
}


read_cell_cycle_blob <- function(blob_f) {
  summary <- fread(blob_f)
  stopifnot(nrow(summary) == 1L)
  ellipse <- data.table(tricycle_pc1=numeric(), tricycle_pc2=numeric())
  label <- '99% robust ellipse unavailable'
  if (summary$status == 'ok') {
    covariance <- matrix(c(summary$cov_11, summary$cov_12,
                           summary$cov_12, summary$cov_22), 2)
    angle <- seq(0, 2*pi, length.out=361)
    points <- sqrt(qchisq(summary$coverage, df=2)) *
      cbind(cos(angle), sin(angle)) %*% chol(covariance)
    points <- sweep(points, 2, c(summary$centre_pc1, summary$centre_pc2), '+')
    ellipse <- data.table(tricycle_pc1=points[, 1], tricycle_pc2=points[, 2])
    label <- sprintf('99%% robust ellipse: %.1f%% of cells outside', 100*summary$prop_outside)
  }
  list(summary=summary, ellipse=ellipse, label=label)
}


plot_cell_cycle_origin <- function(diagnostics, origin, sensitivity, blob, phase_boundaries) {
  expand_limits <- function(values) {
    limits <- range(values, finite=TRUE)
    limits + c(-1, 1) * max(diff(limits) * 0.05, 0.01)
  }
  full_x <- expand_limits(c(diagnostics$tricycle_pc1, origin$tricycle_pc1_origin))
  full_y <- expand_limits(c(diagnostics$tricycle_pc2, origin$tricycle_pc2_origin))
  zoom_x <- if (nrow(blob$ellipse)) {
    expand_limits(c(blob$ellipse$tricycle_pc1, origin$tricycle_pc1_origin))
  } else full_x
  zoom_y <- if (nrow(blob$ellipse)) {
    expand_limits(c(blob$ellipse$tricycle_pc2, origin$tricycle_pc2_origin))
  } else full_y
  g <- ggplot(diagnostics, aes(tricycle_pc1, tricycle_pc2)) +
    geom_hex(aes(colour=after_stat(count), fill=after_scale(colour)),
             bins=180, linewidth=0) +
    geom_path(data=blob$ellipse, colour='white', linewidth=1.1) +
    geom_path(data=blob$ellipse, colour='#D55E00', linewidth=0.6) +
    geom_path(data=sensitivity[order(origin_line_fraction)],
              aes(tricycle_pc1_origin, tricycle_pc2_origin), inherit.aes=FALSE,
              colour='white', linewidth=1.3) +
    geom_path(data=sensitivity[order(origin_line_fraction)],
              aes(tricycle_pc1_origin, tricycle_pc2_origin), inherit.aes=FALSE,
              colour='grey50', linewidth=0.55) +
    geom_point(data=sensitivity[selected == FALSE],
               aes(tricycle_pc1_origin, tricycle_pc2_origin, fill=origin_line_fraction),
               inherit.aes=FALSE, shape=21, size=2.5, colour='grey20') +
    geom_segment(data=phase_boundaries, aes(x, y, xend=xend, yend=yend),
                 inherit.aes=FALSE, linewidth=0.4, linetype='dotted', colour='grey40') +
    geom_point(data=origin, aes(tricycle_pc1_origin, tricycle_pc2_origin),
               inherit.aes=FALSE, shape=4, size=4.5, stroke=1.5) +
    scale_colour_viridis_c(trans='log10', name='cell count',
      guide=guide_colourbar(order=1, direction='horizontal')) +
    scale_fill_gradient(low='grey25', high='grey85', limits=c(0, 1),
      breaks=scales::pretty_breaks(),
      labels=function(x) vapply(x, function(value) {
        format(value, trim=TRUE, scientific=FALSE)
      }, character(1)),
      name='origin line fraction',
      guide=guide_colourbar(order=2, direction='horizontal')) +
    scale_x_continuous(breaks=scales::pretty_breaks()) +
    scale_y_continuous(breaks=scales::pretty_breaks()) +
    theme_bw() + labs(x='tricycle PC1', y='tricycle PC2')
  full <- g + coord_equal(xlim=full_x, ylim=full_y, expand=FALSE) +
    labs(title='All QC-passed cells', subtitle=blob$label,
         caption=cell_cycle_origin_caption(origin))
  zoom <- g + coord_equal(xlim=zoom_x, ylim=zoom_y, expand=FALSE) +
    labs(title='Origin close-up')
  full + zoom + patchwork::plot_layout(guides='collect') &
    theme(legend.position='bottom')
}


# Read marker expression captured from each full run-level count matrix during
# tricycle projection, then retain the QC-passed cells used by the report.
read_cell_cycle_marker_expression <- function(marker_expression_f, scores_qc) {
  marker_expression <- fread(marker_expression_f)
  annotation <- cell_cycle_marker_annotation()
  required <- c("cell_id", as.character(annotation$gene))
  missing <- setdiff(required, names(marker_expression))
  if (length(missing)) {
    stop("marker-expression table is missing: ", paste(missing, collapse = ", "))
  }
  marker_expression <- merge(
    scores_qc[, .(
      cell_id, sample_id, tricycle_pc1, tricycle_pc2, tricycle_theta
    )],
    marker_expression[, ..required],
    by = "cell_id"
  )
  melt(
    marker_expression,
    id.vars = c(
      "cell_id", "sample_id", "tricycle_pc1", "tricycle_pc2",
      "tricycle_theta"
    ),
    measure.vars = as.character(annotation$gene),
    variable.name = "gene", value.name = "log_normalized_expression"
  )
}


cell_cycle_marker_annotation <- function() {
  data.table(
    gene = c(
      "CDKN1A",
      "PCNA", "MCM4", "MCM5", "TYMS",
      "TOP2A", "CCNB1", "CDK1", "TPX2", "SMC2", "MKI67", "UBE2C"
    ),
    marker_phase = factor(
      rep(c("checkpoint", "S", "G2/M"), c(1L, 4L, 7L)),
      levels = c("checkpoint", "S", "G2/M")
    )
  )
}


assign_approximate_cell_cycle_phase <- function(theta) {
  phase <- fifelse(theta < 0.25 * pi | theta >= 1.75 * pi, "G1/G0",
    fifelse(theta < 0.5 * pi, "G1/S",
      fifelse(theta < pi, "S", "G2/M")))
  factor(phase, levels = c("G1/G0", "G1/S", "S", "G2/M"))
}


smooth_periodic_expression <- function(expression_dt, n_bins = 72L, span = 0.5) {
  if (!nrow(expression_dt)) {
    return(list(binned = data.table(), curves = data.table()))
  }

  # Average within angular bins, then repeat a cycle on either side for smoothing.
  binned <- copy(expression_dt)[
    , theta_bin := pmin(n_bins - 1L, floor(tricycle_theta / (2 * pi) * n_bins))
  ][
    , .(expression = mean(log_normalized_expression), n_cells = .N),
    by = .(gene, theta_bin)
  ][
    , theta := (theta_bin + 0.5) * 2 * pi / n_bins
  ]
  curves <- binned[, {
    augmented <- rbindlist(list(
      .SD[, .(theta = theta - 2 * pi, expression, n_cells)],
      .SD[, .(theta, expression, n_cells)],
      .SD[, .(theta = theta + 2 * pi, expression, n_cells)]
    ))
    grid <- seq(0, 2 * pi, length.out = 241L)
    if (uniqueN(.SD$theta) < 4L) {
      fitted <- rep(NA_real_, length(grid))
    } else {
      model <- tryCatch(
        loess(
          expression ~ theta, data = augmented, weights = n_cells,

          # Three copies provide circular neighbours. Divide the requested
          # span by three so it retains its usual one-cycle interpretation.
          span = span / 3, degree = 1, surface = "direct"
        ),
        error = function(e) NULL
      )
      fitted <- if (is.null(model)) {
        rep(NA_real_, length(grid))
      } else {
        # Local fits can overshoot in sparse angular regions. Expression is
        # non-negative, so constrain the diagnostic smoother to that domain.
        pmax(0, as.numeric(predict(model, grid)))
      }
    }
    .(
      theta = grid,
      mean_expression = fitted
    )
  }, by = gene][
    is.finite(mean_expression)
  ][
    , curve_peak := max(mean_expression, na.rm = TRUE),
    by = gene
  ][
    , mean_scaled := fifelse(
      curve_peak > 0, mean_expression / curve_peak, 0
    )
  ]

  list(binned = binned, curves = curves)
}
