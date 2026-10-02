suppressPackageStartupMessages({
  library(data.table)
})


format_cell_cycle_target_prop <- function(x) {
  vapply(x, function(value) {
    paste0(format(100 * value, trim=TRUE, scientific=FALSE, digits=15), '%')
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
