suppressPackageStartupMessages({
  library(data.table)
})


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
