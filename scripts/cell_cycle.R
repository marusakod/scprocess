suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(rhdf5)
})


# Read observed expression for a small marker panel from the project HVG matrix.
# The full sparse arrays are read once, but only the selected marker rows are
# materialized as a dense matrix.
read_cell_cycle_marker_expression <- function(
    counts_h5_f, scores_f, coldata_f, rowdata_f,
    marker_genes = c("TOP2A", "SMC2", "MKI67", "PCNA", "MCM5", "CCNB1", "UBE2C", "TYMS")) {
  features <- as.character(h5read(counts_h5_f, "matrix/features/name"))
  base_ids <- sub("_[SUA]$", "", sub("\\.[0-9]+$", "", features))

  rowdata <- fread(rowdata_f, select = c("ensembl_id", "symbol"))
  rowdata[, ensembl_id := sub("\\.[0-9]+$", "", ensembl_id)]
  rowdata <- rowdata[
    !is.na(ensembl_id) & !is.na(symbol) & nzchar(symbol) & !duplicated(ensembl_id)
  ]
  symbols <- rowdata$symbol[match(base_ids, rowdata$ensembl_id)]
  symbols[is.na(symbols)] <- base_ids[is.na(symbols)]
  marker_keys <- toupper(symbols)
  selected_rows <- which(marker_keys %in% toupper(marker_genes))
  if (!length(selected_rows)) return(data.table())

  data <- h5read(counts_h5_f, "matrix/data")
  indices <- h5read(counts_h5_f, "matrix/indices")
  indptr <- h5read(counts_h5_f, "matrix/indptr")
  shape <- h5read(counts_h5_f, "matrix/shape")
  barcodes <- as.character(h5read(counts_h5_f, "matrix/barcodes"))
  mat <- new(
    "dgCMatrix", i = as.integer(indices), p = as.integer(indptr),
    x = as.numeric(data), Dim = as.integer(shape),
    Dimnames = list(features, barcodes)
  )

  scores <- fread(scores_f, select = c("cell_id", "sample_id", "tricycle_theta"))
  coldata <- fread(coldata_f, select = c("cell_id", "keep", "sum"))
  cells <- merge(scores, coldata[keep == TRUE], by = "cell_id")
  cell_idx <- match(cells$cell_id, colnames(mat))
  cells <- cells[!is.na(cell_idx)]
  cell_idx <- cell_idx[!is.na(cell_idx)]
  if (!nrow(cells)) return(data.table())

  marker_counts <- as(mat[selected_rows, cell_idx, drop = FALSE], "dgCMatrix")
  selected_keys <- marker_keys[selected_rows]
  unique_keys <- unique(selected_keys)
  collapse <- sparseMatrix(
    i = match(selected_keys, unique_keys), j = seq_along(selected_keys), x = 1,
    dims = c(length(unique_keys), length(selected_keys))
  )
  marker_counts <- collapse %*% marker_counts
  lib_size <- cells$sum
  valid <- is.finite(lib_size) & lib_size > 0
  cells <- cells[valid]
  marker_counts <- marker_counts[, valid, drop = FALSE]
  normalized <- marker_counts %*% Diagonal(x = 1e4 / lib_size[valid])
  normalized@x <- log1p(normalized@x)

  expression <- as.data.table(t(as.matrix(normalized)))
  setnames(expression, unique_keys)
  expression[, cell_id := cells$cell_id]
  expression <- melt(
    expression, id.vars = "cell_id", variable.name = "gene",
    value.name = "log_normalized_expression"
  )
  merge(expression, cells[, .(cell_id, sample_id, tricycle_theta)], by = "cell_id")
}


assign_approximate_cell_cycle_phase <- function(theta) {
  phase <- fifelse(theta < 0.25 * pi | theta >= 1.75 * pi, "G1/G0",
    fifelse(theta < 0.5 * pi, "G1/S",
      fifelse(theta < pi, "S", "G2/M")))
  factor(phase, levels = c("G1/G0", "G1/S", "S", "G2/M"))
}


smooth_periodic_expression <- function(expression_dt, n_bins = 72L, span = 0.25) {
  if (!nrow(expression_dt)) return(data.table())
  binned <- copy(expression_dt)[
    , theta_bin := pmin(n_bins - 1L, floor(tricycle_theta / (2 * pi) * n_bins))
  ][
    , .(expression = mean(log_normalized_expression), n_cells = .N),
    by = .(gene, theta_bin)
  ][
    , theta := (theta_bin + 0.5) * 2 * pi / n_bins
  ]
  fits <- binned[, {
    augmented <- rbindlist(list(
      .SD[, .(theta = theta - 2 * pi, expression, n_cells)],
      .SD[, .(theta, expression, n_cells)],
      .SD[, .(theta = theta + 2 * pi, expression, n_cells)]
    ))
    grid <- seq(0, 2 * pi, length.out = 241L)
    model <- loess(expression ~ theta, data = augmented, weights = n_cells,
                   span = span, degree = 2, surface = "direct")
    .(theta = grid, fitted_expression = as.numeric(predict(model, grid)))
  }, by = gene]
  list(binned = binned, fits = fits)
}
