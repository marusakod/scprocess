suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(rhdf5)
})


.tricycle_log <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                  paste0(..., collapse = "")))
  flush.console()
}


.read_10x_csc <- function(path, run) {
  # The current DropletUtils 10x writer stores each sparse array in one large
  # compressed HDF5 chunk. Incremental BPCells reads repeatedly decompress that
  # chunk, so read each array exactly once and construct the in-memory CSC matrix.
  .tricycle_log("Reading CSC arrays from ", path)
  data <- h5read(path, "matrix/data")
  indices <- h5read(path, "matrix/indices")
  indptr <- h5read(path, "matrix/indptr")
  shape <- h5read(path, "matrix/shape")
  barcodes <- as.character(h5read(path, "matrix/barcodes"))
  features <- as.character(h5read(path, "matrix/features/name"))
  if (length(shape) != 2L || length(indptr) != shape[2L] + 1L ||
      length(data) != length(indices)) {
    stop("Invalid 10x CSC matrix structure in ", path)
  }
  prefix <- paste0(run, ":")
  if (!all(startsWith(barcodes, prefix))) {
    barcodes <- paste(run, barcodes, sep = ":")
  }
  mat <- new(
    "dgCMatrix",
    i = as.integer(indices),
    p = as.integer(indptr),
    x = as.numeric(data),
    Dim = as.integer(shape),
    Dimnames = list(features, barcodes)
  )
  if (anyDuplicated(colnames(mat))) stop("Input matrix contains duplicate cell IDs")
  .tricycle_log(
    "Loaded ", nrow(mat), " features x ", ncol(mat), " cells with ",
    length(mat@x), " non-zero entries"
  )
  mat
}


.reference_rotation <- function(species, gene_id_type) {
  data_env <- new.env(parent = emptyenv())
  data("neuroRef", package = "tricycle", envir = data_env)
  reference <- data_env[["neuroRef"]]
  rotation <- as.matrix(reference[, seq_len(2L)])
  if (species == "human") {
    rownames(rotation) <- reference$SYMBOL
  } else if (gene_id_type == "ENSEMBL") {
    rownames(rotation) <- sub("\\.[0-9]+$", "", reference$ensembl)
  } else {
    rownames(rotation) <- reference$symbol
  }
  colnames(rotation) <- c("PC1", "PC2")
  rotation[!is.na(rownames(rotation)) & !duplicated(rownames(rotation)), , drop = FALSE]
}


.feature_keys <- function(features, species, rowdata_f) {
  base_ids <- sub("_[SUA]$", "", features)
  base_ids <- sub("\\.[0-9]+$", "", base_ids)
  ensembl_pattern <- if (species == "mouse") "^ENSMUSG[0-9]+" else "^ENSG[0-9]+"
  gene_id_type <- if (mean(grepl(ensembl_pattern, base_ids)) >= 0.5) "ENSEMBL" else "SYMBOL"

  if (gene_id_type == "SYMBOL") {
    return(list(keys = base_ids, type = gene_id_type))
  }
  if (species == "mouse") {
    return(list(keys = base_ids, type = gene_id_type))
  }

  rowdata <- fread(rowdata_f)
  required <- c("ensembl_id", "symbol")
  missing <- setdiff(required, names(rowdata))
  if (length(missing)) stop("rowdata is missing: ", paste(missing, collapse = ", "))
  rowdata[, ensembl_id := sub("\\.[0-9]+$", "", ensembl_id)]
  rowdata <- rowdata[
    !is.na(ensembl_id) & !is.na(symbol) & nzchar(symbol) & !duplicated(ensembl_id)
  ]
  list(keys = rowdata$symbol[match(base_ids, rowdata$ensembl_id)], type = gene_id_type)
}


.project_selected_cells <- function(
    mat, selected, species, rowdata_f, min_reference_genes) {
  .tricycle_log("Calculating library sizes for ", length(selected), " cells")
  cell_idx <- match(selected, colnames(mat))
  if (anyNA(cell_idx)) stop("Some selected cells are absent from the tricycle matrices")
  # Matrix::colSums traverses the in-memory CSC arrays once. Calculate all
  # columns before indexing to avoid copying most of a run-sized sparse matrix.
  lib_size <- Matrix::colSums(mat)[cell_idx]
  if (any(!is.finite(lib_size)) || any(lib_size <= 0)) {
    stop("Tricycle input contains a non-positive library size")
  }

  .tricycle_log("Matching count features to the tricycle reference")
  feature_info <- .feature_keys(rownames(mat), species, rowdata_f)
  rotation <- .reference_rotation(species, feature_info$type)
  n_reference_total <- nrow(rotation)
  matched_rows <- which(!is.na(feature_info$keys) & feature_info$keys %in% rownames(rotation))
  matched_keys <- feature_info$keys[matched_rows]
  unique_keys <- unique(matched_keys)
  if (length(unique_keys) < min_reference_genes) {
    stop(
      "Only ", length(unique_keys), " tricycle reference genes matched; at least ",
      min_reference_genes, " are required"
    )
  }

  # Materialize only the small fixed-reference subset, then sum S/U/A rows.
  .tricycle_log(
    "Materializing ", length(matched_rows), " matched count rows (",
    length(unique_keys), " unique reference genes)"
  )
  reference_counts <- as(mat[matched_rows, cell_idx], "dgCMatrix")
  collapse <- sparseMatrix(
    i = match(matched_keys, unique_keys),
    j = seq_along(matched_keys),
    x = 1,
    dims = c(length(unique_keys), length(matched_rows))
  )
  reference_counts <- collapse %*% reference_counts
  rownames(reference_counts) <- unique_keys
  colnames(reference_counts) <- selected

  .tricycle_log("Normalizing reference-gene counts")
  normalized <- reference_counts %*% Diagonal(x = 1e4 / lib_size)
  normalized@x <- log1p(normalized@x)
  rotation <- rotation[rownames(normalized), , drop = FALSE]

  # Do not centre within this run. The aggregation rule subtracts one common
  # project mean from these raw coordinates. By linearity, that is equivalent
  # to centring the pooled normalized expression before applying the rotation.
  .tricycle_log("Calculating uncentred tricycle reference projection")
  embedding <- as.matrix(crossprod(normalized, rotation))
  rownames(embedding) <- selected
  colnames(embedding) <- c("PC1", "PC2")
  if (any(!is.finite(embedding))) stop("tricycle returned non-finite coordinates")

  list(embedding = embedding, gene_id_type = feature_info$type,
       n_reference_genes = length(unique_keys),
       n_reference_total = n_reference_total,
       n_reference_missing = n_reference_total - length(unique_keys),
       n_duplicate_feature_mappings = length(matched_keys) - length(unique_keys))
}


#' Estimate uncentred fixed-reference tricycle coordinates for one physical run.
estimate_run_tricycle <- function(
    counts_h5_f, coldata_f, rowdata_f, run_column, run_value, species,
    out_scores_f, out_summary_f,
    min_reference_genes = 100L) {
  .tricycle_log("Reading cell metadata for ", run_column, " ", run_value)
  coldata <- fread(coldata_f)
  required <- c("cell_id", "sample_id", "keep", "scdbl_class", run_column)
  missing <- setdiff(required, names(coldata))
  if (length(missing)) stop("coldata is missing: ", paste(missing, collapse = ", "))
  if (anyDuplicated(coldata$cell_id)) stop("coldata contains duplicate cell IDs")

  known_doublet <- coldata$scdbl_class == "doublet"
  if ("demux_class" %in% names(coldata)) {
    known_doublet <- known_doublet | coldata$demux_class == "doublet"
  }
  known_doublet[is.na(known_doublet)] <- FALSE
  eligible <- coldata$keep == TRUE | known_doublet
  eligible[is.na(eligible)] <- FALSE
  in_group <- !is.na(coldata[[run_column]]) & coldata[[run_column]] == run_value
  selected <- coldata[eligible & in_group, cell_id]
  if (!length(selected)) stop("No eligible cells found for ", run_column, " ", run_value)

  .tricycle_log("Loading one run-level count matrix")
  mat <- .read_10x_csc(counts_h5_f, run_value)
  selected <- selected[selected %in% colnames(mat)]
  if (!length(selected)) stop("No eligible cells were present in the supplied matrices")
  .tricycle_log("Selected ", length(selected), " eligible cells present in the count matrix")
  projection <- .project_selected_cells(
    mat, selected, species, rowdata_f, min_reference_genes
  )
  embedding <- projection$embedding

  scores <- data.table(
    cell_id = rownames(embedding),
    sample_id = coldata$sample_id[match(rownames(embedding), coldata$cell_id)],
    tricycle_projection_group = paste0("run:", run_value),
    tricycle_raw_pc1 = embedding[, 1L],
    tricycle_raw_pc2 = embedding[, 2L]
  )
  .tricycle_log("Writing ", nrow(scores), " uncentred cell projections")
  dir.create(dirname(out_scores_f), recursive = TRUE, showWarnings = FALSE)
  fwrite(scores, out_scores_f)

  summary <- data.table(
    projection_group = run_value,
    projection_group_column = run_column,
    centring_scope = "project_aggregation",
    species = species,
    gene_id_type = projection$gene_id_type,
    n_cells = nrow(embedding),
    n_output_cells = nrow(scores),
    n_reference_genes = projection$n_reference_genes,
    n_reference_total = projection$n_reference_total,
    n_reference_missing = projection$n_reference_missing,
    n_duplicate_feature_mappings = projection$n_duplicate_feature_mappings,
    reference = "tricycle::neuroRef",
    tricycle_version = as.character(packageVersion("tricycle"))
  )
  fwrite(summary, out_summary_f)
  .tricycle_log("Finished ", run_column, " ", run_value)
  invisible(scores)
}
