suppressPackageStartupMessages({
  library(BPCells)
  library(data.table)
  library(Matrix)
})


.tricycle_log <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                  paste0(..., collapse = "")))
  flush.console()
}


.open_and_combine_10x <- function(paths) {
  if (!length(paths)) stop("At least one counts HDF5 is required")
  if (is.null(names(paths)) || any(!nzchar(names(paths)))) {
    stop("counts_h5_fs must be named with the corresponding run IDs")
  }
  mats <- Map(function(path, run) {
    mat <- open_matrix_10x_hdf5(path)
    prefix <- paste0(run, ":")
    if (!all(startsWith(colnames(mat), prefix))) {
      colnames(mat) <- paste(run, colnames(mat), sep = ":")
    }
    mat
  }, paths, names(paths))
  features <- rownames(mats[[1L]])
  for (i in seq_along(mats)[-1L]) {
    if (!identical(features, rownames(mats[[i]]))) {
      stop("All tricycle input matrices must have identical features in identical order")
    }
  }
  mat <- if (length(mats) == 1L) mats[[1L]] else do.call(cbind, mats)
  if (anyDuplicated(colnames(mat))) stop("Input matrices contain duplicate cell IDs")
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
  lib_size <- BPCells::colSums(mat[, cell_idx])
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

  # Algebraically identical to tricycle's documented fixed-reference formula:
  # scale(t(log-expression), center=TRUE, scale=FALSE) %*% rotation.
  .tricycle_log("Projecting cells into the tricycle reference")
  gene_means <- Matrix::rowMeans(normalized)
  offset <- as.numeric(gene_means %*% rotation)
  embedding <- as.matrix(crossprod(normalized, rotation))
  embedding <- sweep(embedding, 2L, offset, "-")
  rownames(embedding) <- selected
  colnames(embedding) <- c("PC1", "PC2")
  if (any(!is.finite(embedding))) stop("tricycle returned non-finite coordinates")

  list(embedding = embedding, gene_id_type = feature_info$type,
       n_reference_genes = length(unique_keys),
       n_reference_total = n_reference_total,
       n_reference_missing = n_reference_total - length(unique_keys),
       n_duplicate_feature_mappings = length(matched_keys) - length(unique_keys))
}


#' Estimate fixed-reference tricycle coordinates for a metadata group.
#'
#' Biological samples use group_column="sample_id". The optional unassigned
#' mode exists only for known HTO/custom doublets that lack a sample assignment;
#' these are centred with eligible cells from their run and clearly recorded.
estimate_tricycle_group <- function(
    counts_h5_fs, coldata_f, rowdata_f, group_column, group_value, species,
    out_scores_f, out_summary_f, output_unassigned_only = FALSE,
    min_reference_genes = 100L) {
  .tricycle_log("Reading cell metadata for ", group_column, " ", group_value)
  coldata <- fread(coldata_f)
  required <- c("cell_id", "sample_id", "keep", "scdbl_class", group_column)
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
  in_group <- !is.na(coldata[[group_column]]) & coldata[[group_column]] == group_value
  selected <- coldata[eligible & in_group, cell_id]
  if (!length(selected)) stop("No eligible cells found for ", group_column, " ", group_value)

  .tricycle_log("Opening ", length(counts_h5_fs), " count matrix file(s)")
  mat <- .open_and_combine_10x(counts_h5_fs)
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
    tricycle_center_group = if (output_unassigned_only) {
      paste0("run:", group_value)
    } else {
      paste0("sample:", group_value)
    },
    tricycle_pc1 = embedding[, 1L],
    tricycle_pc2 = embedding[, 2L]
  )
  if (output_unassigned_only) {
    scores <- scores[is.na(sample_id) | sample_id == ""]
  }
  .tricycle_log("Writing ", nrow(scores), " cell scores")
  dir.create(dirname(out_scores_f), recursive = TRUE, showWarnings = FALSE)
  fwrite(scores, out_scores_f)

  summary <- data.table(
    projection_group = group_value,
    projection_group_column = group_column,
    centring_scope = if (output_unassigned_only) "run_fallback" else "biological_sample",
    output_unassigned_only = output_unassigned_only,
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
  .tricycle_log("Finished ", group_column, " ", group_value)
  invisible(scores)
}


#' Estimate fixed-reference tricycle coordinates for one biological sample.
estimate_sample_tricycle <- function(
    counts_h5_fs, coldata_f, rowdata_f, sample_id, species,
    out_scores_f, out_summary_f,
    min_reference_genes = 100L) {
  estimate_tricycle_group(
    counts_h5_fs = counts_h5_fs,
    coldata_f = coldata_f,
    rowdata_f = rowdata_f,
    group_column = "sample_id",
    group_value = sample_id,
    species = species,
    out_scores_f = out_scores_f,
    out_summary_f = out_summary_f,
    min_reference_genes = min_reference_genes
  )
}
