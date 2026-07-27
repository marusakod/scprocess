fit_markers_edger_bp <- function(counts, coldata, clusters, workers = 1L,
  block_size = 1000L) {
  clusters = as.character(clusters)
  coldata$cluster = as.character(coldata$cluster)
  coldata$sample = as.character(coldata$sample)
  stopifnot(
    is.data.frame(coldata),
    ncol(counts) == nrow(coldata),
    identical(colnames(counts), as.character(coldata$group_id)),
    all(c("cluster", "sample") %in% names(coldata)),
    all(clusters %in% coldata$cluster)
  )

  audit = edger.bp::bp_check_versions(error = TRUE)
  versions = setNames(audit$installed, audit$package)
  message(
    "  edger.bp streamed backend: edger.bp ",
    as.character(utils::packageVersion("edger.bp")),
    "; ",
    paste(names(versions), versions, sep = " ", collapse = "; ")
  )

  group = factor(coldata$cluster, levels = unique(coldata$cluster))
  dispersion_design = stats::model.matrix(~ group)

  message("  preparing disk-backed counts")
  prepared = edger.bp::bp_prepare_dge(
    counts,
    design = dispersion_design,
    group = group,
    norm.method = "TMM",
    block.size = block_size,
    workers = workers,
    verbose = TRUE,
    min.count = 1
  )

  message("  estimating streamed dispersions")
  dispersion_fit = edger.bp::bp_estimate_disp_stream(
    prepared$counts,
    design = dispersion_design,
    group = group,
    lib.size = prepared$effective.lib.size,
    offset = prepared$offset,
    block.size = block_size,
    workers = workers,
    verbose = TRUE
  )
  dispersion = dispersion_fit$tagwise.dispersion
  if (is.null(dispersion)) {
    dispersion = dispersion_fit$trended.dispersion
  }
  if (is.null(dispersion)) {
    dispersion = dispersion_fit$common.dispersion
  }

  # edgeR reuses AveLogCPM calculated with the common dispersion when fitting
  # QL models. Obtain the same row statistic through the public blockwise LRT
  # wrapper, then reuse it for every one-vs-rest Treat fit.
  message("  calculating streamed average logCPM")
  ave_fit = edger.bp::bp_glm_lrt(
    prepared$counts,
    design = dispersion_design,
    dispersion = dispersion_fit$common.dispersion,
    offset = prepared$offset,
    block.size = block_size,
    workers = workers
  )
  ave_log_cpm = ave_fit$table$logCPM
  rm(ave_fit)

  message("  running one-vs-rest streamed Treat tests")
  tables = vector("list", length(clusters))
  for (i in seq_along(clusters)) {
    cluster = clusters[[i]]
    message("    ", cluster)
    selected = group == cluster
    design = stats::model.matrix(~ selected)
    fit = edger.bp::bp_glm_treat(
      prepared$counts,
      design = design,
      dispersion = dispersion,
      offset = prepared$offset,
      coef = "selectedTRUE",
      ave.log.cpm = ave_log_cpm,
      block.size = block_size,
      workers = workers
    )

    out = fit$table
    out$FDR = stats::p.adjust(out$PValue, method = "BH")
    out$gene_id = rownames(out)
    out$cluster = cluster
    rownames(out) = NULL
    tables[[i]] = out
  }

  do.call(rbind, tables)
}
