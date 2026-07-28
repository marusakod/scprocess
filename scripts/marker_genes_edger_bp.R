.bp_row_blocks <- function(n, block_size) {
  starts = seq.int(1L, n, by = block_size)
  lapply(starts, function(start) seq.int(start, min(start + block_size - 1L, n)))
}


.bp_lapply <- function(x, workers, fun) {
  if (workers > 1L && length(x) > 1L) {
    parallel::mclapply(x, fun, mc.cores = workers, mc.preschedule = TRUE)
  } else {
    lapply(x, fun)
  }
}


.h5ad_path <- function(entry) {
  if (is.character(entry)) entry else entry[["path"]]
}


load_gene_biotypes_bpcells <- function(gtf_dt_f) {
  biotypes_dt = data.table::fread(
    gtf_dt_f,
    select = c("gene_id", "symbol", "gene_type")
  )
  if (anyDuplicated(biotypes_dt$gene_id)) {
    stop("Gene metadata contains duplicate gene IDs.", call. = FALSE)
  }
  biotypes_dt
}


build_pseudobulk_cache_bpcells <- function(pb_dir, integration_f, h5ads_yaml_f,
  sel_res, batch_var, min_cl_size, n_cores = 8L, bp_tmpdir) {
  stopifnot(
    is.character(pb_dir), length(pb_dir) == 1L,
    is.character(bp_tmpdir), length(bp_tmpdir) == 1L, dir.exists(bp_tmpdir),
    !dir.exists(pb_dir)
  )

  cl_var = paste0("RNA_snn_res.", sel_res)
  message("    loading integration columns")
  int_dt = data.table::fread(
    integration_f,
    select = c("cell_id", batch_var, cl_var)
  )
  data.table::setnames(int_dt, c(batch_var, cl_var), c("sample", "cluster"))
  int_dt = int_dt[
    !is.na(sample) & nzchar(as.character(sample)) & !is.na(cluster)
  ]
  int_dt[, `:=`(sample = as.character(sample), cluster = as.character(cluster))]

  message("    excluding tiny clusters")
  cl_ns = int_dt[, .N, by = cluster]
  keep_cls = cl_ns[N >= min_cl_size, cluster]
  int_dt = int_dt[cluster %in% keep_cls]
  stopifnot(nrow(int_dt) > 0L, length(keep_cls) > 1L)

  group_map = int_dt[, .(n_cells = .N), by = .(cluster, sample)]
  data.table::setorder(group_map, cluster, sample)
  group_map[, group_id := sprintf("pb%08d", .I)]
  int_dt = group_map[int_dt, on = .(cluster, sample)]
  data.table::setkey(int_dt, sample, cell_id)

  h5ad_paths = yaml::read_yaml(h5ads_yaml_f)
  h5ad_batches = names(h5ad_paths)
  if (!length(h5ad_batches)) stop("The H5AD YAML contains no entries.", call. = FALSE)

  chunk_dir = file.path(bp_tmpdir, "sample_chunks")
  dir.create(chunk_dir, recursive = TRUE)
  message("    pseudobulking ", length(h5ad_batches), " H5AD batches to disk (BPCells, ",
    n_cores, " workers)")

  build_one = function(batch_name) {
    entry = h5ad_paths[[batch_name]]
    h5ad_path = .h5ad_path(entry)
    message("      ", batch_name)
    if (!file.exists(h5ad_path) || file.size(h5ad_path) == 0L) {
      warning("skipping missing or empty H5AD for ", batch_name, ": ", h5ad_path)
      return(NULL)
    }

    mat = BPCells::open_matrix_anndata_hdf5(h5ad_path)
    ok_cells = intersect(colnames(mat), int_dt$cell_id)
    if (!length(ok_cells)) return(NULL)
    batch_int = int_dt[match(ok_cells, cell_id)]
    stopifnot(!anyNA(batch_int$group_id))

    pb_one = BPCells::pseudobulk_matrix(
      mat[, ok_cells, drop = FALSE],
      batch_int$group_id,
      method = "sum"
    )
    pb_one = Matrix::Matrix(pb_one, sparse = TRUE)
    out_dir = file.path(
      chunk_dir,
      sprintf("batch_%08d", match(batch_name, h5ad_batches))
    )
    BPCells::write_matrix_dir(
      BPCells::convert_matrix_type(pb_one, "uint32_t"),
      dir = out_dir
    )
    out_dir
  }

  chunk_paths = .bp_lapply(h5ad_batches, n_cores, build_one)
  chunk_paths = Filter(Negate(is.null), chunk_paths)
  if (!length(chunk_paths)) stop("No pseudobulk chunks were created.", call. = FALSE)

  message("    combining disk-backed sample chunks")
  chunk_mats = lapply(chunk_paths, BPCells::open_matrix_dir)
  genes = rownames(chunk_mats[[1L]])
  same_genes = vapply(chunk_mats, function(x) identical(rownames(x), genes), logical(1))
  if (!all(same_genes)) stop("Gene rows differ between H5AD pseudobulk chunks.", call. = FALSE)
  combined = do.call(cbind, chunk_mats)
  if (anyDuplicated(colnames(combined))) {
    message("    summing pseudobulk groups represented in multiple H5AD batches")
    combined = BPCells::pseudobulk_matrix(
      combined,
      colnames(combined),
      method = "sum"
    )
    combined = Matrix::Matrix(combined, sparse = TRUE)
  }

  coldata = group_map[match(colnames(combined), group_id)]
  stopifnot(!anyNA(coldata$group_id), identical(coldata$group_id, colnames(combined)))

  dir.create(pb_dir, recursive = TRUE)
  counts_col_dir = file.path(pb_dir, "counts_col")
  counts_row_dir = file.path(pb_dir, "counts_row")
  message("    writing persistent column-ordered pseudobulk matrix")
  counts_col = BPCells::write_matrix_dir(
    BPCells::convert_matrix_type(combined, "uint32_t"),
    dir = counts_col_dir
  )
  message("    writing persistent row-ordered pseudobulk matrix")
  counts_row = edger.bp::bp_ensure_storage_order(
    counts_col,
    order = "row",
    outdir = counts_row_dir,
    tmpdir = bp_tmpdir
  )

  data.table::fwrite(
    coldata[, .(group_id, cluster, sample, n_cells)],
    file.path(pb_dir, "coldata.csv.gz")
  )
  data.table::fwrite(
    data.table::data.table(gene_id = rownames(counts_row)),
    file.path(pb_dir, "rowdata.csv.gz")
  )
  list(counts = counts_row, coldata = as.data.frame(coldata))
}


open_pseudobulk_cache_bpcells <- function(pb_dir, storage_order = c("row", "column")) {
  storage_order = match.arg(storage_order)
  matrix_dir = file.path(pb_dir, if (storage_order == "row") "counts_row" else "counts_col")
  coldata_f = file.path(pb_dir, "coldata.csv.gz")
  rowdata_f = file.path(pb_dir, "rowdata.csv.gz")
  stopifnot(dir.exists(matrix_dir), file.exists(coldata_f), file.exists(rowdata_f))

  counts = BPCells::open_matrix_dir(matrix_dir)
  coldata = data.table::fread(coldata_f)
  rowdata = data.table::fread(rowdata_f)
  stopifnot(
    identical(colnames(counts), coldata$group_id),
    identical(rownames(counts), rowdata$gene_id)
  )
  list(counts = counts, coldata = coldata, rowdata = rowdata)
}


filter_pseudobulk_columns_bpcells <- function(counts, coldata, min_cells) {
  stopifnot(ncol(counts) == nrow(coldata), identical(colnames(counts), coldata$group_id))
  lib_size = as.numeric(edger.bp::bp_library_sizes(counts))
  keep = coldata$n_cells >= min_cells & lib_size > 0

  for (cluster in unique(coldata$cluster)) {
    idx = which(coldata$cluster == cluster & keep)
    if (length(idx) > 2L) {
      low = scater::isOutlier(lib_size[idx], log = TRUE, type = "lower", nmads = 3)
      keep[idx[low]] = FALSE
    }
  }
  if (!any(keep)) stop("No pseudobulk columns remain after filtering.", call. = FALSE)
  list(keep = keep, lib_size = lib_size)
}


fit_markers_edger_bp <- function(counts, coldata, clusters, effective_lib_size,
  workers = 1L, block_size = 1000L) {
  clusters = as.character(clusters)
  coldata$cluster = as.character(coldata$cluster)
  coldata$sample = as.character(coldata$sample)
  stopifnot(
    is.data.frame(coldata),
    ncol(counts) == nrow(coldata),
    length(effective_lib_size) == ncol(counts),
    all(is.finite(effective_lib_size)),
    all(effective_lib_size > 0),
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
  offset = log(effective_lib_size)

  message("  estimating streamed dispersions")
  dispersion_fit = edger.bp::bp_estimate_disp_stream(
    counts,
    design = dispersion_design,
    group = group,
    lib.size = effective_lib_size,
    offset = offset,
    block.size = block_size,
    workers = workers,
    verbose = TRUE
  )
  dispersion = dispersion_fit$tagwise.dispersion
  if (is.null(dispersion)) dispersion = dispersion_fit$trended.dispersion
  if (is.null(dispersion)) dispersion = dispersion_fit$common.dispersion

  message("  calculating streamed average logCPM")
  ave_fit = edger.bp::bp_glm_lrt(
    counts,
    design = dispersion_design,
    dispersion = dispersion_fit$common.dispersion,
    offset = offset,
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
      counts,
      design = design,
      dispersion = dispersion,
      offset = offset,
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

  table = do.call(rbind, tables)
  table
}


calc_pseudobulk_variability_bpcells <- function(counts, coldata, effective_lib_size,
  biotypes_dt, pb_hvgs_f, min_hvg_cells = 2L, block_size = 1000L,
  workers = 1L, pseudo_count = 10) {
  use = coldata$n_cells >= min_hvg_cells
  if (sum(use) < 2L) stop("Fewer than two pseudobulks remain for variability ranking.", call. = FALSE)
  y = counts[, use, drop = FALSE]
  lib_size = effective_lib_size[use]
  blocks = .bp_row_blocks(nrow(y), block_size)

  message("  calculating blockwise TMM-logCPM variability")
  block_results = .bp_lapply(blocks, workers, function(rows) {
    block = as.matrix(y[rows, , drop = FALSE])
    logcpm = log(sweep(block, 2L, lib_size, "/") * 1e6 + pseudo_count)
    means = rowMeans(logcpm)
    variances = if (ncol(logcpm) > 1L) {
      rowSums((logcpm - means)^2) / (ncol(logcpm) - 1L)
    } else {
      rep(NA_real_, nrow(logcpm))
    }
    data.table::data.table(
      gene_id = rownames(y)[rows],
      logcpm_var = variances
    )
  })

  hvgs_dt = data.table::rbindlist(block_results)
  hvgs_dt = merge(hvgs_dt, biotypes_dt, by = "gene_id", all.x = TRUE)
  hvgs_dt[, variability_method := "tmm_logcpm_variance"]
  data.table::setorder(hvgs_dt, -logcpm_var)
  data.table::fwrite(hvgs_dt, pb_hvgs_f)
  hvgs_dt
}


calc_marker_context_bpcells <- function(counts, coldata, effective_lib_size,
  block_size = 1000L, workers = 1L, pseudo_count = 10) {
  stopifnot(
    ncol(counts) == nrow(coldata),
    length(effective_lib_size) == ncol(counts)
  )
  clusters = unique(coldata$cluster)
  cluster_cols = split(seq_len(nrow(coldata)), coldata$cluster)
  blocks = .bp_row_blocks(nrow(counts), block_size)

  message("  calculating blockwise marker expression context")
  block_results = .bp_lapply(blocks, workers, function(rows) {
    block = as.matrix(counts[rows, , drop = FALSE])
    gene_id = rownames(counts)[rows]
    data.table::rbindlist(lapply(clusters, function(cluster) {
      cols = cluster_cols[[cluster]]
      selected = block[, cols, drop = FALSE]
      logcpm = log(sweep(selected, 2L, effective_lib_size[cols], "/") * 1e6 + pseudo_count)
      data.table::data.table(
        gene_id = gene_id,
        cluster = cluster,
        logcpm.sel = rowMeans(logcpm),
        n_zero = rowSums(selected == 0),
        n_cl = length(cols)
      )
    }))
  })

  out = data.table::rbindlist(block_results)
  out[, weighted_logcpm_all := sum(n_cl * logcpm.sel), by = gene_id]
  out[, n_all := sum(n_cl), by = gene_id]
  out[, logcpm.other := data.table::fifelse(
    n_all > n_cl,
    (weighted_logcpm_all - n_cl * logcpm.sel) / (n_all - n_cl),
    NA_real_
  )]
  out[, c("weighted_logcpm_all", "n_all") := NULL]
  out[]
}


make_join_pseudobulks_bpcells <- function(integration_f, h5ads_yaml_f, batch_var,
  pb_dir, sel_res, min_cl_size, n_cores = 4L) {
  tmp_base = Sys.getenv("TMPDIR", unset = tempdir())
  bp_tmpdir = tempfile("scprocess-edger-bp-", tmpdir = tmp_base)
  dir.create(bp_tmpdir, recursive = TRUE)
  on.exit(unlink(bp_tmpdir, recursive = TRUE), add = TRUE)

  build_pseudobulk_cache_bpcells(
    pb_dir, integration_f, h5ads_yaml_f, sel_res, batch_var,
    min_cl_size = min_cl_size, n_cores = n_cores, bp_tmpdir = bp_tmpdir
  )
  invisible(pb_dir)
}


prepare_join_pseudobulks_bpcells <- function(pb_dir, prepared_coldata_f,
  prepared_rowdata_f, min_cells, n_cores = 4L, block_size = 1000L) {
  cache = open_pseudobulk_cache_bpcells(pb_dir, storage_order = "row")
  column_filter = filter_pseudobulk_columns_bpcells(
    cache$counts, cache$coldata, min_cells = min_cells
  )
  keep = column_filter$keep
  counts = cache$counts[, keep, drop = FALSE]
  coldata = cache$coldata[keep, , drop = FALSE]
  clusters = unique(coldata$cluster)
  if (length(clusters) < 2L) {
    stop("Fewer than two clusters remain after pseudobulk filtering.", call. = FALSE)
  }

  group = factor(coldata$cluster, levels = clusters)
  dispersion_design = stats::model.matrix(~ group)
  message("  preparing disk-backed counts (", n_cores, " prescheduled workers)")
  prepared = edger.bp::bp_prepare_dge(
    counts,
    design = dispersion_design,
    group = group,
    norm.method = "TMM",
    block.size = block_size,
    workers = n_cores,
    verbose = TRUE,
    min.count = 1
  )

  prepared_coldata = data.table::as.data.table(coldata)
  prepared_coldata[, `:=`(
    raw_lib_size = column_filter$lib_size[keep],
    lib_size = prepared$lib.size,
    norm_factor = prepared$norm.factors,
    effective_lib_size = prepared$effective.lib.size
  )]
  prepared_rowdata = data.table::data.table(gene_id = rownames(prepared$counts))
  data.table::fwrite(prepared_coldata, prepared_coldata_f)
  data.table::fwrite(prepared_rowdata, prepared_rowdata_f)
  invisible(list(coldata = prepared_coldata, rowdata = prepared_rowdata))
}


calculate_join_hvgs_bpcells <- function(pb_dir, prepared_coldata_f, pb_hvgs_f,
  biotypes_dt, n_cores = 4L, block_size = 1000L) {
  cache = open_pseudobulk_cache_bpcells(pb_dir, storage_order = "row")
  coldata = data.table::fread(prepared_coldata_f)
  cols = match(coldata$group_id, colnames(cache$counts))
  if (anyNA(cols)) {
    stop("Prepared pseudobulk columns are absent from the BPCells store.", call. = FALSE)
  }
  counts = cache$counts[, cols, drop = FALSE]
  calc_pseudobulk_variability_bpcells(
    counts = counts,
    coldata = coldata,
    effective_lib_size = coldata$effective_lib_size,
    biotypes_dt = biotypes_dt,
    pb_hvgs_f = pb_hvgs_f,
    block_size = block_size,
    workers = n_cores
  )
  invisible(pb_hvgs_f)
}


calculate_join_marker_genes_bpcells <- function(pb_dir, prepared_coldata_f,
  prepared_rowdata_f, mkrs_f, biotypes_dt, n_cores = 4L, block_size = 1000L) {
  cache = open_pseudobulk_cache_bpcells(pb_dir, storage_order = "row")
  coldata = data.table::fread(prepared_coldata_f)
  rowdata = data.table::fread(prepared_rowdata_f)
  cols = match(coldata$group_id, colnames(cache$counts))
  rows = match(rowdata$gene_id, rownames(cache$counts))
  if (anyNA(cols)) {
    stop("Prepared pseudobulk columns are absent from the BPCells store.", call. = FALSE)
  }
  if (anyNA(rows)) {
    stop("Prepared genes are absent from the BPCells store.", call. = FALSE)
  }
  counts = cache$counts[rows, cols, drop = FALSE]
  clusters = unique(coldata$cluster)

  markers = fit_markers_edger_bp(
    counts, coldata, clusters,
    effective_lib_size = coldata$effective_lib_size,
    workers = n_cores, block_size = block_size
  )

  context = calc_marker_context_bpcells(
    counts,
    coldata,
    coldata$effective_lib_size,
    block_size = block_size,
    workers = n_cores
  )
  markers = data.table::as.data.table(markers)
  markers = merge(markers, context, by = c("gene_id", "cluster"))
  markers = merge(markers, biotypes_dt, by = "gene_id", all.x = TRUE)
  data.table::setorder(markers, cluster, PValue)
  data.table::fwrite(markers, mkrs_f)
  invisible(markers)
}


calculate_marker_genes_bpcells <- function(integration_f, h5ads_yaml_f, batch_var,
  pb_dir, mkrs_f, pb_hvgs_f, biotypes_dt, sel_res, min_cl_size, min_cells,
  n_cores = 4L, block_size = 1000L) {
  prep_dir = tempfile("scprocess-edger-bp-prepared-")
  dir.create(prep_dir)
  on.exit(unlink(prep_dir, recursive = TRUE), add = TRUE)
  prepared_coldata_f = file.path(prep_dir, "coldata.csv.gz")
  prepared_rowdata_f = file.path(prep_dir, "rowdata.csv.gz")

  make_join_pseudobulks_bpcells(
    integration_f, h5ads_yaml_f, batch_var, pb_dir, sel_res, min_cl_size, n_cores
  )
  prepare_join_pseudobulks_bpcells(
    pb_dir, prepared_coldata_f, prepared_rowdata_f, min_cells, n_cores, block_size
  )
  calculate_join_hvgs_bpcells(
    pb_dir, prepared_coldata_f, pb_hvgs_f, biotypes_dt, n_cores, block_size
  )
  calculate_join_marker_genes_bpcells(
    pb_dir, prepared_coldata_f, prepared_rowdata_f, mkrs_f, biotypes_dt,
    n_cores, block_size
  )
  invisible(NULL)
}


make_report_logcpms_bpcells <- function(pb_dir, prepared_coldata_f, mkrs_f,
  pb_hvgs_f, out_f,
  min_cpm_mkr, not_ok_re = "(lincRNA|lncRNA|pseudogene|antisense)",
  pseudo_count = 10) {
  cache = open_pseudobulk_cache_bpcells(pb_dir, storage_order = "row")
  coldata = data.table::fread(prepared_coldata_f)
  counts = cache$counts[, match(coldata$group_id, colnames(cache$counts)), drop = FALSE]

  hvgs = data.table::fread(pb_hvgs_f)
  markers = data.table::fread(mkrs_f)
  top_hvgs = hvgs[!stringr::str_detect(gene_type, not_ok_re)][order(-logcpm_var)][1:100]
  top_markers = markers[
    !stringr::str_detect(gene_type, not_ok_re) &
      logcpm.sel >= log(min_cpm_mkr + 1)
  ] %>% get_top_markers()
  gene_ids = unique(c(top_hvgs$gene_id, top_markers$gene_id))
  rows = match(gene_ids, rownames(counts))
  rows = rows[!is.na(rows)]
  block = as.matrix(counts[rows, , drop = FALSE])
  logcpm = log(sweep(block, 2L, coldata$effective_lib_size, "/") * 1e6 + pseudo_count)

  out = data.table::as.data.table(logcpm, keep.rownames = "gene_id")
  out = data.table::melt(
    out,
    id.vars = "gene_id",
    variable.name = "group_id",
    value.name = "logcpm"
  )
  out = merge(
    out,
    coldata[, .(group_id, cluster, sample_id = sample, n_cells)],
    by = "group_id"
  )
  data.table::fwrite(out, out_f)
  invisible(out)
}
