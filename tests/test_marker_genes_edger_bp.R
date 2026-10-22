library(BPCells)
library(edgeR)
library(Matrix)
library(magrittr)

source("scripts/marker_genes_edger_bp.R")

set.seed(17)

parent_pid = Sys.getpid()
parallel_retry = .bp_lapply(seq_len(8L), 2L, function(i) {
  if (i == 3L && Sys.getpid() != parent_pid) stop("synthetic worker error")
  i
})
stopifnot(identical(unlist(parallel_retry), seq_len(8L)))

cache_test_dir = tempfile("scprocess-edger-bp-cache-")
dir.create(cache_test_dir)
on.exit(unlink(cache_test_dir, recursive = TRUE), add = TRUE)
h5ad_a = file.path(cache_test_dir, "batch_a.h5ad")
h5ad_b = file.path(cache_test_dir, "batch_b.h5ad")
h5ad_empty = file.path(cache_test_dir, "batch_empty.h5ad")
batch_a = matrix(
  c(1, 2, 3, 4, 5, 6),
  nrow = 3L,
  dimnames = list(paste0("g", 1:3), c("c1", "c2"))
)
batch_b = matrix(
  c(10, 20, 30, 40, 50, 60),
  nrow = 3L,
  dimnames = list(paste0("g", 1:3), c("c3", "c4"))
)
write_matrix_anndata_hdf5(
  convert_matrix_type(Matrix(batch_a, sparse = TRUE), "uint32_t"),
  h5ad_a
)
write_matrix_anndata_hdf5(
  convert_matrix_type(Matrix(batch_b[c(3L, 1L, 2L), , drop = FALSE], sparse = TRUE), "uint32_t"),
  h5ad_b
)
rhdf5::h5createFile(h5ad_empty)
integration_f = file.path(cache_test_dir, "integration.csv")
data.table::fwrite(
  data.table::data.table(
    cell_id = c("c1", "c2", "c3", "c4"),
    sample_id = c("s1", "s2", "s1", "s1"),
    project_id = c("P1", "P2", "P1", "P1"),
    RNA_snn_res.0.5 = c("A", "A", "A", "B")
  ),
  integration_f
)
h5ads_yaml_f = file.path(cache_test_dir, "h5ads.yaml")
yaml::write_yaml(
  list(
    sequencing_batch_a = list(path = h5ad_a),
    sequencing_batch_b = list(path = h5ad_b),
    sequencing_batch_empty = list(path = h5ad_empty)
  ),
  h5ads_yaml_f
)
cache_dir = file.path(cache_test_dir, "pseudobulk.bpcells")
cache_tmpdir = file.path(cache_test_dir, "tmp")
dir.create(cache_tmpdir)
gene_order_warnings = character()
cache = withCallingHandlers(build_pseudobulk_cache_bpcells(
  pb_dir = cache_dir,
  integration_f = integration_f,
  h5ads_yaml_f = h5ads_yaml_f,
  sel_res = "0.5",
  batch_var = "sample_id",
  min_cl_size = 1L,
  n_cores = 1L,
  bp_tmpdir = cache_tmpdir
), warning = function(w) {
  if (grepl("Same genes, different order", conditionMessage(w), fixed = TRUE)) {
    gene_order_warnings <<- c(gene_order_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
})
stopifnot(length(gene_order_warnings) == 1L,
  grepl("sequencing_batch_b", gene_order_warnings, fixed = TRUE))
# Real gene-set mismatches and ambiguous identifiers must still fail.
expect_gene_error = function(ids, pattern) {
  bad = batch_b
  rownames(bad) = ids
  error = tryCatch(.bp_align_gene_rows(list(batch_a, bad)), error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), fixed = TRUE))
}
expect_gene_error(c("g1", "g2", "g4"), "Gene sets differ")
expect_gene_error(c("g1", "g1", "g3"), "unique, nonmissing and nonblank")
expect_gene_error(c("g1", NA_character_, "g3"), "unique, nonmissing and nonblank")
expect_gene_error(c("g1", "", "g3"), "unique, nonmissing and nonblank")
stopifnot(identical(.bp_align_gene_rows(list(batch_a, batch_b)), list(batch_a, batch_b)))
cache_dense = .bp_as_matrix(cache$counts)
expected_cache = cbind(
  A_s1 = batch_a[, "c1"] + batch_b[, "c3"],
  A_s2 = batch_a[, "c2"],
  B_s1 = batch_b[, "c4"]
)
colnames(expected_cache) = cache$coldata$group_id
stopifnot(
  dir.exists(file.path(cache_dir, "counts_col")),
  dir.exists(file.path(cache_dir, "counts_row")),
  isTRUE(all.equal(unname(cache_dense), unname(expected_cache))),
  identical(cache$coldata$sample, c("s1", "s2", "s1")),
  identical(cache$coldata$project_id, c("P1", "P2", "P1")),
  identical(cache$coldata$cluster, c("A", "A", "B")),
  identical(cache$coldata$n_cells, c(2L, 1L, 1L))
)

cache_coldata = data.table::as.data.table(cache$coldata)
cache_coldata[, `:=`(
  raw_lib_size = colSums(cache_dense),
  lib_size = colSums(cache_dense),
  norm_factor = 1,
  effective_lib_size = colSums(cache_dense)
)]
cache_prepared_coldata_f = file.path(cache_test_dir, "prepared_coldata.csv.gz")
data.table::fwrite(cache_coldata, cache_prepared_coldata_f)
cache_hvgs_f = file.path(cache_test_dir, "hvgs.csv.gz")
data.table::fwrite(
  data.table::data.table(
    gene_id = rownames(cache_dense),
    symbol = rownames(cache_dense),
    gene_type = "protein_coding",
    logcpm_var = c(3, 2, 1)
  ),
  cache_hvgs_f
)
cache_markers_f = file.path(cache_test_dir, "markers.csv.gz")
data.table::fwrite(
  data.table::data.table(
    gene_id = rownames(cache_dense),
    cluster = "A",
    gene_type = "protein_coding",
    FDR = 0.001,
    logFC = 1,
    n_zero = 0L,
    n_cl = 2L,
    logcpm.sel = log(100)
  ),
  cache_markers_f
)
get_top_markers = function(input_mkrs, ...) input_mkrs
cache_plot_f = file.path(cache_test_dir, "plot_data.csv.gz")
make_report_logcpms_bpcells(
  pb_dir = cache_dir,
  prepared_coldata_f = cache_prepared_coldata_f,
  mkrs_f = cache_markers_f,
  pb_hvgs_f = cache_hvgs_f,
  out_f = cache_plot_f,
  min_cpm_mkr = 0
)
cache_plot = data.table::fread(cache_plot_f)
stopifnot(
  nrow(cache_plot) == nrow(cache_dense) * ncol(cache_dense),
  all(c("gene_id", "group_id", "logcpm", "cluster", "sample_id", "n_cells") %in%
    names(cache_plot))
)

n_genes = 80L
clusters = rep(c("A", "B", "C"), each = 4L)
samples = rep(paste0("s", seq_len(4L)), times = 3L)
projects = rep(c("P1", "P1", "P2", "P2"), times = 3L)
group_ids = sprintf("pb%08d", seq_along(clusters))

means = matrix(12, nrow = n_genes, ncol = length(clusters))
means[seq_len(8L), clusters == "A"] = 45
means[9:16, clusters == "B"] = 45
means[17:24, clusters == "C"] = 45
counts = matrix(
  rnbinom(length(means), mu = as.vector(means), size = 8),
  nrow = n_genes,
  dimnames = list(paste0("gene", seq_len(n_genes)), group_ids)
)
coldata = data.frame(
  group_id = group_ids,
  cluster = clusters,
  sample = samples,
  project_id = projects,
  stringsAsFactors = FALSE
)

bp_dir = tempfile("scprocess-edger-bp-test-")
on.exit(unlink(bp_dir, recursive = TRUE), add = TRUE)
bp_counts = write_matrix_dir(
  convert_matrix_type(Matrix(counts, sparse = TRUE), "uint32_t"),
  dir = bp_dir
)
group = factor(clusters, levels = unique(clusters))
project_id = factor(projects, levels = unique(projects))
for (adjust_project_id in c(FALSE, TRUE)) {
  dispersion_design = if (adjust_project_id) {
    model.matrix(~ project_id + group)
  } else {
    model.matrix(~ group)
  }
  prepared = edger.bp::bp_prepare_dge(
    bp_counts,
    design = dispersion_design,
    group = group,
    workers = 1L,
    block.size = 20L,
    min.count = 1
  )
  observed = fit_markers_edger_bp(
    prepared$counts,
    coldata = coldata,
    clusters = c("C", "A", "B"),
    effective_lib_size = prepared$effective.lib.size,
    workers = 1L,
    block_size = 20L,
    adjust_project_id = adjust_project_id
  )

  dge = DGEList(counts)
  keep = filterByExpr(
    dge,
    design = dispersion_design,
    group = group,
    min.count = 1
  )
  dge = dge[keep, , keep.lib.sizes = FALSE]
  dge = normLibSizes(dge, method = "TMM")
  dge = estimateDisp(dge, design = dispersion_design)

  expected = do.call(rbind, lapply(c("A", "B", "C"), function(cluster) {
    selected = group == cluster
    design = if (adjust_project_id) {
      model.matrix(~ project_id + selected)
    } else {
      model.matrix(~ selected)
    }
    fit = glmQLFit(
      dge,
      design = design,
      dispersion = dge$tagwise.dispersion,
      legacy = FALSE
    )
    table = glmTreat(fit, coef = "selectedTRUE")$table
    table$FDR = p.adjust(table$PValue, method = "BH")
    table$gene_id = rownames(table)
    table$cluster = cluster
    rownames(table) = NULL
    table
  }))

  key = paste(observed$cluster, observed$gene_id)
  expected = expected[match(key, paste(expected$cluster, expected$gene_id)), ]
  stopifnot(
    identical(unique(observed$cluster), c("A", "B", "C")),
    identical(key, paste(expected$cluster, expected$gene_id)),
    isTRUE(all.equal(observed$logFC, expected$logFC, tolerance = 1e-8)),
    isTRUE(all.equal(observed$logCPM, expected$logCPM, tolerance = 1e-8)),
    isTRUE(all.equal(observed$PValue, expected$PValue, tolerance = 1e-8)),
    isTRUE(all.equal(observed$FDR, expected$FDR, tolerance = 1e-8))
  )
}

variability_f = tempfile(fileext = ".csv.gz")
on.exit(unlink(variability_f), add = TRUE)
coldata$n_cells = 10L
effective_lib_size = colSums(counts)
biotypes = data.frame(
  gene_id = rownames(counts),
  symbol = rownames(counts),
  gene_type = "protein_coding"
)
variability = calc_pseudobulk_variability_bpcells(
  counts = bp_counts,
  coldata = coldata,
  effective_lib_size = effective_lib_size,
  biotypes_dt = biotypes,
  pb_hvgs_f = variability_f,
  block_size = 13L,
  workers = 1L
)
expected_logcpm = log(sweep(counts, 2L, effective_lib_size, "/") * 1e6 + 10)
expected_variability = apply(expected_logcpm, 1L, var)
observed_variability = variability$logcpm_var[
  match(names(expected_variability), variability$gene_id)
]
stopifnot(
  file.exists(variability_f),
  identical(unique(variability$variability_method), "tmm_logcpm_variance"),
  isTRUE(all.equal(
    unname(observed_variability),
    unname(expected_variability),
    tolerance = 1e-12
  ))
)

# Exercise the public wrapper with a real PSOCK pool and a persistent fit.
checkpoint_dir = tempfile("scprocess-marker-checkpoint-")
dir.create(checkpoint_dir)
edger.bp::bp_ensure_storage_order(bp_counts, order = "row",
  outdir = file.path(checkpoint_dir, "counts_row"), tmpdir = tempdir())
data.table::fwrite(coldata, file.path(checkpoint_dir, "coldata.csv.gz"))
data.table::fwrite(data.frame(gene_id = rownames(counts)),
  file.path(checkpoint_dir, "rowdata.csv.gz"))
prepared_cols = file.path(checkpoint_dir, "prepared_cols.csv.gz")
prepared_rows = file.path(checkpoint_dir, "prepared_rows.csv.gz")
coldata$effective_lib_size = effective_lib_size
data.table::fwrite(coldata, prepared_cols)
data.table::fwrite(data.frame(gene_id = rownames(counts)), prepared_rows)
markers_f = file.path(checkpoint_dir, "markers.csv.gz")
checkpoint_f = paste0(markers_f, ".fitted.rds")
original_fit = fit_markers_edger_bp
fit_calls = 0L
fit_markers_edger_bp = function(...) {
  fit_calls <<- fit_calls + 1L
  original_fit(...)
}
original_context = calc_marker_context_bpcells
context_calls = 0L
fail_context = FALSE
calc_marker_context_bpcells = function(..., workers, block_apply) {
  stopifnot(Sys.getpid() == parent_pid, workers == 1L,
    identical(block_apply, .bp_lapply))
  context_calls <<- context_calls + 1L
  if (fail_context && context_calls %% 2L == 0L) stop("synthetic summary failure")
  original_context(..., workers = workers, block_apply = block_apply)
}
run_markers = function(block_size = 20L) calculate_join_marker_genes_bpcells(
  checkpoint_dir, prepared_cols, prepared_rows, markers_f, biotypes,
  n_cores = 2L, block_size = block_size)
fail_context = TRUE
failed = tryCatch(run_markers(), error = identity)
stopifnot(inherits(failed, "error"), grepl("synthetic summary", failed$message),
  fit_calls == 1L, file.exists(checkpoint_f), !file.exists(markers_f))
fail_context = FALSE
recovered = run_markers()
stopifnot(fit_calls == 1L, file.exists(markers_f), nrow(recovered) == n_genes * 3L)
again = run_markers()
stopifnot(fit_calls == 1L, isTRUE(all.equal(recovered, again)))
changed = run_markers(block_size = 13L)
stopifnot(fit_calls == 2L, isTRUE(all.equal(recovered, changed, tolerance = 1e-8)))
saved = readRDS(checkpoint_f)
for (field in c("counts_md5", "coldata", "rowdata", "adjust_project_id",
                "versions", "implementations", "package_code", "package_binaries")) {
  stale = saved$signature
  stale[[field]] = NULL
  stopifnot(is.null(.read_marker_fit_checkpoint(checkpoint_f, stale)))
}
# Changing a helper must also change the provenance without changing inputs.
signature_before = .marker_fit_signature(checkpoint_dir, coldata,
  data.frame(gene_id = rownames(counts)), FALSE, 13L, 2L)
original_matrix = .bp_as_matrix
.bp_as_matrix = function(x) { return(.without_bpcells_dense_message(as.matrix(x))) }
signature_after = .marker_fit_signature(checkpoint_dir, coldata,
  data.frame(gene_id = rownames(counts)), FALSE, 13L, 2L)
stopifnot(!identical(signature_before, signature_after))
.bp_as_matrix = original_matrix
coldata$effective_lib_size[1L] = coldata$effective_lib_size[1L] * 1.1
data.table::fwrite(coldata, prepared_cols)
updated_metadata = run_markers(block_size = 13L)
stopifnot(fit_calls == 3L,
  !isTRUE(all.equal(updated_metadata$PValue, recovered$PValue)))
writeBin(charToRaw("interrupted checkpoint write"), checkpoint_f)
run_markers(block_size = 13L)
stopifnot(fit_calls == 4L, is.list(readRDS(checkpoint_f)),
  !length(list.files(checkpoint_dir, pattern = "^\\.marker-fit-", all.files = TRUE)))
fit_markers_edger_bp = original_fit
calc_marker_context_bpcells = original_context
unlink(checkpoint_dir, recursive = TRUE)
