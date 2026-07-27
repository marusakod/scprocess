library(BPCells)
library(edgeR)
library(Matrix)

source("scripts/marker_genes_edger_bp.R")

set.seed(17)
n_genes = 80L
clusters = rep(c("A", "B", "C"), each = 4L)
samples = rep(paste0("s", seq_len(4L)), times = 3L)
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
  stringsAsFactors = FALSE
)

bp_dir = tempfile("scprocess-edger-bp-test-")
on.exit(unlink(bp_dir, recursive = TRUE), add = TRUE)
bp_counts = write_matrix_dir(
  convert_matrix_type(Matrix(counts, sparse = TRUE), "uint32_t"),
  dir = bp_dir
)
observed = fit_markers_edger_bp(
  bp_counts,
  coldata = coldata,
  clusters = c("A", "B", "C"),
  workers = 1L,
  block_size = 20L
)

group = factor(clusters, levels = unique(clusters))
dispersion_design = model.matrix(~ group)
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
  design = model.matrix(~ selected)
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
  identical(key, paste(expected$cluster, expected$gene_id)),
  isTRUE(all.equal(observed$logFC, expected$logFC, tolerance = 1e-8)),
  isTRUE(all.equal(observed$logCPM, expected$logCPM, tolerance = 1e-8)),
  isTRUE(all.equal(observed$PValue, expected$PValue, tolerance = 1e-8)),
  isTRUE(all.equal(observed$FDR, expected$FDR, tolerance = 1e-8))
)
