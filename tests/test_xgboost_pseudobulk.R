# Report logCPMs use the stored all-gene library sizes, keep label names with
# '/', and drop label-batch rows below min_cells.
suppressPackageStartupMessages(library(data.table))
source('scripts/utils.R')
source('scripts/train_xgboost.R')
source('scripts/marker_genes.R')

pb_f = tempfile(fileext = '.h5ad')
counts = Matrix::Matrix(rbind(c(8, 6), c(5, 4), c(8, 7)), sparse = TRUE,
  dimnames = list(paste0('pb', 0:2), c('G2_ENSG2', 'G0_ENSG0')))
anndataR::AnnData(X = counts,
  obs = data.frame(label = c('Ependymal / ciliated', 'b', 'b'), batch = c('s1', 's1', 's2'),
    n_cells = c(2L, 2L, 1L), lib_size = c(32, 27, 24), row.names = rownames(counts)),
  var = data.frame(gene_id = colnames(counts), row.names = colnames(counts)))$write_h5ad(pb_f)

cpms_dt = calc_xgboost_pseudobulk_logcpms(pb_f, min_cells = 2)
stopifnot(nrow(cpms_dt) == 4,
  setequal(cpms_dt$cluster, c('b', 'Ependymal / ciliated')),
  nrow(cpms_dt[batch == 's2']) == 0,
  identical(sort(unique(cpms_dt$symbol)), c('G0', 'G2')),
  isTRUE(all.equal(cpms_dt[cluster == 'b' & symbol == 'G2', logcpm], log(5 / 27 * 1e6 + 10))))
# Depth augmentation sentence for the report, with and without thinned copies
stopifnot(identical(describe_depth_augmentation(list(downsample_copies = 1L,
    downsample_min_umis = 300L, downsample_max_umis = 3000L)),
    'This model was trained with depth augmentation (1 thinned copy per cell, 300 to 3,000 UMIs).'),
  grepl('retrain with depth augmentation', describe_depth_augmentation(list(downsample_copies = 0L))))
message('SUCCESS: classifier pseudobulk logCPMs use stored library sizes and keep labels; augmentation sentence')
