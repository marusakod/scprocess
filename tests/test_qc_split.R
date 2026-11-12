#!/usr/bin/env Rscript
# QC filtering works from the table of cells with doublet calls alone, and a run excluded by
# cellbender (empty table) gives empty outputs

source('scripts/SampleQC.R')
source('scripts/utils.R')

tmp_dir = tempfile('qc-split-')
dir.create(tmp_dir)

# 8 cells of one sample: cell_ok passes; each other cell fails one criterion
n_cells = 8L
coldata_dt = data.table(
  cell_id         = paste0('pool_a:cell_', seq_len(n_cells)),
  pool_id         = 'pool_a',
  sample_id       = 'sample_1',
  demux_class     = 'singlet',
  valid_hto       = TRUE,
  scdbl_class     = 'singlet',
  sum             = 5000L,
  total           = 5000L,
  detected        = 2000L,
  subsets_mito_sum = 250L,
  logit_spliced   = qlogis(0.9))
coldata_dt[2, scdbl_class := 'doublet']
coldata_dt[3, demux_class := 'doublet']
coldata_dt[4, `:=`(sum = 300L, total = 300L)]
coldata_dt[5, subsets_mito_sum := 2000L]
coldata_dt[6, logit_spliced := qlogis(0.99)]
coldata_dt[7, detected := 100L]
coldata_dt[8, valid_hto := FALSE]

cuts_f = file.path(tmp_dir, 'cuts.csv')
fwrite(data.table(sample_id = 'sample_1', min_counts = 500, max_counts = NA_real_, min_feats = 300,
  max_feats = NA_real_, min_mito = 0.003, max_mito = 0.1, min_splice = 0.75, max_splice = 0.97,
  min_cells = 100), cuts_f)

res = .filter_qc(copy(coldata_dt), cuts_f, 'sample_id', 'pool_id', 'hto', 200, 100, 0.5)
keep_by_cell = setNames(res$coldata_dt$keep, res$coldata_dt$cell_id)
stopifnot(
  nrow(res$coldata_dt) == n_cells,
  sum(keep_by_cell) == 1L, keep_by_cell[['pool_a:cell_1']],
  # singlets only enter the QC table: cells 2 (scDblFinder) and 3 (demultiplexing) are excluded
  !(c('pool_a:cell_2', 'pool_a:cell_3') %in% res$qc_all$cell_id) %>% all())

# changing a threshold changes keep without any other input: a looser max_splice keeps cell 6
cuts_dt = fread(cuts_f)[, max_splice := 0.995]
loose_f = file.path(tmp_dir, 'cuts_loose.csv')
fwrite(cuts_dt, loose_f)
res_loose = .filter_qc(copy(coldata_dt), loose_f, 'sample_id', 'pool_id', 'hto', 200, 100, 0.5)
stopifnot(sum(res_loose$coldata_dt$keep) == 2L,
  res_loose$coldata_dt[cell_id == 'pool_a:cell_6', keep])

# main_filter_qc: reads the table from a file and writes the QC and cell tables
coldata_f = file.path(tmp_dir, 'coldata_with_doublets.csv.gz')
fwrite(coldata_dt, coldata_f)
qc_out_f = file.path(tmp_dir, 'qc.csv.gz')
coldata_out_f = file.path(tmp_dir, 'coldata.csv.gz')
main_filter_qc('pool_a', cuts_f, coldata_f, qc_out_f, coldata_out_f, 200, 100, 0.5,
  run_var = 'pool_id', demux_type = 'hto', batch_var = 'sample_id')
stopifnot(
  nrow(fread(coldata_out_f)) == n_cells,
  sum(fread(coldata_out_f)$keep) == 1L,
  # the QC table has the cells that are singlets with a valid hashtag: cells 1, 4, 5, 6 and 7
  nrow(fread(qc_out_f)) == 5L)

# a run excluded by cellbender has an empty table: empty outputs, no error
empty_f = file.path(tmp_dir, 'empty.csv.gz')
file.create(empty_f)
qc_empty_f = file.path(tmp_dir, 'qc_empty.csv.gz')
coldata_empty_f = file.path(tmp_dir, 'coldata_empty.csv.gz')
main_filter_qc('pool_b', cuts_f, empty_f, qc_empty_f, coldata_empty_f, 200, 100, 0.5,
  run_var = 'pool_id', demux_type = 'hto', batch_var = 'sample_id')
stopifnot(file.exists(qc_empty_f), file.size(qc_empty_f) == 0L,
  file.exists(coldata_empty_f), file.size(coldata_empty_f) == 0L)

message('test_qc_split.R passed')
