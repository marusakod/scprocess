# source('scripts/pseudobulk.R')

suppressPackageStartupMessages({
  library('RColorBrewer')
  library("BiocParallel")
  library('circlize')
  library('magrittr')
  library('data.table')
  library('stringr')
  library('assertthat')
  library('viridis')
  library('scales')
  library('ggplot2')
  library('patchwork')
  library('forcats')
  library('readxl')

  library('future')
  library('SingleCellExperiment')

  library('scater')
  library('Seurat')

  library('ComplexHeatmap')
  library('seriation')
  library('purrr')
  library('xgboost')
  library('ggrepel')
  
  library('dplyr')
  library('rhdf5')
  library('Matrix')
  # library('Matrix.utils')
  library('parallel')
  library('edgeR')
  library('yaml')
})

# h5_paths_f = '/pmount/projects/site/pred/mus-brain-analysis/studies/test_multiplexed_project/output/test_hvg/hvg_paths_test_multiplexed_project_2025-01-01.csv'
# qc_stats_f = '/pmount/projects/site/pred/mus-brain-analysis/studies/test_multiplexed_project/output/test_qc/qc_pool_id_statistics_test_multiplexed_project_2025-01-01.csv'
# batch_lu_f = '/pmount/projects/site/pred/mus-brain-analysis/studies/test_multiplexed_project/output/test_pseudobulk/runs_to_batches_test_multiplexed_project_2025-01-01.csv'
# qc_f = '/pmount/projects/site/pred/mus-brain-analysis/studies/test_multiplexed_project/output/test_qc/coldata_dt_all_cells_test_multiplexed_project_2025-01-01.csv.gz'
# batch_var = 'pool_id'
# run_var = 'pool_id'
# subset_f = NULL
# subset_col = NULL 
# subset_str = NULL
# n_cores=2 


make_pb_cells <- function(sel_run, batch_lu_f, qc_stats_f, h5_paths_f, coldata_f, run_var, 
                          batch_var, pb_cells_f, subset_f = NULL, subset_col = NULL, subset_str = NULL) {
  # get batches corresponding to this run
  batch_lu    = fread(batch_lu_f)
  sel_batches = batch_lu[ run_var == sel_run ] %>% .$batch_var
  
  # remove all samples excluded after qc
  qc_stats_dt = fread(qc_stats_f)
  bad_var     = paste0("bad_", batch_var)
  bad_batches = qc_stats_dt[ get(bad_var) == TRUE ] %>% .[[batch_var]]
  
  # filter out bad ones
  ok_batches  = setdiff(sel_batches, bad_batches)
  # write an empty file if nothing lefts
  if (length(ok_batches) == 0) {
    file.create(pb_cells_f)
    return(NULL)
  }
  
  # get what we need
  h5_fs_dt  = fread(h5_paths_f)
  h5_f      = h5_fs_dt[get(run_var) == sel_run]$amb_filt_f %>% unique
  tmp_vars  = c(run_var, batch_var) %>% unique
  keep_dt   = fread(coldata_f) %>% .[ keep == TRUE] %>% .[get(batch_var) %in% ok_batches ] %>%
   .[, c(tmp_vars, "cell_id"), with = FALSE] %>%
    setkey("cell_id")
  keep_ids    = keep_dt %>% .$cell_id
  
  # get full alevin matrix
  message('    loading counts for ', sel_run)
  cells_mat   = .get_h5_mx(h5_f, sel_s = paste0(sel_run, ':')) %>% .sum_SUA
  
  # check is ok
  assert_that( all(keep_ids %in% colnames(cells_mat)) )
  cells_mat   = cells_mat[, keep_ids]
  
  # turn into sce
  sce             = SingleCellExperiment( assays = list(counts = cells_mat) )
  sce$cell_id     = colnames(sce)
  sce[[run_var]]  = sel_run
  if (batch_var != run_var) {
    sce[[batch_var]]  = cell_lu[ colnames(sce) ] %>% .[[ batch_var ]]
  }
  
  # subset if required
  if (!is.null(subset_f)) {
    message('    subsetting sce object')
    # unpack 
    subset_vals = str_split(subset_str, pattern = ',') %>% unlist
    subset_dt   = fread(subset_f)
    assert_that(subset_col %in% colnames(subset_dt))
    assert_that(all(c("cell_id", batch_var) %in% colnames(subset_dt)))
    
    # keep just cells in sel_b with selected labels
    subset_dt   = subset_dt %>%
      .[ get(batch_var) %in% ok_batches ] %>%
      .[ get(subset_col) %in% subset_vals ]
    
    # subset sce object
    sce       = sce[, subset_dt$cell_id]
  }
  
  message('   running aggregateData')
  agg_fn      = "sum"
  pb_mat      = .aggregateData_datatable(sce, by_vars = batch_var, fun = agg_fn)
  # colnames(pb) = sel_b
  
  # store
  message('  save')
  pb          = SingleCellExperiment( assays = list(counts = pb_mat) )
  saveRDS(pb, pb_cells_f)
  
  message('done!')
}


merge_pbs_cells <- function(cells_paths_f, rowdata_f, batch_var, pb_cells_f) {
  # get paths to ok pb cells files
  paths_dt    = fread(cells_paths_f)
  pb_cells    = paths_dt$pb_path %>% lapply(FUN = readRDS) %>% 
    do.call('cbind', .)
  
  # get nice rows
  rows_dt     = fread(rowdata_f) %>% setkey('ensembl_id')
  keep_ids    = rows_dt$ensembl_id
  assert_that(all(keep_ids %in% rownames(pb_cells)))
  pb_cells    = pb_cells[keep_ids, ]
  rows_dt     = rows_dt[ rownames(pb_cells) ]
  assert_that( all(rownames(pb_cells) == rows_dt$ensembl_id) )
  rownames(pb_cells) = rows_dt$gene_id
  
  # make one object with pb_mats as assays
  pb_cells    = SingleCellExperiment( counts(pb_cells), rowData = rows_dt )
  
  # add batch variable
  pb_cells[[ batch_var ]] = colnames(pb_cells)
  
  # store
  message(' save')
  saveRDS(pb_cells, pb_cells_f)
  
  message(' done!')
}



# 
# make_pb_cells_from_h5 <- function(h5_paths_f, qc_stats_f, batch_lu_f, coldata_f, rowdata_f, batch_var, run_var,
#   pb_cells_f, subset_f = NULL, subset_col = NULL, subset_str = NULL, n_cores = 4) {
#   
#   # get runs where any batches are left
#   batch_lu     = fread(batch_lu_f)
#   qc_stats_dt  = fread(qc_stats_f)
#   bad_var      = paste0("bad_", batch_var)
#   keep_batches = qc_stats_dt[ get(bad_var) == FALSE] %>% .[[batch_var]]
#   keep_runs    = batch_lu[batch_var %in% keep_batches, run_var]
#   
#   # get list of h5 files
#   h5_fs_dt    = fread(h5_paths_f) %>%
#     .[get(run_var) %in% keep_runs] %>%
#     unique()
#   h5_fs_ls    = h5_fs_dt$amb_filt_f %>% 
#     setNames(h5_fs_dt[[run_var]])
#   
#   # get all cells kept after qc
#   keep_dt     = fread(coldata_f) %>% .[ keep == TRUE ]
#   
#   # set up parallelization
#   bpparam       = MulticoreParam(workers = n_cores, tasks = length(h5_fs_ls))
#   on.exit(bpstop(bpparam))
#   
#   cell_pbs     = bpmapply( sel_run = names(h5_fs_ls), h5_f = unname(h5_fs_ls),
#    FUN = .get_one_cells_pb, SIMPLIFY = FALSE, BPPARAM = bpparam, 
#    MoreArgs = list(keep_dt, keep_batches, batch_lu, batch_var, run_var, 
#     subset_f = subset_f, subset_col = subset_col, subset_str = subset_str)
#   )
#   
#   # merge sce objects
#   pb_sce = Reduce(function(x, y) {cbind(x, y)}, cell_pbs)
#   
#   # subset rows
#   rows_dt     = fread(rowdata_f) %>% setkey('ensembl_id')
#   keep_ids    = rows_dt$ensembl_id
#   assert_that(all(keep_ids %in% rownames(pb_sce)))
#   pb_sce      = pb_sce[keep_ids, ]
#   rows_dt     = rows_dt[ rownames(pb_sce) ]
#   assert_that( all(rownames(pb_sce) == rows_dt$ensembl_id) )
#   rownames(pb_sce) = rows_dt$gene_id
#   
#   # add batch variable
#   # pb_cells[[ batch_var ]] = colnames(pb_cells)
#   
#   # store
#   message('  save')
#   saveRDS(pb_sce, pb_cells_f)
#   
#   message('done!')
# }
# 
# 
# 
# .get_one_cells_pb <- function(sel_run, h5_f, keep_dt, keep_batches, batch_lu, batch_var, run_var,
#   subset_f = NULL, subset_col = NULL, subset_str = NULL, agg_fn = c("sum", "prop.detected")){
#   
#   browser()
#   
#   agg_fn         = match.arg(agg_fn)
#   
#   # get batches to keep in this run
#   tmp_vars       = c(run_var, batch_var) %>% unique
#   batches_in_run  = batch_lu[run_var == sel_run, batch_var]
#   keep_batches    = intersect(keep_batches, batches_in_run)
#   cell_lu         = keep_dt[get(batch_var) %in% keep_batches,  c(tmp_vars, "cell_id"), with = FALSE] %>%
#     setkey("cell_id")
#   
#   message('    loading matrix from h5 file for run ', sel_run)
#   cells_mat      = .get_h5_mx(h5_f, sel_s = paste0(sel_run, ':'))
#   cells_mat      = .sum_SUA(cells_mat)
#   assert_that(all(cell_lu$cell_id %in% colnames(cells_mat)))
#   cells_mat      = cells_mat[, cell_lu$cell_id]
#   
#   # turn into sce
#   coldata = DataFrame(cell_lu)
#   rownames(coldata) = coldata$cell_id
#   sce  = SingleCellExperiment( assays = list(counts = cells_mat),colData = coldata )
#   
#   rm(cells_mat)
# 
#   # subset if required # this is wrong !!!!!
#   if (!is.null(subset_f)) {
#     message('    subsetting sce object')
#    # unpack
#     subset_vals = str_split(subset_str, pattern = ',') %>% unlist
#     subset_dt   = fread(subset_f)
#     assert_that(subset_col %in% colnames(subset_dt))
#     assert_that(all(c("cell_id", batch_var) %in% colnames(subset_dt)))
#   
#    # keep just cells in sel_b with selected labels
#     subset_dt   = subset_dt %>%
#     .[ get(batch_var) %in% keep_batches ] %>%
#     .[ get(subset_col) %in% subset_vals ]
#   
#   # subset sce object
#   sce       = sce[, subset_dt$cell_id]
#   }
#   
#   message('   running aggregateData')
#   pb_mat      = .aggregateData_datatable(sce, by_vars = batch_var, fun = agg_fn)
#   pb          = SingleCellExperiment( assays = list(counts = pb_mat) )
#   return(pb)
# }



make_pb_empty <- function(sel_run, af_paths_f, pb_empty_f, ambient_method, run_var = 'sample_id') {
  message('create pseudobulk matrix for empty droplets')
  # get files
  af_paths_df   = fread(af_paths_f)
  
  # get bad cellbender samples
  if (ambient_method == 'cellbender') {
    bad_runs      = af_paths_df[ bad_run == TRUE, get(run_var) ]
    if (sel_run %in% bad_runs) {
      # write an empty file
      file.create(pb_empty_f)
      return(NULL)
    }
  }
  
  # do some checks
  af_mat_f     = af_paths_df %>% .[get(run_var) == sel_run, af_mat_f] %>% unique()
  af_knee_f    = af_paths_df %>% .[get(run_var) == sel_run, af_knee_f] %>% unique() 
  assert_that(str_extract(af_mat_f, '\\/af_.*\\/') == str_extract(af_knee_f, '\\/af_.*\\/'))
  
  # get empty pseudobulk
  empty_pb     = .get_one_empty_pb( sample_id = sel_run, af_mat_f = af_mat_f, af_knee_f = af_knee_f)
  
  # save empty pseudobulk for one sample
  message(' saving empty pseudobulk for run ', sel_run)
  saveRDS(empty_pb, file = pb_empty_f)
  
  message(' done!')
}


merge_pbs_empty <- function(af_paths_f, rowdata_f, pb_empty_f, ambient_method) {
  # get paths to empty pb files for for runs that weren't excluded
  af_paths_df = fread(af_paths_f)  
  if (ambient_method == 'cellbender') {
    af_paths_df = af_paths_df %>%
      .[ bad_run == FALSE ]
  }
  
  empty_pb_fs = af_paths_df$pb_tmp_f %>% unique()
  empty_pbs   = empty_pb_fs %>% lapply(FUN = readRDS)
  empty_mat   = do.call('cbind', empty_pbs)
  
  # get nice rows
  rows_dt     = fread(rowdata_f) %>% setkey('ensembl_id')
  keep_ids    = rows_dt$ensembl_id
  assert_that(all(keep_ids %in% rownames(empty_mat)))
  empty_mat   = empty_mat[keep_ids, ]
  
  # add nice rows
  rows_dt     = rows_dt[ rownames(empty_mat) ]
  assert_that( all(rownames(empty_mat) == rows_dt$ensembl_id) )
  rownames(empty_mat)  = rows_dt$gene_id

  # make sce object
  pb_empty    = SingleCellExperiment( assays = list(counts = empty_mat), rowData = as(rows_dt, "DataFrame") )
  
  # store
  message(' save')
  saveRDS(pb_empty, pb_empty_f)
  
  message(' done!')
}


.get_one_empty_pb <- function(sample_id, af_mat_f, af_knee_f) {
  
  # get empty barcodes
  knee_df     = fread(af_knee_f)
  empty_bcs   = knee_df %>%
    .[in_empty_plateau == TRUE, barcode]
  
  rm(knee_df)

  # get full alevin matrix
  empty_mat      = .get_h5_mx(af_mat_f, sel_s = '')
  empty_mat      =  empty_mat[, empty_bcs]
  empty_mat      = .sum_SUA(empty_mat)

  # sum over all empties per samples
  empty_pb    = Matrix::rowSums(empty_mat) %>%
    setNames(rownames(empty_mat)) %>%
    as.matrix()
  colnames(empty_pb) = sample_id
  
  # make sparse
  empty_pb    = as(empty_pb, "dgCMatrix")
  
  return(empty_pb)
}


.aggregateData_datatable <- function(sce, by_vars = c('sample_id'),
  fun = c("sum", "mean", "median", "prop.detected", "num.detected")) {
  fun       = match.arg(fun)
  assay     = "counts"
  # check dimensions; if there are no cells don't return anything
  if (dim(sce)[2] == 0) { return(NULL) }
  
  # check have matrix format we need
  if (!("TsparseMatrix" %in% class(counts(sce))))
    counts(sce)   = counts(sce) %>% as("TsparseMatrix")

  # get counts data
  t_start   = Sys.time()
  mat_dt    = data.table(
    i         = counts(sce)@i + 1,
    j         = counts(sce)@j + 1,
    count     = counts(sce)@x
  ) %>% 
    .[, cell_id := colnames(sce)[j] ] %>% 
    setkey("cell_id")
  t_stop    = Sys.time()
  message('  time to create dt: ', round(difftime(t_stop, t_start, units = 'secs')), ' seconds')

  # make by dt
  by_dt = colData(sce)[, c('cell_id', by_vars)] %>% as.data.table %>% setkey("cell_id")
  
  # join
  t_start   = Sys.time()
  pb_dt     = mat_dt %>% merge(by_dt, by = "cell_id") %>% 
    .[, .(sum = sum(count), .N), by = c("i", by_vars) ]
  t_stop    = Sys.time()
  message('  time to aggregate: ', round(difftime(t_stop, t_start, units = 'secs')), ' seconds')
  
  # add n_cells
  n_cells_dt  = by_dt[, .(n_cells = .N), by = by_vars ]
  rows_dt     = data.table(i = seq.int(nrow(sce)), gene_id = rownames(sce))
  pb_dt       = pb_dt %>% 
    merge(n_cells_dt, by = by_vars) %>% 
    merge(rows_dt, by = "i") %>% 
    .[, i             := NULL ] %>% 
    .[, prop.detected := N / n_cells ]
  
  # (what comes makes only sense when by_vars = 'sample_id')
  genes_ls    = rownames(sce)
  batch_ls    = unique(by_dt[[by_vars]])
  n_genes     = length(genes_ls)
  n_batches   = length(batch_ls)
  
  message('making pseudobulk matrix')
  # do pseudobulk
  cast_frm    = sprintf("gene_id ~ %s", by_vars)
  pb_mat    = pb_dt %>% 
    dcast( as.formula(cast_frm), value.var = fun, fill = 0 ) %>% 
    as.matrix(rownames = "gene_id")
    
  # make full 0 matrix, fill it in to keep order of genes and samples
  mat         = matrix(0, nrow=n_genes, ncol=n_batches) %>%
    set_rownames(genes_ls) %>%
    set_colnames(batch_ls)
  mat[ rownames(pb_mat), colnames(pb_mat) ] = pb_mat
  
  # make sparse
  mat         = as(mat, "sparseMatrix") 
  
  return(mat)
}




plot_empty_plateau <- function(knees, empty_plat_df) {
  # get sample order
  s_ord = names(knees)
  
  # reformat
  knee_data <- lapply(s_ord, function(s) {
    x = knees[[ s ]]
    x %>%
      as_tibble() %>%
      group_by(lib_size = total) %>%
      summarize(n_bc = n()) %>%
      arrange(desc(lib_size)) %>%
      mutate(bc_rank = cumsum(n_bc)) %>%
      mutate(sample_id = s)
  }) %>% bind_rows()
  knee_data$sample_id = factor( knee_data$sample_id, levels = s_ord )
  
  empty_plat_df =  empty_plat_df %>%
    filter(sample_id %in% s_ord) %>%
    mutate(sample_id = factor(sample_id, levels = s_ord))
  
  # plot everything above low count threshold
  p_labels =  c("1", "10", "100", "1k", "10k", "100k", "1M")
  p_breaks =  c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6)

  p = ggplot(knee_data) +
    aes( x = bc_rank, y = lib_size ) +
    geom_line(linewidth = 0.3, color = '#283747') +
    geom_rect(data = empty_plat_df,
              mapping = aes(xmin = empty_start, xmax = empty_end, ymin = 0, ymax = Inf),
              inherit.aes = F, alpha = 0.1, fill = '#7c4b73', linewidth = 0) +
    facet_wrap( ~ sample_id, ncol = 4 ) +
    scale_x_log10(labels = p_labels, breaks = p_breaks) +
    scale_y_log10(labels = p_labels, breaks = p_breaks) +
    theme_classic(base_size = 9) +
    theme( legend.position = 'none' ) +
    labs(x = 'barcode rank', y = 'library size')
    
  return(p)
}


calc_empty_genes <- function(pb_cells_f, pb_empty_f, fdr_thr, logfc_thr, empty_gs_f) {
  message('calculate empty genes')
  
  message(' load pseudobulk matrices')
  cells_mat   = readRDS(pb_cells_f) %>% assays %>% .[[1]]
  empties_mat = readRDS(pb_empty_f) %>% assays %>% .[[1]]
  
  common_gs   = intersect(rownames(cells_mat), rownames(empties_mat))
  cells_mat   = cells_mat[ common_gs, ]
  empties_mat = empties_mat[ common_gs, ]

  message('  run edgeR')
  empty_genes = .run_edger_empties(cells_mat, empties_mat)
  
  # label empty genes
  empty_genes = empty_genes %>%
    .[, is_ambient := fcase(logFC > 0 & FDR < fdr_thr, TRUE, default = FALSE)]
  
  # store
  message('save results')
  fwrite( empty_genes, file = empty_gs_f )
  
  message('done!')
}


.run_edger_empties <- function(cells_mat, empties_mat) {
  assert_that( nrow(cells_mat) == nrow(empties_mat) )
  assert_that( all(rownames(cells_mat) == rownames(empties_mat)) )
  # set up DGE objects
  message(" remove outlier samples")
  name_ls   = c("cells", "empties")
  mat_ls    = list(cells_mat, empties_mat) %>% setNames(name_ls)
  x_ls      = name_ls %>% lapply(function(nn) {
    # get matrix
    mat       = mat_ls[[ nn ]]
    
    # remove samples with zero counts
    lib_sizes = colSums(mat)
    nzero_idx = lib_sizes > 0
    mat       = mat[, nzero_idx ]
    lib_sizes = lib_sizes[ nzero_idx ]
    
    # remove samples with outlier library sizes
    outliers  = scater::isOutlier(lib_sizes, log = TRUE, type = "lower", nmads = 3)
    mat       = mat[, !outliers, drop = FALSE] %>% 
      set_colnames( paste0( nn, "-", colnames(.) ) )
    
    return(mat)
  }) %>% setNames(c("cells", "empties"))
  
  # remove genes with zero counts, give better names
  message(" remove zero genes")
  row_sums    = lapply(x_ls, rowSums) %>% Reduce("+", .)
  keep_rows   = row_sums > 0

  # make big matrix
  message(" make DGE objects for each celltype")
  dge_ls    = name_ls %>% lapply(function(nn) {
    this_x    = x_ls[[ nn ]]
    dge_tmp   = DGEList(this_x, remove.zeros = FALSE)
    dge_tmp   = normLibSizes(dge_tmp, method = "TMMwsp")
  }) %>% setNames(name_ls)
  
  # join together
  dge       = do.call(cbind, dge_ls)
  
  # make design matrix
  message(" set up design")
  col_ns    = colnames(dge)
  des_all   = data.table(
    celltype  = str_extract(col_ns, "^[^-]+") %>% fct_relevel("empties", after = Inf), 
    sample_id = str_extract(col_ns, "[^-]+$")
  )
  
  # filter out tiny genes
  message(" filter gene matrix")
  mm_all    = model.matrix( ~ celltype, data = des_all )
  keep_gs   = filterByExpr(dge, group = des_all$celltype, min.count = 1)
  dge       = dge[ keep_gs, ]
  
  # estimate dispersion
  message(" estimate dispersion")
  if (length(unique(des_all$sample_id)) > 500) {
    message("  estimating dispersion without design matrix (because super slow when many samples)")
    dge       = estimateDisp(dge)
  } else {
    message("  estimating dispersion with design matrix")
    dge       = estimateDisp(dge, design = mm_all)
  }
  
  # fit for this celltype
  message(" fit edgeR")
  fit         = glmQLFit(dge, design = mm_all)
  fit         = glmTreat(fit, coef = "celltypeempties")
  exclude_dt  = topTags(fit, n = Inf, sort.by = "none") %>%
    as.data.frame %>% as.data.table(keep.rownames = 'gene_id')
  
  # calculate logcpms
  message(' calculate mean expression in all clusters')
  mean_cpms = name_ls %>% lapply(function(cl) {
    this_x    = x_ls[[ cl ]]
    lib_sizes = dge_ls[[ cl ]] %>% getNormLibSizes
    this_cpms = t(t(this_x) / lib_sizes) * 1e6
    means_tmp = data.table(
      col_lab    = paste0("mean_logcpm.", cl),
      gene_id     = rownames(this_cpms),
      mean_logcpm = rowMeans(log(this_cpms + 1))
    )
    return(means_tmp)
  }) %>% rbindlist %>% 
    dcast.data.table( gene_id ~ col_lab, value.var = "mean_logcpm" )
  
  # join
  message(' tidy up')
  exclude_dt  = exclude_dt %>% merge(mean_cpms, by = "gene_id") %>%
    .[ order(log(PValue) * sign(logFC), -logFC) ]
  
  message(' done!')
  
  return(exclude_dt)
}


.ignore_this <- function() {
  min_cpm = 100
  empties_dt %>% 
    .[ (mean_logcpm.cells >= log(min_cpm + 1)) | (mean_logcpm.empties >= log(min_cpm + 1)) ] %>% 
    .[, .(gene_id, log2fc = logFC %>% round(1), padj = FDR %>% signif(2), 
      CPM.cells = (exp(mean_logcpm.cells) - 1) %>% signif(2) %>% round,
      CPM.empties = (exp(mean_logcpm.empties) - 1) %>% signif(2) %>% round
    ) ] %>% print(topn = 20)
}
