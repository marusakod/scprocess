

source('/home/kodermam/packages/scprocess/scripts/utils.R')

make_clean_sces(
sel_b         = 'GSM8612262_SRR31256771', 
sel_run       = 'GSM8612262_SRR31256771', 
integration_f = '/pstore/data/mus-brain-analysis/studies/test_project/output/test_integration/integrated_dt_test_project_2025-01-01.csv.gz', 
h5_paths_f    = '/pstore/data/mus-brain-analysis/studies/test_project/output/test_hvg/hvg_paths_test_project_2025-01-01.csv', 
coldata_f     = '/pstore/data/mus-brain-analysis/studies/test_project/output/test_qc/coldata_dt_all_cells_test_project_2025-01-01.csv.gz',
run_var       = 'sample_id',
batch_var     = 'sample_id',
clean_sce_f   = '/pstore/data/mus-brain-analysis/studies/test_project/output/test_integration/test_sce_cells_clean_GSM8612262_SRR31256771_test_project_2025-01-01.rds',
rowdata_f     = '/pstore/data/mus-brain-analysis/studies/test_project/output/test_qc/rowdata_dt_test_project_2025-01-01.csv.gz')


make_clean_sces <- function(sel_b, sel_run, integration_f, h5_paths_f, 
  coldata_f, rowdata_f, run_var, batch_var, clean_sce_f) {
  # load, exclude doublets
  int_dt        = fread(integration_f) %>%
    .[ is_dbl == FALSE & in_dbl_cl == FALSE ]
  
  # check if ok
  ok_batches    = int_dt[[batch_var]] %>% unique
  if (!(sel_b %in% ok_batches)) {
    message('excluded ', sel_b, '; creating empty sce file.')
    file.create(clean_sce_f)
    return(NULL)
  }
  message('creating clean sce file for ', sel_b)

  # get some subsets
  message('  getting values specific to ', sel_b)
  batch_int   = int_dt %>% .[ get(batch_var) == sel_b ]
  batch_ids   = batch_int$cell_id
  batch_cols  = fread(coldata_f) %>% setkey("cell_id") %>% .[ batch_ids ]

  # get sce object
  message('  loading counts into sce')
  h5_paths    = fread(h5_paths_f)
  filtered_f  = h5_paths[ get(run_var) == sel_run ]$amb_filt_f %>% unique
  assert_that(file.exists(filtered_f))
  sce         = .get_sce_clean(filtered_f, sel_run, run_var, rowdata_f,  subset_cells = batch_ids)
  
  # add things to colData
  message('  adding coldata')
  sce         = .add_coldata(sce, batch_cols)
  message('  adding integration outputs')
  sce         = .add_int_variables(sce, batch_int)

  # save
  message('  saving sce')
  saveRDS(sce, clean_sce_f, compress = FALSE)

  message('done!')
}


.get_sce_clean <- function(filtered_f, sel_run, run_var, rowdata_f, subset_cells = batch_ids){
  # read matrix
  mat         = .get_h5_mx(filtered_f, paste0(sel_run, ':')) %>% .sum_SUA
  if (!is.null(subset_cells)) {
    message('    subsetting sce to specified cells')
    assert_that( all(subset_cells %in% colnames(mat)) )
    mat         = mat[, subset_cells]
  }

  # read rowdata
  rows_dt       = fread(rowdata_f) %>% setkey('ensembl_id')
  keep_ids      = rows_dt$ensembl_id
  assert_that(all(keep_ids %in% rownames(mat)))
  mat           = mat[keep_ids, ]
  rows_dt       = rows_dt[ rownames(mat) ]
  assert_that( identical(rownames(mat), rows_dt$ensembl_id) )
  rownames(mat) = rows_dt$gene_id
  
  # make sce object
  sce_clean           = SingleCellExperiment(assays = list(counts = mat), rowData = rows_dt)
  sce_clean$cell_id   = colnames(mat)
  rownames(sce_clean) = rowData(sce_clean)$gene_id
  
  # convert to TsparseMatrix
  counts(sce_clean) = counts(sce_clean) %>% as("TsparseMatrix")

  return(sce_clean)
}

.add_coldata <- function(sce, batch_cols) {
  # some checks
  assert_that( colnames(colData(sce)) == "cell_id" )
  assert_that( all(colnames(sce) == batch_cols$cell_id) )

  # add columns
  cols_df       = batch_cols %>% as("DataFrame") %>% set_rownames(colnames(sce))
  colData(sce)  = cols_df

  return(sce)
}

.add_int_variables <- function(sce_clean, batch_int) {
  assert_that( all(colnames(sce_clean) == batch_int$cell_id) )
  # get useful integration variables, add to sce object
  int_vs      = c('UMAP1', 'UMAP2', str_subset(names(batch_int), "RNA_snn_res"))
  for (v in int_vs) {
    if (str_detect(v, "RNA_snn_res")) {
      colData(sce_clean)[[ v ]] = batch_int[[ v ]] %>% factor
    } else {
      colData(sce_clean)[[ v ]] = batch_int[[ v ]]
    }
  }

  return(sce_clean)
}