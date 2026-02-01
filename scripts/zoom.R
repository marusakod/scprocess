
suppressPackageStartupMessages({
  library("data.table")
  library("magrittr")
  library("assertthat")
  library("SingleCellExperiment")
  library("stringr")
  library("zellkonverter")
})


make_subset_objects <- function(sel_b, batch_var, smpl_stats_f, h5ads_yaml_f, 
  subset_f, subset_col, subset_str, integration_f, save_sce, subset_sce_f, save_adata, subset_h5ad_f) {
  
  # get all good samples
  zoom_stats_dt = fread(smpl_stats_f)
  bad_var       = paste0("bad_", batch_var)
  ok_batches    = zoom_stats_dt[ get(bad_var) == FALSE ] %>% .[[ batch_var ]]
  
  # write empty file if sample was excluded
  if (!(sel_b %in% ok_batches)) {
    message('sample ', sel_b, ' was excluded. Creating empty files')
    if(save_sce) file.create(subset_sce_f)
    if(save_anndata) file.create(subset_h5ad_f)
    return(NULL)
  }

  # get list of input files
  message('Creating sce file for sample ', sel_b)
  paths_ls    = yaml::read_yaml(h5ads_yaml_f) %>% unlist()

  # get input sce file
  in_h5ad_f    = paths_ls[[sel_b]]
  assert_that( file.exists(in_h5ad_f) )
  in_sce      = readH5AD(in_h5ad_f)
  
  # load subset data.table
  subset_dt   = fread(subset_f)
  assert_that( subset_col %in% colnames(subset_dt) )
  assert_that( all(c("cell_id", batch_var) %in% colnames(subset_dt)) )
  
  # get cells matching subset
  subset_vals = str_split(subset_str, pattern = ',') %>% unlist
  subset_dt   = subset_dt %>%
    .[ get(batch_var) == sel_b ] %>%
    .[ get(subset_col) %in% subset_vals ] %>%
    .[, c(batch_var, 'cell_id', subset_col), with = FALSE]
  assert_that( all(subset_dt$cell_id %in% colnames(in_sce)) )

  # subset to cells we need
  sce_zoom    = in_sce[, subset_dt$cell_id]
  
  # get integration results
  int_dt      = fread(integration_f) 
  batch_int   = int_dt %>% .[ get(batch_var) == sel_b ]
  batch_cells = batch_int$cell_id

  # remove umap and clustering cols from before and add new ones, add label
  rm_cols     = c('UMAP1', 'UMAP2', str_subset(names(colData(sce_zoom)), "RNA_snn_res"))
  new_coldata = colData(sce_zoom) %>% as.data.table %>%
    .[, (rm_cols) := NULL ] %>%
    merge(subset_dt, by = c('cell_id', batch_var)) %>%
    setnames(new = 'label', old = subset_col) %>%
    setkey("cell_id") %>% 
    .[ batch_cells ]
  
  # tidy up sce
  sce_zoom    = sce_zoom[, batch_cells]
  colData(sce_zoom) = as(new_coldata, "DataFrame")
  
  # get useful integration variables
  int_vs   = c('UMAP1', 'UMAP2', str_subset(names(batch_int), "RNA_snn_res"))

  # add these to sce object
  for (v in int_vs) {
    if (str_detect(v, "RNA_snn_res")) {
      colData(sce_zoom)[[ v ]] = batch_int[[ v ]] %>% factor
    } else {
      colData(sce_zoom)[[ v ]] = batch_int[[ v ]]
    }
  }
  
  # save  
  if(save_sce) saveRDS(sce_zoom, subset_sce_f, compress = FALSE)
  if(save_adata) writeH5AD(sce_zoom, subset_h5ad_f)
  message('done!')
}
