
suppressPackageStartupMessages({
  library("assertthat")
  library("magrittr")
  library("forcats")
  library("stringr")

  library("BiocParallel")
  RhpcBLASctl::omp_set_num_threads(1L)
  library("rhdf5")
  library("data.table")
  library("fishpond")

  library("SingleCellExperiment")
  library("scattermore")
  library("ggridges")
  library("MetBrewer")
  library("scales")
  library("Matrix")
  library("scuttle")
})

save_hto_sce <- function(sce_df_f, sce_hto_f, n_cores){
  # unpack some inputs
  samples_dt  = fread(sce_df_f)
  samples     = samples_dt$pool_id
  bcs_ls      = samples_dt$bcs_csv
  hto_mats_ls = samples_dt$hto_f
  trans_ls    = samples_dt$wl_trans_f

  # demultiplex each pool and get sce hto objects
  bpparam     = MulticoreParam(workers = n_cores, tasks = length(samples))
  on.exit(bpstop(bpparam))

  message("  demultiplexing with hto counts")
  sce_ls      = bplapply(seq_along(samples), function(i) {
    # get sample and file
    sel_s       = samples[[ i ]]
    message(sel_s)

    bcs_f     = bcs_ls[[ i ]]
    hto_mat_f = hto_mats_ls[[ i ]]
    trans_f   = trans_ls[[ i ]]

    hto_sce = .get_one_hto_sce(sel_s, bcs_f, hto_mat_f, trans_f)
    return(hto_sce)
  }, BPPARAM = bpparam)

  # check no surprises
  assert_that( length(unique(sapply(sce_ls, nrow))) == 1 )

  # concatenate counts matrices
  message('  joining many hto matrices (takes a while)')
  counts_mat  = lapply(sce_ls, counts) %>% .join_spmats

  # get annotations for cells
  message('  joining colData info')
  cells_dt    = sce_ls %>%
    lapply(function(s) colData(s) %>% as.data.frame %>% as.data.table) %>%
    rbindlist
  assert_that( all.equal(colnames(counts_mat), cells_dt$cell_id) )

  rm(sce_ls); gc()

  # put into one big file
  message('  making sce object')
  sce = SingleCellExperiment(
    list(counts = counts_mat),
    colData = cells_dt
    )

  message(' saving hto sce object')
  saveRDS(sce, file = sce_hto_f, compress = FALSE)
  message('done!')

}


.get_one_hto_sce <- function(sel_s, bcs_f, hto_mat_f, trans_f){
  # get file for barcode translation
  bc_dict = trans_f %>% fread(header = FALSE) %>%
    set_colnames(c("bc_rna", "bc_hto"))

  # get all barcodes called as cells
  cell_bcs = fread(bcs_f, header = FALSE) %>%
    set_colnames("cell_bc")

  # get hto counts
  hto_counts = .get_alevin_mx(hto_mat_f, sel_s = '')

  # translate hto bcs to match rna barcodes
  hto_true_bcs = bc_dict[bc_hto %chin% colnames(hto_counts)] %>%
    .[order(match(bc_hto, colnames(hto_counts))), bc_rna]

  colnames(hto_counts) = hto_true_bcs

  # keep only cell barcodes
  keep_bcs = cell_bcs %>%
  .[cell_bc %chin% colnames(hto_counts), cell_bc]

  hto_counts = hto_counts[, keep_bcs]
  colnames(hto_counts) = paste(sel_s, colnames(hto_counts), sep = ':')

  # create a seurat object
  hto_seu = CreateSeuratObject(counts = hto_counts, assay = 'HTO')
  hto_seu = NormalizeData(hto_seu, assay = "HTO", normalization.method = "CLR")

  message("  demultiplexing sample ", sel_s)
  hto_seu = HTODemux(hto_seu, assay = "HTO", positive.quantile = 0.99)

  # get demultiplexing metadata
  demux_dt  = hto_seu[[]] %>% as.data.table(keep.rownames = "cell_id") %>%
    .[, .(cell_id, HTO_classification, HTO_classification.global, hash.ID)] %>%
    .[, guess := hash.ID %>% str_replace_all("-", "_") ] %>%
    .[, hash.ID := NULL] %>%
    .[, pool_id := sel_s] %>%
    .[, HTO_classification.global := tolower(HTO_classification.global)]

  # get sce object
  hto_sce = SingleCellExperiment(list(counts = hto_counts),
                       colData = demux_dt)

  return(hto_sce)

}



# functions for multiplexing QC

get_hto_dt <- function(pool, hto_sce){

  pool_cells = colData(hto_sce)$pool_id == pool
  pool_sce = hto_sce[, pool_cells]
  pool_counts = counts(pool_sce)

  pool_seu = CreateSeuratObject(counts = pool_counts, assay = 'HTO')
  pool_seu = NormalizeData(pool_seu, assay = "HTO", normalization.method = "CLR")

  pool_meta = colData(pool_sce) %>% as.data.table() %>%
    .[, hto_total := counts(pool_sce) %>% colSums ] %>%
    .[, n_htos    := nrow(pool_sce)]

  pool_norm_counts = GetAssayData(pool_seu, assay = "HTO", layer = "data") %>%
    t() %>%
    as.data.table(keep.rownames = 'cell_id') %>%
    melt(id.vars = 'cell_id', variable.name = 'hto_id', value.name = 'norm_count')

  pool_counts = GetAssayData(pool_seu, assay = 'HTO', layer = 'counts') %>%
    t() %>%
    as.data.table(keep.rownames = 'cell_id') %>%
    melt(id.vars = 'cell_id', variable.name = 'hto_id', value.name = 'count')

  pool_all = pool_norm_counts %>%
    merge(pool_counts, by = c('cell_id', 'hto_id')) %>%
    merge(pool_meta, by = 'cell_id') %>%
    .[, prop        := count / sum(count), by = .(pool_id, cell_id) ]  %>%
    .[, hto_id      := factor(hto_id)]

  return(pool_all)

}


hto_ridges <- function(sel_sample, proj_meta, hto_dt_ls){

  # get the right pool dt
  smpl_meta = proj_meta %>%
    .[sample_id == sel_sample]

  pool = smpl_meta$pool_id %>% unique()
  hto_id_smpl = smpl_meta$hto_id %>% unique()

  # get all htos in pool
  pool_htos = proj_meta %>%
    .[pool_id == pool, hto_id]

  pool_dt = hto_norm_counts_ls[[pool]] %>%
    .[guess == hto_id_smpl] %>%

    cols = met.brewer(name="Manet", n= length(pool_htos), type="discrete")

  p = ggplot(pool_dt, aes(x = norm_count, y = hto_id, fill = hto_id)) +
    geom_density_ridges(scale = 1, alpha = 0.8) +
    theme_minimal() +
    labs(
      title = paste(sel_sample, hto_id_smpl, sep = ', '),
      x = "hto_counts",
      y = NULL
    ) +
    theme(legend.position = "none") +
    scale_fill_manual(values = cols)

  return(p)
}


# make pairwise plots
hto_pairwise <- function(pool_dt, var = c("prop", "norm_count")){

  p_var = match.arg(var)

  plot_dt     = merge(
    pool_dt[, c('cell_id', 'guess', 'hto_id', p_var), with = FALSE],
    pool_dt[, c('cell_id', 'guess', 'hto_id', p_var), with = FALSE],
    by = c("cell_id", "guess"), allow.cartesian = TRUE, suffixes = c(".x", ".y")
  ) %>%
    .[ as.integer(hto_id.x) > as.integer(hto_id.y) ]

  hto_vals    = plot_dt$guess %>% unique() %>% setdiff(c('Doublet', 'Negative')) %>% sort
  cols_tmp    = MetBrewer::met.brewer( name = 'Johnson', n = length(hto_vals),
                                       type = 'discrete' ) %>% setNames(hto_vals)
  hto_cols    = c(cols_tmp, Negative = "grey20", Doublet = "grey80")

  x_var = paste0(p_var, '.x')
  y_var = paste0(p_var, '.y')

  g = ggplot(plot_dt) +
    aes( x = get(x_var), y = get(y_var), colour = guess ) +
    geom_point( size = 0.2 ) +
    scale_x_continuous( breaks = pretty_breaks() ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    scale_colour_manual( values = hto_cols, breaks = names(hto_cols) ) +
    guides( colour = guide_legend(override.aes = list(size = 3) ) ) +
    facet_grid( hto_id.y ~ hto_id.x ) +
    theme_classic() +
    labs(color = 'HTO guess')

  return(g)
}


.get_alevin_mx <- function(af_mat_f, sel_s){
  # get this file
  h5_filt   = H5Fopen(af_mat_f, flags = "H5F_ACC_RDONLY" )

  # get indices of barcodes
  mat       = sparseMatrix(
    i = as.vector(h5_filt$matrix$indices +1),
    p = as.vector(h5_filt$matrix$indptr),
    x = as.vector(h5_filt$matrix$data),
    repr = "C",
    dims = h5_filt$matrix$shape
  ) %>% as("TsparseMatrix")

  # add names
  bcs           = h5_filt$matrix$barcodes
  colnames(mat) = paste0(sel_s, bcs)
  rownames(mat) = as.character(h5_filt$matrix$features$name)

  return(mat)
}


