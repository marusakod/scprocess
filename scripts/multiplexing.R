
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
  library("yaml")

  library("SingleCellExperiment")
  library("scattermore")
  library("ggridges")
  library("MetBrewer")
  library("scales")
  library("Matrix")
  library("scuttle")
})

get_one_hto_sce <- function(sel_pool, sample_stats_f, amb_yaml_f, hto_mat_f, trans_f, hto_sce_f, ambient_method, seurat_quantile) {
  # if ambient method is cellbender exclude bad samples
  smpl_status = FALSE
  
  if(ambient_method == 'cellbender') {
    # loading file with bad bender samples
    message(' loading cellbender sample stats file')
    sample_stats_df = fread(sample_stats_f)
    smpl_status     = unique(sample_stats_df[ pool_id == sel_pool, bad_run])
    
    if(smpl_status) {
      message('  sample ', sel_pool, ' has been excluded. Saving empty results file')
      file.create(hto_sce_f)
      message('done!')
      
      return(NULL)
    }
  }
  
  # get file with all barcodes called as cells
  yaml_data       = yaml.load_file(amb_yaml_f)
  bcs_f = yaml_data$bcs_f
  
  # get file for barcode translation
  bc_dict = trans_f %>% fread(header = FALSE) %>%
    set_colnames(c("bc_rna", "bc_hto"))
  
  # get all barcodes called as cells
  cell_bcs = bcs_f %>% fread(header = FALSE) %>%
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
  colnames(hto_counts) = paste(sel_pool, colnames(hto_counts), sep = ':')
  
  # create a seurat object
  hto_seu   = CreateSeuratObject(counts = hto_counts, assay = 'HTO')
  hto_seu   = NormalizeData(hto_seu, assay = "HTO", normalization.method = "CLR")
  
  message("  demultiplexing sample ", sel_pool)
  hto_seu   = HTODemux(hto_seu, assay = "HTO", positive.quantile = seurat_quantile)
  
  # get demultiplexing metadata
  demux_dt  = hto_seu[[]] %>% as.data.table(keep.rownames = "cell_id") %>%
    .[, .(cell_id, guess = hash.ID, pool_id = sel_pool, HTO_classification, 
      HTO_classification.global = tolower(HTO_classification.global))]
  
  # get sce object
  hto_sce   = SingleCellExperiment(list(counts = hto_counts), colData = demux_dt)
  
  # create output file
  message(" saving results")
  saveRDS(hto_sce, file = hto_sce_f, compress = FALSE)
  
  message(" done!")
  return(NULL)
}


# functions for multiplexing QC
get_hto_dt <- function(pool_sce) {
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

  pool_counts_long = GetAssayData(pool_seu, assay = "HTO", layer = "counts") %>%
    t() %>%
    as.data.table(keep.rownames = 'cell_id') %>%
    melt(id.vars = 'cell_id', variable.name = 'hto_id', value.name = 'count')

  pool_all = pool_norm_counts %>%
    merge(pool_counts_long, by = c('cell_id', 'hto_id')) %>%
    merge(pool_meta, by = 'cell_id') %>%
    .[, prop        := count / sum(count), by = .(pool_id, cell_id) ]  %>%
    .[, hto_id      := factor(hto_id)]

  return(pool_all)
}


hto_ridges <- function(sel_sample, proj_meta, hto_dt_ls) {
  # get the right pool dt
  smpl_meta = proj_meta %>%
    .[ sample_id == sel_sample ]
  pool = smpl_meta$pool_id %>% unique()
  hto_id_smpl = smpl_meta$hto_id %>% unique()

  # get all htos in pool
  pool_htos = proj_meta %>%
    .[pool_id == pool, hto_id]
  pool_dt = hto_dt_ls[[pool]] %>%
    .[guess == hto_id_smpl] 

  cols  = MetBrewer::met.brewer( name = 'Renoir', n = length(pool_htos),
    type = 'discrete' ) %>% setNames(sort(pool_htos))
  pool_dt$hto_id = factor(pool_dt$hto_id, levels = names(cols))

  g = ggplot(pool_dt, aes(x = norm_count, y = hto_id %>% fct_rev, fill = hto_id)) +
    geom_density_ridges(scale = 1, alpha = 0.8) +
    theme_minimal() +
    labs(
      title = paste(sel_sample, hto_id_smpl, sep = ', '),
      x = "Normalized HTO counts",
      y = NULL
    ) +
    theme(legend.position = "none") +
    scale_fill_manual(values = cols)

  return(g)
}


# make pairwise plots
hto_pairwise <- function(pool_dt, var = c("prop", "norm_count")) {

  p_var = match.arg(var)

  plot_dt     = merge(
    pool_dt[, c('cell_id', 'guess', 'hto_id', p_var), with = FALSE],
    pool_dt[, c('cell_id', 'guess', 'hto_id', p_var), with = FALSE],
    by = c("cell_id", "guess"), allow.cartesian = TRUE, suffixes = c(".x", ".y")
  ) %>%
    .[ as.integer(hto_id.x) > as.integer(hto_id.y) ]
  # colour options w 8: "Homer1", "Redon", "Renoir", "Signac"
  hto_vals    = plot_dt$guess %>% unique() %>% setdiff(c('Doublet', 'Negative')) %>% sort
  cols_tmp    = MetBrewer::met.brewer( name = 'Renoir', n = length(hto_vals),
                                       type = 'discrete' ) %>% setNames(hto_vals)
  hto_cols    = c(cols_tmp, Negative = "grey20", Doublet = "grey80")

  x_var = paste0(p_var, '.x')
  y_var = paste0(p_var, '.y')

  if(var == 'prop') {
    x_title = 'HTO proportion 1'
    y_title = 'HTO_proportion 2'
  }else{
    x_title = 'Normalized HTO counts 1'
    y_title = 'Normalized HTO counts 2'
  }

  g = ggplot(plot_dt) +
    aes( x = get(x_var), y = get(y_var), colour = guess ) +
    geom_point( size = 0.2 ) +
    scale_x_continuous( breaks = pretty_breaks() ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    scale_colour_manual( values = hto_cols, breaks = names(hto_cols) ) +
    guides( colour = guide_legend(override.aes = list(size = 3) ) ) +
    facet_grid( hto_id.y ~ hto_id.x ) +
    theme_classic() +
    labs(color = 'HTO guess', x = x_title, y = y_title)

  return(g)
}


hto_barplot <- function(hto_dt_ls) {
  
  plot_dt = rbindlist(hto_dt_ls) %>%
    .[, .(pool_id, guess)] %>%
    .[, .(count = .N), by = .(pool_id, guess) ] %>%
    .[, pct := count / sum(count) * 100, by = pool_id ]
  
  hto_vals    = plot_dt$guess %>% unique() %>% setdiff(c('Doublet', 'Negative')) %>% sort
  cols_tmp    = MetBrewer::met.brewer( name = 'Renoir', n = length(hto_vals),
                                       type = 'discrete' ) %>% setNames(hto_vals)
  
  hto_cols    = c(cols_tmp, Negative = "grey20", Doublet = "grey80")
  plot_dt$guess = factor(plot_dt$guess, levels = names(hto_cols) )
  
  g = ggplot(plot_dt) +
    aes( y = pool_id, x = pct, fill = guess ) +
    geom_col( position="fill" ) +
    scale_fill_manual( values = hto_cols, breaks = names(hto_cols), guide = guide_legend(ncol = 2) ) +
    scale_x_continuous( breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20) ) +
    theme_classic() + theme( legend.position = "bottom" ) +
    labs(fill = 'HTO guess', x = 'pct. of barcodes', y = '')
  
  return(g)
}
