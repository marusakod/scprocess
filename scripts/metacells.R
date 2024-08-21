# metacells.R
suppressPackageStartupMessages({
  library('SuperCell')
})

apply_supercell <- function(sce_in_f, sce_out_f, mc_map_f, max_cells, n_cores = 1) {
  message('running SuperCell to make metacells')

  # set up cluster
  message('  setting up cluster')
  bpparam     = MulticoreParam( workers = n_cores )
  register(bpparam)
  on.exit(bpstop(bpparam))

  # load big sce file
  message('  loading sce file')
  sce_in      = readRDS(sce_in_f)

  # get samples
  message('  filtering samples')
  sample_ls   = sce_in$sample_id %>% unique %>% sort

  # aggregate each sample
  message('  running supercell on each sample')
  super_ls    = bplapply(sample_ls, function(s)
      .apply_one_supercell( sce_in[, sce_in$sample_id == s], max_cells ),
    BPPARAM = bpparam)

  # concatenate together
  message('  joining them together')
  assert_that( all(sapply(super_ls, function(s) 
    all( rownames(s) == rownames(super_ls[[1]]) ) 
    )) )
  sce_out     = super_ls %>% lapply( function(l) l$sce_out ) %>% do.call(cbind, .)
  mc_map_dt   = super_ls %>% lapply(function(l) l$mc_map) %>% rbindlist

  # add rows
  message('    adding row and column data')
  rowData(sce_out) = rowData(sce_in)

  # save sce_out
  message('    saving')
  saveRDS( sce_out, file = sce_out_f, compress = FALSE )
  fwrite( mc_map_dt, file = mc_map_f )

  message('done!')
}


.apply_one_supercell <- function(sce_tmp, max_cells, k_knn = 5) {
  # specify QC variables
  qc_vars     = c('sum', 'subsets_mito_sum', 'subsets_mito_percent', 'total',
    'total_spliced', 'total_unspliced', 'logit_spliced', 'total_w_ribo',
    'total_rrna', 'total_mt_rrna')

  # get sample-level variables
  cols_in     = colData(sce_tmp) %>% as.data.table
  var_counts  = cols_in %>% sapply(function(x) x %>% unique %>% length)
  meta_vars   = names(which(var_counts == 1)) %>%
    setdiff('sample_id') %>%
    setdiff( c(qc_vars, 'prob_cell') ) %>%
    str_subset( 'sample_(tube_)?id', negate = TRUE ) %>%
    str_subset( 'RNA_snn_res', negate = TRUE ) %>%
    str_subset( '^bender', negate = TRUE )
  meta_dt     = cols_in[, c('sample_id', meta_vars), with = FALSE] %>% unique

  # define all variables
  all_vars    = c('sample_id', 'cell_id', qc_vars, meta_vars, 'n_cells')

  # do log counts
  sce_tmp     = sce_tmp %>% logNormCounts
  n_sce_tmp   = ncol(sce_tmp)

  # check if easy
  if (n_sce_tmp <= max_cells) {
    # fiddle with columns
    cols_df     = colData(sce_tmp) %>%
      .[, c('sample_id', 'cell_id', qc_vars, meta_vars) ]
    cols_df$max_cells = 1
    assert_that( all(colnames(cols_df) == all_vars) )

    # copy sce
    sce_out     = SingleCellExperiment(
      assays = list(counts = counts(sce_tmp) ),
      colData = cols_df )

    # add cell map
    mc_map      = cols_in[, .(cell_id, meta_cell_id = cell_id)]

    # store
    return(sce_out)
  }

  # if not easy, set gamma
  gamma       = n_sce_tmp / max_cells

  # run SuperCell
  sc_obj      = logcounts(sce_tmp) %>%
    SCimplify( k.knn = k_knn, gamma = gamma, n.var.genes = 1000 )
  cell_mbrshp = sc_obj$membership
  super_mat   = counts(sce_tmp) %>%
    supercell_GE( cell_mbrshp, mode = 'sum' )

  # aggregate QC metrics
  qc_tmp      = cols_in[, c('sample_id', 'cell_id', qc_vars), with = FALSE] %>%
    .[, meta_cell := cell_mbrshp ] %>%
    .[, lapply(.SD, sum),
      by = .(sample_id, meta_cell), .SDcols = qc_vars ]
  ids_tmp     = cols_in[, .(cell_id, meta_cell = cell_mbrshp)] %>%
    .[, .(cell_id = .SD[1]$cell_id, n_cells = .N), by = meta_cell ]
  qc_tmp      = merge(ids_tmp, qc_tmp, by = 'meta_cell') %>%
    .[, subsets_mito_percent  := subsets_mito_sum / sum * 100 ] %>%
    .[, logit_spliced         := (total_spliced + 1) %>%
      `/`(total_unspliced + total_spliced + 2) %>% qlogis ]

  # save metacell map
  mc_map      = cols_in[, .(cell_id, meta_cell = cell_mbrshp)] %>%
    merge( ids_tmp[, .(meta_cell, meta_cell_id = cell_id)], by = "meta_cell" ) %>%
    .[, .(cell_id, meta_cell_id)]

  # join this to metadata
  cols_out    = qc_tmp %>% merge(meta_dt, by = 'sample_id') %>%
    .[ order(meta_cell) ] %>%
    .[, meta_cell := NULL ]

  # put in nice order
  cols_out    = cols_out[, all_vars, with = FALSE ]

  # check that we counted the number of cells correctly
  assert_that( all(cols_out$n_cells == table(cell_mbrshp)) )

  # turn into sce object
  colnames(super_mat) = cols_out$cell_id
  cols_df     = cols_out %>% as('DataFrame') %>% set_rownames( .$cell_id )
  sce_out     = SingleCellExperiment(
    assays = list( counts = super_mat ),
    colData = cols_df )

  return(list(sce_out = sce_out, mc_map = mc_map))
}
