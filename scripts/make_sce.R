# see Marusa's file:
# /projects/site/pred/neurogenomics/users/kodermam/Siletti_2022/code/Siletti_SampleQC.R

# make_sce.R
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
  library("rtracklayer")
  library("Matrix")
  library("scuttle")
})


source('scripts/ambient.R') # to get the function that reads a matrix from .h5 file


save_cellbender_as_sce <- function(sce_df_f, metadata_f, gtf_dt_f, mito_str, sce_f, demux_f,
  bender_prob = 0.5, n_cores = 8, demux_type = "", keep_smpls_str) {

  # get name of run column
  if(demux_type == ""){
    run = 'sample_id'
  }else{
    run = 'pool_id'
  }

  # unpack some inputs
  samples_dt  = fread(sce_df_f)
  samples     = samples_dt[, get(run)]
  cb_full_ls  = samples_dt$cb_full
  cb_filt_ls  = samples_dt$cb_filt

  # check some inputs
  assert_that(
    length(samples) == length(cb_full_ls),
    length(samples) == length(cb_filt_ls)
  )
  assert_that(
    all(str_detect(cb_full_ls, samples)),
    all(str_detect(cb_filt_ls, samples))
  )

  # get some necessary inputs
  message('loading cellbender outputs into sce')
  message('  loading gene annotations')
  gene_annots = .get_gene_annots(gtf_dt_f)

  # check h5 files are ok
  h5closeAll()

  # run this in parallel
  bpparam     = MulticoreParam(workers = n_cores, tasks = length(samples))
  on.exit(bpstop(bpparam))
  sce_ls      = bplapply(seq_along(samples), function(i) {
    # get sample and file
    sel_s       = samples[[ i ]]
    message(sel_s)
    cb_full_f   = cb_full_ls[[ i ]]
    cb_filt_f   = cb_filt_ls[[ i ]]

    # get matrix
    bender_ls   = .get_one_bender_outputs(cb_full_f, cb_filt_f, sel_s, bender_prob)

    # turn into sce
    sce         = .make_one_sce(bender_ls, sel_s, run, gene_annots, mito_str)

    return(sce)
    }, BPPARAM = bpparam)

  # check no surprises
  assert_that( length(unique(sapply(sce_ls, nrow))) == 1 )

  # concatenate counts matrices
  message('  joining many matrices (takes a while)')
  counts_mat  = lapply(sce_ls, counts) %>% .join_spmats

  # remove weird genes
  weird_gs    = str_detect(rownames(counts_mat), "unassigned_gene")
  if ( any( weird_gs ) ) {
    warnings("removing some weird genes with names like 'unassigned_gene_1'")
    counts_mat  = counts_mat[ !weird_gs, ]
  }

  # get annotations for cells
  message('  joining colData info')
  cells_dt    = sce_ls %>%
    lapply(function(s) colData(s) %>% as.data.frame %>% as.data.table) %>%
    rbindlist
  assert_that( all.equal(colnames(counts_mat), cells_dt$cell_id) )
  rm(sce_ls); gc()

  # put into one big file
  message('  making sce object')
  sce         = SingleCellExperiment(list(counts = counts_mat),
    colData = cells_dt)

  # get annotations for rows
  message('  adding gene annotations')
  sce         = .add_gene_annots(sce, gene_annots)
  mt_gs       = str_detect(rownames(sce), mito_str)
  message(sprintf("  the following %d genes were selected as mitochondrial: ", sum(mt_gs)))
  message("    ", rowData(sce)$symbol[ mt_gs ] %>% paste(collapse = ", "))

  # add metadata
  message('  adding metadata')
  if(demux_type != ""){
    sce = sce %>%.add_demux_metadata(metadata_f, demux_f, demux_type)

    # remove any unwanted sample_ids from the object
    keep_smpls = str_split(keep_smpls_str, pattern = ',') %>% unlist()
    sce = sce[, colData(sce)$sample_id %in% keep_smpls]

  }else{
    sce = sce %>% .add_metadata(metadata_f)
  }

  message('  saving file')
  saveRDS(sce, file = sce_f, compress = FALSE)
  message('done!')
}



save_noncb_as_sce <- function(sce_df_f, ambient_method, metadata_f, gtf_dt_f, mito_str, sce_f, demux_f,
                              min_counts = 100, n_cores = 8, demux_type = "", keep_smpls_str ) {
  # get name of run column
  if(demux_type == ""){
    run = 'sample_id'
  }else{
    run = 'pool_id'
  }

  # unpack some inputs
  samples_dt  = fread(sce_df_f)
  samples     = samples_dt[, get(run)]

  # get a list of all matrices
  if(ambient_method == 'decontx'){
    all_cell_mat_fs = samples_dt$dcx_filt

  }else{
    all_cell_mat_fs = samples_dt$bcs_filt
  }

  # get gene annotations
  message('  loading gene annotations')
  gene_annots = .get_gene_annots(gtf_dt_f)

  # get sce objects
  message(' making sce objects for all samples')
  bpparam     = MulticoreParam(workers = n_cores, tasks = length(samples))
  on.exit(bpstop(bpparam))
  sce_ls      = bplapply(seq_along(samples), function(i) {
    # get sample and file
    sel_s       = samples[[ i ]]
    message(sel_s)
    mat_f     = all_cell_mat_fs[[ i ]]

    # turn into sce
    sce         = .get_one_nonbender_sce(mat_f, sel_s, mito_str, gene_annots, min_counts)

    return(sce)
  }, BPPARAM = bpparam)


  # check no surprises
  assert_that( length(unique(sapply(sce_ls, nrow))) == 1 )

  # concatenate counts matrices
  message('  joining many matrices (takes a while)')
  counts_mat  = lapply(sce_ls, counts) %>% .join_spmats

  # double-check for weird genes
  weird_gs = str_detect(rownames(counts_mat), "unassigned_gene")
  assert_that( all(!weird_gs) )

  # get annotations for cells
  message('  joining colData info')
  cells_dt    = sce_ls %>%
    lapply(function(s) colData(s) %>% as.data.frame %>% as.data.table) %>%
    rbindlist
  assert_that( all.equal(colnames(counts_mat), cells_dt$cell_id) )

  # put into one big file
  message('  making sce object')
  sce         = SingleCellExperiment(list(counts = counts_mat),
                                     colData = cells_dt)

  # get annotations for rows
  message('  adding gene annotations')
  sce         = .add_gene_annots(sce, gene_annots)
  rm(sce_ls); gc()

  # add metadata
  message('  adding metadata')
  if(demux_type != ""){
    sce = sce %>%.add_demux_metadata(metadata_f, demux_f, demux_type)
    
    keep_smpls = str_split(keep_smpls_str, pattern = ',') %>% unlist()
    # remove any unwanted sample_ids from the object
    sce = sce[, colData(sce)$sample_id %in% keep_smpls]

  }else{
    sce = sce %>% .add_metadata(metadata_f)
  }

  message('  saving file')
  saveRDS(sce, file = sce_f, compress = FALSE)
  message('done!')
}




.get_one_bender_outputs <- function(cb_full_f, cb_filt_f, sel_s, bender_prob) {
  # get this file
  h5_filt   = H5Fopen(cb_filt_f, flags = "H5F_ACC_RDONLY" )

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
  colnames(mat) = paste0(sel_s, ":", bcs)
  rownames(mat) = h5_filt$matrix$features$name

  # # find latent cell probabilities
  cell_probs  = h5_filt$droplet_latents$cell_probability
  if (is.null(cell_probs)) {
    cell_probs  = h5_filt$matrix$latent_cell_probability
  }
  bc_sums     = colSums(mat)

  # assemble dataframe
  bender_dt   = data.table(
    cell_id     = paste0(sel_s, ":", bcs),
    barcode     = bcs,
    prob_cell   = cell_probs,
    bc_count    = bc_sums
    ) %>%
    .[, bender_n_used   := .N ] %>%
    .[, bender_n_ok     := sum(prob_cell > 0.5) ] %>%
    .[, bender_prop_ok  := bender_n_ok / bender_n_used ] %>%
    .[, bender_logit_ok := qlogis( (bender_n_ok + 1) / (bender_n_used + 2) ) ]
  assert_that( all(bender_dt$cell_id == colnames(mat)) )

  # close h5 object
  H5Fclose(h5_filt)


  return(list(mat = mat, bender_dt = bender_dt))
}


.make_one_sce <- function(bender_ls, sel_s, run_var, gene_annots, mito_str) {
  # unpack inputs
  mat         = bender_ls$mat
  bender_dt   = bender_ls$bender_dt

  # split rownames into S / U / A
  splice_ns   = c("U", "S", "A")
  usa_ls      = splice_ns %>% lapply(function(l) {
    regex_str   = sprintf("_%s$", l)
    sel_gs      = str_subset(rownames(mat), regex_str)
    return(sel_gs)
    }) %>% setNames(splice_ns)

  # couple of checks that the gene names are sensible
  assert_that( length(table(sapply(usa_ls, length))) == 1 )
  g_ns_chk    = lapply(usa_ls, function(gs) str_replace(gs, "_[USA]$", ""))
  assert_that( all(sapply(g_ns_chk, function(l) all(l == g_ns_chk[[ 1 ]]))) )
  proper_gs   = g_ns_chk[[ 1 ]]

  # calculate spliced values
  total_spliced   = mat[ usa_ls[[ "S" ]], ] %>% Matrix::colSums(.)
  total_unspliced = mat[ usa_ls[[ "U" ]], ] %>% Matrix::colSums(.)
  usa_mat_ls      = lapply(usa_ls, function(gs) mat[ gs, ] )
  counts_mat      = Reduce("+", usa_mat_ls, accumulate = FALSE) %>%
    set_rownames( proper_gs )
  assert_that( all(colnames(counts_mat) == colnames(mat)) )

  # check for any missing genes
  missing_gs    = setdiff(rownames(counts_mat), gene_annots$ensembl_id)
  assert_that( (length(missing_gs) == 0) | (all(str_detect(missing_gs, "unassigned_gene"))) )

  # get totals with rRNA genes included
  total_raw     = colSums(counts_mat)

  # get counts for rRNA genes, exclude them
  rrna_ens_ids  = gene_annots[ gene_type == 'rRNA' ]$ensembl_id
  rrna_gs       = rownames(counts_mat) %in% rrna_ens_ids
  rrna_sum      = colSums(counts_mat[ rrna_gs, ] )
  mt_rrna_ens_ids   = gene_annots[ gene_type == 'Mt_rRNA' ]$ensembl_id
  mt_rrna_gs        = rownames(counts_mat) %in% mt_rrna_ens_ids
  mt_rrna_sum       = colSums(counts_mat[ mt_rrna_gs, ] )
  keep_gs       = and(!rrna_gs, !mt_rrna_gs)
  counts_mat    = counts_mat[ keep_gs, ]

  # get counts for mitochondrial genes
  mt_ensembls   = gene_annots[ str_detect(gene_annots$symbol, mito_str) ]$ensembl_id
  mt_gs         = rownames(counts_mat) %in% mt_ensembls
  mito_mat      = counts_mat[ mt_gs, ]
  mito_sum      = colSums(mito_mat)
  mito_detected = colSums(mito_mat > 0 )

  # make sce object
  sce_tmp             = SingleCellExperiment( assays = list(counts = counts_mat) )
  sce_tmp[[run_var]]  = sel_s
  sce_tmp$cell_id     = colnames(counts_mat)

  # add to sce object
  sce_tmp$sum         = colSums(counts_mat)
  sce_tmp$detected    = colSums(counts_mat > 0)
  sce_tmp$subsets_mito_sum = mito_sum
  sce_tmp$subsets_mito_detected = mito_detected
  sce_tmp$subsets_mito_percent = mito_sum / sce_tmp$sum * 100
  sce_tmp$total       = sce_tmp$sum

  # add splicing
  sce_tmp$total_spliced   = total_spliced
  sce_tmp$total_unspliced = total_unspliced
  sce_tmp$logit_spliced   = qlogis((total_spliced + 1) / (total_spliced + total_unspliced + 2))

  # add rrna
  sce_tmp$total_w_ribo  = total_raw
  sce_tmp$total_rrna    = rrna_sum
  sce_tmp$total_mt_rrna = mt_rrna_sum

  # add cellbender info
  assert_that( all(colnames(sce_tmp) == bender_dt$cell_id) )
  bender_vars   = c("prob_cell", "bc_count", "bender_n_used", "bender_n_ok",
    "bender_prop_ok", "bender_logit_ok", "mean_ok")
  for (bv in bender_vars) {
    sce_tmp[[ bv ]] = bender_dt[[ bv ]]
  }

  # exclude barcodes with 0 counts
  bc_totals     = Matrix::colSums(counts_mat)
  sce_tmp       = sce_tmp[, bc_totals > 0 ]

  return(sce_tmp)
}

.join_spmats <- function(mat_ls) {
  message("trying to join multiple sparse matrices")

  message("  removing any empty matrices")
  # remove any empty matrices
  mat_ls    = mat_ls[ sapply(mat_ls, function(m) length(m@x) > 0) ]

  message("  checking we can actually make the large matrix")
  # get how many non-zeros per matrix
  n_nnz     = vapply(mat_ls, function(m) length(m@x), FUN.VALUE = integer(1))
  cum_nz    = c(0, cumsum(n_nnz))
  n_cols    = vapply(mat_ls, ncol, FUN.VALUE = integer(1))
  cum_cols  = c(0, cumsum(n_cols))
  n_rows    = nrow(mat_ls[[1]])

  # initialize matrix
  n_total   = sum(n_nnz)
  assert_that( n_total <= 2^31,
    msg = "too many non-zero entries! can't make sparse matrix" )
  message("  initializing")
  mat_all   = new("dgTMatrix",
    i    = integer(n_total),
    j    = integer(n_total),
    x    = numeric(n_total),
    Dim  = c(n_rows, sum(n_cols))
    )

  # add each matrix to this
  message("  adding ", length(mat_ls), " samples:\n  ", appendLF = FALSE)
  for (i in seq_along(mat_ls)) {
    message(".", appendLF = FALSE)
    # get this matrix
    m         = mat_ls[[i]]

    # how to adjust locations?
    m_idx     = seq(cum_nz[i] + 1, cum_nz[i + 1], by = 1)
    j_delta   = as.integer(cum_cols[i])

    # put values in appropriate places
    mat_all@i[ m_idx ] = m@i
    mat_all@j[ m_idx ] = m@j + j_delta
    mat_all@x[ m_idx ] = m@x
  }
  message()

  # get all colnames
  message("  adding column and row names")
  all_cols  = lapply(mat_ls, colnames) %>% do.call(c, .)
  assert_that( length(all_cols) == ncol(mat_all) )
  colnames(mat_all) = all_cols
  rownames(mat_all) = rownames(mat_ls[[1]])
  message("done!")

  return(mat_all)
}

.get_gene_annots <- function(gtf_dt_f) {

  gene_annots   = fread(gtf_dt_f) %>%
    .[, chromosome := chromosome %>% fct_infreq ]
  if ("NC_007605.1" %in% levels(gene_annots$chromosome))
    gene_annots[, chromosome := chromosome %>% fct_relevel("NC_007605.1", after = Inf) ]

  assert_that( all(table(gene_annots$gene_id) == 1) )

  return(gene_annots)
}

.add_gene_annots <- function(sce_in, gene_annots) {
  # get current rows
  setkey(gene_annots, 'ensembl_id')
  assert_that( all(rownames(sce_in) %in% gene_annots$ensembl_id) )

  # add better annotations
  annots_dt     = gene_annots[ rownames(sce_in) ]

  # get nice ordering of genes
  annots_dt     = annots_dt[ order(chromosome, start, end) ]
  nice_order    = annots_dt$ensembl_id
  annots_df     = annots_dt[, .(gene_id, ensembl_id, symbol, gene_type)] %>%
    as('DataFrame') %>% set_rownames(.$gene_id)
  counts_mat    = counts(sce_in)[nice_order, ] %>%
    set_rownames(annots_df$gene_id)

  # put in order of chromosome and start
  sce_out       = SingleCellExperiment(
    assays = list(counts = counts_mat),
    colData = colData(sce_in), rowData = annots_df)
  rownames(sce_out) = rowData(sce_out)$gene_id

  return(sce_out)
}


.add_qc_metrics <- function(sce, mito_str, bpparam) {
  mt_gs     = str_detect(rownames(sce), mito_str)
  message(sprintf("  the following %d genes are selected as mitochondrial: ", length(mt_gs)))
  message("    ", rowData(sce)$symbol[ mt_gs ] %>% paste(collapse = ", "))
  sce       = sce %>%
    addPerCellQCMetrics( subsets = list( mito = mt_gs ), BPPARAM = bpparam )

  return(sce)
}


.add_metadata <- function(sce, metadata_f) {
  # get all metadata
  metadata_all  = fread(metadata_f)
  assert_that( "sample_id" %in% names(metadata_all) )
  assert_that( all(unique(sce$sample_id) %in% metadata_all$sample_id) )

  # join to coldata
  coldata_in    = colData(sce) %>% as.data.frame %>% as.data.table
  coldata_out   = merge(coldata_in, metadata_all, by = "sample_id") %>%
    setkey("cell_id")
  coldata_out   = coldata_out[ colnames(sce) ]
  coldata_df    = coldata_out %>% as('DataFrame') %>% set_rownames(.$cell_id)
  assert_that( all(colnames(sce) == coldata_df$cell_id) )

  # put this back
  colData(sce)  = coldata_df
  assert_that( !is.null(colnames(sce)) )

  return(sce)
}


.add_demux_metadata <- function(sce, metadata_f, demux_f, demux_type){
  metadata_all = fread(metadata_f)
  assert_that( all(unique(sce$pool_id) %in% metadata_all$pool_id))

  coldata_in = colData(sce) %>% as.data.frame() %>% as.data.table()

  if(demux_type == 'af'){
    hto_sce = readRDS(demux_f)

    # get demultiplexing metadata
    hto_coldata = colData(hto_sce) %>% as.data.frame() %>% as.data.table() %>%
     # fix labels of doublets from the same sample
    .[, c("hto1", "hto2") := tstrsplit(HTO_classification, "_", fixed = TRUE)] %>%
    .[, hto_id := fifelse(hto1 == hto2, hto1, HTO_classification)] %>%
    .[, c('hto1', 'hto2'): NULL] %>%
     setnames("HTO_classification.global", "demux_class")

    coldata_out = hto_coldata %>%
    # merge with sample metadata
     merge(metadata_all, by = c("hto_id", "pool_id"), all_x = TRUE) %>%
    # merge with rest of sce metadata
     merge(coldata_in, by = c("cell_id", "pool_id"))

  }else if(demux_type == 'custom'){
    demux_out = fread(demux_f)  %>%
      .[, cell_id := paste(pool_id, cell_id, sep = ":" )]

    common_bcs = length(intersect(demux_out$cell_id), coldata_in$cell_id)
    assert_that(common_bcs > 0)

    message(common_bcs, " matching between custom demultiplexing file and input sce")
    # discard all cells in demux_output but not in sce
    trash_h = length(setdiff(demux_out$cell_id, coldata_in$cell_id))
    message(trash_n, " cells from custom demultiplexing file discared")

    coldata_out = coldata_in %>%
      merge(demux_out, by = 'cell_id', all.x = TRUE, all.y = FALSE) %>%
      merge(metadata_all, by =c('pool_id', 'sample_id'))

    # check if column global class exists and if ok
    if('class' %in% colnames(demux_out)){
      class_vals == unique(demux_out$class)
      assert_that(all(class_vals %in% c("singlet", "negative", "doublet")))

      # label cells in sce but not in demux_output as "Negative"
      coldata_out = coldata_out %>%
      .[is.na(class()), class := "negative"] %>%
      setnames('class', 'demux_class')
    }else{
      coldata_out = coldata_out %>%
      .[, demux_class:= fifelse(is.na(sample_id), "negative", "singlet")]
    }
  }

  coldata_df    = coldata_out %>% as('DataFrame') %>% set_rownames(.$cell_id)
  assert_that( all(colnames(sce) == coldata_df$cell_id) )

  # put updated metadata back to original sce
  colData(sce)  = coldata_df
  assert_that( !is.null(colnames(sce)) )

  return(sce)

}


calc_gene_totals <- function(sce_input) {
  gene_totals_dt  = rowData(sce_input) %>% as.data.frame %>%
    as.data.table %>%
    .[, gene_total    := rowSums(counts(sce_input)) ] %>%
    .[, n_expressed   := rowSums(counts(sce_input) > 0) ] %>%
    .[, prop_exp      := rowMeans(counts(sce_input) > 0) ] %>%
    .[, gene_cpm      := gene_total / sum(gene_total) * 1e6 ]

  return(gene_totals_dt)
}


.get_one_nonbender_sce <- function(mat_f, sel_s, mito_str, gene_annots, min_counts) {

  # read matrix
  counts = .get_alevin_mx(af_mat_f = mat_f, sel_s = paste0(sel_s, ':'))


  # exclude barcodes with v low counts (maybe remove?)
  keep_idx    = colSums(counts) >= min_counts
  if ( sum(keep_idx) == 0 ) {
    warning("no barcodes kept for ", sel_s, "; skipping")
    return(NULL)
  }
  mat    = counts[, keep_idx]

  # split rownames into S / U / A
  splice_ns   = c("U", "S", "A")
  usa_ls      = splice_ns %>% lapply(function(l) {
    regex_str   = sprintf("_%s$", l)
    sel_gs      = str_subset(rownames(mat), regex_str)
    return(sel_gs)
  }) %>% setNames(splice_ns)

  # couple of checks that the gene names are sensible
  assert_that( length(table(sapply(usa_ls, length))) == 1 )
  g_ns_chk    = lapply(usa_ls, function(gs) str_replace(gs, "_[USA]$", ""))
  assert_that( all(sapply(g_ns_chk, function(l) all(l == g_ns_chk[[ 1 ]]))) )
  proper_gs   = g_ns_chk[[ 1 ]]

  # calculate spliced values
  total_spliced   = mat[ usa_ls[[ "S" ]], ] %>% Matrix::colSums(.)
  total_unspliced = mat[ usa_ls[[ "U" ]], ] %>% Matrix::colSums(.)
  usa_mat_ls      = lapply(usa_ls, function(gs) mat[ gs, ] )
  counts_mat      = Reduce("+", usa_mat_ls, accumulate = FALSE) %>%
    set_rownames( proper_gs )
  assert_that( all(colnames(counts_mat) == colnames(mat)) )

  # check for any missing genes
  missing_gs    = setdiff(rownames(counts_mat), gene_annots$ensembl_id)
  assert_that( (length(missing_gs) == 0) | (all(str_detect(missing_gs, "unassigned_gene"))) )

  # get totals with rRNA genes included
  total_raw     = colSums(counts_mat)

  # get counts for rRNA genes, exclude them
  rrna_ens_ids  = gene_annots[ gene_type == 'rRNA' ]$ensembl_id
  rrna_gs       = rownames(counts_mat) %in% rrna_ens_ids
  rrna_sum      = colSums(counts_mat[ rrna_gs, ] )
  mt_rrna_ens_ids   = gene_annots[ gene_type == 'Mt_rRNA' ]$ensembl_id
  mt_rrna_gs        = rownames(counts_mat) %in% mt_rrna_ens_ids
  mt_rrna_sum       = colSums(counts_mat[ mt_rrna_gs, ] )
  keep_gs       = and(!rrna_gs, !mt_rrna_gs)
  counts_mat    = counts_mat[ keep_gs, ]

  # get counts for mitochondrial genes
  mt_ensembls   = gene_annots[ str_detect(gene_annots$symbol, mito_str) ]$ensembl_id
  mt_gs         = rownames(counts_mat) %in% mt_ensembls
  mito_mat      = counts_mat[ mt_gs, ]
  mito_sum      = colSums(mito_mat)
  mito_detected = colSums(mito_mat > 0 )

  # make sce object
  sce_tmp             = SingleCellExperiment( assays = list(counts = counts_mat) )
  sce_tmp$sample_id   = sel_s
  sce_tmp$cell_id     = colnames(counts_mat)

  # add to sce object
  sce_tmp$sum         = colSums(counts_mat)
  sce_tmp$detected    = colSums(counts_mat > 0)
  sce_tmp$subsets_mito_sum = mito_sum
  sce_tmp$subsets_mito_detected = mito_detected
  sce_tmp$subsets_mito_percent = mito_sum / sce_tmp$sum * 100
  sce_tmp$total       = sce_tmp$sum

  # add splicing
  sce_tmp$total_spliced   = total_spliced
  sce_tmp$total_unspliced = total_unspliced
  sce_tmp$logit_spliced   = qlogis((total_spliced + 1) / (total_spliced + total_unspliced + 2))

  # add rrna
  sce_tmp$total_w_ribo  = total_raw
  sce_tmp$total_rrna    = rrna_sum
  sce_tmp$total_mt_rrna = mt_rrna_sum

  # convert to TsparseMatrix
  counts(sce_tmp) = counts(sce_tmp) %>% as("TsparseMatrix")

  return(sce_tmp)
}


save_hto_sce <- function(sce_df_f, sce_hto_f){
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

    hto_sce = .get_one_hto_sce(sel_s, bcs_f, hto_mat_fm, trans_f)
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

  messafe(' saving hto sce object')
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
  hto_counts = hto_counts[, hto_true_bcs]
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





# maybe useful for later...
#
# get_bcs_dt <- function(fry_dirs, name_regex, n_cores) {
#   # get subdirectories
#   fry_subs  = lapply(fry_dirs, function(fry_dir)
#     .get_fry_sub_dirs(fry_dir, name_regex, what = 'qc')) %>%
#     do.call(c, .)
#
#   # load featureDump.txt from each
#   message('loading barcodes:\n  ', appendLF = FALSE)
#   bpparam   = MulticoreParam(workers = n_cores, tasks = length(fry_subs))
#   on.exit(bpstop(bpparam))
#   bcs_dt    = bplapply(names(fry_subs), function(fry_name) {
#     message(fry_name, " ", appendLF = FALSE)
#     # get date
#     fry_sub   = fry_subs[[fry_name]]
#
#     # get date
#     fs        = list.files(fry_sub, recursive = FALSE )
#     date_str  = fs %>% str_subset("outs_alevin_[0-9]{4}-[0-9]{2}-[0-9]{2}") %>%
#       str_extract("[0-9]{4}-[0-9]{2}-[0-9]{2}") %>% unique %>%
#       sort(decreasing = TRUE) %>% .[ 1 ]
#     fry_dir   = paste0('outs_fry_', date_str)
#     alvn_dir  = paste0('outs_alevin_', date_str)
#
#     # get data
#     tmp_dt    = readAlevinFryQC(
#         mapDir    = file.path(fry_sub, alvn_dir),
#         permitDir = file.path(fry_sub, fry_dir),
#         quantDir  = file.path(fry_sub, fry_dir)
#       ) %>% use_series('cbTable') %>% as.data.table %>%
#       .[, sample_id := fry_name]
#     }, BPPARAM = bpparam) %>% rbindlist %>%
#     setcolorder('sample_id') %>% janitor::clean_names(.)
#   message()
#
#   return(bcs_dt)
# }

# make_alevinQC_reports <- function(fry_dirs, save_dir, overwrite = FALSE) {
#   # get all directories
#   fry_subs  = lapply(fry_dirs, function(fry_dir)
#     .get_fry_sub_dirs(fry_dir, name_regex, what = 'qc')) %>% do.call(c, .)
#
#   # for each one, make report
#   for (fry_name in names(fry_subs)) {
#     # get relevant subdirectory
#     fry_sub   = fry_subs[[fry_name]]
#     # check if already done
#     report_f  = sprintf('alevinQC_%s.html', fry_name)
#     if (file.exists(file.path(save_dir, report_f)) & overwrite == FALSE) {
#       message(fry_name, ' report already done; skipping')
#       message()
#       next
#     }
#
#     message()
#     message('making alevinQC report for ', fry_name)
#
#     # else make report
#     alevinQC::alevinFryQCReport(
#       mapDir    = file.path(fry_sub, 'outs_alevin'),
#       permitDir = file.path(fry_sub, 'outs_fry'),
#       quantDir  = file.path(fry_sub, 'outs_fry'),
#       sampleId = fry_name, outputFile = report_f, outputDir = save_dir,
#       forceOverwrite = overwrite
#       )
#   }
# }

# plot_knees <- function(bcs_dt, what = c('barcodes', 'genes'), do_facet = TRUE) {
#   what        = match.arg(what)
#
#   if (what == 'barcodes') {
#     # what to summarise by?
#     agg_var     = 'original_freq'
#     permit_ls   = c(TRUE, FALSE)
#
#     # define breaks
#     x_brks      = 10^seq(0, 7) %>% log10
#     x_labs      = c("1", "10", "100", "1k", "10k", "100k", "1M", "10M")
#     y_brks      = x_brks
#     y_labs      = x_labs
#
#     # define labels
#     x_lab       = "Cell barcode rank"
#     y_lab       = "Cell barcode frequency"
#
#     # what to do with x axis
#     .trans_fn   = log10
#
#   } else {
#     # what to summarise by?
#     agg_var     = 'nbr_genes_above_zero'
#     permit_ls   = c(TRUE)
#
#     # define breaks
#     x_brks      = seq(0, 2e4, 5e3)
#     x_labs      = as.character(x_brks)
#     y_brks      = c(0, 10^seq(0, 7)) %>% `+`(1) %>% log10
#     y_labs      = c("0", "1", "10", "100", "1k", "10k", "100k", "1M", "10M")
#
#     # define labels
#     x_lab       = "Cell barcode rank"
#     y_lab       = "No. genes observed"
#
#     # what to do with x axis
#     .trans_fn   = identity
#   }
#   assert_that( length(x_brks) == length(x_labs) )
#   assert_that( length(y_brks) == length(y_labs) )
#
#   # aggregate
#   plot_dt     = bcs_dt[ in_permit_list %in% permit_ls ] %>%
#     setnames(agg_var, 'agg_var', skip_absent = TRUE) %>%
#     .[, .(n_bc = .N), by = c("sample_id", "agg_var", "in_permit_list") ] %>%
#     .[ order(sample_id, -agg_var) ] %>%
#     .[, cum_bcs   := cumsum(n_bc), by = sample_id ] %>%
#     .[, n_ok      := max(.SD[ in_permit_list == TRUE ]$cum_bcs), by = sample_id ] %>%
#     .[, sample_id := fct_reorder(sample_id, n_ok) ]
#
#   # make plot
#   g = ggplot(plot_dt) +
#     aes(x = .trans_fn(cum_bcs), y = log10(agg_var), group = sample_id ) +
#     geom_step( aes(color = in_permit_list),
#       alpha = ifelse(do_facet, 1, 0.5), direction = 'vh' ) +
#     # scale_x_continuous( breaks = log_brks, labels = log_labs ) +
#     scale_x_continuous( breaks = x_brks, labels = x_labs ) +
#     scale_y_continuous( breaks = y_brks, labels = y_labs ) +
#     scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "grey")) +
#     theme_bw() +
#     theme(
#       panel.grid.minor = element_blank(),
#       panel.grid.major.x = element_blank(),
#       axis.title = element_text(size = 12)) +
#     labs(x = x_lab, y = y_lab)
#
#   if (do_facet)
#     g = g + facet_wrap( ~ sample_id )
#
#   return(g)
# }



