suppressPackageStartupMessages({
  library('magrittr')
  library('data.table')
  library('forcats')
  library('assertthat')
  library('stringr')
  library('ggplot2')
  library('scales')
  library('patchwork')
  library('ComplexHeatmap')
  library('ggbeeswarm')
  library('ggh4x')
  library('Matrix')
  library('SingleCellExperiment')
  library('yaml')
  library('BiocParallel')
  library('scater')
  library('tidyverse')
  library('scDblFinder')
  library('DelayedArray')
  library('HDF5Array')
  library('rhdf5')
  library('ggplot.multistats')
  library('viridis')
  library('patchwork')
  library('UpSetR')
})

# define some breaks
n_brks      = c(1e1, 2e1, 5e1, 1e2, 2e2, 5e2, 
  1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5) %>% log10
n_labs      = c("10", "20", "50", "100", "200", "500",
  "1k", "2k", "5k", "10k", "20k", "50k", "100k")
log_brks    = c(1e1, 3e1, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5) %>%
  log10
log_labs    = c("10", "30", "100", "300", "1k", "3k", "10k", "30k", "100k", "300k")
mito_brks   = c(1e-4, 1e-3, 1e-2, 0.10, 0.50, 0.90, 0.99, 0.999) %>% qlogis
mito_labs   = c("0.01%", "0.1%", "1%", "10%", "50%", "90%", "99%", "99.9%")
splice_brks = c(1e-4, 1e-3, 1e-2, 0.10, 0.50, 0.90, 0.99, 0.999) %>% qlogis
splice_labs = c("0.01%", "0.1%", "1%", "10%", "50%", "90%", "99%", "99.9%")

main_qc <- function(run_name, metadata_f, cuts_f, amb_yaml_f, run_stats_f, demux_f, 
  gtf_dt_f, ambient_method, sce_fs_str, all_batches_str, rowdata_f, dbl_f, qc_f,
  coldata_f, mito_str, exclude_mito, hard_min_counts, hard_min_feats, hard_max_mito,
  run_var = "sample_id", demux_type = "none", batch_var = "sample_id", dbl_min_feats = 100, 
  dbl_min_cells = 100) {
  # check inputs
  exclude_mito  = as.logical(exclude_mito)

  # split output files and check if ok
  all_batches = str_split(all_batches_str, pattern = ',') %>% unlist()
  sce_fs_ls   = str_split(sce_fs_str, pattern = ',') %>% unlist()
  sce_fs_dirs = lapply(sce_fs_ls, FUN = dirname)
  assert_that(all(sapply(sce_fs_dirs, dir.exists)))
  assert_that(length(all_batches) == length(sce_fs_ls))
  assert_that(all(str_detect(sce_fs_ls, all_batches)))
  sce_fs_ls   = sce_fs_ls %>% setNames(all_batches)
  
  # check whether to exclude bc of cellbender
  exc_by_cb   = .was_excluded_by_cellbender(ambient_method, run_stats_f, run_var, run_name,
    sce_fs_ls, dbl_f, qc_f, rowdata_f, coldata_f)
  if (exc_by_cb)
    return(NULL)

  # get filtered ambient file
  yaml_data     = yaml.load_file(amb_yaml_f)
  filt_counts_f = yaml_data$filt_counts_f

  # get gene annotations
  gene_annots = .get_gene_annots(gtf_dt_f)

  # calculate qc metrics
  message('  loading counts into sce')
  sce         = .get_sce(filt_counts_f, run_name, mito_str, exclude_mito, gene_annots, run_var)

  # add annotations for rows
  message('  adding gene annotations')
  sce         = .add_gene_annots(sce, gene_annots)

  # add metadata
  message('  adding metadata')
  if (demux_type == "none") {
    sce         = sce %>% .add_metadata(metadata_f)
  } else {
    sce         = sce %>% .add_demux_metadata(metadata_f, demux_f, demux_type)
  }

  # do doublet calcs
  message('  starting doublet detection')
  dbl_dt      = .run_scdblfinder(sce, run_name, run_var, ambient_method, dbl_f)
  
  # add doublet info to coldata
  message('  adding doublet info to column data')
  sce         = .add_dbl_info(sce, dbl_dt, run_var, demux_type)

  # save rowdata
  rd          = rowData(sce) %>% as.data.table()
  
  message('  saving row data')
  fwrite(rd, file = rowdata_f)
  
  # do qc filtering, save table with qc for all singlets, coldata for all cells, get sce object with only singlets that pass qc
  message('  filtering cells')
  tmp_ls      = .filter_qc(sce, cuts_f, batch_var, run_var, demux_type, hard_min_counts, hard_min_feats, hard_max_mito)

  message('  saving column and qc data')
  fwrite(tmp_ls$coldata_dt, file = coldata_f)
  fwrite(tmp_ls$qc_all, file = qc_f)

 
  message('done!')
   
  return(NULL)
}

.was_excluded_by_cellbender <- function(ambient_method, run_stats_f, run_var, run_name,
  sce_fs_ls, dbl_f, qc_f, rowdata_f, coldata_f) {
  exc_by_cb   = FALSE
  if (ambient_method == 'cellbender') {
    # loading file with bad bender samples
    message('  loading cellbender sample stats file')
    run_stats_df  = fread(run_stats_f) 
    exc_by_cb     = unique(run_stats_df[get(run_var) == run_name, bad_run])
    if (exc_by_cb) {
      message('  run ', run_name, ' has been excluded. Saving empty results file')
      lapply(sce_fs_ls, file.create)
      # file.create(dimred_f)
      file.create(dbl_f)
      file.create(qc_f)
      file.create(rowdata_f)
      file.create(coldata_f)
      message('done!')
    }
  }
  return(exc_by_cb)
}

.get_gene_annots <- function(gtf_dt_f) {
  gene_annots   = fread(gtf_dt_f) %>%
    .[, chromosome := chromosome %>% fct_infreq ]
  if ("NC_007605.1" %in% levels(gene_annots$chromosome))
    gene_annots[, chromosome := chromosome %>% fct_relevel("NC_007605.1", after = Inf) ]

  assert_that( all(table(gene_annots$gene_id) == 1) )

  return(gene_annots)
}

.get_sce <- function(mat_f, sel_run, mito_str, exclude_mito, gene_annots, run_var, subset_cells = NULL) {
  # read matrix
  mat         = .get_h5_mx(mat_f, paste0(sel_run, ':'))
  if (!is.null(subset_cells)) {
    message('    subsetting sce to specified cells')
    assert_that( all(subset_cells %in% colnames(mat)) )
    mat         = mat[, subset_cells]
  }

  # heck for weird genes
  weird_gs    = str_detect(rownames(mat), "unassigned_gene")
  assert_that( all(!weird_gs) )

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
  total_raw     = Matrix::colSums(counts_mat)

  # get counts for rRNA genes, exclude them
  rrna_ens_ids  = gene_annots[ gene_type == 'rRNA' ]$ensembl_id
  rrna_gs       = rownames(counts_mat) %in% rrna_ens_ids
  rrna_sum      = Matrix::colSums(counts_mat[ rrna_gs, ] )
  mt_rrna_ens_ids   = gene_annots[ gene_type == 'Mt_rRNA' ]$ensembl_id
  mt_rrna_gs        = rownames(counts_mat) %in% mt_rrna_ens_ids
  mt_rrna_sum       = Matrix::colSums(counts_mat[ mt_rrna_gs, ] )
  keep_gs       = and(!rrna_gs, !mt_rrna_gs)
  counts_mat    = counts_mat[ keep_gs, ]

  # get counts for mitochondrial genes
  mt_ensembls   = gene_annots[ str_detect(gene_annots$symbol, mito_str) ]$ensembl_id
  mt_gs         = rownames(counts_mat) %in% mt_ensembls
  mito_mat      = counts_mat[ mt_gs, ]
  mito_sum      = Matrix::colSums(mito_mat)
  mito_detected = Matrix::colSums(mito_mat > 0 )
  
  # make sce object
  sce_tmp             = SingleCellExperiment( assays = list(counts = counts_mat) )
  sce_tmp[[run_var]]  = sel_run
  sce_tmp$cell_id     = colnames(counts_mat)

  # add to sce object
  sce_tmp$sum         = Matrix::colSums(counts_mat)
  sce_tmp$detected    = Matrix::colSums(counts_mat > 0)
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
  
  # remove mitochondrial genes
  if(exclude_mito == TRUE) {
    sce_tmp = sce_tmp[!mt_gs, ]
  }
  
  # convert to TsparseMatrix
  counts(sce_tmp) = counts(sce_tmp) %>% as("TsparseMatrix")

  # exclude any zeros
  nzero_idx   = colSums(counts(sce_tmp)) > 0
  sce_tmp     = sce_tmp[, nzero_idx]

  return(sce_tmp)
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

.add_metadata <- function(sce, meta_f) {
  # get all metadata
  metadata_all  = fread(meta_f)
  assert_that( "sample_id" %in% names(metadata_all) )
  assert_that( all(unique(sce$sample_id) %in% metadata_all$sample_id) )

  # join to coldata
  coldata_in    = colData(sce) %>% as.data.frame() %>% as.data.table
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

.add_demux_metadata <- function(sce, meta_f, demux_f, demux_type) {
  # set up
  metadata_all  = fread(meta_f)
  assert_that( all(unique(sce$pool_id) %in% metadata_all$pool_id))
  coldata_in    = colData(sce) %>% as.data.frame() %>% as.data.table()

  # do stuff for hto
  if (demux_type == 'hto') {
    hto_sce       = readRDS(demux_f)
    # convert all underscores in hto_id column to match seurat output
    metadata_all$hto_id = gsub('_', '-', metadata_all$hto_id)

    # get demultiplexing metadata
    hto_coldata   = colData(hto_sce) %>%
      as.data.frame() %>%
      as.data.table() %>%
      setnames(
        old = c("HTO_classification.global", "HTO_classification"),
        new = c("demux_class", "hto_id")
      )

    coldata_out = hto_coldata %>%
      # merge with sample metadata
      merge(metadata_all, by = c("hto_id", "pool_id"), all.x = TRUE) %>%
      # merge with rest of sce metadata
      merge(coldata_in, by = c("cell_id", "pool_id"), all.y = TRUE) %>%
      # label all cells missing from hto mat as negative
      .[is.na(demux_class), demux_class := 'negative'] %>%
      .[is.na(guess),   guess   := 'Negative'] %>%
      .[is.na(hto_id),  hto_id  := 'Negative']

  # do stuff for custom
  } else if (demux_type == 'custom') {
    # load data, check
    demux_out   = fread(demux_f)
    if (!("cell_id" %in% colnames(demux_out)))
      demux_out   = demux_out %>% .[, cell_id := paste(pool_id, barcode, sep = ":" )]
    demux_out   = demux_out[, .(cell_id, sample_id)]
    common_bcs  = length(intersect(demux_out$cell_id, coldata_in$cell_id))
    assert_that(common_bcs > 0)

    # discard all cells in demux_output but not in sce
    message(common_bcs, " matching between custom demultiplexing file and input sce")
    demux_out   = demux_out[ cell_id %in% common_bcs ]

    # merge together
    coldata_out = coldata_in %>%
      merge(demux_out, by = c('cell_id'), all.x = TRUE, all.y = FALSE) %>%
      merge(metadata_all, by = c('pool_id', 'sample_id'), all.x = TRUE)

    # check if column global class exists and if ok
    if ('class' %in% colnames(demux_out)) {
      class_vals = unique(demux_out$class)
      assert_that(all(class_vals %in% c("singlet", "negative", "doublet")))

      # label cells in sce but not in demux_output as "Negative"
      coldata_out = coldata_out %>%
        .[is.na(class), class := "negative"] %>%
        setnames('class', 'demux_class')

    # otherwise assume everything is fine!
    } else {
      coldata_out = coldata_out %>%
        .[, demux_class := fifelse(is.na(sample_id), "negative", "singlet") ]
    }
  }

  # prep for sce colData
  coldata_df    = coldata_out %>% as('DataFrame') %>% set_rownames(.$cell_id)
  coldata_df    = coldata_df[colnames(sce), ]
  assert_that( all(colnames(sce) == coldata_df$cell_id) )

  # put updated metadata back to original sce
  colData(sce)  = coldata_df
  assert_that( !is.null(colnames(sce)) )

  return(sce)
}

.run_scdblfinder <- function(sce, run_name, run_var = 'sample_id', ambient_method,
  dbl_f, min_feats = 100, min_cells = 100) {
  # subset to just this sample
  run_idx     = sce[[run_var]] == run_name
  sce         = sce[, run_idx ]
    
  # exclude tiny cells
  message('  filtering out cells with low counts')
  # assert_that( 'detected' %in% colnames(colData(sce)) )
  detected    = colSums(counts(sce) > 0)
  keep_idx    = sce$detected >= min_feats
  assert_that( sum(keep_idx) >= min_cells, msg = "insufficient cells to run scDblFinder :(")
  message(sprintf('    keeping %d / %d cells (%.0f%%)',
    sum(keep_idx), length(keep_idx), 100 * sum(keep_idx) / length(keep_idx)))
  sce         = sce[, keep_idx]
    
  # run scDblFinder
  message('  running scDblFinder')
  dbl_dt      = scDblFinder(sce, returnType = 'table', multiSampleMode = 'singleModel', 
      verbose = FALSE ) %>% as.data.table(keep.rownames = 'cell_id') %>%
    .[, (run_var) := run_name ] %>% setcolorder(run_var) %>%
    .[ type == 'real' ] %>% setnames("class", "scdbl_class")
  
  # check if class is available from demultiplexing
  if ('demux_class' %in% colnames(colData(sce))) {
    # extract demux_class
    demux_dt = colData(sce) %>% as.data.table(keep.rownames = 'cell_id') %>%
      .[, c("cell_id", run_var, "demux_class"), with = FALSE]

    # join to dbl data
    dbl_dt = dbl_dt %>%
      merge(demux_dt, by = c('cell_id', run_var), all.x = TRUE, all.y = FALSE)
  }
  
  # put in nice order    
  setkeyv(dbl_dt, "cell_id")
  dbl_dt      = dbl_dt[ colnames(sce) ]

  # save
  assert_that( all(dbl_dt$cell_id == colnames(sce)) )
  message('  saving results')
  fwrite(dbl_dt, file = dbl_f)   

  message('done!')
    
  return(dbl_dt)
}

.add_dbl_info <- function(sce, dbl_dt, sample_var, demux_type = "none") {
  # get input
  coldata_in    = colData(sce) %>% as.data.frame() %>% as.data.table
  missing_cells = setdiff(coldata_in$cell_id, dbl_dt$cell_id)
  
  # define columns we need  
  if (demux_type == "none") {
    keep_cols     = c("cell_id", sample_var, "scdbl_class")
    by_cols       = c("cell_id", sample_var)
  } else {
    keep_cols     = c("cell_id", sample_var, "scdbl_class", "demux_class")
    by_cols       = c("cell_id", sample_var, "demux_class")
  }
  
  # merge data together
  dbl_data      = dbl_dt %>%
    .[, ..keep_cols, with = FALSE ]  
  coldata_out   = merge(coldata_in, dbl_data, by = by_cols, all.x = TRUE)

  # set dbl_class to singlet or replace with demux_class if available
  # these seem to be extremely poor quality cells, so ok to set as singlet here and have them removed at qc
  if (length(missing_cells) != 0) {
    coldata_out   = coldata_out %>%
      .[is.na(scdbl_class), scdbl_class := 'singlet']
  }
  
  # put back into colData in sce  
  coldata_out   = coldata_out %>%
    as('DataFrame') %>% set_rownames(.$cell_id)
  coldata_out   = coldata_out[colnames(sce), ]
  assert_that( identical(colnames(sce), coldata_out$cell_id))
  colData(sce)  = coldata_out
  
  return(sce)
}

.calc_one_dimred <- function(sce, sel_sample) {
  # run PCA on this sce
  sce       = sce %>% logNormCounts %>% runPCA
  pca_dt    = reducedDim(sce, "PCA") %>%
    as.data.table %>%
    .[,1:2] %>% set_colnames(c('pc1', 'pc2'))
  dimred_dt = data.table(
    cell_id   = colnames(sce),
    sample_id = sel_sample
  ) %>% cbind(pca_dt)
  
  return(dimred_dt)
}

.filter_qc <- function(sce, cuts_f, batch_var, run_var, demux_type, hard_min_counts, hard_min_feats, hard_max_mito) {
  # store initial coldata
  coldata_in  = colData(sce) %>% as.data.table()
  cuts_dt     = fread(cuts_f)
  
  # restrict to singlets
  if (demux_type == "none") {
    keep_idx    = coldata_in$scdbl_class == 'singlet'
  } else {
    if (batch_var == "sample_id") {
      keep_idx    = (coldata_in$scdbl_class == 'singlet') & (coldata_in$demux_class == 'singlet')
    } else if (batch_var == "pool_id") {
      keep_idx    = coldata_in$scdbl_class == 'singlet'
    }
  }

  # make qc dt
  suppressMessages({
    qc_all      = coldata_in[keep_idx] %>% 
      make_qc_dt(batch_var = batch_var, run_var = run_var,
        qc_names = c('log_counts', 'log_feats', 'logit_mito', 'logit_spliced') )    
  })

  # baseline filtering
  qc_dt       = qc_all %>%
    .[ log_counts     >= log10(hard_min_counts) ] %>%
    .[ log_feats      >= log10(hard_min_feats) ] %>%
    .[ logit_mito     <  qlogis(hard_max_mito) ]
  assert_that( all(unique(qc_dt[[batch_var]]) %in% cuts_dt[[batch_var]]))
  
  # apply additional filtering criteria
  keep_dt     = qc_dt %>%
    merge(cuts_dt, by = batch_var) %>% 
    .[ log_counts     >= log10(min_counts)] %>%
    .[ log_feats      >= log10(min_feats)] %>%
    .[ logit_mito     > qlogis(min_mito)] %>%
    .[ logit_mito     < qlogis(max_mito)] %>%
    .[ logit_spliced  > qlogis(min_splice)] %>%
    .[ logit_spliced  < qlogis(max_splice)]
  
  # record which kept  
  qc_all      = qc_all %>%
    .[, keep_hard := cell_id %in% qc_dt$cell_id ] %>%
    .[, keep      := cell_id %in% keep_dt$cell_id]
  
  # add keep column to original coldata
  message('  subsetting column and qc data')
  keep_ids    = qc_all[ keep == TRUE ]$cell_id
  coldata_out = coldata_in %>%
    .[, keep := fifelse(cell_id %chin% keep_ids, TRUE, FALSE)]
  
  return(list( coldata_dt = coldata_out, qc_all = qc_all ))
}

.save_or_save_empty <- function(sce_obj, f, label) {
  if (ncol(sce_obj) == 0) {
    message("No cells passed qc for sample ", label, ". Saving empty file")
    file.create(f)
  } else {
    assert_that(length(f) == 1)
    saveRDS(sce_obj, file = f, compress = FALSE)
  }    
}

## Helper to merge QC and coldata files using data.table::fread/fwrite
merge_qc_files <- function(qc_files, col_files, qc_out, col_out) {
  # helper to select only existing, non-empty files
  sel_existing = function(fs) {
    fs[ file.exists(fs) & (file.info(fs)$size > 0) ]
  }

  qc_files_ok  = sel_existing(qc_files)
  col_files_ok = sel_existing(col_files)

  # If no qc files available, create empty outputs and exit
  if (length(qc_files_ok) == 0) {
    # create empty gz files to match expected outputs
    file.create(qc_out)
    file.create(col_out)
    return(invisible(NULL))
  }

  # read and rbind qc tables
  qc_dt_list = lapply(qc_files_ok, function(f) fread(f))
  qc_dt_all  = rbindlist(qc_dt_list, use.names = TRUE, fill = TRUE)
  # write gzipped output
  fwrite(qc_dt_all, file = qc_out, compress = 'gzip')

  # handle coldata merging
  if (length(col_files_ok) == 0) {
    # create only qc output if no coldata
    file.create(col_out)
    return(invisible(NULL))
  }
  col_dt_list <- lapply(col_files_ok, function(f) fread(f))
  col_dt_all  <- rbindlist(col_dt_list, use.names = TRUE, fill = TRUE)
  fwrite(col_dt_all, file = col_out, compress = 'gzip')

  return(invisible(NULL))
}

make_cuts_dt <- function(cuts_f, batch_var, b_lvls) {
  cuts_tmp  = fread(cuts_f) %>% 
    .[, .(
      batch_var         = get(batch_var), 
      log_counts_min    = log10(min_counts), 
      log_feats_min     = log10(min_feats), 
      logit_mito_min    = qlogis(min_mito), 
      logit_mito_max    = qlogis(max_mito), 
      logit_spliced_min = qlogis(min_splice), 
      logit_spliced_max = qlogis(max_splice), 
      n_cells_min       = log10(min_cells)
    )] %>% 
    melt( id = "batch_var", variable.name = "cut_var", value.name = "cut_point") %>% 
    .[, qc_var  := str_extract(cut_var, "^(.+)(?=_(min|max))") ] %>% 
    .[, minmax  := str_extract(cut_var, "(min|max)") ]

  # possibly subset
  assert_that( all(b_lvls %in% cuts_tmp$batch_var) )
  cuts_tmp  = cuts_tmp[ batch_var %in% b_lvls ]

  # make wide
  cuts_dt   = cuts_tmp %>%
    dcast(batch_var + qc_var ~ minmax, value.var = "cut_point") %>% 
    .[ is.na(min), min := -Inf ] %>% 
    .[ is.na(max), max := Inf ]
  assert_that( all(!is.na(cuts_dt$batch_var)) )

  return(cuts_dt)
}

plot_qc_ranges_marginals <- function(qc_input, b_lvls, qc_names, qc_lu, cuts_dt, batch_var) {
  # melt, add names
  tmp_names = intersect(qc_names, colnames(qc_input))
  qc_melt   = copy(qc_input) %>%
    melt(measure = tmp_names, val = 'qc_val', var = 'qc_var') %>%
    .[, qc_full   := qc_lu[ as.character(qc_var) ] ] %>%
    .[, qc_var    := factor(qc_var, levels = tmp_names) ] %>%
    .[, qc_full   := fct_reorder(qc_full, as.integer(qc_var)) ]
  hlines_dt = cuts_dt %>% copy %>%
    .[, qc_full   := qc_lu[ qc_var ] ] %>%
    .[, qc_var    := factor(qc_var, levels = qc_names) ] %>%
    .[, qc_full   := fct_reorder(qc_full, as.integer(qc_var)) ]

  # calculate medians etc
  qc_meds = qc_melt %>%
    .[, .(
      log10_N = log10(.N),
      q50     = median(qc_val, na.rm = TRUE),
      q10     = quantile(qc_val, 0.1, na.rm = TRUE),
      q90     = quantile(qc_val, 0.9, na.rm = TRUE),
      q025    = quantile(qc_val, 0.025, na.rm = TRUE),
      q975    = quantile(qc_val, 0.975, na.rm = TRUE)
      ), by = c('batch_var', 'qc_var', 'qc_full') ]

  # bar width
  bar_w     = 0.4
  n_dt      = qc_meds[, .(batch_var, `no. of cells` = log10_N) ] %>% unique %>%
    melt.data.table( id = "batch_var", var = "var", val = "value")

  # put in nice order
  n_dt      = n_dt %>%
    .[, batch_var := factor(batch_var, levels = rev(b_lvls)) ]
  n_lims    = c( min(n_dt$value) + log10(0.5), max(n_dt$value) + log10(2) )
  qc_melt   = qc_melt %>%
    .[, batch_var := factor(batch_var, levels = rev(b_lvls)) ]
  qc_meds   = qc_meds %>%
    .[, batch_var := factor(batch_var, levels = rev(b_lvls)) ]
  cuts_dt   = cuts_dt %>% copy %>% 
    .[, batch_var := factor(batch_var, levels = rev(b_lvls)) ]
  hlines_dt = hlines_dt %>%
    .[, batch_var := factor(batch_var, levels = rev(b_lvls)) ]

  # make plot of n_cells by sample
  g_n = ggplot() +
    geom_rect( data = hlines_dt[ qc_var == "n_cells" ],
      aes(
        xmin = as.integer(batch_var) - 0.4, xmax = as.integer(batch_var) + 0.4, 
        ymin = min, ymax = max
      ),
      fill = 'grey80', colour = NA, alpha = 0.2 ) +
    geom_point( data = n_dt, aes( y = value, x = as.integer(batch_var) ),
      size = 4, shape = 21, fill = 'grey40') +
    scale_x_continuous( breaks = seq.int(length(b_lvls)), labels = rev(b_lvls) ) +
    facet_grid( . ~ var, scales = 'free', space = 'free_y' ) +
    scale_y_continuous(breaks = n_brks, labels = n_labs) +
    expand_limits( y = n_lims ) +
    coord_flip( xlim = c(0.5, length(b_lvls) + 0.5), expand = FALSE ) +
    theme_classic() +
    theme(
      axis.text.x       = element_text( angle = 90, hjust = 1, vjust = 0.5 ),
      strip.text.y      = element_blank()
      ) +
    labs( y = NULL, x = batch_var )

  # make plots of qc metrics by sample
  g_violin = ggplot() +
    geom_rect( data = hlines_dt[ qc_var != "n_cells" ], 
      aes(
        xmin = as.integer(batch_var) - 0.4, xmax = as.integer(batch_var) + 0.4, 
        ymin = min, ymax = max
      ),
      fill = 'grey80', colour = NA, alpha = 0.2 ) +
    geom_violin( data = qc_melt[ !is.na(qc_val) ],
      aes( x = batch_var, y = qc_val ), colour = NA, fill = 'grey40',
      kernel = 'rectangular', adjust = 0.1, scale = 'width', width = 0.8) +
    # geom_hline( data = hlines_dt, aes( yintercept = cut_point ),
    #   colour = 'black', linetype = 'dashed', linewidth = 0.5, alpha = 0.5 ) +
    facet_grid( . ~ qc_full, scales = 'free', space = 'free_y' ) +
    scale_x_discrete( breaks = levels(qc_meds$batch_var), drop = FALSE ) +
    facetted_pos_scales(
      y = list(
        qc_full == "no. of UMIs"     ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full == "no. of genes"  ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full == "mito. pct."         ~
          scale_y_continuous(breaks = mito_brks, labels = mito_labs),
        qc_full == "spliced pct."      ~
          scale_y_continuous(breaks = splice_brks, labels = splice_labs)
        )
      ) +
    coord_flip() +
    theme_classic() +
    theme(
      axis.ticks.y      = element_blank(),
      axis.text.y       = element_blank(),
      axis.text.x       = element_text( angle = 90, hjust = 1, vjust = 0.5 )
      ) +
    labs( x = NULL, y = NULL )

  g = g_n + g_violin + plot_layout(widths = c(1, 5))

  return(g)
}

plot_qc_ranges_pairwise <- function(qc_input, qc_names, qc_lu, thrshlds_dt) {
 
  # calc medians etc
  qc_meds   = qc_input %>%
    melt(measure = qc_names, val = 'qc_val', var = 'qc_var') %>%
    .[, .(
      log10_N   = log10(.N),
      q50       = median(qc_val, na.rm = TRUE),
      q10       = quantile(qc_val, 0.1, na.rm = TRUE),
      q90       = quantile(qc_val, 0.9, na.rm = TRUE),
      q025      = quantile(qc_val, 0.025, na.rm = TRUE),
      q975      = quantile(qc_val, 0.975, na.rm = TRUE)
      ),
      by = c('sample_id', 'qc_var')] %>%
    .[, qc_full := qc_lu[ qc_var ] ]

  # make pairwise plot
  pairs_dt  = merge(qc_meds, qc_meds,
    by = c('sample_id', 'log10_N'), allow.cartesian = TRUE) %>%
    .[, qc_var.x  := factor(qc_var.x, levels = qc_names) ] %>%
    .[, qc_full.x := fct_reorder(qc_full.x, as.integer(qc_var.x)) ] %>%
    .[, qc_var.y  := factor(qc_var.y, levels = qc_names) ] %>%
    .[, qc_full.y := fct_reorder(qc_full.y, as.integer(qc_var.y)) ] %>%
    .[ as.integer(qc_var.x) > as.integer(qc_var.y) ]

  scales_x_ls = list(
        qc_full.x == "no. of UMIs"    ~
          scale_x_continuous(breaks = log_brks, labels = log_labs),
        qc_full.x == "no. of genes" ~
          scale_x_continuous(breaks = log_brks, labels = log_labs),
        qc_full.x == "mito. pct."        ~
          scale_x_continuous(breaks = mito_brks, labels = mito_labs),
        qc_full.x == "spliced pct."     ~
          scale_x_continuous(breaks = splice_brks, labels = splice_labs)
  )

  scales_y_ls = list(
        qc_full.y == "no. of UMIs"    ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full.y == "no. of genes" ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_full.y == "mito. pct."        ~
          scale_y_continuous(breaks = mito_brks, labels = mito_labs),
        qc_full.y == "spliced pct."     ~
          scale_y_continuous(breaks = splice_brks, labels = splice_labs)
  )

  # make plot
  g = ggplot(pairs_dt) +
    geom_point(
      aes( x = q50.x, y = q50.y, size = log10_N ),
      colour = 'black', shape = 21, alpha = 0.8
      ) +
    scale_x_continuous( breaks = pretty_breaks() ) +
    scale_y_continuous( breaks = pretty_breaks() ) +
    scale_size( range = c(1, 4), breaks = log_brks, labels = log_labs ) +
    guides(fill = guide_legend(override.aes = list(size = 3, shape = 21) ) ) +
    facet_grid( qc_full.y ~ qc_full.x, scales = 'free' ) +
    facetted_pos_scales(
      x = scales_x_ls,
      y = scales_y_ls
      ) +
    theme_bw() +
    theme(
      axis.text.x       = element_text( angle = 90, hjust = 1, vjust = 0.5 ),
      panel.grid        = element_blank(),
      strip.background  = element_rect( fill = 'white')
      ) +
    labs( x = NULL, y = NULL, size = 'no. of cells\nin sample' )

  return(g)
}

make_qc_dt_file <- function(sce_f, qc_f, overwrite=FALSE) {
  # check if already done
  if (file.exists(qc_f) & overwrite == FALSE)
    return(fread(qc_f))

  # load sce file
  sce         = sce_f %>% readRDS

  # calculate reads
  all_counts  = sce %>% counts %>% Matrix::colSums(.)

  # calculate mito reads
  mt_idx      = rowData(sce)$symbol %>%
    str_detect('^MT-')
  assert_that( sum(mt_idx) == 13 )
  mt_counts   = sce[mt_idx, ] %>% counts %>% Matrix::colSums(.)

  # calculate no. features
  all_feats   = sce %>% counts %>% `>`(0) %>% Matrix::colSums(.)

  # assemble
  qc_dt = data.table(
    cell_id     = colnames(sce),
    all_counts  = all_counts,
    all_feats   = all_feats,
    mito_counts = mt_counts
    )
  fwrite(qc_dt, file=qc_f)

  return(qc_dt)
}

get_cols_dt <- function(sce_f, cols_f, overwrite=FALSE) {
  # check if already done
  if (file.exists(cols_f) & overwrite == FALSE)
    return(fread(cols_f))

  # load sce file
  sce     = sce_f %>% readRDS

  # assemble
  cols_dt   = colData(sce) %>% as.data.frame() %>%
    as.data.table(keep.rownames = 'cell_id')
  fwrite(cols_dt, file=cols_f)

  return(cols_dt)
}

plot_qc_metric_scatter <- function(dt, qc_names, qc_lu, cuts_one, name) {
  qc_names = intersect(qc_names, colnames(dt))
  melt_dt = dt %>%
    melt( measure = qc_names, value.name='qc_val', variable.name='qc_var') %>%
    .[, qc_full := qc_lu[ as.character(qc_var) ] ] %>%
    .[, qc_var  := factor(qc_var, levels = qc_names) ] %>%
    .[, qc_full := fct_reorder(qc_full, as.integer(qc_var)) ]

  # get what to plot
  plot_dt = merge(
    melt_dt[, .(cell_id, qc_x = qc_full, val_x = qc_val)],
    melt_dt[, .(cell_id, qc_y = qc_full, val_y = qc_val)],
    by = 'cell_id', allow.cartesian = TRUE
    ) %>% .[ as.integer(qc_x) > as.integer(qc_y) ]

  # make some lines
  rects_dt  = merge(
    cuts_one[ qc_var %in% qc_names ] %>% 
      .[, .(dummy = "dummy", qc_var.y = qc_var, min.y = min, max.y = max)],
    cuts_one[ qc_var %in% qc_names ] %>% 
      .[, .(dummy = "dummy", qc_var.x = qc_var, min.x = min, max.x = max)],
    by = "dummy", allow.cartesian = TRUE) %>% 
    .[, qc_x      := qc_lu[ qc_var.x ] ] %>%
    .[, qc_var.x  := factor(qc_var.x, levels = qc_names) ] %>%
    .[, qc_x      := fct_reorder(qc_x, as.integer(qc_var.x)) ] %>%
    .[, qc_y      := qc_lu[ qc_var.y ] ] %>%
    .[, qc_var.y  := factor(qc_var.y, levels = qc_names) ] %>%
    .[, qc_y      := fct_reorder(qc_y, as.integer(qc_var.y)) ] %>%
    .[ as.integer(qc_x) > as.integer(qc_y)]

  # get sample name
  g = ggplot() +
    geom_rect( data = rects_dt, 
      aes( xmin = min.x, xmax = max.x, ymin = min.y, ymax = max.y ),
      fill = 'grey80', colour = NA, alpha = 0.4 ) +
    geom_bin2d(data = plot_dt, aes( x=val_x, y=val_y ), bins = 40) + 
    scale_fill_distiller( palette='RdBu', trans='log10' ) +
    facet_grid( qc_y ~ qc_x, scales='free' ) +
    theme_classic() +
    theme(
      panel.grid        = element_blank(),
      strip.background  = element_rect(fill = 'white')
    ) +
    facetted_pos_scales(
      x = list(
        qc_x == "no. of UMIs"    ~
          scale_x_continuous(breaks = log_brks, labels = log_labs),
        qc_x == "no. of genes" ~
          scale_x_continuous(breaks = log_brks, labels = log_labs),
        qc_x == "mito. pct."        ~
          scale_x_continuous(breaks = mito_brks, labels = mito_labs),
        qc_x == "spliced pct."     ~
          scale_x_continuous(breaks = splice_brks, labels = splice_labs)
        ),
      y = list(
        qc_y == "no. of UMIs"    ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_y == "no. of genes" ~
          scale_y_continuous(breaks = log_brks, labels = log_labs),
        qc_y == "mito. pct."        ~
          scale_y_continuous(breaks = mito_brks, labels = mito_labs),
        qc_y == "spliced pct."     ~
          scale_y_continuous(breaks = splice_brks, labels = splice_labs)
        )
      ) +
    labs(
      x     = 'QC metric 1',
      y     = 'QC metric 2',
      fill  = 'no. of cells',
      title = name
      )

  return(g)
}

plot_totals_split_by_meta <- function(pre_dt, post_dt, meta_dt) {
  # define age splits
  age_breaks  = c(0, 40, 50, 60, 70, 80, 100)
  age_labels  = paste0('<=', age_breaks[-1])
  yrs_breaks  = c(0, 10, 20, 30, 40, 60)
  yrs_labels  = paste0(yrs_breaks[-length(yrs_breaks)], ' to ', yrs_breaks[-1])

  # load metadata
  meta_dt     = copy(meta_dt) %>%
    .[, lesion_type := lesion_type %>%
      fct_recode(`WM (ctrl)` = "WM", `GM (ctrl)` = "GM")] %>%
    .[, age_cat     := cut(age_at_death, breaks = age_breaks,
      labels = age_labels) %>% factor(levels = age_labels) ] %>%
    .[, yrs_w_ms    := cut(years_w_ms, breaks = yrs_breaks,
      labels = yrs_labels) %>% factor(levels = yrs_labels) ] %>%
    .[, .(sample_id, subject_id, matter, lesion_type, diagnosis,
      sex, age_cat, yrs_w_ms, pmi_cat, brain_bank = sample_source, seq_pool)]

  # join to keep totals
  pre_n_dt    = pre_dt[, .(pre_n = .N), by = sample_id]
  post_n_dt   = post_dt[, .(post_n = .N), by = sample_id]
  keep_n_dt   = merge(pre_n_dt, post_n_dt, by = 'sample_id', all.x = TRUE) %>%
    .[ is.na(post_n), post_n := 0 ] %>%
    merge(meta_dt, by = 'sample_id')

  # pick which to show, get their values
  meta_vars   = c('matter', 'lesion_type', 'diagnosis', 'sex', 'age_cat',
    'yrs_w_ms', 'pmi_cat', 'brain_bank', 'seq_pool')
  level_ord   = meta_vars %>%
    lapply(function(v) levels(factor(keep_n_dt[[v]]))) %>%
    do.call(c, .)

  # melt by variable
  tmp_dt      = keep_n_dt %>%
    melt(id = c('sample_id', 'subject_id', 'pre_n', 'post_n'),
      measure = meta_vars, variable.name = 'meta_var', value.name = 'meta_val') %>%
    .[, meta_val := factor(meta_val, levels=level_ord)] %>%
    melt(measure = c('pre_n', 'post_n'),
      value.name = 'n_cells', variable.name = 'qc_status') %>%
    .[n_cells > 0, .(
      n_cells   = round( sum(n_cells) / 1e3 ) %>% as.integer,
      n_samples = .N,
      n_donors  = length(unique(subject_id))
      ), by = .(meta_var, meta_val, qc_status)] %>%
    melt(measure = c('n_cells', 'n_samples', 'n_donors'),
      value.name = 'count_val', variable.name = 'count_var')

  # calculate lost cells etc
  lost_dt     = merge(
    tmp_dt[qc_status == 'pre_n', .(meta_var, meta_val, count_var, pre_val=count_val)],
    tmp_dt[qc_status == 'post_n', .(meta_var, meta_val, count_var, post_val=count_val)],
    by = c('meta_var', 'meta_val', 'count_var')
    ) %>%
    .[, lost_val := pre_val - post_val]
  all_dt      = rbind(
    tmp_dt[qc_status == 'post_n', .(meta_var, meta_val, qc_status = 'kept', count_var, count_val)],
    lost_dt[, .(meta_var, meta_val, qc_status = 'excluded', count_var, count_val=lost_val)]
    ) %>% .[, qc_status := qc_status %>% factor(levels=c('excluded', 'kept')) ] %>%
    .[, count_var := count_var %>% fct_recode(`#k cells` = 'n_cells',
      `# samples` = 'n_samples', `# donors` = 'n_donors') ]

  # plot
  g = ggplot(all_dt[count_val > 0]) +
    aes( y = fct_rev(meta_val), x = count_val,
      label = count_val, fill = qc_status ) +
    geom_col(colour = NA) +
    geom_text(colour = 'black', position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_manual(values = c(kept = 'grey50', excluded = 'grey80'),
      breaks = c('kept', 'excluded')) +
    scale_x_continuous(breaks = pretty_breaks()) +
    facet_grid( meta_var ~ count_var, scales='free', space='free_y') +
    theme_classic() +
    labs(y = NULL, x = NULL, fill = 'QC status',
      title = 'Summary of cells, samples and donors retained')

  return(g)
}

plot_qc_summary_heatmap <- function(qc_stats, meta_input) {
  # make matrix of z-scores
  stats_tmp = copy(qc_stats) %>%
    .[ sample_id %in% meta_input$sample_id ] %>%
    .[, z := scale(med_val), by = qc_var ] %>%
    .[ qc_var %in% c('logit_mito', 'logit_splice'), z := z * -1 ]

  # make matrix of z-scores
  z_wide    = stats_tmp %>%
    dcast.data.table( sample_id ~ qc_var, value.var = 'z')
  z_mat     = z_wide[, -'sample_id', with = FALSE] %>%
    as.matrix %>% set_rownames(z_wide$sample_id)

  # make matrix of text labels
  lab_wide  = stats_tmp %>%
    .[ qc_var == "log10_N",
      lab := sprintf("%0.1fk", 10^med_val / 1e3) ] %>%
    .[ qc_var == "log_counts",
      lab := sprintf("%0.1fk", 10^med_val / 1e3) ] %>%
    .[ qc_var == "log_feats",
      lab := sprintf("%0.1fk", 10^med_val / 1e3) ] %>%
    .[ qc_var == "logit_mito",
      lab := sprintf("%0.1f%%", plogis(med_val) * 1e2) ] %>%
    .[ qc_var == "logit_splice",
      lab := sprintf("%0.0f%%", 1/(1 + 2^med_val) * 1e2) ] %>%
    dcast.data.table( sample_id ~ qc_var, value.var = 'lab')
  lab_mat   = lab_wide[, -'sample_id', with = FALSE] %>%
    as.matrix %>% set_rownames(lab_wide$sample_id)

  # lots of colours
  z_cols    = cols_fn(seq(-2, 2, 0.2), 0.2, 'RdBu',
    pal_dir = -1, range = 'symmetric')

  # define function for labelling
  labelling_fn <- function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%s", lab_mat[i, j]), x, y,
      gp = gpar(fontsize = 6))
  }

  # make sample annotations
  rows_dt   = copy(meta_input) %>% setkey('sample_id') %>%
    .[rownames(z_mat)]

  # do colours for patients
  pats_dt     = meta_input[, .(subject_n = .N), by = subject_id] %>%
    .[ subject_n > 1, .(subject_id, subject_n) ] %>% .[ order(-subject_n)]
  pat_list    = pats_dt$subject_id
  pat_cols    = nice_cols[seq.int(length(pat_list))] %>% setNames(pat_list)
  donor_vals  = rows_dt$subject_id
  donor_vals[ !(donor_vals %in% names(pat_cols)) ] = NA
  les_vals    = rows_dt$lesion_type %>% fct_drop %>% levels
  les_cols    = brewer.pal(length(les_vals), 'BrBG') %>%
    setNames(les_vals)

  # make col annotations
  row_annots  = rowAnnotation(
    lesion    = rows_dt$lesion_type,
    diagnosis = rows_dt$diagnosis,
    donor     = donor_vals,
    sex       = rows_dt$sex,
    age       = rows_dt$age_at_death,
    pmi_cat   = rows_dt$pmi_cat,
    col     = list(
      lesion    = les_cols,
      diagnosis = disease_cols,
      donor     = pat_cols,
      sex       = c(M = 'black', F = 'white'),
      age       = cols_fn(rows_dt$age_at_death, 10, 'Greys',
        pal_dir = 1, range='natural'),
      pmi_cat   = c(up_to_6H = '#f7fcf5', `6H_to_12H` = '#41ab5d',
        over_12H = '#005a32')
    ),
    annotation_name_side='top', show_legend = c(donor = FALSE)
  )

  # do column titles
  var_lookup  = c(
    log10_N       = "no. of cells",
    log_counts    = "no. of UMIs",
    log_feats     = "no. features",
    logit_mito    = "mito. pct.",
    logit_splice  = "pct. spliced"
  )

  # do heatmap for these genes
  name_str  = "median\nQC value\n(z-scored,\nsign flipped:\n+ve = good)"
  hm_obj    = Heatmap(
    matrix = z_mat, col = z_cols, cell_fun = labelling_fn,
    column_labels = var_lookup[ colnames(z_mat) ],
    cluster_rows = TRUE, cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 7),
    name = name_str, heatmap_legend_param = list(title = name_str),
    left_annotation = row_annots,
    show_row_names = FALSE, show_column_names = FALSE,
    # row_split = row_split, cluster_row_slices = FALSE,
    row_names_side = "left", column_names_side = "top")

  return(hm_obj)
}

calc_qc_summary <- function(qc_dt, kept_dt, cuts_dt, qc_lu, batch_var) {
  # get totals by sample
  tbl_tmp     = merge(
      qc_dt[, .(n_pre_QC = .N), by = .(batch_var)],
      kept_dt[, .(n_post_QC = .N), by = .(batch_var)],
      by = "batch_var", all = TRUE) %>%
    .[ is.na(n_post_QC), n_post_QC := 0 ] %>%
    .[, n_excluded      := n_pre_QC - n_post_QC ] %>%
    .[, pct_excluded    := round(100*(1 - n_post_QC / n_pre_QC), 1) ] %>%
    .[, batch_excluded := n_post_QC == 0 ]

  # get exclusions split by qc_var
  exclude_dt  = .calc_exclusions_by_qc_var(qc_dt, cuts_dt, qc_lu)

  # join together and make pretty
  tbl_pretty  = tbl_tmp %>% 
    .[, .(batch_var, `excluded` = batch_excluded, 
      `N pre-QC` = n_pre_QC, `N post-QC` = n_post_QC, 
      `N excluded` = n_excluded, `pct. excluded` = pct_excluded
    )] %>% 
    merge(exclude_dt, by = "batch_var") %>% 
    setnames("batch_var", batch_var) %>% 
    .[ order(-`pct. excluded`, -`N post-QC`) ]

  return(tbl_pretty)
}

# calc proportions excluded for each threshold
.calc_exclusions_by_qc_var <- function(qc_dt, cuts_dt, qc_lu) {
  # get relevant vars
  cut_vars    = intersect(cuts_dt$qc_var, colnames(qc_dt))
  qc_tmp      = qc_dt %>% 
    .[, c("batch_var", "cell_id", cut_vars), with = FALSE]

  # calc exclusion for each
  all_batches = unique(qc_tmp$batch_var)
  exclude_tmp = lapply(cut_vars, function(vv) {
    tmp_tmp     = lapply(all_batches, function(bb) {
      # get spec
      cuts_tmp    = cuts_dt[ (batch_var == bb) & (qc_var == vv) ]
      spec        = c(cuts_tmp$min, cuts_tmp$max)

      # find how many meet this
      tmp_dt  = qc_tmp[ batch_var == bb ] %>% 
        .[, .(.N, N_keep = sum(get(vv) %between% spec)), by = .(batch_var) ] %>% 
        .[, batch_var := bb ] %>% 
        .[, qc_var := vv ] %>% 
        .[, col_title := sprintf("pct. excluded by %s", qc_lu[[ vv ]]) ]
      }) %>% rbindlist
    }) %>% rbindlist %>% 
    .[, pct_exc   := round((N - N_keep) / N * 100, 1) ] %>% 
    .[, qc_var    := factor(qc_var, levels = names(qc_lu)) ] %>% 
    .[, col_title := fct_reorder(col_title, as.integer(qc_var)) ]

  # make nice
  exclude_dt  = exclude_tmp %>% 
    dcast( batch_var ~ col_title, value.var = "pct_exc" )

  return(exclude_dt)
}

############## for saving sce files. put here bc uses other functions


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

############## FUNCTIONS FROM SampleQC package

# from make_qc_dt.R
make_qc_dt <- function(qc_df, batch_var = 'sample_id', run_var = "sample_id", 
  qc_names = c('log_counts', 'log_feats', 'logit_mito'), annot_vars = NULL) {

  # some checks
  if ( 'DFrame' %in% class(qc_df) )
    qc_df      = as.data.frame(qc_df)
  assert_that( is.data.frame(qc_df), msg = "qc_df must be a data.frame" )

  assert_that( batch_var %in% colnames(qc_df),
    msg = sprintf("%s is listed as variable for samples but is not in data.frame", batch_var))
  if (batch_var != run_var)
    assert_that( run_var %in% colnames(qc_df),
      msg = sprintf("%s is listed as variable for samples but is not in data.frame", run_var))

  reserved_ns  = c(batch_var, 'group_id', 'cell_id')
  assert_that( length(intersect(annot_vars, reserved_ns)) == 0,
    msg = paste0("The following variable names are reserved and cannot be used as annot_vars:\n",
      paste(reserved_ns, collapse = ", ")))

  assert_that( all(annot_vars %in% names(qc_df)),
    msg = sprintf("the following variables are listed in annot_vars but not in qc_df:\n%s",
      paste(setdiff(annot_vars, names(qc_df)), collapse = ", ")))

  # set up qc_dt
  qc_dt   = .init_qc_dt(qc_df, batch_var, run_var)

  # add known metrics
  if ('log_counts' %in% qc_names) {
    qc_dt   = .add_log_counts(qc_dt, qc_df)
  }
  if ('log_feats' %in% qc_names) {
    qc_dt   = .add_log_feats(qc_dt, qc_df)
  }
  if ('logit_mito' %in% qc_names) {
    qc_dt   = .add_logit_mito(qc_dt, qc_df)
  }
  if ('splice_ratio' %in% qc_names) {
    qc_dt   = .add_splice_ratio(qc_dt, qc_df)
  }

  # add unknown metrics
  qc_dt   = .add_unknown_metrics(qc_dt, qc_df, qc_names)

  # add some useful annotations
  qc_dt   = .add_qc_annots(qc_dt, batch_var)

  # add specified annotation variables
  qc_dt   = .add_annot_vars(qc_dt, qc_df, annot_vars)

  # put in nice order
  tmp_vars  = c(batch_var, run_var) %>% unique
  setcolorder(qc_dt, c('cell_id', tmp_vars, qc_names))

  # double-check everything is ok
  .check_qc_dt(qc_dt, qc_names, annot_vars)

  return(qc_dt)
}


.init_qc_dt <- function(qc_df, batch_var, run_var) {
  # add cell identifiers
  if ('cell_id' %in% colnames(qc_df)) {
    qc_dt   = data.table(cell_id = qc_df$cell_id)
  } else if ( !is.null(rownames(qc_df)) ) {
    qc_dt   = data.table(cell_id = rownames(qc_df))
  } else {
    stop("input data.frame must have either rownames or 'cell_id' as a column")
  }
  assert_that( length(unique(qc_dt$cell_id)) == nrow(qc_dt),
               msg = "cell identifiers are not unique")

  # add sample identifiers
  qc_dt[, (batch_var) := qc_df[[batch_var]] ]
  if (batch_var != run_var)
    qc_dt[, (run_var) := qc_df[[run_var]] ]

  # check no missing values or NAs
  assert_that( all(!is.na(qc_dt$cell_id)), msg = "missing values in cell_id")
  assert_that( all(!is.na(qc_dt[[batch_var]])), msg = sprintf("missing values in %s", batch_var))

  return(qc_dt)
}


.add_log_counts <- function(qc_dt, qc_df) {
  # what names do we have, and want?
  df_names  = colnames(qc_df)
  valid_ns  = c('log_counts', 'total', 'sum', 'nCount_RNA')

  # check which are present
  here_ns   = vapply(valid_ns, function(v) v %in% df_names, logical(1))
  assert_that( sum(here_ns) >= 1,
               msg = paste0(
                 "no valid column present for log_counts\n",
                 paste0("valid columns are: ", paste(valid_ns, collapse = ", "))
               ))
  to_use    = valid_ns[here_ns][[1]]

  # add values
  if (to_use %in% 'log_counts') {
    qc_dt[, log_counts := qc_df[[ to_use ]] ]

  } else if (to_use %in% c('total', 'sum', 'nCount_RNA')) {
    assert_that( all(qc_df[[ to_use ]] > 0) )
    qc_dt[, log_counts := log10(qc_df[[ to_use ]]) ]

  } else {
    stop("log_counts requested but required variables not present")

  }

  # do some checks
  assert_that( "log_counts" %in% names(qc_dt) )
  assert_that( !any(is.na(qc_dt$log_counts)),
               msg = "some log_counts values are NA")
  assert_that( !any(is.infinite(qc_dt$log_counts)),
               msg = "some log_counts values are infinite")
  assert_that( all(qc_dt$log_counts >= 0),
               msg = "some log_counts values are <= 0")

  return(qc_dt)
}


.add_log_feats <- function(qc_dt, qc_df) {
  # what names do we have, and want?
  df_names  = colnames(qc_df)
  valid_ns  = c('log_feats', 'detected', 'nFeature_RNA')

  # check which are present
  here_ns   = vapply(valid_ns, function(v) v %in% df_names, logical(1))
  assert_that( sum(here_ns) >= 1,
               msg = paste0(
                 "no valid column present for log_feats\n",
                 paste0("valid columns are: ", paste(valid_ns, collapse = ", "))
               ))
  to_use    = valid_ns[here_ns][[1]]

  # add values
  if (to_use %in% 'log_feats') {
    qc_dt[, log_feats := qc_df[[ to_use ]] ]

  } else if (to_use %in% c('detected', 'nFeature_RNA')) {
    assert_that( all(qc_df[[ to_use ]] > 0) )
    qc_dt[, log_feats := log10(qc_df[[ to_use ]]) ]

  } else {
    stop("log_feats requested but required variables not present")

  }

  # do some checks
  assert_that( "log_feats" %in% names(qc_dt) )
  assert_that( !any(is.na(qc_dt$log_feats)),
               msg = "some log_feats values are NA")
  assert_that( !any(is.infinite(qc_dt$log_feats)),
               msg = "some log_feats values are infinite")
  assert_that( all(qc_dt$log_feats >= 0),
               msg = "some log_feats values are <= 0")

  return(qc_dt)
}


.add_logit_mito <- function(qc_dt, qc_df) {
  # what names do we have, and want?
  df_names  = colnames(qc_df)

  # add logit-transformed mitochondrial proportion to qc_dt
  if ('logit_mito' %in% df_names) {
    qc_dt[, logit_mito  := qc_df$logit_mito ]

  } else if ( ('subsets_mito_sum' %in% df_names) & ('total' %in% df_names) ) {
    qc_dt[, logit_mito  := qlogis( (qc_df$subsets_mito_sum + 1) / (qc_df$total + 2) ) ]

  } else if ( ('subsets_mt_sum' %in% df_names) & ('total' %in% df_names) ) {
    qc_dt[, logit_mito  := qlogis( (qc_df$subsets_mt_sum + 1) / (qc_df$total + 2) ) ]

  } else if ( ('percent.mt' %in% df_names) & ('nCount_RNA' %in% df_names) ) {
    total_counts  = qc_df$nCount_RNA
    mt_counts     = qc_df$nCount_RNA * qc_df$percent.mt / 100
    assert_that( all(abs(mt_counts - round(mt_counts, 0)) < 1e-10) )
    qc_dt[, logit_mito  := qlogis( (mt_counts + 1) / (total_counts + 2) ) ]

  } else if ( ('mito_prop' %in% df_names) & ('log_counts' %in% df_names) ) {
    total_counts  = 10^qc_df$log_counts
    mt_counts     = qc_df$mito_prop * total_counts
    assert_that( all(abs(mt_counts - round(mt_counts, 0)) < 1e-8) )
    qc_dt[, logit_mito  := qlogis( (mt_counts + 1) / (total_counts + 2) ) ]

  } else {
    stop("logit_mito requested but required variables not present")
  }

  # do some checks
  assert_that( "logit_mito" %in% names(qc_dt) )
  assert_that( !any(is.na(qc_dt$logit_mito)),
               msg = "some logit_mito values are NA")
  assert_that( !any(is.infinite(qc_dt$logit_mito)),
               msg = "some logit_mito values are infinite")

  return(qc_dt)
}


.add_splice_ratio <- function(qc_dt, qc_df) {
  # what names do we have, and want?
  df_names  = colnames(qc_df)

  # add logit-transformed mitochondrial proportion to qc_dt
  if ('splice_ratio' %in% df_names) {
    qc_dt[, splice_ratio  := qc_df$splice_ratio ]

  } else if ( ('total_spliced' %in% df_names) & ('total_unspliced' %in% df_names) ) {
    qc_dt[, splice_ratio  := log2( (qc_df$total_spliced + 1) / (qc_df$total_unspliced + 1) ) ]

  } else {
    stop("logit_mito requested but required variables not present")

  }

  # do some checks
  assert_that( "splice_ratio" %in% names(qc_dt) )
  assert_that( !any(is.na(qc_dt$logit_mito)),
               msg = "some logit_mito values are NA")
  assert_that( !any(is.infinite(qc_dt$logit_mito)),
               msg = "some logit_mito values are infinite")

  return(qc_dt)
}


list_known_metrics <- function() {
  return(c('log_counts', 'log_feats', 'logit_mito', 'splice_ratio'))
}


.add_unknown_metrics <- function(qc_dt, qc_df, qc_names)  {
  # anything to add?
  to_add  = setdiff(qc_names, list_known_metrics())
  if ( length(to_add) == 0 )
    return(qc_dt)

  # add them
  message("adding the following metrics that are not known to `SampleQC`:")
  message(paste(to_add, collapse = ", "))
  for (v in to_add) {
    assert_that( v %in% names(qc_df), msg = paste0(v, " missing from qc_df"))
    set(qc_dt, i = NULL, v, qc_df[[v]])
    assert_that( !any(is.na(qc_dt$v)), msg = paste0("NA values for ", v))
    assert_that( !any(is.infinite(qc_dt$v)), msg = paste0("infinite values for ", v))
  }

  return(qc_dt)
}


.add_qc_annots <- function(qc_dt, batch_var) {
  # add annotations for sample size
  qc_dt[, log_N  := log10(.N), by=batch_var]

  # and factor version
  N_cuts    = c(1,100,200,400,1000,2000,4000,10000,20000,40000,Inf)
  N_labs    = paste0('<=', N_cuts[-1])
  qc_dt[, N_cat  := factor(
    cut(10^log_N, breaks = N_cuts, labels = N_labs),
    levels = N_labs), by = batch_var]

  # add annotations relating to no. of UMIss
  if ('log_counts' %in% names(qc_dt) ) {
    # add median log counts per sample
    qc_dt[, med_counts  := median(log_counts), by=batch_var]

    # put mito level into categories
    counts_cuts = c(1,100,300,1000,3000,10000,30000, Inf)
    counts_labs = paste0('<=', counts_cuts[-1])
    qc_dt[, counts_cat  := factor(
      cut(10^med_counts, breaks = counts_cuts, labels = counts_labs),
      levels = counts_labs), by = batch_var]
  }

  # add annotations relating to features
  if ('log_feats' %in% names(qc_dt) ) {
    # add median log feats per sample
    qc_dt[, med_feats   := median(log_feats), by=batch_var]

    # put mito level into categories
    feats_cuts = c(1,100,300,1000,3000,10000,30000, Inf)
    feats_labs = paste0('<=', feats_cuts[-1])
    qc_dt[, feats_cat  := factor(
      cut(10^med_feats, breaks = feats_cuts, labels = feats_labs),
      levels = feats_labs), by = batch_var]
  }

  # add annotations relating to mitochondrial proportions
  if ('logit_mito' %in% names(qc_dt) ) {
    # add median mito proportion
    qc_dt[, med_mito  := median(plogis(logit_mito)), by=batch_var]

    # put mito level into categories
    mito_cuts   = c(0,0.01,0.05,0.1,0.2,0.5,1)
    mito_labs   = paste0('<=', mito_cuts[-1])
    qc_dt[, mito_cat  := factor(
      cut(med_mito, breaks = mito_cuts, labels = mito_labs),
      levels = mito_labs), by = batch_var]
  }

  # add annotations relating to mitochondrial proportions
  if ('splice_ratio' %in% names(qc_dt) ) {
    # add median mito proportion
    qc_dt[, med_splice  := median(plogis(splice_ratio)), by=batch_var]

    # put mito level into categories
    splice_cuts = c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 1)
    splice_labs = paste0('<=', splice_cuts[-1])
    qc_dt[, splice_cat  := factor(
      cut(med_splice, breaks = splice_cuts, labels = splice_labs),
      levels = splice_labs), by = batch_var]
  }

  return(qc_dt)
}


.add_annot_vars <- function(qc_dt, qc_df, annot_vars) {
  # check all present
  assert_that( all(annot_vars %in% names(qc_df)) )
  # add them
  for (v in annot_vars)
    qc_dt[[v]] = qc_df[[v]]

  # check that they all sample level
  for (v in annot_vars) {
    check_dt  = qc_dt[, c('sample_id', v), with = FALSE] %>%
      .[, .N, by = c('sample_id', v) ]
    assert_that( nrow(check_dt) == length(unique(check_dt$sample_id)),
                 msg = paste0("annotation variable ", v, " has more than one value per\n",
                              "sample (should be sample-level only)"))
  }

  return(qc_dt)
}


.check_qc_dt <- function(qc_dt, qc_names, annot_vars) {
  # unpack
  col_names   = colnames(qc_dt)

  # check specific names
  if ('log_counts' %in% col_names)
    assert_that( all(qc_dt$log_counts >= 0) )
  if ('log_feats' %in% col_names)
    assert_that( all(qc_dt$log_feats >= 0) )
  if ('logit_mito' %in% col_names)
    assert_that( all(is.finite(qc_dt$logit_mito)) )
  if ('splice_ratio' %in% col_names)
    assert_that( all(is.finite(qc_dt$splice_ratio)) )

  # check qc metrics and annotations for NAs
  for (n in qc_names) {
    assert_that( all(!is.na(qc_dt[[n]])) )
  }
  annots_auto   = c(
    "med_counts", "counts_cat",
    "med_feats", "feats_cat",
    "med_mito", "mito_cat",
    "med_splice", "splice_cat",
    "log_N", "N_cat")
  for (n in c(annots_auto, annot_vars)) {
    if ( n %in% names(qc_dt) )
      assert_that( all(!is.na(qc_dt[[n]])),
                   msg = paste0('NA present in an annotation variable, ', n) )
  }
}

plot_upset_of_exclusions <- function(qc_tmp, qc_names, qc_lu, cuts_dt) {
  # make list for upsets
  qc_tmp    = qc_tmp[ keep_hard == TRUE ]
  var_lu    = c(
    log_counts    = "umis", 
    log_feats     = "genes", 
    logit_mito    = "mito", 
    logit_spliced = "spliced"
  )
  eps       = 1e-10

  # get relevant vars
  cuts_tmp  = cuts_dt[ batch_var == unique(qc_tmp$batch_var) ] %>% 
    .[ qc_var != "n_cells" ]
  cut_vars  = unique(cuts_tmp$qc_var)
  tmp_ls    = lapply(cut_vars, function(cc) {
    # find spec
    tmp_dt    = cuts_tmp[ qc_var == cc ]
    cut_spec  = c(cuts_tmp[ qc_var == cc ]$min, cuts_tmp[ qc_var == cc ]$max)

    # find excluded cells
    exc_cells   = list(
      low   = qc_tmp[ get(cc) + eps < cut_spec[1] ]$cell_id,
      high  = qc_tmp[ get(cc) - eps > cut_spec[2] ]$cell_id
    )
    # make names nice
    names(exc_cells)  = names(exc_cells) %>% paste(var_lu[[ cc ]], sep = "_")

    # remove anything that is blank
    n_cells   = sapply(exc_cells, length)
    exc_cells = exc_cells[ n_cells > 0 ]

    return(exc_cells)
    }) %>% do.call(c, .)

  # check no overlaps with good cells
  ok_cells  = qc_tmp[ keep == TRUE ]$cell_id
  assert_that( all(sapply(tmp_ls, function(l) length(intersect(l, ok_cells))) == 0),
    msg = "some overlap between passed and excluded barcodes" )

  # turn into full list
  upset_ls  = tmp_ls %>% c( list(passed_qc = ok_cells) )
  upset_dt  = names(upset_ls) %>% lapply(function(nn) 
      data.table(set = nn, cell_id = upset_ls[[nn]])) %>% 
    rbindlist %>% .[, dummy := 1 ] %>% 
    dcast( cell_id ~ set, value.var = "dummy", fill = 0 )

  # make ratios nice
  n_cols    = ncol(upset_dt)
  mat_prop  = 0.4

  # do nicer colours for up / down
  row_ord   = upset_dt[, -c('cell_id')] %>% as.matrix %>% colSums %>%
    sort(decreasing = TRUE) %>% names
  row_cols  = rep("#FB8072", length(row_ord)) %>% setNames(row_ord)
  row_cols[ "passed_qc" ] = "#7BAFDE"

  # plot upset
  return(upset(upset_dt, sets = colnames(upset_dt)[-1], order.by = 'freq',
     mb.ratio = c(1 - mat_prop, mat_prop), sets.bar.color = row_cols))   
}
