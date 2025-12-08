suppressPackageStartupMessages({
  library("Matrix")
  library("SingleCellExperiment")
  library("edgeR")
  library("DESeq2")
  library("scater")
  library("BiocParallel")
  RhpcBLASctl::omp_set_num_threads(1L)
})


calculate_marker_genes <- function(integration_f, sces_yaml_f, pb_f, mkrs_f, pb_hvgs_f,
  fgsea_go_bp_f = "", fgsea_go_cc_f = "", fgsea_go_mf_f = "", fgsea_paths_f = "", fgsea_hlmk_f = "",
  do_gsea, gtf_dt_f, gsea_dir, sel_res, min_cl_size, min_cells, zoom = FALSE, 
  n_cores = 4) {
  # define some variables
  batch_var   = "sample_id"
  zoom        = as.logical(zoom)
  
  # check some inputs
  assert_that(
    is.character(sces_yaml_f),
    is.character(sel_res),
    is.character(pb_f),
    is.character(mkrs_f),
    is.character(pb_hvgs_f),
    is.character(gtf_dt_f)
  )
  assert_that(
    file.exists(sces_yaml_f),
    file.exists(gtf_dt_f)
  )
  assert_that(
    is.numeric(min_cl_size),
    is.numeric(min_cells),
    is.numeric(n_cores)
  )
  setDTthreads(n_cores)

  # load gene biotypes
  message("calculating marker genes")
  message("  loading gene biotypes")
  biotypes_dt = load_gene_biotypes(gtf_dt_f)

  # make_pb_object
  message("  making pseudobulk object")
  pb          = make_pseudobulk_object(pb_f, integration_f, sces_yaml_f, sel_res, batch_var,
    min_cl_size = min_cl_size, agg_fn = "sum", zoom = zoom, n_cores = n_cores)

  # calc cpms
  message("  calculating logCPMs")
  cpms_dt     = pb %>%
    make_logcpms_all(batch_var, lib_size_method = "edger",
      min_cells = min_cells, n_cores = n_cores) %>%
    merge(biotypes_dt[, .(gene_id, gene_type)], by = "gene_id")
  rm(pb); gc()

  # calc hvgs
  message("  calculating highly variable genes")
  calc_hvgs_pseudobulk(pb_hvgs_f, cpms_dt, batch_var, n_cores = n_cores)

  # calculate markers
  message("  calculating markers")
  mkrs_dt     = calc_find_markers_pseudobulk(mkrs_f, cpms_dt, biotypes_dt, batch_var,
    n_cores = n_cores)

  message("done!")
}

make_pseudobulk_object <- function(pb_f, integration_f, sces_yaml_f, sel_res, batch_var,
  min_cl_size = 1e2, agg_fn = c("sum", "prop.detected"), zoom = FALSE, n_cores = 8) {
  # check inputs
  agg_fn      = match.arg(agg_fn)

  # load up clusters
  message('    loading integration output')
  int_dt     = fread(integration_f)
  if (zoom == FALSE) {
    #exclude doublets
    int_dt = int_dt %>%
      .[ (is_dbl == FALSE) & (in_dbl_cl == FALSE) ]
  }

  # exclude tiny clusters
  message('    excluding tiny clusters')  
  cl_var      = paste0("RNA_snn_res.", sel_res)
  assert_that( cl_var %in% names(int_dt))
  cl_ns       = int_dt[[ cl_var ]] %>% table
  keep_cls    = names(cl_ns)[ cl_ns >= min_cl_size ]
  int_dt     = int_dt %>%
    .[ get(cl_var) %in% keep_cls ]
  
  # make pseudobulks for selected clusters for each sample
  message('    making pseudobulk counts for individual samples')
  batches     = int_dt[[ batch_var ]] %>% unique()
  sce_paths   = yaml::read_yaml(sces_yaml_f)
  
  # set up cluster
  bpparam     = MulticoreParam(workers = n_cores, tasks = length(batches))  
  if (zoom) {
    pb_ls       = bplapply(batches, FUN = .make_one_zoom_pseudobulk, BPPARAM = bpparam, 
      sce_paths = sce_paths, int_dt = int_dt, batch_var = batch_var, cl_var = cl_var, keep_cls = keep_cls, agg_fn = agg_fn)
  } else {
    pb_ls       = lapply(batches, FUN = .make_one_pseudobulk, sce_paths = sce_paths, 
      batch_var = batch_var, cl_var = cl_var, keep_cls = keep_cls, agg_fn = agg_fn)
    # pb_ls       = bplapply(batches, FUN = .make_one_pseudobulk, BPPARAM = bpparam, 
    #   sce_paths = sce_paths, batch_var = batch_var, cl_var = cl_var, keep_cls = keep_cls, agg_fn = agg_fn)
  }
  
  message('    merging pseudobulk counts')
  assay_ls    = lapply(keep_cls, function(cl){
    assays      = lapply(pb_ls, function(pb) assay(pb, cl))
    assay_mat   = Reduce(cbind, assays)  
    return(assay_mat)
    }) %>% setNames(keep_cls)
  
  n_cells_ls  = sapply(pb_ls, function(pb) int_colData(pb)$n_cells)
  
  pb          = SingleCellExperiment(assays = assay_ls)
  pb@metadata$agg_pars = list(
    assay = "counts", 
    by    = c("cluster", batch_var),
    fun   = 'sum'
  )
  int_colData(pb)$n_cells = n_cells_ls
  rowData(pb) = rowData(pb_ls[[1]])
  
  # save results
  message('   saving outputs')
  saveRDS(pb, file = pb_f, compress = FALSE)

  return(pb)
}

.make_one_pseudobulk <- function(sel_b, sce_paths, batch_var, cl_var, keep_cls, agg_fn) {
  message(sel_b)
  sce_f     = sce_paths[[sel_b]]
  tmp_sce   = readRDS(sce_f)
  
  # add clusters to sce
  colData(tmp_sce)[["cluster"]] = colData(tmp_sce)[[cl_var]] 

  # filter sce
  tmp_sce   = tmp_sce[, colData(tmp_sce)$cluster %in% keep_cls]
  pb        = aggregateData_datatable(tmp_sce, batch_var, by_vars = c("cluster", batch_var), 
    fun = agg_fn, all_cls = keep_cls)
  
  # make sure all assays are in the object, if not, add columns with zeros
  missing_assays = setdiff(keep_cls, assayNames(pb))  
  if ( length(missing_assays) > 0 ) {
    message('  adding assays with zero counts')
    for(assay in missing_assays){
      missing_counts  = Matrix(0, nrow = nrow(pb), ncol = 1, sparse = FALSE, 
        dimnames = list(rownames(pb), sel_b))
      assay(pb, assay) = missing_counts
    }
  }
  
  message('done!')
  return(pb)
}

.make_one_zoom_pseudobulk <- function(sel_b, sce_paths, int_dt, batch_var, cl_var, keep_cls, agg_fn) {
  message(sel_b)
  sce_f       = sce_paths[[sel_b]]
  tmp_sce     = readRDS(sce_f)
  smpl_int_dt = copy(int_dt) %>% .[ get(batch_var) == sel_b ]
  
  # remove umap and clustering cols from before and add new ones
  rm_cols = c('UMAP1', 'UMAP2', str_subset(names(colData(tmp_sce)), "RNA_snn_res"))
  new_coldata = colData(tmp_sce) %>% as.data.table %>%
    .[ , (rm_cols) := NULL] %>%
    as.data.frame() %>%
    set_rownames(.$cell_id)
  
  colData(tmp_sce) = DataFrame(new_coldata)
  
  assert_that(all(smpl_int_dt$cell_id %in% colnames(tmp_sce)))
  # reorder cells in sce
  tmp_sce = tmp_sce[, smpl_int_dt$cell_id]
  
  # add clusters to sce
  colData(tmp_sce)[["cluster"]] = smpl_int_dt[[cl_var]] 
  
  pb  = aggregateData_datatable(tmp_sce, batch_var, by_vars = c("cluster", batch_var), 
                                fun = agg_fn, all_cls = keep_cls)
  
  # make sure all assays are in the object, if not, add columns with zeros
  missing_assays = setdiff(keep_cls, assayNames(pb))  
  if ( length(missing_assays) > 0 ) {
    message('  adding assays with zero counts')
    for(assay in missing_assays){
      missing_counts  = Matrix(0, nrow = nrow(pb), ncol = 1, 
                               sparse = FALSE, dimnames = list(rownames(pb), sel_b))
      assay(pb, assay) = missing_counts
    }
  }
  
  message('done!')
  return(pb)
}

aggregateData_datatable <- function(sce, batch_var, by_vars = c("cluster", "sample_id"),
  fun = c("sum", "mean", "median", "prop.detected", "num.detected"), all_cls) {
  fun       = match.arg(fun)
  assay     = "counts"
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
  by_1      = by_vars[[1]]
  by_2      = by_vars[[2]]
  by_dt     = data.table(
    cell_id   = colnames(sce),
    cluster   = factor(sce[[ by_1 ]]),
    batch_var = factor(sce[[ by_2 ]])
    ) %>% setkey("cell_id") %>% setnames("batch_var", batch_var)

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

  # define various things
  genes_ls    = rownames(sce)
  cl_ls       = levels(by_dt$cluster)
  batches     = levels(by_dt[[ batch_var ]])
  n_genes     = length(genes_ls)
  n_batches   = length(batches)

  # do pseudobulk
  message('  assembling list of matrices')
  cast_frm    = sprintf("gene_id ~ %s", batch_var)
  mat_ls      = lapply(cl_ls, function(cl) {
    # get matrix for this cluster
    outs_mat    = pb_dt[ cluster == cl ] %>%
      dcast.data.table( as.formula(cast_frm), value.var = fun, fill = 0 ) %>%
      as.matrix(rownames = "gene_id")

    # make full 0 matrix, fill it in
    mat         = matrix(0, nrow=n_genes, ncol=n_batches) %>%
        set_rownames(genes_ls) %>%
        set_colnames(batches)
    mat[ rownames(outs_mat), colnames(outs_mat) ] = outs_mat

    return(mat)
    }) %>% setNames(cl_ls)

  message('  making pseudobulk object')
  md          = metadata(sce)
  md$agg_pars = list(assay = assay, by = by_vars, fun = fun, scale = scale)
  pb          = SingleCellExperiment(mat_ls,
    rowData = rowData(sce), metadata = md)
  cd          = data.frame(colData(sce)[, by_vars])
  # for (i in names(cd)) if (is.factor(cd[[i]]))
  #   cd[[i]]     = droplevels(cd[[i]])
  cd$cluster  = factor(cd$cluster, levels = all_cls)
  ns          = table(cd)
  if (length(by_vars) == 2) {
    ns    = asplit(ns, 2)
    ns    = purrr::map(ns, ~c(unclass(.)))
  } else {
    ns     = c(unclass(ns))
  }
  int_colData(pb)$n_cells = ns
  if (length(by_vars) == 2) {
    cd      = colData(sce)
    ids     = colnames(pb)
    counts  = vapply(ids, function(u) {
      m     = as.logical(match(cd[, by_vars[2]], u, nomatch = 0))
      vapply(cd[m, ], function(u) length(unique(u)), numeric(1))
    }, numeric(ncol(colData(sce))))
    cd_keep     = apply(counts, 1, function(u) all(u == 1))
    cd_keep     = setdiff(names(which(cd_keep)), by_vars)
    if (length(cd_keep) != 0) {
      m     = match(ids, cd[, by_vars[2]], nomatch = 0)
      cd    = cd[m, cd_keep, drop = FALSE]
      rownames(cd)    = ids
      colData(pb)     = cd
    }
  }

  # add helpful row data
  assert_that( all(rownames(pb) == rownames(sce)) )
  return(pb)
}

make_logcpms_all_rmd <- function(pb, batch_var, lib_size_method = c("edger", "raw", "pearson",
  "vst", "rlog"), exc_regex = NULL, min_cells = 10, n_cores = 4) {  
  # check inputs
  lib_size_method   = match.arg(lib_size_method)
  
  # set up cluster
  cl_ls       = assayNames(pb)
 
  # exclude genes if requested
  if (!is.null(exc_regex)) {
    exc_idx     = rownames(pb) %>% str_detect(exc_regex)
    exc_gs_str  = rowData(pb)$symbol[ exc_idx ] %>% paste0(collapse = " ")
    sprintf("    excluding %d genes: %s", sum(exc_idx), exc_gs_str) %>% message
    pb          = pb[ !exc_idx, ]
  }
  
  # calculate logcpms
  logcpms_all = lapply(cl_ls, function(sel_cl) {
    message(sel_cl, appendLF = FALSE)
    # message(sel_cl)
    tmp_dt    = .get_logcpm_dt_one_cl(pb, batch_var, cl = sel_cl, min_cells = min_cells, 
      lib_size_method = lib_size_method)
    if (!is.null(tmp_dt))
      tmp_dt   = tmp_dt[, cluster := sel_cl ]
    return(tmp_dt)
  }) %>% rbindlist
  
  # add # cells
  ncells_dt   = muscat_n_cells(pb) %>%
    as.data.table %>% set_colnames(c("cluster", batch_var, "n_cells"))
  logcpms_all = merge( logcpms_all, ncells_dt, by = c("cluster", batch_var) )
  
  return(logcpms_all)
}

make_logcpms_all <- function(pb, batch_var, lib_size_method = c("edger", "raw", "pearson",
  "vst", "rlog"), exc_regex = NULL, min_cells = 10, n_cores = 4) {  
  # check inputs
  lib_size_method   = match.arg(lib_size_method)

  # set up cluster
  cl_ls       = assayNames(pb)
  bpparam     = MulticoreParam(workers = n_cores, tasks = length(cl_ls))
  register(bpparam)
  on.exit(bpstop(bpparam))

  # exclude genes if requested
  if (!is.null(exc_regex)) {
    exc_idx     = rownames(pb) %>% str_detect(exc_regex)
    exc_gs_str  = rowData(pb)$symbol[ exc_idx ] %>% paste0(collapse = " ")
    sprintf("    excluding %d genes: %s", sum(exc_idx), exc_gs_str) %>% message
    pb          = pb[ !exc_idx, ]
  }

  # calculate logcpms
  logcpms_all = bplapply(cl_ls, function(sel_cl) {
      message(sel_cl, appendLF = FALSE)
      # message(sel_cl)
      tmp_dt    = .get_logcpm_dt_one_cl(pb, batch_var, cl = sel_cl,
        min_cells = min_cells, lib_size_method = lib_size_method)
      if (!is.null(tmp_dt))
       tmp_dt   = tmp_dt[, cluster := sel_cl ]
      return(tmp_dt)
    }, BPPARAM = bpparam) %>% rbindlist

  # add # cells
  ncells_dt   = muscat_n_cells(pb) %>%
    as.data.table %>% set_colnames(c("cluster", batch_var, "n_cells"))
  logcpms_all = merge( logcpms_all, ncells_dt, by = c("cluster", batch_var) )

  return(logcpms_all)
}

.get_logcpm_dt_one_cl <- function(pb, batch_var, cl, min_cells = 10, pseudo_count = 10,
  lib_size_method = c("raw", "edger", "pearson", "vst", "rlog")) {

  # check inputs
  lib_size_method   = match.arg(lib_size_method)

  # extract raw counts for non-tiny samples
  use_idx     = muscat_n_cells(pb)[cl, ] >= min_cells
  if (sum(use_idx) == 0) {
    message("no samples with sufficient cells for ", cl, "; skipping")
    return(NULL)
  }
  x           = assay(pb[, use_idx], cl)

  # exclude tiny samples
  ls          = colSums(x)
  out_idx     = scater::isOutlier(ls, log = TRUE, type = "lower", nmads = 3)
  x           = x[, !out_idx, drop = FALSE]

  # do DESeq2 options
  if (lib_size_method %in% c("rlog", "vst")) {
    # do DESeq2
    dds   = DESeq2::DESeqDataSetFromMatrix(countData = x,
      colData = data.frame(dummy = rep(1, ncol(x))), design = ~ 1)
    if (lib_size_method == "vst") {
      mat   = assay(DESeq2::vst(dds, blind = TRUE))
    } else if (lib_size_method == "rlog") {
      mat   = assay(DESeq2::rlog(dds, blind = TRUE))
    }

    # convert this to output data.table
    logcpm_dt   = mat %>%
      as.data.table(keep.rownames = "gene_id") %>%
      melt.data.table(id = "gene_id", value.name = "logcpm", variable.name = batch_var)

    return(logcpm_dt)
  }

  # calculate library sizes
  if (lib_size_method %in% c("edger")) {
    # make DGE object, get library sizes
    suppressMessages({
      dge_obj     = DGEList(x, remove.zeros = TRUE) %>% normLibSizes
    })
    lib_sizes   = getNormLibSizes(dge_obj)
  } else if (lib_size_method == "raw") {
    # just take raw column totals
    lib_sizes   = colSums(x)
  }

  # calc library sizes
  libsizes_dt = data.table( batch_var = colnames(x), lib_size = lib_sizes ) %>% 
    setnames("batch_var", batch_var)

  # calculate logcpms
  logcpm_dt   = as.matrix(x) %>%
    as.data.table(keep.rownames = "gene_id") %>%
    melt.data.table(id = "gene_id", value.name = "count", variable.name = batch_var) %>%
    merge(libsizes_dt, by = batch_var) %>%
    .[, logcpm  := log(count / lib_size * 1e6 + pseudo_count) ] %>%
    .[, symbol  := str_extract(gene_id, "^[^_]+") ]

  return(logcpm_dt)
}

make_props_dt <- function(pb_prop, batch_var, exc_regex = NULL, min_cells = 10, n_cores = 4) {
  # set up cluster
  cl_ls       = assayNames(pb_prop)
  bpparam     = MulticoreParam(workers = n_cores, tasks = length(cl_ls), progressbar = TRUE)
  register(bpparam)
  on.exit(bpstop(bpparam))

  # exclude genes if requested
  if (!is.null(exc_regex)) {
    exc_idx     = rownames(pb_prop) %>% str_detect(exc_regex)
    exc_gs_str  = rowData(pb_prop)$symbol[ exc_idx ] %>% paste0(collapse = " ")
    sprintf("    excluding %d genes: %s", sum(exc_idx), exc_gs_str) %>% message
    pb_prop     = pb_prop[ !exc_idx, ]
  }

  # calculate logcpms
  props_dt    = bplapply(cl_ls, function(sel_cl) {
      # extract raw counts for non-tiny samples
      use_idx     = muscat_n_cells(pb_prop)[sel_cl, ] >= min_cells
      if (sum(use_idx) == 0) {
        message("no samples with sufficient cells for ", sel_cl, "; skipping")
        return(NULL)
      }
      x           = assay(pb_prop[, use_idx], sel_cl)

      # make data.table
      props_tmp   = x %>%
        as.data.table(keep.rownames = "gene_id") %>%
        melt.data.table(id = "gene_id", value.name = "prop", variable.name = batch_var) %>%
        .[, cluster := sel_cl ] %>%
        .[, symbol  := str_extract(gene_id, "^[^_]+") ]

      return(props_tmp)
    }, BPPARAM = bpparam) %>% rbindlist

  # add # cells
  ncells_dt   = muscat_n_cells(pb_prop) %>%
    as.data.table %>% set_colnames(c("cluster", batch_var, "n_cells"))
  props_dt    = merge( props_dt, ncells_dt, by = c("cluster", batch_var) )

  return(props_dt)
}

calc_hvgs_pseudobulk <- function(pb_hvgs_f, cpms_dt, batch_var, n_cores = 8) {
  # remove outliers
  message('calculating highly variable genes')
  cl_ls     = cpms_dt$cluster %>% unique
  message('  start cluster')
  bpparam   = MulticoreParam(workers = n_cores, tasks = length(cl_ls))
  register(bpparam)
  on.exit(bpstop(bpparam))

  message('  remove outlier samples')
  x_ls      = cl_ls %>% bplapply(function(cl) {
    # make matrix
    this_x    = cpms_dt[ cluster == cl ] %>%
      .[, col_lab := sprintf("%s-%s", cluster, get(batch_var)) ] %>%
      dcast.data.table( gene_id ~ col_lab, value.var = 'count' ) %>%
      as.matrix(rownames = 'gene_id')

    # remove samples with outlier library sizes
    ls_dt     = cpms_dt[ cluster == cl ] %>% 
      .[, c("cluster", batch_var, "lib_size"), with = FALSE ] %>%
      unique %>% .[, col_lab := sprintf("%s-%s", cluster, get(batch_var)) ] %>%
      .[ order(col_lab) ]
    ls        = ls_dt$lib_size
    ol        = scater::isOutlier(ls, log = TRUE, type = "lower", nmads = 3)
    this_x    = this_x[, !ol, drop = FALSE]

    return(this_x)
  }, BPPARAM = bpparam) %>% setNames(cl_ls)

  # make big matrix
  message('  make big matrix')
  bulk_mat  = do.call(cbind, x_ls)

  message('  run VST')
  suppressMessages({
    dds       = DESeqDataSetFromMatrix(countData = bulk_mat,
      colData = data.frame(dummy = rep(1, ncol(bulk_mat))), design = ~ 1)
  })
  # sometimes it doesn't work
  vst_obj     = tryCatch({
    vst(dds, blind = TRUE)
    }, error = function(cond) {
    varianceStabilizingTransformation(dds, blind = TRUE)
    })

  # get variances, add biotypes
  biotypes_dt = cpms_dt[, .(gene_id, symbol, gene_type)] %>% unique
  hvgs_dt     = data.table(
    gene_id     = rownames(vst_obj),
    vst_var     = rowVars(assay(vst_obj))
    ) %>%
    merge(biotypes_dt, by = "gene_id") %>% .[ order(-vst_var) ]

  # save results
  fwrite(hvgs_dt, file = pb_hvgs_f)
  message('done!')

  return(hvgs_dt)
}

calc_find_markers_pseudobulk <- function(mkrs_pb_f, logcpms_all, rows_dt, batch_var,
  method = c("edger", "voom"), n_cores = n_cores) {
  method    = match.arg(method)

  cl_ls     = logcpms_all$cluster %>% unique
  message('  start cluster')
  bpparam   = MulticoreParam(workers = n_cores, tasks = length(cl_ls))
  register(bpparam)
  on.exit(bpstop(bpparam))
  message('  remove outlier samples')
  # x_ls      = cl_ls %>% bplapply(function(cl) {
  x_ls      = cl_ls %>% lapply(function(cl) {
    # make matrix
    this_x    = logcpms_all[ cluster == cl ] %>%
      .[, col_lab := sprintf("%s-%s", cluster, get(batch_var)) ] %>%
      dcast.data.table( gene_id ~ col_lab, value.var = 'count' ) %>%
      as.matrix(rownames = 'gene_id')

    # remove samples with outlier library sizes
    ls_dt     = logcpms_all[ cluster == cl ] %>% 
      .[, c("cluster", batch_var, "lib_size"), with = FALSE ] %>%
      unique %>% .[, col_lab := sprintf("%s-%s", cluster, get(batch_var)) ] %>%
      .[ order(col_lab) ]
    ls        = ls_dt$lib_size
    ol        = scater::isOutlier(ls, log = TRUE, type = "lower", nmads = 3)
    this_x    = this_x[, !ol, drop = FALSE]

    return(this_x)
  }) %>% setNames(cl_ls)
  # }, BPPARAM = bpparam) %>% setNames(cl_ls)
  # bpstop(bpparam)

  # make big matrix
  message('  make big DGE object')
  x_full    = do.call(cbind, x_ls)
  dge_all   = DGEList(x_full, remove.zeros = TRUE)

  # make design matrix
  message('  set up design')
  col_ns    = colnames(dge_all)
  des_all   = data.table(
    cluster   = str_extract(col_ns, "^[^-]+"),
    batch_var = str_extract(col_ns, "[^-]+$")
  )

  # filter out tiny genes
  message('  filter genes and estimate library sizes')
  mm_all    = model.matrix( ~ cluster, data = des_all )
  keep_gs   = filterByExpr(dge_all, group = des_all$cluster, min.count = 1)
  dge       = dge_all[ keep_gs, ]

  # calculate library sizes
  dge_all   = normLibSizes(dge_all, method = "TMMwsp")

  # run DE
  if (method == "edger") {
    # calculate dispersion
    message("  estimating dispersion without design matrix")
    dge       = estimateDisp(dge)

    # fit model to each cluster
    message('  run edgeR on each cluster')
    message('    ', appendLF = FALSE)
    # mkrs_pb_dt  = cl_ls %>% bplapply( function(sel_cl) {
    mkrs_pb_dt  = cl_ls %>% lapply( function(sel_cl) {
      message(sel_cl, " ", appendLF = FALSE)
      # make design matrix for this celltype
      this_d    = copy(des_all) %>% .[, is_cluster := cluster == sel_cl ] %>%
        model.matrix( ~ is_cluster, data = . )

      # fit for this celltype
      fit       = glmQLFit(dge, design = this_d)
      fit       = glmTreat(fit, coef = "is_clusterTRUE")
      tmp_dt    = topTags(fit, n = Inf, sort.by = "none") %>%
        as.data.frame %>% as.data.table(keep.rownames = 'gene_id') %>%
        .[, cluster := sel_cl ]

      return(tmp_dt)
    }) %>% rbindlist
    message()
    # }, BPPARAM = bpparam) %>% rbindlist
    # bpstop(bpparam)
  } else if (method == "voom") {
    # estimate replicate correlations
    do_rnd    = FALSE
    if (do_rnd) {
      message('  calculate intra-sample dependencies')
      corfit    = duplicateCorrelation(voom(dge_all, mm_all),
        mm_all, block = des_all$sample_id)
    }

    # fit model to each cluster
    message('  run voom on each cluster')
    mkrs_pb_dt  = cl_ls %>% bplapply( function(sel_cl) {
      # make design matrix for this celltype
      mm_cl       = copy(des_all) %>% .[, is_cluster := cluster == sel_cl ] %>%
        model.matrix( ~ is_cluster, data = . )

      # # remove tiny genes
      # dge_cl      = dge[ filterByExpr(dge, mm_cl), ]

      if (do_rnd) {
        voom_obj  = voom(dge, mm_cl,
          block = des_all$sample_id, correlation = corfit$consensus.correlation)
        fit       = lmFit(voom_obj, mm_cl,
          block = des_all$sample_id, correlation = corfit$consensus.correlation)
      } else {
        # run voom
        voom_obj  = voom(dge, mm_cl, plot = FALSE, save.plot = TRUE)
        fit       = lmFit(voom_obj, mm_cl)
      }

      # fit for this celltype
      suppressMessages({
        # tmp_dt    = eBayes(fit) %>% topTable(n = Inf, sort.by = "none") %>%
        tmp_dt    = treat(fit) %>% topTable(n = Inf, sort.by = "none") %>%
          as.data.frame %>% as.data.table(keep.rownames = 'gene_id') %>%
          .[, cluster := sel_cl ]
      })

      return(tmp_dt)
    }, BPPARAM = bpparam) %>% rbindlist

    # edit names
    mkrs_pb_dt  = mkrs_pb_dt %>%
      .[, .(gene_id, cluster, logFC, logCPM = AveExpr,
        PValue = P.Value, FDR = adj.P.Val, t)]
  }

  # calculate logcpms
  message('  calculate mean expression in all clusters')
  cpms_ok   = cl_ls %>% bplapply(function(cl) {
    ok_ls     = x_ls[[ cl ]] %>% colnames %>% str_extract("(?<=-).+")
    cpms_ok   = logcpms_all[ (cluster == cl) & (get(batch_var) %in% ok_ls) ]
    return(cpms_ok)
    }, BPPARAM = bpparam) %>% rbindlist
  cpm_means = cpms_ok[, .(mean_logcpm = mean(logcpm), n_batches = .N),
    by = .(cluster, gene_id)]

  # calculate relative cpms
  message('  calculate means in selected and other clusters')
  rel_cpms  = cl_ls %>% bplapply(function(cl) {
    # make matrix
    cpms_sel    = cpm_means[ cluster == cl ] %>%
      .[, .(logcpm = mean_logcpm, cluster, gene_id)]
    cpms_other  = cpm_means[ cluster != cl ] %>%
      .[, .(logcpm = sum(n_batches * mean_logcpm) / sum(n_batches)), by = gene_id ]

    # join together
    rel_sel     = merge(cpms_sel, cpms_other, by = "gene_id",
      suffixes = c(".sel", ".other"))

    return(rel_sel)
  }, BPPARAM = bpparam) %>% rbindlist

  message('  adding no. samples with zero counts')
  zeros_dt    = cl_ls %>% bplapply(function(cl) {
    # make matrix
    this_x    = x_ls[[ cl ]]
    tmp_dt    = data.table(
      cluster   = cl,
      gene_id   = rownames(this_x),
      n_zero    = rowSums(this_x == 0),
      n_cl      = ncol(this_x)
    )

    return(tmp_dt)
  }, BPPARAM = bpparam) %>% rbindlist
  mkrs_pb_dt  = merge(mkrs_pb_dt, zeros_dt, by = c("cluster", "gene_id"))

  message('  adding logcpm data')
  mkrs_pb_dt  = merge(mkrs_pb_dt, rel_cpms, by = c("cluster", "gene_id")) %>%
    .[ order(cluster, PValue) ]

  message('  adding row data from sce')
  mkrs_pb_dt  = merge(mkrs_pb_dt, rows_dt, by = 'gene_id') %>%
    .[ order(cluster, PValue) ]

  # save results
  fwrite(mkrs_pb_dt, file = mkrs_pb_f)

  return(mkrs_pb_dt)
}

get_top_markers <- function(input_mkrs, fdr_cut = 0.01, n_top = 10, max_zero_p = 0.5) {
  # order and filter
  top_mkrs_tmp  = input_mkrs %>% .[ order(cluster, FDR, -abs(logFC)) ]

  # take some top genes
  top_mkrs_dt   = rbind(
    top_mkrs_tmp[ (FDR < fdr_cut) & (logFC > 0) & (n_zero / n_cl <= max_zero_p) ],
    top_mkrs_tmp[ (FDR < fdr_cut) & (logFC < 0) ],
    top_mkrs_tmp[ FDR >= fdr_cut & ((n_zero / n_cl <= max_zero_p) | (logFC < 0)) ]
    ) %>%
    .[, .SD[ 1:n_top ], by = cluster] %>%
    .[, prop_zero := n_zero / n_cl ]

  return(top_mkrs_dt)
}

plot_dotplot <- function(pb_exp, pb_prop, known_dt, sel_types,
  cl_order = NULL, min_prop = 0.01) {
  # check inputs
  if (is.null(cl_order)) {
    cl_order  = assayNames(pb)
  } else {
    assert_that(
      all(sort(assayNames(pb)) == sort(cl_order)) ,
      all(sort(assayNames(pb_prop)) == sort(cl_order))
    )
  }
  # define pseudocount
  p_count   = 1
  # check objects are similar
  assert_that(all(rownames(pb_exp) == rownames(pb_prop)),
  all(colnames(pb_exp) == colnames(pb_prop)))
  gs_idx    = rowData(pb_exp)$symbol %in% known_dt$symbol

  # get expression (at celltype level)
  exp_all   = lapply(assayNames(pb_exp), function(n) {
      assay(pb_exp, n) %>% .[gs_idx, , drop = FALSE] %>%
        as.data.table(keep.rownames = 'gene_id') %>%
        melt.data.table(id = 'gene_id', variable.name = 'sample_id',
          value.name = 'count') %>% .[, cluster := n]
    }) %>% rbindlist
  sizes_dt  = sapply(assays(pb_exp), colSums) %>%
    as.data.table(keep.rownames = 'sample_id') %>%
    melt.data.table(id = 'sample_id', variable.name = 'cluster', value.name = 'lib_size')
  exp_dt  = merge(exp_all, sizes_dt, by = c('cluster', 'sample_id')) %>%
    .[, .(count = sum(count), lib_size = sum(lib_size)),
      by = c('cluster', 'gene_id')] %>%
    .[, CPM   := ifelse(lib_size == 0, 0, count / lib_size)] %>%
    .[, logcpm  := log10(CPM*1e6 + p_count)]

  # get proportions (at celltype level)
  prop_all  = lapply(assayNames(pb_prop), function(n) {
      assay(pb_prop, n) %>% .[gs_idx, , drop = FALSE] %>%
        as.data.table(keep.rownames = 'gene_id') %>%
        melt.data.table(id = 'gene_id', variable.name = 'sample_id',
           value.name = 'prop_detected') %>% .[, cluster := n]
    }) %>% rbindlist
  assert_that(
    !any(prop_all$prop_detected > 1, na.rm = TRUE),
    !any(prop_all$prop_detected < 0, na.rm = TRUE)
    )
  cells_dt  = int_colData(pb_prop)$n_cells %>% do.call(rbind, .) %>%
    as.data.table(keep.rownames = 'sample_id') %>%
    melt.data.table(id = 'sample_id', variable.name = 'cluster', value.name = 'n_cells')
  props_dt  = merge(prop_all, cells_dt, by = c('cluster', 'sample_id')) %>%
    .[, cells_detected := prop_detected * n_cells] %>%
    .[, .(cells_detected = sum(cells_detected), n_cells = sum(n_cells)),
      by = c('cluster', 'gene_id')] %>%
    .[, prop_detected := cells_detected / n_cells]
  assert_that(
    !any(props_dt$prop_detected > 1, na.rm = TRUE),
    !any(props_dt$prop_detected < 0, na.rm = TRUE)
    )
  assert_that(
    all(props_dt$cluster == exp_dt$cluster),
    all(props_dt$gene_id == exp_dt$gene_id)
    )

  # add labels to genes
  rows_dt   = rowData(pb_exp) %>% as.data.frame %>%
    as.data.table(keep.rownames = 'gene_id') %>% .[, .(gene_id, symbol)]

  # add labels to celltypes
  symbol_ord  = known_dt$symbol %>% as.character %>% fct_inorder %>% levels
  dots_dt     = merge(exp_dt, props_dt, by = c('cluster', 'gene_id')) %>%
    merge(rows_dt, by = 'gene_id') %>%
    merge(known_dt, by = 'symbol', allow.cartesian = TRUE) %>%
    .[, symbol  := factor(symbol, levels = symbol_ord) ] %>%
    .[, cluster := factor(cluster, levels = cl_order) ]

  # define how to split out markers
  title_str = sprintf('Expression of known marker genes across\n%s clusters',
    sel_types %>% tolower %>% paste(collapse = " + "))

  # dotplot
  brk_vals  = log10(p_count + c(0, 10^seq(0, 4)))
  brk_labs  = c('0', '1', '10', '100', '1k', '10k')
  plot_dt   = dots_dt[count > 0] %>%
    .[ prop_detected > min_prop ] %>% setorder( CPM )

  g = ggplot(plot_dt) +
    aes(x = cluster, y = fct_rev(symbol), fill = logcpm, size = prop_detected) +
    geom_blank( data = dots_dt, aes(x = cluster, y = fct_rev(symbol)) ) +
    geom_point(shape = 21) +
    scale_fill_viridis(breaks = brk_vals, labels = brk_labs, option = 'magma') +
    scale_size(range = c(0, 4), breaks = pretty_breaks()) +
    expand_limits(size = 0, fill = log10(p_count)) +
    facet_grid( type_broad ~ ., space = 'free', scales = 'free') +
    theme_bw() +
    theme(
      axis.text.x       = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid        = element_blank(),
      strip.text        = element_text(size = 8),
      strip.background  = element_rect( fill = "white" )
    ) +
    labs(x = 'cluster', y = 'marker gene',
      fill = 'CPM', size = 'propn. cells\nexpressing',
      title = title_str)

  return(g)
}

plot_selected_genes <- function(sel_dt, cpms_dt, cl_order = NULL, pseudo_count = 10,
  ncol = 2, nrow = NULL) {
  # check cl_order ok
  if (!is.null(cl_order)) {
    assert_that( all(sort(unique(cpms_dt$cluster)) == sort(cl_order)) )
  }
  # infer number of rows
  if (is.null(nrow)) {
    nrow    = max(5, ceiling( nrow(sel_dt) / 2 ))
  }
  # which genes?
  if ("vst_var" %in% colnames(sel_dt)) {
    sel_mkrs  = copy(sel_dt) %>%
      .[, symbol_lab := sprintf("%s (variance = %.2f)", symbol, vst_var) %>% fct_inorder ] %>%
      .[, .(gene_id, symbol_lab) ]
  } else {
    sel_mkrs  = copy(sel_dt) %>%
      .[, symbol_lab := symbol %>% fct_inorder ] %>%
      .[, .(gene_id, symbol_lab) ]
  }

  # get data
  plot_dt   = cpms_dt %>% merge( sel_mkrs, by = 'gene_id' )
  if (!is.null(cl_order)) {
    plot_dt   = plot_dt[, cluster := factor(cluster, levels = cl_order) ]
  }

  # get nice colours
  cl_ls     = plot_dt$cluster %>% unique
  cl_cols   = nice_cols[ seq_along(cl_ls) ] %>% setNames(cl_ls)

  # plot!
  log_brks  = c(0, 10, 20, 50, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4) %>%
    `+`(pseudo_count) %>% log
  log_labs  = c('0', '10', '20', '50', '100', '200', '500',
    '1k', '2k', '5k', '10k', '20k', '50k')
  g = ggplot(plot_dt) +
    aes( x = cluster, y = logcpm, fill = cluster, size = n_cells ) +
    geom_quasirandom( colour = 'black', shape = 21 ) +
    scale_y_continuous( breaks = log_brks, labels = log_labs ) +
    expand_limits( y = (c(0, 100) + 10) %>% log ) +
    scale_fill_manual( values = cl_cols, guide = "none" ) +
    scale_size( range = c(0, 6) ) + expand_limits( size = 0 ) +
    facet_wrap( ~ symbol_lab, scales = 'free_y', ncol = ncol, nrow = nrow ) +
    theme_classic() +
    theme(
      axis.text.x       = element_text( angle = 90, hjust = 1, vjust = 0.5 )
    ) +
    labs(
      x = "cluster", y = "counts per million", size = "# cells"
    )

  return(g)
}

plot_top_marker_genes <- function(sel_cl, top_mkrs_dt, logcpms_all,
  cl_order = NULL, pseudo_count = 10) {
  # check cl_order ok
  if (!is.null(cl_order)) {
    assert_that( all(sort(unique(cpms_dt$cluster)) == sort(cl_order)) )
  }

  # which genes?
  fdr_show  = 0.01
  sel_mkrs  = top_mkrs_dt[ cluster == sel_cl ] %>%
    .[, symbol_lab := ifelse(FDR < fdr_show,
        sprintf("%s (FDR = %.1e)", symbol, FDR),
        sprintf("%s (FDR = %.3f)", symbol, FDR)
      ) %>% fct_inorder ] %>%
    .[, .(gene_id, symbol_lab) ]

  # get data
  plot_dt   = logcpms_all %>% merge( sel_mkrs, by = 'gene_id' ) %>%
    .[, is_sel := ifelse(cluster == sel_cl, "test cluster", "other") ]
  if (is.null(cl_order)) {
    plot_dt   = plot_dt[, cluster := fct_relevel(cluster, sel_cl) ]
  } else {
    plot_dt   = plot_dt[, cluster := factor(cluster, levels = cl_order) ]
  }

  # plot!
  log_brks  = c(0, 10, 20, 50, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4) %>%
    `+`(pseudo_count) %>% log
  log_labs  = c('0', '10', '20', '50', '100', '200', '500',
    '1k', '2k', '5k', '10k', '20k', '50k')
  g = ggplot(plot_dt) +
    aes( x = cluster, y = logcpm, fill = is_sel, size = n_cells ) +
    geom_quasirandom( colour = 'black', shape = 21 ) +
    scale_y_continuous( breaks = log_brks, labels = log_labs ) +
    expand_limits( y = (c(0, 100) + pseudo_count) %>% log ) +
    scale_fill_manual( values = c(`test cluster` = 'grey10', other = 'grey80') ) +
    scale_size( range = c(0, 6) ) + expand_limits( size = 0 ) +
    facet_wrap( ~ symbol_lab, scales = 'free_y', ncol = 2, nrow = 5 ) +
    theme_classic() +
    theme(
      axis.text.x       = element_text( angle = 90, hjust = 1, vjust = 0.5 )
    ) +
    labs(
      x = "cluster", y = "counts per million", fill = 'cluster status', size = "# cells",
      title = sprintf("Top %d marker genes for %s", nrow(sel_mkrs), sel_cl)
    )

  return(g)
}

plot_clusters_by_metadata <- function(meta_dt, clusts_dt, meta_vars = NULL,
  cl_order = NULL, min_cells = 10) {
  # check inputs
  assert_that( all(meta_vars %in% colnames(meta_dt)) )

  # join metadata and clusters
  melt_dt   = clusts_dt %>%
    .[, .N, by = .(sample_id, cluster) ] %>%
    dcast.data.table( sample_id ~ cluster, value.var = "N", fill = 0) %>% 
    melt.data.table( id = "sample_id", variable.name = "cluster", value.name = "N") %>% 
    merge(meta_dt[, c('sample_id', meta_vars), with = FALSE ], by = "sample_id") %>%
    melt.data.table( measure = meta_vars, var = "meta_var", val = "meta_val" ) %>%
    .[, .(N = sum(N)), by = .(meta_var, meta_val, cluster) ] %>%
    .[, prop := N / sum(N), by = .(meta_var, cluster) ]

  # put clusters in order
  if (!is.null(cl_order)) {
    # in case there are missing clusters:
    cl_all    = unique(melt_dt$cluster) %>% sort
    cl_tmp    = c(cl_order, setdiff(cl_all, cl_order))
    melt_dt = melt_dt[ , cluster := factor(cluster, levels = cl_tmp) ]
  }
  # define palettes
  pal_ls      = c("Greys", "Reds", "Blues", "Greens", "Purples", "Oranges", "YlOrBr")

  # temp function for plotting
  .plot_fn <- function( ii ) {
    this_var  = meta_vars[[ ii ]]
    # plot quasirandom dots facetted by meta_var
    plot_dt   = melt_dt %>% .[ meta_var == this_var ]
    g = ggplot(plot_dt) +
      aes( x = cluster, y = 100 * prop, fill = meta_val ) +
      geom_col( colour = 'black' ) +
      # geom_quasirandom( colour = 'black', shape = 21, size = 3 ) +
      scale_y_continuous( breaks = pretty_breaks() ) +
      scale_fill_brewer( palette = pal_ls[[ ii ]], direction = -1 ) +
      # scale_size( range = c(0, 6), breaks = sqrt_brks, labels = sqrt_labs ) +
      expand_limits( size = 0, y = 0 ) +
      theme_classic() +
      theme(
        axis.text.x       = element_text( angle = 90, hjust = 1, vjust = 0.5 )
      ) +
      labs(
        y = "percentage", size = "# cells", x = NULL, fill = this_var
      )
    }
  g     = lapply(seq_along(meta_vars), .plot_fn) %>% wrap_plots( ncol = 1 )

  return(g)
}

plot_marker_dotplot_pb <- function(cpms_dt, props_dt, markers_dt,
  nice_order = FALSE, min_p = 0.01) {
  # calc expression per cluster
  p_count   = 1
  dots_all  = merge(
    cpms_dt[, .(cluster, sample_id, gene_id, symbol, count, lib_size)],
    props_dt[, .(cluster, sample_id, gene_id, n_detected = prop * n_cells, n_cells)],
    by = c("sample_id", "cluster", "gene_id")) %>%
    .[, .(
      count         = sum(count),
      lib_size      = sum(lib_size),
      n_detected    = sum(n_detected),
      n_cells       = sum(n_cells)
      ), by = c("cluster", "gene_id", "symbol")] %>%
    .[, CPM       := ifelse(lib_size == 0, 0, count / lib_size * 1e6) ] %>%
    .[, log10cpm  := log10( CPM + p_count ) ] %>%
    .[, p_detect  := n_detected / n_cells ]

  # restrict to desired markers
  plot_dt   = merge(dots_all, markers_dt, by = 'symbol') %>%
    .[, symbol    := symbol %>% factor(levels = markers_dt$symbol) ]

  if (nice_order) {
    # put matrix in nice order
    order_dt    = plot_dt %>%
      dcast(symbol ~ cluster, value.var = "log10cpm", fill = 0) %>%
      melt( id = "symbol", var = "cluster" ) %>%
      .[, symbol := factor(symbol, levels = markers_dt$symbol) ] %>%
      .[ order(cluster, symbol) ] %>%
      .[, smth_score := ksmooth(as.numeric(symbol), value,
        kernel = "normal", x.points = as.numeric(symbol))$y, by = cluster ] %>%
      .[, .SD[ which.max(smth_score) ], by = cluster ] %>%
      .[ order(symbol, -smth_score) ]
    assert_that( all( sort(order_dt$cluster) == sort(unique(plot_dt$cluster)) ) )

    # add this order to dots
    plot_dt     = plot_dt[, cluster := factor(cluster, levels = order_dt$cluster) ]
  }

  # dotplot
  brk_vals  = log10(p_count + c(0, 10^seq(0, 4)))
  brk_labs  = c('0', '1', '10', '100', '1k', '10k')
  plot_dt   = plot_dt[ p_detect >= min_p ] %>% setorder( CPM )
  g = ggplot(plot_dt) +
    aes(x = cluster, y = fct_rev(symbol), fill = log10cpm, size = p_detect) +
    geom_point(shape = 21) +
    scale_fill_viridis(breaks = brk_vals, labels = brk_labs, option = 'magma') +
    scale_size(range = c(0, 4), breaks = pretty_breaks()) +
    expand_limits(size = 0, fill = log10(p_count)) +
    facet_grid( marker_cat ~ ., scales = 'free_y', space = 'free_y' ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid  = element_blank(),
      panel.spacing     = unit(0.2, "lines"),
      strip.text.y      = element_text(size = 8),
      strip.background  = element_rect( fill = 'white')
    ) +
    labs(
      x = "cluster", y = "marker gene",
      fill = "CPM", size = "propn. cells\nexpressing"
      )

  return(g)
}

calc_confuse_dt <- function(cl1_dt, cl2_dt, cl1, cl2) {
  assert_that( cl1 %in% names(cl1_dt) )
  assert_that( cl2 %in% names(cl2_dt) )
  assert_that( !is.null(levels(cl1_dt[[ cl1 ]])) )
  assert_that( !is.null(levels(cl2_dt[[ cl2 ]])) )
  assert_that( cl2 %in% names(cl2_dt) )
  assert_that(
    "cell_id" %in% names(cl1_dt),
    "cell_id" %in% names(cl2_dt)
  )

  confuse_dt = merge(
    cl1_dt[, .(cell_id, cl1 = get(cl1)) ],
    cl2_dt[, .(cell_id, cl2 = get(cl2)) ],
    by = "cell_id") %>%
    .[, .N, by = .(cl1, cl2) ]
  match_dt    = expand.grid(
    cl1  = unique(confuse_dt$cl1),
    cl2  = unique(confuse_dt$cl2)
    )
  confuse_dt  = confuse_dt %>%
    merge( match_dt, by = c("cl1", "cl2"), all = TRUE ) %>%
    .[ is.na(N), N := 0 ] %>%
    .[, N0        := N + 0.1 ] %>%
    .[, log_N     := log(N0) ] %>%
    .[, p_cl1     := N / sum(N), by = cl1 ] %>%
    .[, log_p_cl1 := (N0 / sum(N0)) %>% log, by = cl1 ] %>%
    .[, p_cl2     := N / sum(N), by = cl2 ] %>%
    .[, log_p_cl2 := (N0 / sum(N0)) %>% log, by = cl2 ]

  # # sort factor levels
  # lvls_cl1    = confuse_dt$cl1 %>% levels %>% sort
  # confuse_dt  = confuse_dt[, cl1 := factor(cl1, levels = lvls_cl1)]
  # lvls_cl2    = confuse_dt$cl2 %>% levels %>% sort
  # confuse_dt  = confuse_dt[, cl2 := factor(cl2, levels = lvls_cl2)]

  return(confuse_dt)
}

plot_cluster_comparison_heatmap <- function(confuse_dt, cl1_lab, cl2_lab,
  plot_var = c("log_p_cl1", "p_cl1", "log_p_cl2", "p_cl2", "N", "log_N"),
  do_sort = c("no", "hclust", "seriate"), order_cl1 = NULL, order_cl2 = NULL,
  lbl_threshold = 0.05) {
  # check inputs
  plot_var    = match.arg(plot_var)
  if (!is.null(order_cl1))
    assert_that( all(unique(confuse_dt$cl1) %in% order_cl1) )
  if (!is.null(order_cl2))
    assert_that( all(unique(confuse_dt$cl2) %in% order_cl2) )

  # don't make any changes!
  copy_dt     = copy(confuse_dt)
  if (!is.null(order_cl1))
    copy_dt     = copy_dt[, cl1 := factor(cl1, levels = order_cl1)]
  if (!is.null(order_cl2))
    copy_dt     = copy_dt[, cl2 := factor(cl2, levels = order_cl2)]

  # decide what to plot
  if (plot_var == "N") {
    data_wide   = dcast(copy_dt, cl1 ~ cl2, value.var = "N", fill = 0)
    value_name  = "no. cells"

    # define colours
    max_val     = copy_dt$N %>% max
    chunk_opts  = c(5e2, 1e3, 5e3, 1e4)
    best_chunk  = ((max_val / 10) / chunk_opts) %>% `<`(1) %>% which %>% min %>% chunk_opts[ . ]
    n_brks      = seq(0, ceiling(max_val / best_chunk) * best_chunk, by = best_chunk)
    n_labs      = n_brks
    mat_cols  = cols_fn(seq(min(n_brks), max(n_brks), best_chunk),
      res = best_chunk, pal = "viridis", pal_dir = 1, range = "natural")
    lgd       = list(title = "no. cells\nin sample",
      at = n_brks, labels = n_labs)

  } else if (plot_var == "log_N") {
    data_wide   = dcast(copy_dt, cl1 ~ cl2, value.var = "log_N")

    # define log breaks
    log_brks  = c(0, 3, 10, 30, 1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5) %>%
      `+`(1) %>% log
    log_labs  = c("0", "3", "10", "30", "100", "300",
      "1k", "3k", "10k", "30k", "100k", "300k")
    assert_that( length(log_brks) == length(log_labs) )

    # truncate to observed range
    max_val   = copy_dt$log_N %>% max
    assert_that( max_val <= max(log_brks) )
    max_brk   = (log_brks <= max_val) %>% which %>% max
    log_brks  = log_brks[seq.int(max_brk)]
    log_labs  = log_labs[seq.int(max_brk)]

    # get colours
    res       = 0.1
    mat_cols  = cols_fn(seq(min(log_brks), max_val, res),
      res = res, pal = "viridis", pal_dir = 1, range = "natural")
    lgd       = list(title = "no. cells",
      at = log_brks, labels = log_labs)

  } else if (plot_var == "p_cl1") {
    data_wide   = dcast(copy_dt, cl1 ~ cl2, value.var = "p_cl1")

    # define percentage breaks
    pct_brks  = seq(0, 1, 0.2)
    pct_labs  = sprintf("%.0f%%", 100 * pct_brks)
    res       = 0.1
    mat_cols  = cols_fn(seq(min(pct_brks), max(pct_brks), res),
      res = res, pal = "viridis", pal_dir = 1, range = "natural")
    lgd       = list(title = "cluster\nproportion\n(rows sum to 1)",
      at = pct_brks, labels = pct_labs)

  } else if (plot_var == "log_p_cl1") {
    data_wide   = dcast(copy_dt, cl1 ~ cl2, value.var = "log_p_cl1")

    # define colours
    pct_brks  = c(0.001, 0.002, 0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1) %>% log
    pct_labs  = c("0.1%", "0.2%", "0.4%", "1%", "2%", "4%", "10%", "20%", "40%", "100%")
    res       = 0.1
    mat_cols  = cols_fn(seq(min(pct_brks), max(pct_brks), res),
      res = res, pal = "viridis", pal_dir = 1, range = "natural")
    lgd       = list(title = "cluster\nproportion\n(rows sum to 1)",
      at = pct_brks, labels = pct_labs)

  } else if (plot_var == "p_cl2") {
    data_wide   = dcast(copy_dt, cl1 ~ cl2, value.var = "p_cl2")

    # define percentage breaks
    pct_brks  = seq(0, 1, 0.2)
    pct_labs  = sprintf("%.0f%%", 100 * pct_brks)
    res       = 0.1
    mat_cols  = cols_fn(seq(min(pct_brks), max(pct_brks), res),
      res = res, pal = "viridis", pal_dir = 1, range = "natural")
    lgd       = list(title = "original\nproportion\n(cols sum to 1)",
      at = pct_brks, labels = pct_labs)

  } else if (plot_var == "log_p_cl2") {
    data_wide   = dcast( copy_dt, cl1 ~ cl2, value.var = "log_p_cl2" )

    # define colours
    pct_brks  = c(0.001, 0.002, 0.004, 0.01, 0.02, 0.04, 0.1, 0.2, 0.4, 1) %>% log
    pct_labs  = c("0.1%", "0.2%", "0.4%", "1%", "2%", "4%", "10%", "20%", "40%", "100%")
    res       = 0.1
    mat_cols  = cols_fn(seq(min(pct_brks), max(pct_brks), res),
      res = res, pal = "viridis", pal_dir = 1, range = "natural")
    lgd       = list(title = "original\nproportion\n(cols sum to 1)",
      at = pct_brks, labels = pct_labs)

  }

  # add text annotations
  if ( plot_var %in% c("p_cl1", "log_p_cl1", "p_cl2", "log_p_cl2")) {
    sel_cl    = str_extract(plot_var, "cl[0-9]")
    # define annotations
    txt_mat   = dcast( copy_dt, cl1 ~ cl2, value.var = paste0("p_", sel_cl) ) %>%
      as.matrix(rownames = "cl1")
    .lbl_fn <- function(j, i, x, y, width, height, fill) {
      val     = txt_mat[i, j]
      if (val < lbl_threshold) {
        s       = ""
      } else {
        # n_dec   = ifelse( log10(val) > 1, 0, 1)
        # s       = paste0("%.", n_dec, "f%%") %>% sprintf(val)
        s       = sprintf("%.0f%%", 100 * val)
      }
      return(grid.text(s, x, y, gp = gpar(fontsize = 8)))
    }

  } else if ( plot_var %in% c("N", "log_N")) {
    # define annotations
    txt_mat   = dcast( copy_dt, cl1 ~ cl2, value.var = "N" ) %>%
      as.matrix(rownames = "cl1")
    .lbl_fn <- function(j, i, x, y, width, height, fill)
      grid.text(sprintf("%s", signif(txt_mat[i, j], 2)),
        x, y, gp = gpar(fontsize = 8))
  }

  # turn into matrix
  data_mat    = data_wide %>% as.matrix( rownames = "cl1" )

  # make data for annotations
  rows_dt     = copy_dt[, .(total_cl1 = sum(N)), by = cl1 ] %>%
    .[, log_total_cl1 := log(total_cl1) ] %>%
    setkey("cl1") %>% .[ rownames(data_mat) ]
  cols_dt     = copy_dt[, .(total_cl2 = sum(N)), by = cl2 ] %>%
    .[, log_total_cl2 := log(total_cl2) ] %>%
    setkey("cl2") %>% .[ colnames(data_mat) ]
  assert_that(
    nrow(data_mat) == nrow(rows_dt),
    ncol(data_mat) == nrow(cols_dt)
    )

  # do annotations
  # if (plot_var == "log_N") {
    log_brks  = c(0, 1, 10, 1e2, 1e3, 1e4, 1e5) %>%
      `+`(1) %>% log
    log_labs  = c("0", "1", "10", "100", "1k", "10k", "100k")
    res       = 0.1
    log_cols  = cols_fn(seq(min(log_brks), max(log_brks), res),
      res = res, pal = "magma", pal_dir = 1, range = "natural")

    # label with broad types
    row_annots  = rowAnnotation(
      `cl1 total`  = rows_dt$log_total_cl1,
      col = list(`cl1 total` = log_cols),
      annotation_name_side = "top",
      annotation_legend_param = list(
        `cl1 total` = list(at = log_brks, labels = log_labs)
        )
      )
    col_annots  = HeatmapAnnotation(
      `cl2 total`  = cols_dt$log_total_cl2,
      col = list(`cl2 total` = log_cols),
      annotation_name_side = "right",
      annotation_legend_param = list(
        `cl2 total` = list(at = log_brks, labels = log_labs)
        )
      )

  if (do_sort == "no") {
    rows_flag   = FALSE
    cols_flag   = FALSE
    row_order   = NULL
    col_order   = NULL

  } else if (do_sort == "hclust") {
    # define vars
    rows_flag   = TRUE
    cols_flag   = TRUE
    row_order   = NULL
    col_order   = NULL

  } else if (do_sort == "seriate") {
    # do seriate
    data_min    = data_mat %>% as.vector %>% min(na.rm = TRUE)
    data_mat[is.na(data_mat)] = data_mat
    seriate_obj = seriate(data_mat - data_min, method = "BEA")

    # define vars
    rows_flag   = FALSE
    cols_flag   = FALSE
    row_order   = get_order(seriate_obj, 1)
    col_order   = get_order(seriate_obj, 2)
  }

  # heatmap
  hm_obj      = Heatmap(
    matrix = data_mat, col = mat_cols,
    cell_fun = .lbl_fn,
    cluster_rows = rows_flag, cluster_columns = cols_flag,
    row_order = row_order, column_order = col_order,
    # row_names_gp = gpar(fontsize  =  8), column_names_gp = gpar(fontsize  =  8),
    row_title = cl1_lab, column_title = cl2_lab,
    left_annotation = row_annots, top_annotation = col_annots,
    heatmap_legend_param = lgd,
    row_names_side = "left", column_names_side = "top",
    na_col = "grey"
    )

  return(hm_obj)
}

plot_metadata_over_umap <- function(meta_dt, int_dt, meta_var) {
  assert_that( meta_var %in% colnames(meta_dt) )

  plot_dt   = int_dt[, .(sample_id, cell_id, UMAP1, UMAP2)] %>%
    merge(meta_dt[, .(sample_id, meta_var = get(meta_var))], by = "sample_id")

  g = ggplot(plot_dt) +
    aes( x = UMAP1, y = UMAP2 ) +
    geom_bin2d( bins = 50, aes(fill = after_stat( density ) %>% pmin(0.01) %>% pmax(0.0001)) ) +
    scale_fill_distiller( palette = "RdBu", trans = "log10", limits = c(0.0001, 0.01) ) +
    facet_wrap( ~ meta_var ) +
    theme_classic() +
    theme( aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank() ) +
    labs( fill = sprintf("%s\ndensity", meta_var) )

  return(g)
}

plot_clusters_annotated_by_densities = function(int_dt, v, plot_ratio = sqrt(2)) {
  # check input
  assert_that( v %in% names(int_dt) )

  # define rows and cols
  n_vals    = unique(int_dt[[v]]) %>% length
  n_rows    = sqrt(n_vals / plot_ratio) %>% floor
  n_cols    = ceiling(n_vals / n_rows)

  # do plot
  g = ggplot(int_dt) +
    aes( x = UMAP1, y = UMAP2 ) +
    geom_bin2d( bins = 50, 
      aes(fill = after_stat( density ) %>% multiply_by(100) %>% pmin(1) %>% pmax(0.01)) ) +
    scale_fill_distiller( palette = "RdBu", trans = "log10", limits = c(0.01, 1) ) +
    facet_wrap( sprintf("~ %s", v), nrow = n_rows, ncol = n_cols ) +
    theme_classic() +
    theme( axis.text = element_blank(), axis.ticks = element_blank(), 
      panel.grid = element_blank(), aspect.ratio = 1 ) +
    labs( fill = "pct. of\nbin in UMAP" )

  return(g)
}

load_cell_expression <- function(sce_f, sel_dt, int_dt) {
  sce         = sce_f %>% readRDS %>% logNormCounts

  assert_that( all( sel_dt$gene_id %in% rownames(sce) ) )
  sel_idx     = rownames(sce) %in% sel_dt$gene_id
  # assert_that( sum( sel_idx ) == nrow(sel_dt) )

  cell_exp_dt = logcounts(sce)[ sel_idx, ] %>%
    as.data.table( keep.rownames = "gene_id" ) %>%
    melt.data.table( id = "gene_id", var = "cell_id", val = "logcount" ) %>%
    merge( sel_dt, by = "gene_id", allow.cartesian = TRUE ) %>%
    merge( int_dt[, .(cell_id, UMAP1, UMAP2)], by = "cell_id" ) %>%
    .[, max_val   := quantile(logcount, 0.99), by = gene_id] %>%
    .[ is.na(max_val) | max_val == 0, max_val := 1, by = gene_id ] %>%
    .[, norm_val  := logcount %>% `/`(max_val) %>% pmin(1)] %>%
    .[, symbol    := symbol %>% factor( levels = unique(sel_dt$symbol) ) ]

  return(cell_exp_dt)
}

plot_selected_genes_umap <- function(sel_dt, cols_to_rows = 1.25) {
  # guess at nice numbers for rows and cols
  n_genes = length(unique(sel_dt$symbol))
  n_row   = sqrt(n_genes / cols_to_rows) %>% ceiling
  n_col   = ceiling(n_genes / n_row)

  # make plot
  g = ggplot(sel_dt) +
    aes( x = UMAP1, y = UMAP2, colour = norm_val ) +
    geom_point( size = 0.1 ) +
    scale_colour_viridis( breaks = pretty_breaks(),
      guide = guide_legend( override.aes = list(size = 3) )) +
    facet_wrap( ~ symbol, ncol = n_col, nrow = n_row ) +
    theme_bw() +
    theme( axis.text = element_blank(), axis.ticks = element_blank(), 
      panel.grid = element_blank(), strip.background = element_rect( fill = "white" ), 
      aspect.ratio = 1 ) +
    labs( colour = "scaled log\nexpression\n(max val. = 1)" )
}

process_custom_markers <- function(custom_f, biotypes_dt, hvgs_dt, marker_calcs_dt) {
  cust_mkrs_dt  = fread(custom_f)
  label_lvls    = cust_mkrs_dt$label %>% unique 
  
  # Ensure symbol is present; if not, try to retrieve it using ensembl_id
  if (!"symbol" %in% colnames(cust_mkrs_dt)) {
    # Check if all ensembl_id values are valid
    assert_that(
      all(cust_mkrs_dt$ensembl_id %in% biotypes_dt$ensembl_id),
      msg = sprintf("Warning: Some 'ensembl_id' values in '%s' do not match any in the GTF file.", custom_f)
    )
    cust_mkrs_dt = merge(cust_mkrs_dt, biotypes_dt[, .(symbol, ensembl_id)], by = "ensembl_id", all.x = TRUE, all.y = FALSE)
  } else {
    # Check if all symbol values are valid
    assert_that(
      all(cust_mkrs_dt$symbol %in% biotypes_dt$symbol), 
      msg = sprintf("Warning: Some 'symbol' values in '%s' do not match any in the GTF file.", custom_f)
    )
    assert_that(
      all(cust_mkrs_dt[, .N, by = symbol]$N == 1), 
      msg = sprintf("Some genes in the following custom markers file are repeated across groups:\n%s", custom_f)
    )
  }
  
  # Check that at least two symbols are in marker_calcs_dt
  ok_symbols = intersect(marker_calcs_dt$symbol, cust_mkrs_dt$symbol)
  
  if (length(ok_symbols) < 2) {
    message(sprintf(
      "Skipping '%s': Less than 2 symbols found in the marker genes file. Heatmap for these markers will not be shown.", 
      custom_f
    ))
    return(NULL)
  }
   
  # add gene_id
  if (!('gene_id' %in% names(cust_mkrs_dt))) {
    gene_lu = hvgs_dt[ symbol %in% cust_mkrs_dt$symbol ] %>% 
      .[, .SD[ which.max(vst_var) ], by = symbol ] %>% 
      .[, .(symbol, gene_id) ]
    cust_mkrs_dt = cust_mkrs_dt %>% merge(gene_lu, by = "symbol")
  }
  # put in nice order
  cust_mkrs_dt = cust_mkrs_dt[, label := factor(label, levels = label_lvls) ] %>% 
    .[order(label)]
  
  return(cust_mkrs_dt)
}

plot_heatmap_of_selected_genes <- function(mkrs_dt, panel_dt, max_fc = 3, min_cpm = 10, 
  pseudocount = 10, annotate_genes = TRUE, cluster_rows = TRUE) {
  # unpack
  sel_gs      = panel_dt$gene_id
  
  # get marker values
  missing_gs  = setdiff( sel_gs, unique(mkrs_dt$gene_id))
  if ( length(missing_gs) > 0 )
    message('the following genes are not in the marker genes file:\n  ', 
      paste(missing_gs, collapse = ','))
  mkrs_sel    = mkrs_dt[ (gene_id %in% sel_gs) ]
  max_cpms    = mkrs_sel[, .(max_logcpm = max(logcpm.sel)), by = symbol ]
  keep_gs     = max_cpms[ max_logcpm >= log(min_cpm + pseudocount) ]$symbol
  mkrs_sel    = mkrs_sel[ symbol %in% keep_gs ]

  # make nice matrices
  log2fc_mat  = mkrs_sel %>% 
    .[, .(symbol, cluster, log2fc = logFC) ] %>% 
    dcast( symbol ~ cluster, value.var = 'log2fc' ) %>% 
    as.matrix(rownames = 'symbol')

  # annotate with adjusted p-values
  fdr_mat     = mkrs_sel %>% 
    .[, .(symbol, cluster, 
      signif = ifelse( logcpm.sel < log(min_cpm + 1), '',
        ifelse( (FDR > 0.05) | (logFC < 0), '', 
          ifelse(FDR > 0.01, '*', ifelse(FDR > 0.001, '**', '***'))))) ] %>%
    dcast( symbol ~ cluster, value.var = 'signif', fill = '' ) %>% 
    as.matrix(rownames = 'symbol') %>% t

  # make colours
  fc_cols     = cols_fn(seq(-max_fc, max_fc, 1), res = 0.1, 
    pal = "RdBu", pal_dir = -1, range = "natural")

  # define annotations
  .lbl_fn <- function(j, i, x, y, width, height, fill)
    grid.text(sprintf("%s", fdr_mat[i, j]), 
      x, y, gp = gpar(fontsize = 6))
  lgd         = list(title = "log2fc in\ncluster", at = c(-max_fc, 0, max_fc), 
    direction = "vertical")

  # maybe do annotations
  n_genesets  = length(unique(panel_dt$label))
  if (annotate_genes & (n_genesets > 1)) {
    # label with buckets
    col_lvls    = panel_dt[ symbol %in% keep_gs ]$label %>% fct_inorder %>% levels
    cols_dt     = panel_dt[, .(symbol,  label = label %>% factor(levels = col_lvls))] %>%
      setkey('symbol') %>% .[ rownames(log2fc_mat) ]
    col_cols    = MetBrewer::met.brewer( name = 'Signac', n = length(col_lvls)) %>% setNames(col_lvls)
    col_annots  = HeatmapAnnotation(
      label    = cols_dt$label,
      col      = list(label = col_cols),
      annotation_name_side = "right", 
      annotation_legend_param = list(
        label = list(ncol = 3, labels_gp = gpar(fontsize = 8))
      ) )
    col_split   = cols_dt$label
  } else {
    col_annots  = NULL
    col_split   = NULL
  }

  # heatmap
  hm_obj      = Heatmap(
    matrix = t(log2fc_mat), col = fc_cols, 
    cell_fun = .lbl_fn,
    cluster_rows = cluster_rows, cluster_columns = TRUE,
    column_split = col_split, cluster_column_slices = FALSE, column_title = NULL,
    column_names_gp = gpar( fontsize = 8 ),
    top_annotation = col_annots,
    heatmap_legend_param = lgd,
    row_names_side = "left", column_names_side = "top",
    na_col = "grey"
    )

  return(hm_obj)
}

muscat_n_cells <- function(x) {
  y = int_colData(x)$n_cells
  if (is.null(y)) 
    return(NULL)
  if (length(metadata(x)$agg_pars$by) == 2) 
    y = as.matrix(data.frame(y, check.names = FALSE))
  return(as.table(y))
}
