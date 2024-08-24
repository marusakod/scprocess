# metacells.R
suppressPackageStartupMessages({
  # all SuperCell dependancies 
  library("igraph")
  library("RANN")
  library("WeightedCluster")
  library("corpcor")
  library("weights")
  library("Hmisc") 
  library("Matrix")
  library("patchwork")
  library("plyr")
  library("irlba")
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


# functions from SuperCell package (https://github.com/GfellerLab/SuperCell)

supercell_GE <- function(ge, groups, mode = c("average", "sum"), weights = NULL, do.median.norm = FALSE){
  if(ncol(ge) != length(groups)){
    stop("Length of the vector groups has to be equal to the number of cols in matrix ge")
  }

  mode <- mode[1]
  if(!(mode %in% c("average", "sum"))){
    stop(paste("mode", mode, "is unknown. Available values are 'average' and 'sum'."))
  }

  N.SC <- max(groups)
  supercell_size <- as.vector(table(groups))
  j <- rep(1:N.SC, supercell_size) # column indices of matrix M.AV that, whene GE.SC <- ge %M.AV%

  goups.idx  <- plyr::split_indices(groups)
  i <- unlist(goups.idx) # row indices of matrix M.AV that, whene GE.SC <- ge %M.AV%

  if(is.null(weights)){
    M.AV <- Matrix::sparseMatrix(i = i, j = j)
    GE.SC <- ge %*% M.AV

    if(mode == "average"){
      GE.SC <- sweep(GE.SC, 2, supercell_size, "/")
    }
  } else {

    if(length(weights) != length(groups)){
      stop("weights must be the same length as groups or NULL in case of unweighted averaging")
    }

    if(mode != "average"){
      stop(paste("weighted averaging is supposted only for mode = 'average', not for", mode))
    }

    M.AV <- Matrix::sparseMatrix(i = i, j = j, x = weights[i])
    GE.SC <- ge %*% M.AV

    weighted_supercell_size <- unlist(lapply(goups.idx, FUN = function(x){sum(weights[x])}))
    GE.SC <- sweep(GE.SC, 2, weighted_supercell_size, "/")
  }

  if(do.median.norm){
      GE.SC <- (GE.SC+0.01)/apply(GE.SC+0.01, 1, stats::median)
  }

  return(GE.SC)
}


SCimplify <- function(X,
                      genes.use = NULL,
                      genes.exclude = NULL,
                      cell.annotation = NULL,
                      cell.split.condition = NULL,
                      n.var.genes = min(1000, nrow(X)),
                      gamma = 10,
                      k.knn = 5,
                      do.scale = TRUE,
                      n.pc = 10,
                      fast.pca = TRUE,
                      do.approx = FALSE,
                      approx.N = 20000,
                      block.size = 10000,
                      seed = 12345,
                      igraph.clustering = c("walktrap", "louvain"),
                      return.singlecell.NW = TRUE,
                      return.hierarchical.structure = TRUE,
                      ...){

  N.c <- ncol(X)

  if(gamma > 100 & N.c < 100000){
    warning(paste0("Graining level (gamma = ", gamma, ") seems to be very large! Please, consider using smaller gamma, the suggested range is 10-50."))
  }

  if(is.null(rownames(X))){
    if(!(is.null(genes.use) | is.null(genes.exclude))){
      stop("rownames(X) is Null \nGene expression matrix X is expected to have genes as rownames")
    } else {
      warning("colnames(X) is Null, \nGene expression matrix X is expected to have genes as rownames! \ngenes will be created automatically in a form 'gene_i' ")
      rownames(X) <- paste("gene", 1:nrow(X), sep = "_")
    }
  }

  if(is.null(colnames(X))){
    warning("colnames(X) is Null, \nGene expression matrix X is expected to have cellIDs as colnames! \nCellIDs will be created automatically in a form 'cell_i' ")
    colnames(X) <- paste("cell", 1:N.c, sep = "_")
  }

  cell.ids <- colnames(X)

  keep.genes    <- setdiff(rownames(X), genes.exclude)
  X             <- X[keep.genes,]


  if(is.null(genes.use)){
    n.var.genes <- min(n.var.genes, nrow(X))
    if(N.c > 50000){
      set.seed(seed)
      idx         <- sample(N.c, 50000)
      gene.var    <- apply(X[,idx], 1, stats::var)
    } else {
      gene.var    <- apply(X, 1, stats::var)
    }

    genes.use   <- names(sort(gene.var, decreasing = TRUE))[1:n.var.genes]
  }

  if(length(intersect(genes.use, genes.exclude)) > 0){
    stop("Sets of genes.use and genes.exclude have non-empty intersection")
  }

  genes.use <- genes.use[genes.use %in% rownames(X)]
  X <- X[genes.use,]

  if(do.approx & approx.N >= N.c){
    do.approx <- FALSE
    warning("approx.N is larger or equal to the number of single cells, thus, an exact simplification will be performed")
  }

  if(do.approx & (approx.N < round(N.c/gamma))){
    approx.N <- round(N.c/gamma)
    warning(paste("approx.N is set to N.SC", approx.N))
  }

  if(do.approx & ((N.c/gamma) > (approx.N/3))){
    warning("approx.N is not much larger than desired number of super-cells, so an approximate simplification may take londer than an exact one!")
  }

  if(do.approx){
    set.seed(seed)
    approx.N            <- min(approx.N, N.c)
    presample           <- sample(1:N.c, size = approx.N, replace = FALSE)
    presampled.cell.ids <- cell.ids[sort(presample)]
    rest.cell.ids       <- setdiff(cell.ids, presampled.cell.ids)
  } else {
    presampled.cell.ids <- cell.ids
    rest.cell.ids       <- c()
  }

  X.for.pca            <- Matrix::t(X[genes.use, presampled.cell.ids])
  if(do.scale){ X.for.pca            <- scale(X.for.pca) }
  X.for.pca[is.na(X.for.pca)] <- 0

  if(is.null(n.pc[1]) | min(n.pc) < 1){stop("Please, provide a range or a number of components to use: n.pc")}
  if(length(n.pc)==1) n.pc <- 1:n.pc

  if(fast.pca & (N.c < 1000)){
    warning("Normal pca is computed because number of cell is low for irlba::irlba()")
    fast.pca <- FALSE
  }

  if(!fast.pca){
      PCA.presampled          <- stats::prcomp(X.for.pca, rank. = max(n.pc), scale. = F, center = F)
  } else {
    set.seed(seed)
    PCA.presampled          <- irlba::irlba(X.for.pca, nv = max(n.pc, 25))
    PCA.presampled$x        <- PCA.presampled$u %*% diag(PCA.presampled$d)
    PCA.presampled$rotation <- PCA.presampled$v
  }


  sc.nw <- build_knn_graph(
    X = PCA.presampled$x[,n.pc],
    k = k.knn, from = "coordinates",
    #use.nn2 = use.nn2,
    dist_method = "euclidean",
    #directed = directed,
    #DoSNN = DoSNN,
    #pruning = pruning,
    #which.snn = which.snn,
    #kmin = kmin,
    ...
  )


  #simplify

  k   <- round(N.c/gamma)

  if(igraph.clustering[1] == "walktrap"){
    g.s              <- igraph::cluster_walktrap(sc.nw$graph.knn)
    g.s$membership   <- igraph::cut_at(g.s, k)

  } else if(igraph.clustering[1] == "louvain") {
    warning(paste("igraph.clustering =", igraph.clustering, ", gamma is ignored"))
    g.s    <- igraph::cluster_louvain(sc.nw$graph.knn)

  } else {
    stop(paste("Unknown clustering method (", igraph.clustering, "), please use louvain or walkrtap"))
  }

  membership.presampled        <- g.s$membership
  names(membership.presampled) <- presampled.cell.ids

  ## Split super-cells containing cells from different annotations or conditions
  if(!is.null(cell.annotation) | !is.null(cell.split.condition)){
    if(is.null(cell.annotation)) cell.annotation <- rep("a", N.c)
    if(is.null(cell.split.condition)) cell.split.condition <- rep("s", N.c)
    names(cell.annotation) <- names(cell.split.condition) <- cell.ids

    split.cells <- interaction(cell.annotation[presampled.cell.ids], cell.split.condition[presampled.cell.ids], drop = TRUE)

    membership.presampled.intr <- interaction(membership.presampled, split.cells, drop = TRUE)
    membership.presampled <- as.numeric(membership.presampled.intr)
    names(membership.presampled) <- presampled.cell.ids
  }



  SC.NW                        <- igraph::contract(sc.nw$graph.knn, membership.presampled)
  if(!do.approx){
    SC.NW                        <- igraph::simplify(SC.NW, remove.loops = T, edge.attr.comb="sum")
  }


  if(do.approx){

    PCA.averaged.SC      <- as.matrix(Matrix::t(supercell_GE(t(PCA.presampled$x[,n.pc]), groups = membership.presampled)))
    X.for.roration       <- Matrix::t(X[genes.use, rest.cell.ids])



    if(do.scale){ X.for.roration <- scale(X.for.roration) }
    X.for.roration[is.na(X.for.roration)] <- 0


    membership.omitted   <- c()
    if(is.null(block.size) | is.na(block.size)) block.size <- 10000

    N.blocks <- length(rest.cell.ids)%/%block.size
    if(length(rest.cell.ids)%%block.size > 0) N.blocks <- N.blocks+1


    if(N.blocks>0){
      for(i in 1:N.blocks){ # compute knn by blocks
        idx.begin <- (i-1)*block.size + 1
        idx.end   <- min(i*block.size,  length(rest.cell.ids))

        cur.rest.cell.ids    <- rest.cell.ids[idx.begin:idx.end]

        PCA.ommited          <- X.for.roration[cur.rest.cell.ids,] %*% PCA.presampled$rotation[, n.pc] ###

        D.omitted.subsampled <- proxy::dist(PCA.ommited, PCA.averaged.SC) ###

        membership.omitted.cur        <- apply(D.omitted.subsampled, 1, which.min) ###
        names(membership.omitted.cur) <- cur.rest.cell.ids ###

        membership.omitted   <- c(membership.omitted, membership.omitted.cur)
      }
    }

    membership.all_       <- c(membership.presampled, membership.omitted)
    membership.all        <- membership.all_


    names_membership.all <- names(membership.all_)
    ## again split super-cells containing cells from different annotation or split conditions
    if(!is.null(cell.annotation) | !is.null(cell.split.condition)){

      split.cells <- interaction(cell.annotation[names_membership.all],
                                 cell.split.condition[names_membership.all], drop = TRUE)


      membership.all.intr <- interaction(membership.all_, split.cells, drop = TRUE)

      membership.all      <- as.numeric(membership.all.intr)

    }


    SC.NW                        <- igraph::simplify(SC.NW, remove.loops = T, edge.attr.comb="sum")
    names(membership.all) <- names_membership.all
    membership.all <- membership.all[cell.ids]

  } else {
    membership.all       <- membership.presampled[cell.ids]
  }
  membership       <- membership.all

  supercell_size   <- as.vector(table(membership))

  igraph::E(SC.NW)$width         <- sqrt(igraph::E(SC.NW)$weight/10)

  if(igraph::vcount(SC.NW) == length(supercell_size)){
    igraph::V(SC.NW)$size          <- supercell_size
    igraph::V(SC.NW)$sizesqrt      <- sqrt(igraph::V(SC.NW)$size)
  } else {
    igraph::V(SC.NW)$size          <- as.vector(table(membership.all_))
    igraph::V(SC.NW)$sizesqrt      <- sqrt(igraph::V(SC.NW)$size)
    warning("Supercell graph was not splitted")
  }



  res <- list(graph.supercells = SC.NW,
              gamma = gamma,
              N.SC = length(unique(membership)),
              membership = membership,
              supercell_size = supercell_size,
              genes.use = genes.use,
              simplification.algo = igraph.clustering[1],
              do.approx = do.approx,
              n.pc = n.pc,
              k.knn = k.knn,
              sc.cell.annotation. = cell.annotation,
              sc.cell.split.condition. = cell.split.condition
  )

  if(return.singlecell.NW){res$graph.singlecell <- sc.nw$graph.knn}
  if(!is.null(cell.annotation) | !is.null(cell.split.condition)){
    res$SC.cell.annotation. <- supercell_assign(cell.annotation, res$membership)
    res$SC.cell.split.condition. <- supercell_assign(cell.split.condition, res$membership)
  }

  if(igraph.clustering[1] == "walktrap" & return.hierarchical.structure)  res$h_membership <- g.s

  return(res)
}



