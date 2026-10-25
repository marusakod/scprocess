suppressPackageStartupMessages({library(BPCells);library(Matrix)})
expr <- parse('scripts/shiny.R')
for(name in c('.open_shiny_h5ad','.align_shiny_gene_rows')) {
 hit <- Filter(function(x)is.call(x)&&identical(x[[1L]],as.name('<-'))&&
 identical(x[[2L]],as.name(name)),as.list(expr))
 stopifnot(length(hit)==1L);eval(hit[[1L]])
}
run <- function() {
 tmp <- tempfile('shiny-readonly-');dir.create(tmp)
 on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
 raw <- matrix(c(1,0,3,4,0,6,7,8,0,2,1,9),nrow=3,
 dimnames=list(c('g1','g2','g3'),paste0('c',1:4)))
 paths <- file.path(tmp,c('a.h5ad','b.h5ad'))
 for(i in 1:2) {
  # A minimal AnnData CSR fixture (cells by genes), independent of writer API.
  mx <- as(Matrix(raw[,(2*i-1):(2*i)],sparse=TRUE),'dgCMatrix')
  f <- hdf5r::H5File$new(paths[i],mode='w')
  group <- f$create_group('X')
  group[['data']] <- mx@x;group[['indices']] <- mx@i;group[['indptr']] <- mx@p
  hdf5r::h5attr(group,'encoding-type') <- 'csr_matrix'
  hdf5r::h5attr(group,'shape') <- c(ncol(mx),nrow(mx))
  for(name in c('obs','var')) {
   g <- f$create_group(name)
   g[['_index']] <- if(name=='obs')colnames(mx) else rownames(mx)
   hdf5r::h5attr(g,'_index') <- '_index'
  }
  f$close_all()
  Sys.setFileTime(paths[i],as.POSIXct('2020-01-01',tz='UTC'))
 }
 before_time <- file.info(paths)$mtime;before_hash <- tools::md5sum(paths)
 opened <- lapply(paths,.open_shiny_h5ad)
 on.exit(for(x in opened)x$handle$close_all(),add=TRUE)
 mx <- do.call(cbind,.align_shiny_gene_rows(lapply(opened,`[[`,'matrix')))
 stopifnot(isTRUE(all.equal(unname(as.matrix(mx)),unname(raw))))
 # Force lazy row/column transforms and pseudobulk reads while handles live.
 stopifnot(isTRUE(all.equal(as.numeric(colSums(mx)),as.numeric(colSums(raw)))))
 pb <- cbind(rowSums(mx[,1:2]),rowSums(mx[,3:4]))
 stopifnot(isTRUE(all.equal(unname(as.matrix(pb)),unname(cbind(
  rowSums(raw[,1:2]),rowSums(raw[,3:4]))))))
 for(x in opened)x$handle$close_all()
 stopifnot(identical(before_time,file.info(paths)$mtime),
 identical(before_hash,tools::md5sum(paths)),all(!vapply(opened,function(x)x$handle$is_valid,logical(1))))
 opened <- list()
}
run()
message('SUCCESS: Shiny lazy count/pseudobulk reads preserve source bytes and mtimes and close handles')
