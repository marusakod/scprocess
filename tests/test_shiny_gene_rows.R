# Actual count-preservation regression, including disk-backed BPCells matrices.
suppressPackageStartupMessages({library(Matrix);library(BPCells)})
expressions <- parse("scripts/shiny.R")
definition <- Filter(function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
  identical(x[[2L]], as.name(".align_shiny_gene_rows")), as.list(expressions))
stopifnot(length(definition)==1L)
eval(definition[[1L]])
a <- matrix(c(1,2,3,4,5,6),3,dimnames=list(c("g1","g2","g3"),c("a1","a2")))
b <- matrix(c(10,20,30,40,50,60),3,dimnames=list(c("g1","g2","g3"),c("b1","b2")))
convert <- function(x) as(Matrix(x,sparse=TRUE),"IterableMatrix")
matrices <- list(a=convert(a),b=convert(b[c(3,1,2),,drop=FALSE]))
warning <- character()
aligned <- withCallingHandlers(.align_shiny_gene_rows(matrices),warning=function(w){
  warning <<- c(warning,conditionMessage(w));invokeRestart("muffleWarning")
})
actual <- as.matrix(do.call(cbind,aligned))
stopifnot(identical(rownames(actual),rownames(a)),
  identical(colnames(actual),c(colnames(a),colnames(b))),
  isTRUE(all.equal(unname(actual),unname(cbind(a,b)))),
  length(warning)==1L,grepl("Same genes, different order",warning))
unchanged <- list(a=convert(a),b=convert(b))
stopifnot(identical(.align_shiny_gene_rows(unchanged),unchanged))
for (ids in list(c("g1","g2","other"),c("g1","g1","g3"),c("g1","","g3"),c("g1",NA,"g3"))) {
  bad <- b;rownames(bad)<-ids
  stopifnot(inherits(tryCatch(.align_shiny_gene_rows(list(a,bad)),error=identity),"error"))
}
stopifnot(inherits(tryCatch(.align_shiny_gene_rows(list()),error=identity),"error"))
message("SUCCESS: Shiny gene-row alignment preserves counts and rejects invalid/different genes")
