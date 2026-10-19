# Run from the repository root. Optional args: existing quant directory, output txt.
source('scripts/mapping.R')

# Exercise the fitted-knee branch deterministically without fitting a toy curve.
fit_env <- new.env(parent=globalenv())
fit_env$barcodeRanks <- function(...) list()
fit_env$metadata <- function(x) list(inflection=258, knee=fit_env$fitted_knee)
get_second <- .get_knee_and_shin_2
environment(get_second) <- fit_env

ranks <- data.table(barcode=paste0('bc',seq_len(7)),
                    rank=as.numeric(seq_len(7)),total=c(20000,1431,1430,1428,1000,258,5))
mat <- matrix(1,nrow=1,ncol=nrow(ranks),dimnames=list('gene',ranks$barcode))
run_second <- function(knee, table=ranks, matrix=mat) {
 fit_env$fitted_knee <- knee
 get_second(matrix,table,rank_empty_plateau=0,min_umis_empty=5,
            shin1_x=1,knee2=NA,shin2=NA)
}

# Exact matches preserve the historical rank.
exact <- run_second(1430)
stopifnot(identical(exact$knee2_x,ranks[total==1430,unique(rank)]))
# Regression: 1429 is not observed. Equidistant totals select the first ranked row.
gap <- run_second(1429)
stopifnot(!any(ranks$total==1429),length(gap$knee2_x)==1L,gap$knee2_x==3L,
          unname(gap$sel_knee['knee'])==1429L)
# Downstream interval arithmetic now succeeds without changing the fitted knee.
params <- .get_params_ls(list(ranks_dt=ranks,shin1_x=1),gap,'shin2')
stopifnot(is.finite(params$total_included),params$expected_cells==1,params$lc==258)
# Duplicate observed totals share a rank and still yield one bound.
tied <- copy(ranks);tied[2:3,`:=`(total=1430,rank=2.5)]
stopifnot(identical(run_second(1430,tied)$knee2_x,2.5))
# Out-of-range finite fitted values use the nearest endpoint.
stopifnot(run_second(30000)$knee2_x==1L,run_second(1)$knee2_x==7L)
# A missing fit is a separate failure, not silently replaced with a bound.
error <- tryCatch({run_second(NA_real_);NULL},error=conditionMessage)
stopifnot(is.character(error),grepl('fitted knee is not finite',error,fixed=TRUE))
cat('PASS: exact, missing, tied, endpoint and invalid fitted-knee tests\n')

args <- commandArgs(trailingOnly=TRUE)
if (length(args)) {
 stopifnot(length(args)==2L)
 quant <- args[1];out <- args[2]
 inputs <- c(file.path(quant,'quant.json'),list.files(file.path(quant,'alevin'),full.names=TRUE))
 before <- file.info(inputs)[,c('size','mtime')];hashes <- tools::md5sum(inputs)
 sce <- loadFry(quant,outputFormat=list(S=c('S'),U=c('U'),A=c('A')))
 counts <- do.call(rbind,lapply(assayNames(sce),function(n) assay(sce,n)))
 counts <- counts[,colSums(counts)>0,drop=FALSE]
 rm(sce);gc(FALSE)
 # Use the normal automatic path; do not apply custom sample overrides here.
 result <- calc_ambient_params(counts,run='10X357-2')
 stopifnot(nrow(result)==728182L,
   all(is.finite(result$total_droplets_included)),
   all(is.finite(result$expected_cells)))
 columns <- c('knee1','shin1','knee2','shin2','expected_cells','total_droplets_included','low_count_threshold')
 capture.output({cat('PASS: complete automatic ambient-parameter calculation for existing 10X357-2 counts\n');
   print(unique(result[,..columns]));cat('Input contents, sizes and modification times unchanged.\n')},file=out)
 stopifnot(identical(before,file.info(inputs)[,c('size','mtime')]),identical(hashes,tools::md5sum(inputs)))
 cat('PASS: real-sample parameter calculation; inputs unchanged\n')
}
