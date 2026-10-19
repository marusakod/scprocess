# Run from the scprocess repository root.
source('scripts/utils.R')
source('scripts/mapping.R')
ranks <- data.table(sample_id='test',total=c(100,80,60,60,20,5),rank=1:6,
 knee1=90,knee2=60,shin1=70,shin2=5)
knee <- .get_knee_coordinates_from_ranks(ranks,'sample_id')
shin <- .get_knee_coordinates_from_ranks(ranks,'sample_id','shin')
stopifnot(knee$knee1_rank==1,knee$knee2_rank==3,shin$shin1_rank==2,shin$shin2_rank==6,
 abs(shin$shin1_log10_umis-log10(70))<1e-12)
stopifnot(nrow(.landmark_mcd_ellipse(1:4,1:4))==0,
 nrow(.landmark_mcd_ellipse(rep(1,10),1:10))==0)
set.seed(42)
x <- rnorm(38);y <- .7*x+rnorm(38,sd=.5)
contour <- .landmark_mcd_ellipse(x,y)
fit <- robustbase::covMcd(cbind(x,y),alpha=.5,nsamp='deterministic')
stopifnot(nrow(contour)==181,max(abs(mahalanobis(as.matrix(contour),fit$center,fit$cov)-qchisq(.95,2)))<1e-6)
custom <- tempfile(fileext='.yaml')
yaml::write_yaml(list(sample_id=list(a=list(mapping=list(knee1=100)),b=list(qc=list(min_counts=100)),
 c=list(mapping=list(knee1=NULL))),pool_id=list(pool_a=list(mapping=list(shin2=5)))),custom)
stopifnot(identical(get_mapping_custom_runs(custom),'a'),
 identical(get_mapping_custom_runs(custom,'pool_id'),'pool_a'),length(get_mapping_custom_runs(''))==0)
unlink(custom)
params <- merge(knee,shin,by='sample_id')
for(landmark in c('knee','shin')) {
 g <- plot_knee_coordinates_pairwise(params,landmark=landmark,custom_runs='test')
 b <- ggplot_build(g)
 stopifnot(nrow(g$data)==6,all(g$data$parameter_source=='Custom'),nrow(b$data[[1]])==6,
 all(is.finite(b$data[[1]]$x)),all(is.finite(b$data[[1]]$y)),g$theme$aspect.ratio==1)
}
cat('PASS: knee/shin ranks, covariance contours, insufficient-data fallback, custom-run parsing and plots\n')
