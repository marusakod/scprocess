expr <- parse('scripts/shiny.R')
definition <- Filter(function(x)is.call(x)&&identical(x[[1L]],as.name('<-'))&&
 identical(x[[2L]],as.name('.get_density')),as.list(expr))
stopifnot(length(definition)==1L)
eval(definition[[1L]])
set.seed(42)
x <- rnorm(10003);y <- .5*x+rnorm(length(x))
for (grid in list(25,c(31,27))) {
 ref <- MASS::kde2d(x,y,n=grid)
 expected <- ref$z[cbind(findInterval(x,ref$x),findInterval(y,ref$y))]
 for(chunk in c(113,8192,20000)) {
  actual <- .get_density(x,y,n=grid,chunk_size=chunk)
  stopifnot(isTRUE(all.equal(actual,expected,tolerance=1e-12)))
 }
}
h <- c(.8,1.2);lims <- c(-5,5,-6,6)
ref <- MASS::kde2d(x,y,n=41,h=h,lims=lims)
stopifnot(isTRUE(all.equal(.get_density(x,y,n=41,h=h,lims=lims,chunk_size=257),
 ref$z[cbind(findInterval(x,ref$x),findInterval(y,ref$y))],tolerance=1e-12)))
message('SUCCESS: chunked Shiny density matches MASS::kde2d across chunk sizes, grids, bandwidths and limits')
