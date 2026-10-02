source('scripts/cell_cycle.R')
stopifnot(identical(
  format_cell_cycle_target_prop(c(0.005, 0.01, 0.9, 1, 0.000005)),
  c('0.5%', '1%', '90%', '100%', '0.0005%')))
setDTthreads(1)
set.seed(41)
weights <- c(0.001, 0.01, 0.1, 0.2, 0.2)
original_weights <- weights
original_seed <- .Random.seed
display <- cell_cycle_display_sample(weights, seed=17)
stopifnot(identical(weights, original_weights), identical(.Random.seed, original_seed),
  identical(display, cell_cycle_display_sample(weights, seed=17)),
  all(display$retained[weights == max(weights)]),
  isTRUE(all.equal(display$probabilities, weights / max(weights))),
  isTRUE(all.equal(display$expected_cells, sum(weights) / max(weights))))
for (weights in list(rep(0.005, 10), rep(1, 10))) {
  stopifnot(all(cell_cycle_display_sample(weights, seed=17)$retained))
}
saturated <- c(0.1, 0.7, 1, 1)
display <- cell_cycle_display_sample(saturated, seed=17)
stopifnot(identical(display$probabilities, saturated), all(display$retained[3:4]))
for (bad in list(numeric(), c(0, 0.1), c(NA, 0.1), c(Inf, 0.1), c(-1, 0.1), c(1.1, 0.1))) {
  stopifnot(inherits(try(cell_cycle_display_sample(bad, seed=17), silent=TRUE), 'try-error'))
}
cat('Display rescaling, guaranteed inclusion, reproducibility and unchanged inputs passed\n')
points <- rbind(cbind(rnorm(400), rnorm(400)), cbind(rnorm(20, 8), rnorm(20, 8)))
input <- tempfile(fileext='.csv')
output <- tempfile(fileext='.csv')
fwrite(data.table(tricycle_pc1=points[, 1], tricycle_pc2=points[, 2]), input)
checksum <- tools::md5sum(input)
result <- estimate_cell_cycle_blob(input, output)
stopifnot(result$status == 'ok', identical(checksum, tools::md5sum(input)))
reference <- robustbase::covMcd(points, alpha=0.5, nsamp='deterministic',
                               raw.only=FALSE, use.correction=TRUE)
expected <- mahalanobis(points, reference$center, reference$cov) > qchisq(0.99, 2)
stopifnot(result$n_outside == sum(expected), result$prop_outside == mean(expected))
blob <- read_cell_cycle_blob(output)
stopifnot(nrow(blob$ellipse) == 361L, grepl('99% robust ellipse:', blob$label, fixed=TRUE))
boundary <- as.matrix(blob$ellipse)
stopifnot(max(abs(mahalanobis(boundary, reference$center, reference$cov) - qchisq(0.99, 2))) < 1e-8)
for (bad in list(matrix(1, 10, 2), matrix(1:6, 3, 2), cbind(1:10, 1:10))) {
  fwrite(data.table(tricycle_pc1=bad[, 1], tricycle_pc2=bad[, 2]), input)
  suppressWarnings(estimate_cell_cycle_blob(input, output))
  unavailable <- read_cell_cycle_blob(output)
  stopifnot(unavailable$summary$status == 'unavailable', nrow(unavailable$ellipse) == 0L,
            nzchar(unavailable$summary$reason))
}
cat('MCD reference counts, ellipse, input preservation and degenerate-input tests passed\n')
