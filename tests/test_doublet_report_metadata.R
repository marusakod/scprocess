# Evaluate only the reader: no integration calculation or project data is loaded.
definitions = parse('scripts/integration.R')
for (definition in definitions) {
  if (is.call(definition) && length(definition) >= 3L &&
      identical(as.character(definition[[2L]]), 'read_doublet_cluster_threshold')) {
    eval(definition)
    break
  }
}

manual = data.frame(dbl_threshold = 0.3, dbl_threshold_method = 'manual')
observed = read_doublet_cluster_threshold(manual)
stopifnot(observed$threshold == 0.3, observed$method == 'manual')
adaptive = manual
adaptive$dbl_threshold_method = 'adaptive'
for (key in c('n_low', 'n_high', 'n_equal_50', 'median_logit_low', 'median_logit_high',
    'scaled_mad_logit_low', 'scaled_mad_logit_high', 'weight_low', 'weight_high')) {
  adaptive[[paste0('dbl_threshold_', key)]] = 1
}
stopifnot(read_doublet_cluster_threshold(adaptive)$method == 'adaptive')

# Older reports must not infer the cutoff from a current configuration file.
missing_threshold = tryCatch(read_doublet_cluster_threshold(data.frame(cell_id = 'c1')), error = identity)
stopifnot(inherits(missing_threshold, 'error'),
  grepl('rerun integration', conditionMessage(missing_threshold)))
for (method in c('fixed', 'mean_mad', 'split_mad_mixture')) {
  old = manual
  old$dbl_threshold_method = method
  result = tryCatch(read_doublet_cluster_threshold(old), error = identity)
  stopifnot(inherits(result, 'error'), grepl('Unsupported recorded', conditionMessage(result)))
}
cat('PASS: current manual/adaptive metadata; older incomplete outputs rejected\n')
