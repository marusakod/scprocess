# Presentation helpers for join.Rmd.template; source utils.R first.
prepare_join_label_comparison <- function(predictions, integration, column_var,
                                          min_cl2_p = 0.05) {
  stopifnot(all(c('cell_id', 'predicted_label_naive') %in% names(predictions)),
            all(c('cell_id', column_var) %in% names(integration)),
            !anyDuplicated(predictions$cell_id), !anyDuplicated(integration$cell_id))
  labels <- predictions[, .(cell_id, label = as.character(predicted_label_naive))]
  columns <- integration[, .(cell_id, comparison = as.character(get(column_var)))]
  labelled <- merge(labels, columns, by = 'cell_id')
  labelled <- labelled[!is.na(label) & nzchar(label) &
                         !is.na(comparison) & nzchar(comparison)]
  if (!nrow(labelled)) stop('No matched naive labels for ', column_var)
  labelled[, label := factor(label)]
  labelled[, comparison := factor(comparison)]
  confuse <- calc_confuse_dt(labelled, labelled, 'label', 'comparison',
                            min_cl2_p = min_cl2_p, other_cl1_label = 'Other')
  # Selection uses raw column fractions, never log-color pseudocounts.
  stopifnot(sum(confuse$N) == nrow(labelled),
            all(abs(confuse[, sum(p_cl2), by = cl2]$V1 - 1) < 1e-8),
            all(abs(confuse[, sum(p_cl1), by = cl1]$V1 - 1) < 1e-8),
            all(confuse[as.character(cl1) != 'Other', max(p_cl2), by = cl1]$V1 >= min_cl2_p))
  confuse
}

plot_join_label_comparison <- function(confuse, model, column_title,
                                       normalization = c('column', 'row'),
                                       lbl_threshold = 0.05) {
  normalization <- match.arg(normalization)
  plot_cluster_comparison_heatmap(confuse, paste(model, 'naive labels'), column_title,
    plot_var = if (normalization == 'column') 'log_p_cl2' else 'log_p_cl1',
    do_sort = 'hclust', lbl_threshold = lbl_threshold)
}
