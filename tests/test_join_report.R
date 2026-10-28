source("scripts/utils.R")
source("scripts/integration.R")

ids <- paste0("cell", seq_len(60L))
predictions <- data.table(cell_id = ids,
  predicted_label_naive = c(rep("A", 25L), rep("B", 4L), "rare",
                           rep("A", 5L), rep("B", 24L), "rare"))
metadata <- data.table(cell_id = ids,
  cluster = rep(c("cl1", "cl2"), each = 30L))
comparison <- prepare_join_label_comparison(predictions, metadata, "cluster")
stopifnot(sum(comparison$N) == 60L,
  setequal(as.character(comparison$cl1), c("A", "B", "Other")),
  comparison[cl1 == "Other" & cl2 == "cl1", N] == 1L,
  abs(comparison[cl1 == "A" & cl2 == "cl1", p_cl2] - 25 / 30) < 1e-12,
  all(abs(comparison[, sum(p_cl1), by = cl1]$V1 - 1) < 1e-12),
  all(abs(comparison[, sum(p_cl2), by = cl2]$V1 - 1) < 1e-12))

# A comparison column literally named cl1 must retain numeric heatmap values.
plot_file <- tempfile(fileext = ".pdf")
pdf(plot_file)
for (normalization in c("row", "column")) {
  plot <- plot_join_label_comparison(comparison, "test", "Clusters", normalization)
  draw(plot)
}
dev.off()
stopifnot(file.info(plot_file)$size > 0)
unlink(plot_file)

predictions[1, predicted_label_naive := NA_character_]
without_missing <- prepare_join_label_comparison(predictions, metadata, "cluster")
stopifnot(sum(without_missing$N) == 59L)
error <- tryCatch(prepare_join_label_comparison(rbind(predictions, predictions[1]),
  metadata, "cluster"), error = identity)
stopifnot(inherits(error, "error"))
message("PASS: join label proportions, Other aggregation, both heatmaps and ID validation")
