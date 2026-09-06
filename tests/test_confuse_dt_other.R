suppressPackageStartupMessages({
  library(assertthat)
  library(data.table)
  library(forcats)
  library(magrittr)
})

source("scripts/utils.R")

labels = data.table(
  cell_id = as.character(seq_len(20L)),
  prediction = factor(c(
    rep("kept", 4L), paste0("low", seq_len(6L)),
    rep("kept", 4L), paste0("low", seq_len(4L)), "low7", "low8"
  ))
)
clusters = data.table(
  cell_id = as.character(seq_len(20L)),
  zoom_cluster = factor(rep(c("cl01", "cl02"), each = 10L))
)

confusion = calc_confuse_dt(
  labels,
  clusters,
  "prediction",
  "zoom_cluster",
  min_cl2_p = 0.3,
  other_cl1_label = "Other"
)

stopifnot(
  setequal(as.character(unique(confusion$cl1)), c("kept", "Other")),
  all(abs(confusion[, sum(p_cl2), by = cl2]$V1 - 1) < 1e-12),
  confusion[cl1 == "Other" & cl2 == "cl01", N] == 6L,
  confusion[cl1 == "Other" & cl2 == "cl02", N] == 6L
)
