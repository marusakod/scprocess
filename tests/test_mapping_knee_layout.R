# Run from the scprocess repository root in the workflow R environment.
source("scripts/mapping.R")
library(ggrepel)
folder = tempfile("knee-layout-")
dir.create(folder)
pdf(file.path(folder, "layout-device.pdf"))
on.exit_cleanup = function() {
  dev.off()
  unlink(folder, recursive = TRUE)
}
tryCatch({
  samples = paste0("sample", seq_len(12))
  files = setNames(file.path(folder, paste0(samples, ".csv")), samples)
  for (f in files) fwrite(data.table(total = c(10000, 1000, 100, 10), rank = 1:4), f)
  params = data.table(sample_id = samples, knee1 = 10000, shin1 = 1000,
    knee2 = 100, shin2 = 10, total_droplets_included = 4, expected_cells = 2)
  for (n in c(1, 2, 5, 8, 11, 12)) {
    p = plot_barcode_ranks_w_params(files[seq_len(n)], params,
      "sample_id", show_lines = FALSE, n_slots = 12)
    stopifnot(p$facet$params$nrow == 3, p$facet$params$ncol == 4)
    grob = ggplotGrob(p)
    # Empty slots have no panel grob, but must retain their layout units.
    stopifnot(sum(grid::unitType(grob$heights) == "null") == 3,
      sum(grid::unitType(grob$widths) == "null") == 4)
    ggsave(file.path(folder, paste0("layout-", n, ".png")), p, width = 9, height = 6)
  }
  cat("PASS: 1, 2, 5, 8, 11 and 12 samples preserve a four-column, three-row grid\n")
}, finally = on.exit_cleanup())
