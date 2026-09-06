library(data.table)

source("scripts/shiny.R")

# Four coincident centroids compete for the same nearest cell-free position.
# The allocator must remain deterministic while reserving distinct footprints.
all_meta = CJ(
  UMAP_1 = seq(-0.5, 0.5, length.out = 50L),
  UMAP_2 = seq(-0.5, 0.5, length.out = 50L)
)
centroids = data.table(
  cluster = paste0("cl", seq_len(4L)),
  UMAP_1 = 0,
  UMAP_2 = 0
)

grid_res = 80L
r_label = 4L
label_gap = 1L
padding_frac = 0.05

positions_a = .compute_repel_positions(
  all_meta,
  centroids,
  grid_res = grid_res,
  max_radius_frac = 0.2,
  r_clear = 2L,
  r_label = r_label,
  label_gap = label_gap,
  padding_frac = padding_frac
)
positions_b = .compute_repel_positions(
  all_meta,
  centroids,
  grid_res = grid_res,
  max_radius_frac = 0.2,
  r_clear = 2L,
  r_label = r_label,
  label_gap = label_gap,
  padding_frac = padding_frac
)

stopifnot(
  identical(positions_a, positions_b),
  identical(positions_a$cluster, centroids$cluster)
)

to_grid = function(values, source_values) {
  source_range = range(source_values)
  source_span = diff(source_range)
  grid_min = source_range[1L] - source_span * padding_frac
  grid_step = source_span * (1 + 2 * padding_frac) / grid_res
  round((values - grid_min) / grid_step + 0.5)
}

label_grid = cbind(
  to_grid(positions_a$repel_1, all_meta$UMAP_1),
  to_grid(positions_a$repel_2, all_meta$UMAP_2)
)
label_distances = as.matrix(dist(label_grid, method = "maximum"))
minimum_separation = 2L * r_label + label_gap

stopifnot(
  all(label_distances[upper.tri(label_distances)] >= minimum_separation)
)

# Degenerate embeddings cannot define a two-dimensional placement grid.
degenerate = data.table(UMAP_1 = 0, UMAP_2 = seq_len(5L))
error_message = tryCatch(
  {
    .compute_repel_positions(degenerate, centroids[1L])
    NA_character_
  },
  error = function(e) conditionMessage(e)
)
stopifnot(grepl("must span both dimensions", error_message, fixed = TRUE))
