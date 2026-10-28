"""Aggregate tricycle projections and estimate a pooled project origin."""

import argparse
import gzip
from pathlib import Path

import numpy as np
import polars as pl
from scipy.ndimage import gaussian_filter


TRICYCLE_COLUMNS = [
  "cell_id", "sample_id", "tricycle_projection_group",
  "tricycle_raw_pc1", "tricycle_raw_pc2",
]
CELL_CYCLE_MARKER_GENES = [
  "CDKN1A",
  "PCNA", "MCM4", "MCM5", "TYMS",
  "TOP2A", "CCNB1", "CDK1", "TPX2", "SMC2", "MKI67", "UBE2C",
]
MARKER_EXPRESSION_COLUMNS = ["cell_id", *CELL_CYCLE_MARKER_GENES]

# Fixed numerical/display methods, rather than additional analysis options.
DENSITY_WEIGHT_ITERATIONS = 80
DIAGNOSTIC_DRAWS = 100


def aggregate_tricycle(
  score_fs, marker_expression_fs, summary_fs, coldata_f, required_h5_fs,
  bandwidth_multiplier, kde_grid_size, target_prop, target_prop_grid, seed,
  out_scores_f, out_marker_expression_f, out_origin_f, out_sensitivity_f,
  out_diagnostics_f, manual_origin = None
):
  """Combine run projections, estimate the project origin and save cell phases."""
  scores_df, marker_expression_df, summaries_df = _load_tricycle_inputs(
    score_fs, marker_expression_fs, summary_fs)

  # Centre all projected cells together, then select QC-passed origin candidates.
  scores_df, projection_center  = _center_project_projection(scores_df)
  candidate_scores_df           = _load_origin_candidates(scores_df, coldata_f, required_h5_fs)
  points                        = candidate_scores_df.select("tricycle_pc1", "tricycle_pc2").to_numpy()
  density, bandwidth, n_bins    = _binned_kde_density(
    points, bandwidth_multiplier, kde_grid_size)

  # Expected inclusion probabilities determine the origin; no cells are drawn here.
  sensitivity_df, probabilities, estimated_origin = _estimate_origin_sensitivity(
    points, density, target_prop, target_prop_grid)
  origin_coordinates = estimated_origin.copy()
  if manual_origin is not None:
    manual_origin = np.asarray(manual_origin, dtype = float)
    if manual_origin.shape != (2,) or not np.isfinite(manual_origin).all():
      raise ValueError("manual tricycle origin must contain two finite coordinates")
    origin_coordinates = manual_origin
  origin_method = "manual" if manual_origin is not None else "kde_density_equalized"

  # Assign an angle around the chosen origin to every projected cell.
  all_points  = scores_df.select("tricycle_pc1", "tricycle_pc2").to_numpy()
  theta       = np.mod(np.arctan2(
    all_points[:, 1] - origin_coordinates[1],
    all_points[:, 0] - origin_coordinates[0]), 2 * np.pi)
  scores_df   = scores_df.with_columns(pl.Series("tricycle_theta", theta))

  # This seeded draw is only for diagnostic plots, not origin estimation.
  display_retained, display_draw_index = _diagnostic_downsample(probabilities, seed)
  diagnostics_df  = candidate_scores_df.with_columns(
    pl.Series("kde_density", density),
    pl.Series("inclusion_probability", probabilities),
    pl.Series("diagnostic_retained", display_retained))
  origin_df       = _make_origin_summary(
    origin_method = origin_method, origin_coordinates = origin_coordinates,
    estimated_origin = estimated_origin, projection_center = projection_center,
    n_projection_cells = scores_df.height, n_candidate_cells = candidate_scores_df.height,
    target_prop = float(target_prop), manual_origin = manual_origin,
    probabilities = probabilities, bandwidth = bandwidth,
    bandwidth_multiplier = bandwidth_multiplier, n_bins = n_bins, seed = seed,
    display_draw_index = display_draw_index, display_retained = display_retained,
    summaries_df = summaries_df)

  # Save per-cell projections/markers/diagnostics and the two summary tables.
  for output_f in (out_scores_f, out_marker_expression_f, out_origin_f,
                   out_sensitivity_f, out_diagnostics_f):
    Path(output_f).parent.mkdir(parents = True, exist_ok = True)
  with gzip.open(out_scores_f, "wb") as handle:
    scores_df.write_csv(handle)
  with gzip.open(out_marker_expression_f, "wb") as handle:
    marker_expression_df.write_csv(handle)
  origin_df.write_csv(out_origin_f)
  sensitivity_df.write_csv(out_sensitivity_f)
  with gzip.open(out_diagnostics_f, "wb") as handle:
    diagnostics_df.write_csv(handle)
  return scores_df, marker_expression_df, origin_df, sensitivity_df, diagnostics_df


def _load_tricycle_inputs(score_fs, marker_expression_fs, summary_fs):
  """Load per-cell projections/markers and verify their run-level summaries."""
  # One row per projected cell; each input covers one physical run.
  score_dfs = [pl.read_csv(path) for path in score_fs]
  if not score_dfs:
    raise ValueError("At least one tricycle score table is required")
  for run_df in score_dfs:
    cols_missing = set(TRICYCLE_COLUMNS) - set(run_df.columns)
    if cols_missing:
      raise KeyError(f"tricycle score table is missing: {', '.join(sorted(cols_missing))}")
  scores_df = pl.concat([run_df.select(TRICYCLE_COLUMNS) for run_df in score_dfs])
  if scores_df["cell_id"].n_unique() != scores_df.height:
    raise ValueError("tricycle score tables contain duplicate cell IDs")
  if scores_df.height == 0:
    raise ValueError("tricycle score tables are empty")

  # Sparse marker columns can contain only integer-looking zeros throughout
  # Polars' schema-inference window and non-zero floating-point values later.
  # Their type is known, so declare it rather than relying on sampled values.
  marker_schema = {gene: pl.Float64 for gene in CELL_CYCLE_MARKER_GENES}
  marker_dfs = [
    pl.read_csv(path, schema_overrides = marker_schema)
    for path in marker_expression_fs
  ]
  if len(marker_dfs) != len(score_dfs):
    raise ValueError("Each tricycle score table requires one marker-expression table")
  for run_df in marker_dfs:
    cols_missing = set(MARKER_EXPRESSION_COLUMNS) - set(run_df.columns)
    if cols_missing:
      raise KeyError(
        f"tricycle marker-expression table is missing: {', '.join(sorted(cols_missing))}"
      )
  marker_expression_df = pl.concat([
    run_df.select(MARKER_EXPRESSION_COLUMNS) for run_df in marker_dfs
  ])
  if marker_expression_df["cell_id"].n_unique() != marker_expression_df.height:
    raise ValueError("tricycle marker-expression tables contain duplicate cell IDs")
  if set(marker_expression_df["cell_id"].to_list()) != set(scores_df["cell_id"].to_list()):
    raise ValueError("tricycle scores and marker expression contain different cell IDs")

  # Projection summaries have one row per run and must use the same reference.
  summaries_df = pl.concat([pl.read_csv(path) for path in summary_fs])
  if summaries_df.height != len(summary_fs) or summaries_df.height == 0:
    raise ValueError("Each tricycle score table requires one projection summary")
  for column in ("reference", "tricycle_version", "species"):
    if summaries_df[column].n_unique() != 1:
      raise ValueError(f"Projection summaries disagree on {column}")
  return scores_df, marker_expression_df, summaries_df


def _center_project_projection(scores_df):
  """Subtract one common project mean, including projected known doublets."""
  # Centre all run projections once using their pooled project mean.
  raw_points = scores_df.select("tricycle_raw_pc1", "tricycle_raw_pc2").to_numpy()
  if not np.isfinite(raw_points).all():
    raise ValueError("tricycle raw projections contain non-finite coordinates")
  projection_center = raw_points.mean(axis = 0)
  scores_df = scores_df.with_columns(
    (pl.col("tricycle_raw_pc1") - projection_center[0]).alias("tricycle_pc1"),
    (pl.col("tricycle_raw_pc2") - projection_center[1]).alias("tricycle_pc2"),
  ).drop("tricycle_raw_pc1", "tricycle_raw_pc2")

  return scores_df, projection_center


def _load_origin_candidates(scores_df, coldata_f, required_h5_fs):
  """Return one row per QC-passed cell used to estimate the project origin."""
  # Every cell entering preliminary PCA must have a projection.
  required_ids = []
  for path in required_h5_fs:
    required_ids.extend(_read_barcodes(path))
  missing_required = set(required_ids) - set(scores_df["cell_id"].to_list())
  if missing_required:
    preview = ", ".join(sorted(missing_required)[:20])
    raise ValueError(f"cells entering preliminary PCA lack tricycle scores: {preview}")

  coldata_df = pl.read_csv(coldata_f, ignore_errors = True)
  required_metadata = {"cell_id", "keep"}
  cols_missing = required_metadata - set(coldata_df.columns)
  if cols_missing:
    raise KeyError(f"coldata is missing: {', '.join(sorted(cols_missing))}")
  if coldata_df["cell_id"].n_unique() != coldata_df.height:
    raise ValueError("coldata contains duplicate cell IDs")

  candidate_cells_df = coldata_df.filter(pl.col("keep") == True).select("cell_id")
  candidate_scores_df = candidate_cells_df.join(scores_df, on = "cell_id", how = "inner")
  if candidate_scores_df.height != candidate_cells_df.height:
    raise ValueError("some QC-passed origin candidates lack tricycle scores")
  return candidate_scores_df


def _read_barcodes(path):
  import h5py

  with h5py.File(path, "r") as handle:
    return [x.decode("utf-8") for x in handle["matrix/barcodes"][:]]


def _binned_kde_density(points, bandwidth_multiplier = 1.0, grid_size = 500):
  """Evaluate a scalable Gaussian KDE approximation at the input points."""
  points = np.asarray(points, dtype = float)
  if points.ndim != 2 or points.shape[1] != 2:
    raise ValueError("KDE input must have two columns")
  if points.shape[0] < 3:
    raise ValueError("At least three cells are required to estimate the origin")
  if not np.isfinite(points).all():
    raise ValueError("KDE input contains non-finite coordinates")

  n_bins = int(grid_size)
  if n_bins < 20:
    raise ValueError("kde_grid_size must be at least 20")
  lower = points.min(axis = 0)
  upper = points.max(axis = 0)
  coordinate_span = upper - lower
  if np.any(coordinate_span <= 0):
    raise ValueError("Tricycle coordinates are degenerate")
  # Bin cells on a regular two-dimensional grid before Gaussian smoothing.
  edges = [
    np.linspace(lower[i], upper[i], n_bins + 1)
    for i in range(2)
  ]
  binned_counts, _, _ = np.histogram2d(points[:, 0], points[:, 1], bins = edges)

  # Convert coordinate bandwidths to bin units for the Gaussian filter.
  bandwidth = np.array([
    _bandwidth_nrd(points[:, 0]), _bandwidth_nrd(points[:, 1])
  ]) * float(bandwidth_multiplier)
  bin_width = np.array([edge[1] - edge[0] for edge in edges])
  grid_density = gaussian_filter(
    binned_counts, sigma = bandwidth / bin_width, mode = "constant", truncate = 4.0
  )
  grid_density /= points.shape[0] * np.prod(bin_width)

  # Read each cell's density from the grid bin containing its projection.
  cell_bins = []
  for axis in range(2):
    bin_index = np.searchsorted(edges[axis], points[:, axis], side = "right") - 1
    cell_bins.append(np.clip(bin_index, 0, n_bins - 1))
  density = grid_density[cell_bins[0], cell_bins[1]]
  if np.any(~np.isfinite(density)) or np.any(density <= 0):
    raise ValueError("KDE returned non-positive or non-finite densities")
  return density, bandwidth, n_bins


def _bandwidth_nrd(values):
  """Match MASS::kde2d's effective univariate normal-reference bandwidth."""
  values = np.asarray(values, dtype = float)
  iqr_scale = (np.quantile(values, 0.75) - np.quantile(values, 0.25)) / 1.34
  scale = min(float(np.std(values, ddof = 1)), float(iqr_scale))
  if not np.isfinite(scale) or scale <= 0:
    raise ValueError("KDE coordinates must have positive finite spread")
  return 1.06 * scale * values.size ** (-0.2)


def _estimate_origin_sensitivity(points, density, target_prop, target_prop_grid):
  """Compare weighted origins at the requested expected retained fractions."""
  sensitivity_targets = sorted(set(float(value) for value in [target_prop, *target_prop_grid]))
  if any(not np.isfinite(value) or not 0 < value <= 1 for value in sensitivity_targets):
    raise ValueError("target_prop and target_prop_grid must be in (0, 1]")
  selected_target = float(target_prop)
  sensitivity_rows = []
  selected_probabilities = None
  center = None
  for sensitivity_target in sensitivity_targets:
    probabilities = _density_equalized_probabilities(density, sensitivity_target * len(points))
    candidate_center = np.average(points, axis = 0, weights = probabilities)
    sensitivity_rows.append({
      "target_prop": sensitivity_target,
      "expected_retained_cells": float(probabilities.sum()),
      "tricycle_pc1_origin": candidate_center[0],
      "tricycle_pc2_origin": candidate_center[1],
      "selected": sensitivity_target == selected_target,
    })
    if sensitivity_target == selected_target:
      selected_probabilities = probabilities
      center = candidate_center
  sensitivity_df = pl.DataFrame(sensitivity_rows)
  return sensitivity_df, selected_probabilities, center


def _density_equalized_probabilities(density, target_cells):
  """Return capped inverse-density inclusion probabilities with fixed sum."""
  density = np.asarray(density, dtype = float)
  if np.any(~np.isfinite(density)) or np.any(density <= 0):
    raise ValueError("Density values must be finite and positive")
  expected_target = min(float(target_cells), density.size)
  if not np.isfinite(expected_target) or expected_target <= 0:
    raise ValueError("target_cells must be positive")
  if expected_target == density.size:
    return np.ones(density.size)

  # Lower-density cells receive larger probabilities, capped at one.
  lower_scale = 0.0
  upper_scale = float(density.max())

  # Solve sum(min(1, scale / density)) = expected_target by bisection.
  # The target is an expected number of cells and need not be an integer.
  for _ in range(DENSITY_WEIGHT_ITERATIONS):
    trial_scale = (lower_scale + upper_scale) / 2
    expected_cells = np.minimum(1.0, trial_scale / density).sum()
    if expected_cells < expected_target:
      lower_scale = trial_scale
    else:
      upper_scale = trial_scale
  probabilities = np.minimum(1.0, upper_scale / density)
  probabilities *= expected_target / probabilities.sum()
  return np.minimum(1.0, probabilities)


def _diagnostic_downsample(probabilities, seed, n_draws = DIAGNOSTIC_DRAWS):
  """Choose the closest-sized Bernoulli draw for display; never use it to fit theta."""
  rng = np.random.default_rng(int(seed))
  expected_cells = float(np.sum(probabilities))
  best_retained = None
  closest_size_error = np.inf
  best_draw_index = None
  for index in range(n_draws):
    draw = rng.random(len(probabilities)) < probabilities
    size_error = abs(int(draw.sum()) - expected_cells)
    if size_error < closest_size_error:
      best_retained = draw
      closest_size_error = size_error
      best_draw_index = index
  return best_retained, best_draw_index


def _make_origin_summary(*, origin_method, origin_coordinates, estimated_origin,
  projection_center, n_projection_cells, n_candidate_cells, target_prop,
  manual_origin, probabilities, bandwidth, bandwidth_multiplier, n_bins,
  seed, display_draw_index, display_retained, summaries_df):
  """Record the fitted origin and provenance as a single-row table."""
  origin_df = pl.DataFrame({
    "origin_method": [origin_method],
    "projection_centring_method": ["pooled_project_mean"],
    "tricycle_raw_pc1_project_mean": [projection_center[0]],
    "tricycle_raw_pc2_project_mean": [projection_center[1]],
    "n_projection_center_cells": [n_projection_cells],
    "tricycle_pc1_origin": [origin_coordinates[0]],
    "tricycle_pc2_origin": [origin_coordinates[1]],
    "selected_kde_tricycle_pc1_origin": [estimated_origin[0]],
    "selected_kde_tricycle_pc2_origin": [estimated_origin[1]],
    "target_prop_used_for_origin": [manual_origin is None],
    "n_candidate_cells": [n_candidate_cells],
    "target_prop": [target_prop],
    "expected_retained_cells": [float(probabilities.sum())],
    "bandwidth_pc1": [bandwidth[0]],
    "bandwidth_pc2": [bandwidth[1]],
    "bandwidth_multiplier": [float(bandwidth_multiplier)],
    "kde_grid_bins": [n_bins],
    "seed": [int(seed)],
    "diagnostic_draw": [display_draw_index],
    "diagnostic_draws_tried": [DIAGNOSTIC_DRAWS],
    "n_diagnostic_retained": [int(display_retained.sum())],
    "projection_groups": [summaries_df.height],
    "reference": [summaries_df["reference"][0]],
    "tricycle_version": [summaries_df["tricycle_version"][0]],
    "species": [summaries_df["species"][0]],
    "reference_genes_matched_min": [summaries_df["n_reference_genes"].min()],
    "reference_genes_matched_max": [summaries_df["n_reference_genes"].max()],
    "duplicate_feature_mappings_total": [summaries_df["n_duplicate_feature_mappings"].sum()]
  })

  return origin_df


def main():
  parser = argparse.ArgumentParser()

  # File paths and repeated lists of run-level input tables.
  for argument in ('coldata_f', 'out_scores_f', 'out_marker_expression_f',
                   'out_origin_f', 'out_sensitivity_f', 'out_diagnostics_f'):
    parser.add_argument('--' + argument, required = True)
  for argument in ('score_f', 'marker_expression_f', 'summary_f'):
    parser.add_argument('--' + argument, action = 'append', required = True)
  parser.add_argument('--required_h5_f', action = 'append', default = [])

  # Numerical settings retain the existing CLI interface and defaults.
  parser.add_argument('--kde_grid_size', type = int, default = 500)
  parser.add_argument('--seed', type = int, default = 20230308)
  parser.add_argument('--bandwidth_multiplier', type = float, default = 1.0)
  parser.add_argument('--target_prop', type = float, required = True)
  parser.add_argument('--target_prop_grid', action = 'append', type = float, default = [])
  parser.add_argument("--origin_manual", nargs = 2, type = float)
  args = parser.parse_args()

  aggregate_tricycle(
    score_fs = args.score_f,
    marker_expression_fs = args.marker_expression_f,
    summary_fs = args.summary_f,
    coldata_f = args.coldata_f,
    required_h5_fs = args.required_h5_f,
    bandwidth_multiplier = args.bandwidth_multiplier,
    kde_grid_size = args.kde_grid_size,
    target_prop = args.target_prop,
    target_prop_grid = args.target_prop_grid,
    seed = args.seed,
    out_scores_f = args.out_scores_f,
    out_marker_expression_f = args.out_marker_expression_f,
    out_origin_f = args.out_origin_f,
    out_sensitivity_f = args.out_sensitivity_f,
    out_diagnostics_f = args.out_diagnostics_f,
    manual_origin = args.origin_manual)


if __name__ == "__main__":
  main()
