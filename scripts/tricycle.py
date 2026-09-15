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


def _read_barcodes(path):
  import h5py

  with h5py.File(path, "r") as handle:
    return [x.decode("utf-8") for x in handle["matrix/barcodes"][:]]


def _bandwidth_nrd(values):
  """Match MASS::kde2d's effective univariate normal-reference bandwidth."""
  values = np.asarray(values, dtype=float)
  iqr_scale = (np.quantile(values, 0.75) - np.quantile(values, 0.25)) / 1.34
  scale = min(float(np.std(values, ddof=1)), float(iqr_scale))
  if not np.isfinite(scale) or scale <= 0:
    raise ValueError("KDE coordinates must have positive finite spread")
  return 1.06 * scale * values.size ** (-0.2)


def _binned_kde_density(points, bandwidth_multiplier=1.0, grid_size=500):
  """Evaluate a scalable Gaussian KDE approximation at the input points."""
  points = np.asarray(points, dtype=float)
  if points.ndim != 2 or points.shape[1] != 2:
    raise ValueError("KDE input must have two columns")
  if points.shape[0] < 3:
    raise ValueError("At least three cells are required to estimate the origin")
  if not np.isfinite(points).all():
    raise ValueError("KDE input contains non-finite coordinates")

  n_bins = int(grid_size)
  if n_bins < 20:
    raise ValueError("kde_grid_size must be at least 20")
  lower = points.min(axis=0)
  upper = points.max(axis=0)
  span = upper - lower
  if np.any(span <= 0):
    raise ValueError("Tricycle coordinates are degenerate")
  edges = [
    np.linspace(lower[i], upper[i], n_bins + 1)
    for i in range(2)
  ]
  hist, _, _ = np.histogram2d(points[:, 0], points[:, 1], bins=edges)

  bandwidth = np.array([
    _bandwidth_nrd(points[:, 0]), _bandwidth_nrd(points[:, 1])
  ]) * float(bandwidth_multiplier)
  bin_width = np.array([edge[1] - edge[0] for edge in edges])
  smoothed = gaussian_filter(
    hist, sigma=bandwidth / bin_width, mode="constant", truncate=4.0
  )
  smoothed /= points.shape[0] * np.prod(bin_width)

  indices = []
  for axis in range(2):
    idx = np.searchsorted(edges[axis], points[:, axis], side="right") - 1
    indices.append(np.clip(idx, 0, n_bins - 1))
  density = smoothed[indices[0], indices[1]]
  if np.any(~np.isfinite(density)) or np.any(density <= 0):
    raise ValueError("KDE returned non-positive or non-finite densities")
  return density, bandwidth, n_bins


def _density_equalized_probabilities(density, target_cells):
  """Return capped inverse-density inclusion probabilities with fixed sum."""
  density = np.asarray(density, dtype=float)
  if np.any(~np.isfinite(density)) or np.any(density <= 0):
    raise ValueError("Density values must be finite and positive")
  target = min(int(target_cells), density.size)
  if target < 1:
    raise ValueError("target_cells must be positive")
  if target == density.size:
    return np.ones(density.size)

  lo = 0.0
  hi = float(density.max())
  for _ in range(80):
    mid = (lo + hi) / 2
    total = np.minimum(1.0, mid / density).sum()
    if total < target:
      lo = mid
    else:
      hi = mid
  probabilities = np.minimum(1.0, hi / density)
  probabilities *= target / probabilities.sum()
  return np.minimum(1.0, probabilities)


def aggregate_tricycle(
  score_fs, marker_expression_fs, summary_fs, coldata_f, required_h5_fs,
  bandwidth_multiplier, kde_grid_size, target_cells, target_cells_grid, seed,
  out_scores_f, out_marker_expression_f, out_origin_f, out_sensitivity_f,
  out_diagnostics_f, manual_origin=None
):
  score_tables = [pl.read_csv(path) for path in score_fs]
  if not score_tables:
    raise ValueError("At least one tricycle score table is required")
  for table in score_tables:
    missing = set(TRICYCLE_COLUMNS) - set(table.columns)
    if missing:
      raise KeyError(f"tricycle score table is missing: {', '.join(sorted(missing))}")
  scores = pl.concat([table.select(TRICYCLE_COLUMNS) for table in score_tables])
  if scores["cell_id"].n_unique() != scores.height:
    raise ValueError("tricycle score tables contain duplicate cell IDs")
  if scores.height == 0:
    raise ValueError("tricycle score tables are empty")

  # Sparse marker columns can contain only integer-looking zeros throughout
  # Polars' schema-inference window and non-zero floating-point values later.
  # Their type is known, so declare it rather than relying on sampled values.
  marker_schema = {gene: pl.Float64 for gene in CELL_CYCLE_MARKER_GENES}
  marker_tables = [
    pl.read_csv(path, schema_overrides=marker_schema)
    for path in marker_expression_fs
  ]
  if len(marker_tables) != len(score_tables):
    raise ValueError("Each tricycle score table requires one marker-expression table")
  for table in marker_tables:
    missing = set(MARKER_EXPRESSION_COLUMNS) - set(table.columns)
    if missing:
      raise KeyError(
        f"tricycle marker-expression table is missing: {', '.join(sorted(missing))}"
      )
  marker_expression = pl.concat([
    table.select(MARKER_EXPRESSION_COLUMNS) for table in marker_tables
  ])
  if marker_expression["cell_id"].n_unique() != marker_expression.height:
    raise ValueError("tricycle marker-expression tables contain duplicate cell IDs")
  if set(marker_expression["cell_id"].to_list()) != set(scores["cell_id"].to_list()):
    raise ValueError("tricycle scores and marker expression contain different cell IDs")

  summaries = pl.concat([pl.read_csv(path) for path in summary_fs])
  if summaries.height != len(summary_fs) or summaries.height == 0:
    raise ValueError("Each tricycle score table requires one projection summary")
  for column in ("reference", "tricycle_version", "species"):
    if summaries[column].n_unique() != 1:
      raise ValueError(f"Projection summaries disagree on {column}")

  required_ids = []
  for path in required_h5_fs:
    required_ids.extend(_read_barcodes(path))
  missing_required = set(required_ids) - set(scores["cell_id"].to_list())
  if missing_required:
    preview = ", ".join(sorted(missing_required)[:20])
    raise ValueError(f"cells entering preliminary PCA lack tricycle scores: {preview}")

  coldata = pl.read_csv(coldata_f, ignore_errors=True)
  required_metadata = {"cell_id", "keep"}
  missing = required_metadata - set(coldata.columns)
  if missing:
    raise KeyError(f"coldata is missing: {', '.join(sorted(missing))}")
  if coldata["cell_id"].n_unique() != coldata.height:
    raise ValueError("coldata contains duplicate cell IDs")

  raw_points = scores.select("tricycle_raw_pc1", "tricycle_raw_pc2").to_numpy()
  if not np.isfinite(raw_points).all():
    raise ValueError("tricycle raw projections contain non-finite coordinates")
  projection_center = raw_points.mean(axis=0)
  scores = scores.with_columns(
    (pl.col("tricycle_raw_pc1") - projection_center[0]).alias("tricycle_pc1"),
    (pl.col("tricycle_raw_pc2") - projection_center[1]).alias("tricycle_pc2"),
  ).drop("tricycle_raw_pc1", "tricycle_raw_pc2")

  candidates = coldata.filter(pl.col("keep") == True).select("cell_id")
  candidate_scores = candidates.join(scores, on="cell_id", how="inner")
  if candidate_scores.height != candidates.height:
    raise ValueError("some QC-passed origin candidates lack tricycle scores")
  points = candidate_scores.select("tricycle_pc1", "tricycle_pc2").to_numpy()
  density, bandwidth, n_bins = _binned_kde_density(
    points, bandwidth_multiplier, kde_grid_size
  )
  selected_target = min(int(target_cells), candidate_scores.height)
  sensitivity_targets = sorted({
    min(int(value), candidate_scores.height)
    for value in [target_cells, *target_cells_grid]
  })
  sensitivity_rows = []
  selected_probabilities = None
  center = None
  for sensitivity_target in sensitivity_targets:
    probabilities = _density_equalized_probabilities(density, sensitivity_target)
    candidate_center = np.average(points, axis=0, weights=probabilities)
    sensitivity_rows.append({
      "target_cells": sensitivity_target,
      "expected_retained_cells": float(probabilities.sum()),
      "tricycle_pc1_origin": candidate_center[0],
      "tricycle_pc2_origin": candidate_center[1],
      "selected": sensitivity_target == selected_target,
    })
    if sensitivity_target == selected_target:
      selected_probabilities = probabilities
      center = candidate_center
  sensitivity = pl.DataFrame(sensitivity_rows)
  probabilities = selected_probabilities
  automatic_center = center.copy()
  if manual_origin is not None:
    manual_origin = np.asarray(manual_origin, dtype=float)
    if manual_origin.shape != (2,) or not np.isfinite(manual_origin).all():
      raise ValueError("manual tricycle origin must contain two finite coordinates")
    center = manual_origin
  origin_method = "manual" if manual_origin is not None else "kde_density_equalized"

  all_points = scores.select("tricycle_pc1", "tricycle_pc2").to_numpy()
  theta = np.mod(np.arctan2(all_points[:, 1] - center[1],
                            all_points[:, 0] - center[0]), 2 * np.pi)
  scores = scores.with_columns(pl.Series("tricycle_theta", theta))

  rng = np.random.default_rng(int(seed))
  diagnostics = candidate_scores.with_columns(
    pl.Series("kde_density", density),
    pl.Series("inclusion_probability", probabilities),
    pl.Series("diagnostic_retained", rng.random(probabilities.size) < probabilities)
  )
  origin = pl.DataFrame({
    "origin_method": [origin_method],
    "projection_centring_method": ["pooled_project_mean"],
    "tricycle_raw_pc1_project_mean": [projection_center[0]],
    "tricycle_raw_pc2_project_mean": [projection_center[1]],
    "n_projection_center_cells": [scores.height],
    "tricycle_pc1_origin": [center[0]],
    "tricycle_pc2_origin": [center[1]],
    "selected_kde_tricycle_pc1_origin": [automatic_center[0]],
    "selected_kde_tricycle_pc2_origin": [automatic_center[1]],
    "target_cells_used_for_origin": [manual_origin is None],
    "n_candidate_cells": [candidate_scores.height],
    "target_cells": [selected_target],
    "expected_retained_cells": [float(probabilities.sum())],
    "bandwidth_pc1": [bandwidth[0]],
    "bandwidth_pc2": [bandwidth[1]],
    "bandwidth_multiplier": [float(bandwidth_multiplier)],
    "kde_grid_bins": [n_bins],
    "seed": [int(seed)],
    "projection_groups": [summaries.height],
    "reference": [summaries["reference"][0]],
    "tricycle_version": [summaries["tricycle_version"][0]],
    "species": [summaries["species"][0]],
    "reference_genes_matched_min": [summaries["n_reference_genes"].min()],
    "reference_genes_matched_max": [summaries["n_reference_genes"].max()],
    "duplicate_feature_mappings_total": [summaries["n_duplicate_feature_mappings"].sum()]
  })

  for path in (
    out_scores_f, out_marker_expression_f, out_origin_f, out_sensitivity_f,
    out_diagnostics_f,
  ):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
  with gzip.open(out_scores_f, "wb") as handle:
    scores.write_csv(handle)
  with gzip.open(out_marker_expression_f, "wb") as handle:
    marker_expression.write_csv(handle)
  origin.write_csv(out_origin_f)
  sensitivity.write_csv(out_sensitivity_f)
  with gzip.open(out_diagnostics_f, "wb") as handle:
    diagnostics.write_csv(handle)
  return scores, marker_expression, origin, sensitivity, diagnostics


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--score_f", action="append", required=True)
  parser.add_argument("--marker_expression_f", action="append", required=True)
  parser.add_argument("--summary_f", action="append", required=True)
  parser.add_argument("--coldata_f", required=True)
  parser.add_argument("--required_h5_f", action="append", default=[])
  parser.add_argument("--bandwidth_multiplier", type=float, default=1.0)
  parser.add_argument("--kde_grid_size", type=int, default=500)
  parser.add_argument("--target_cells", type=int, default=5000)
  parser.add_argument("--target_cells_grid", action="append", type=int, default=[])
  parser.add_argument("--origin", nargs=2, type=float)
  parser.add_argument("--seed", type=int, default=20230308)
  parser.add_argument("--out_scores_f", required=True)
  parser.add_argument("--out_marker_expression_f", required=True)
  parser.add_argument("--out_origin_f", required=True)
  parser.add_argument("--out_sensitivity_f", required=True)
  parser.add_argument("--out_diagnostics_f", required=True)
  args = parser.parse_args()
  aggregate_tricycle(
    args.score_f, args.marker_expression_f, args.summary_f, args.coldata_f,
    args.required_h5_f, args.bandwidth_multiplier, args.kde_grid_size,
    args.target_cells, args.target_cells_grid, args.seed, args.out_scores_f,
    args.out_marker_expression_f, args.out_origin_f, args.out_sensitivity_f,
    args.out_diagnostics_f, manual_origin=args.origin
  )


if __name__ == "__main__":
  main()
