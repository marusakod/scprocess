"""Aggregate tricycle projections and estimate a pooled project origin."""

import argparse
import gzip
from pathlib import Path

import numpy as np
import polars as pl
from scipy.ndimage import gaussian_filter


TRICYCLE_COLUMNS = [
  "cell_id", "sample_id", "tricycle_center_group", "tricycle_pc1", "tricycle_pc2"
]


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
  score_fs, summary_fs, coldata_f, required_h5_fs, bandwidth_multiplier, kde_grid_size,
  target_cells, seed, out_scores_f, out_origin_f, out_diagnostics_f
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

  candidates = coldata.filter(pl.col("keep") == True).select("cell_id")
  candidate_scores = candidates.join(scores, on="cell_id", how="inner")
  if candidate_scores.height != candidates.height:
    raise ValueError("some QC-passed origin candidates lack tricycle scores")
  points = candidate_scores.select("tricycle_pc1", "tricycle_pc2").to_numpy()
  density, bandwidth, n_bins = _binned_kde_density(
    points, bandwidth_multiplier, kde_grid_size
  )
  probabilities = _density_equalized_probabilities(density, target_cells)
  center = np.average(points, axis=0, weights=probabilities)

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
    "origin_method": ["kde_density_equalized"],
    "tricycle_pc1_origin": [center[0]],
    "tricycle_pc2_origin": [center[1]],
    "n_candidate_cells": [candidate_scores.height],
    "target_cells": [min(int(target_cells), candidate_scores.height)],
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

  for path in (out_scores_f, out_origin_f, out_diagnostics_f):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
  with gzip.open(out_scores_f, "wb") as handle:
    scores.write_csv(handle)
  origin.write_csv(out_origin_f)
  with gzip.open(out_diagnostics_f, "wb") as handle:
    diagnostics.write_csv(handle)
  return scores, origin, diagnostics


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--score_f", action="append", required=True)
  parser.add_argument("--summary_f", action="append", required=True)
  parser.add_argument("--coldata_f", required=True)
  parser.add_argument("--required_h5_f", action="append", default=[])
  parser.add_argument("--bandwidth_multiplier", type=float, default=1.0)
  parser.add_argument("--kde_grid_size", type=int, default=500)
  parser.add_argument("--target_cells", type=int, default=5000)
  parser.add_argument("--seed", type=int, default=20230308)
  parser.add_argument("--out_scores_f", required=True)
  parser.add_argument("--out_origin_f", required=True)
  parser.add_argument("--out_diagnostics_f", required=True)
  args = parser.parse_args()
  aggregate_tricycle(
    args.score_f, args.summary_f, args.coldata_f, args.required_h5_f,
    args.bandwidth_multiplier, args.kde_grid_size, args.target_cells, args.seed,
    args.out_scores_f, args.out_origin_f, args.out_diagnostics_f
  )


if __name__ == "__main__":
  main()
