import gzip
import importlib.util
import json
import math
from pathlib import Path

import polars as pl
import pytest


ROOT = Path(__file__).resolve().parents[1]


def _load_module(name, path):
  spec = importlib.util.spec_from_file_location(name, path)
  module = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(module)
  return module


utils = _load_module("scprocess_utils", ROOT / "scripts" / "scprocess_utils.py")
zoom = _load_module("zoom", ROOT / "scripts" / "zoom.py")


def test_qc_max_schema_defaults_are_disabled():
  with open(ROOT / "resources" / "schemas" / "config.schema.json") as handle:
    qc_schema = json.load(handle)["properties"]["qc"]["properties"]

  assert qc_schema["qc_max_counts"]["default"] is None
  assert qc_schema["qc_max_feats"]["default"] is None


@pytest.mark.parametrize(
  "qc, message",
  [
    ({"qc_min_counts": 500, "qc_max_counts": 499}, "qc_max_counts"),
    ({"qc_min_feats": 300, "qc_max_feats": 299}, "qc_max_feats"),
  ],
)
def test_qc_max_must_not_be_below_min(qc, message):
  with pytest.raises(ValueError, match=message):
    utils._validate_qc_bounds(qc)


def test_qc_max_equal_to_min_is_valid():
  utils._validate_qc_bounds({
    "qc_min_counts": 500,
    "qc_max_counts": 500,
    "qc_min_feats": 300,
    "qc_max_feats": 300,
  })


def test_custom_batch_max_overrides_are_propagated_and_validated():
  config = {"qc": {
    "qc_min_counts": 500,
    "qc_max_counts": None,
    "qc_min_feats": 300,
    "qc_max_feats": None,
    "qc_min_mito": 0,
    "qc_max_mito": 0.1,
    "qc_min_splice": 0,
    "qc_max_splice": 1,
    "qc_min_cells": 100,
  }}
  custom = {"sample_a": {"qc": {"qc_max_counts": 5000, "qc_max_feats": 2000}}}

  result = utils._get_batch_parameters_one_batch("sample_a", config, custom)
  assert result["qc"]["qc_max_counts"] == 5000
  assert result["qc"]["qc_max_feats"] == 2000

  custom["sample_a"]["qc"]["qc_max_counts"] = 499
  with pytest.raises(ValueError, match="qc for sample_a"):
    utils._get_batch_parameters_one_batch("sample_a", config, custom)


def test_zoom_max_thresholds_are_inclusive(tmp_path):
  labels_f = tmp_path / "labels.csv"
  qc_f = tmp_path / "qc.csv.gz"
  output_f = tmp_path / "filtered.csv.gz"

  pl.DataFrame({
    "cell_id": ["inside", "counts_boundary", "feats_boundary", "too_many_counts", "too_many_feats"],
    "label": ["selected"] * 5,
  }).write_csv(labels_f)
  qc_df = pl.DataFrame({
    "cell_id": ["inside", "counts_boundary", "feats_boundary", "too_many_counts", "too_many_feats"],
    "keep": [True] * 5,
    "log_counts": [math.log10(v) for v in [3000, 5000, 3000, 5001, 3000]],
    "log_feats": [math.log10(v) for v in [1000, 1000, 2000, 1000, 2001]],
    "logit_mito": [0.0] * 5,
    "logit_spliced": [0.0] * 5,
  })
  with gzip.open(qc_f, "wb") as handle:
    qc_df.write_csv(handle)

  zoom.filter_zoom_labels_by_qc(
    labels_f=labels_f,
    qc_all_f=qc_f,
    output_f=output_f,
    labels_col="label",
    sel_labels=["selected"],
    qc_max_counts=5000,
    qc_max_feats=2000,
  )

  kept = pl.read_csv(output_f)["cell_id"].to_list()
  assert kept == ["inside", "counts_boundary", "feats_boundary"]


def test_zoom_sample_exclusions_apply_without_qc_thresholds(tmp_path):
  labels_f = tmp_path / "labels.csv"
  qc_f = tmp_path / "qc.csv.gz"
  output_f = tmp_path / "filtered.csv.gz"

  pl.DataFrame({
    "cell_id": ["a1", "a2", "b1"],
    "sample_id": ["sample_a", "sample_a", "sample_b"],
    "label": ["selected"] * 3,
  }).write_csv(labels_f)
  with gzip.open(qc_f, "wb") as handle:
    pl.DataFrame({"cell_id": [], "keep": []}).write_csv(handle)

  zoom.filter_zoom_labels_by_qc(
    labels_f=labels_f,
    qc_all_f=qc_f,
    output_f=output_f,
    labels_col="label",
    sel_labels=["selected"],
    batch_var="sample_id",
    exclude_batches=["sample_a"],
  )

  kept = pl.read_csv(output_f)["cell_id"].to_list()
  assert kept == ["b1"]
