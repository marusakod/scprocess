from types import SimpleNamespace

import polars as pl

from scripts.scprocess_utils import get_resources


def _resource_params(user_vals=None):
  return {
    "defaults": {},
    "user_vals": user_vals or {},
    "lm_df": pl.read_csv(
      "resources/snakemake/resources_lm_params_2026-07-07.csv"
    ),
    "R1_sizes": {},
    "n_batches": 100,
    "n_run_mapping": 8,
    "run_chemistries": {},
    "RUNS_TO_LIBS": {},
    "chem_stats_paths": {},
    "_resolved_chemistries": {},
  }


def _rule_context(counts_h5_f):
  rules = SimpleNamespace(estimate_run_tricycle=object())
  inputs = SimpleNamespace(counts_h5_f=str(counts_h5_f))
  return rules, inputs


def test_tricycle_memory_model_has_two_gb_floor(tmp_path):
  counts_h5_f = tmp_path / "small.h5"
  counts_h5_f.write_bytes(b"small")
  rules, inputs = _rule_context(counts_h5_f)

  memory = get_resources(
    _resource_params(), rules, inputs,
    "estimate_run_tricycle", "memory", attempt=1,
  )

  assert memory == 2048


def test_tricycle_memory_scales_with_its_h5_input_and_retry(tmp_path):
  counts_h5_f = tmp_path / "large.h5"
  with counts_h5_f.open("wb") as handle:
    handle.truncate(128 * 1024**2)
  rules, inputs = _rule_context(counts_h5_f)
  params = _resource_params()

  first_attempt = get_resources(
    params, rules, inputs,
    "estimate_run_tricycle", "memory", attempt=1,
  )
  second_attempt = get_resources(
    params, rules, inputs,
    "estimate_run_tricycle", "memory", attempt=2,
  )

  assert first_attempt == 3372
  assert second_attempt == 5058


def test_tricycle_runtime_and_user_overrides(tmp_path):
  counts_h5_f = tmp_path / "counts.h5"
  counts_h5_f.write_bytes(b"small")
  rules, inputs = _rule_context(counts_h5_f)

  assert get_resources(
    _resource_params(), rules, inputs,
    "estimate_run_tricycle", "time", attempt=1,
  ) == 10

  overrides = {
    "gb_estimate_run_tricycle": 4,
    "mins_estimate_run_tricycle": 20,
  }
  params = _resource_params(overrides)
  assert get_resources(
    params, rules, inputs,
    "estimate_run_tricycle", "memory", attempt=1,
  ) == 4096
  assert get_resources(
    params, rules, inputs,
    "estimate_run_tricycle", "time", attempt=1,
  ) == 20
