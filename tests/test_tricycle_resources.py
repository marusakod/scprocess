from types import SimpleNamespace
from pathlib import Path

import polars as pl
import pytest

from scripts.scprocess_utils import get_resources, prep_resource_params


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
  rules = SimpleNamespace(calculate_tricycle_projection_one_run=object())
  inputs = SimpleNamespace(counts_h5_f=str(counts_h5_f))
  return rules, inputs


def test_tricycle_memory_model_has_two_gb_floor(tmp_path):
  counts_h5_f = tmp_path / "small.h5"
  counts_h5_f.write_bytes(b"small")
  rules, inputs = _rule_context(counts_h5_f)

  memory = get_resources(
    _resource_params(), rules, inputs,
    "calculate_tricycle_projection_one_run", "memory", attempt=1,
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
    "calculate_tricycle_projection_one_run", "memory", attempt=1,
  )
  second_attempt = get_resources(
    params, rules, inputs,
    "calculate_tricycle_projection_one_run", "memory", attempt=2,
  )

  assert first_attempt == 3372
  assert second_attempt == 5058


def test_tricycle_runtime_and_user_overrides(tmp_path):
  counts_h5_f = tmp_path / "counts.h5"
  counts_h5_f.write_bytes(b"small")
  rules, inputs = _rule_context(counts_h5_f)

  assert get_resources(
    _resource_params(), rules, inputs,
    "calculate_tricycle_projection_one_run", "time", attempt=1,
  ) == 10

  overrides = {
    "gb_calculate_tricycle_projection_one_run": 4,
    "mins_calculate_tricycle_projection_one_run": 20,
  }
  params = _resource_params(overrides)
  assert get_resources(
    params, rules, inputs,
    "calculate_tricycle_projection_one_run", "memory", attempt=1,
  ) == 4096
  assert get_resources(
    params, rules, inputs,
    "calculate_tricycle_projection_one_run", "time", attempt=1,
  ) == 20


@pytest.mark.parametrize('rule_name', [
  'aggregate_tricycle', 'estimate_tricycle_reference_ellipse'])
def test_pooled_tricycle_resources_use_defaults_overrides_and_retries(rule_name):
  root = Path(__file__).resolve().parents[1]
  schema_f = root / 'resources/schemas/config.schema.json'
  rules = SimpleNamespace(**{rule_name: object()})
  inputs = SimpleNamespace()
  defaults = prep_resource_params({}, schema_f, root)
  assert get_resources(defaults, rules, inputs, rule_name, 'memory', 1) == 8192
  assert get_resources(defaults, rules, inputs, rule_name, 'memory', 2) == 12288
  assert get_resources(defaults, rules, inputs, rule_name, 'time', 1) == 30
  overrides = prep_resource_params({'resources': {
    f'gb_{rule_name}': 2, f'mins_{rule_name}': 5}}, schema_f, root)
  assert get_resources(overrides, rules, inputs, rule_name, 'memory', 1) == 2048
  assert get_resources(overrides, rules, inputs, rule_name, 'time', 1) == 5
