"""Canonical inputs remain supported; retired interfaces are rejected."""
from pathlib import Path
import importlib.util
import json
import sys

import polars as pl
import pytest
from jsonschema import Draft202012Validator

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import train_xgboost


@pytest.mark.parametrize('name', ['config', 'join', 'zoom'])
def test_training_rejects_legacy_mapping_input(name):
  schema = json.loads((ROOT / f'resources/schemas/{name}.schema.json').read_text())
  validator = Draft202012Validator(schema['properties']['train_xgboost'])
  assert validator.is_valid({'annots_f': 'cells.csv', 'label_specs_f': 'labels.csv'})
  assert not validator.is_valid({'annots_f': 'cells.csv', 'label_map_f': 'labels.csv'})


@pytest.mark.parametrize('name', ['config', 'join', 'zoom'])
def test_shiny_rejects_old_display_name_array(name):
  schema = json.loads((ROOT / f'resources/schemas/{name}.schema.json').read_text())
  validator = Draft202012Validator(schema['properties']['shiny'])
  assert validator.is_valid({'metadata_labels': {'donor_id': 'Donor'}})
  assert not validator.is_valid({'var_names': ['Donor']})


def test_setup_requires_current_profile_name():
  schema = json.loads((ROOT / 'resources/schemas/setup.schema.json').read_text())
  validator = Draft202012Validator(schema['properties']['user'])
  assert validator.is_valid({'profile_template': 'slurm_shpc'})
  assert validator.is_valid({'local_cores': 2})
  assert not validator.is_valid({'profile': 'slurm_shpc'})
  assert not validator.is_valid({'profile_template': 'slurm_shpc', 'local_cores': 2})


@pytest.mark.parametrize('name,retired,current', [
  ('config', 'gb_estimate_run_tricycle', 'gb_calculate_tricycle_projection_one_run'),
  ('config', 'mins_estimate_run_tricycle', 'mins_calculate_tricycle_projection_one_run'),
  ('join', 'gb_join_zoom_hvg_stats', 'gb_join_zoom_calculate_raw_gene_count_statistics'),
  ('join', 'gb_join_zoom_hvg_trend', 'gb_join_zoom_fit_expected_gene_variances'),
  ('join', 'gb_join_zoom_hvg_clip', 'gb_join_zoom_calculate_clipped_gene_count_statistics'),
  ('join', 'gb_join_zoom_hvgs', 'gb_join_zoom_rank_and_select_highly_variable_genes'),
])
def test_resource_aliases_are_rejected(name, retired, current):
  schema = json.loads((ROOT / f'resources/schemas/{name}.schema.json').read_text())
  validator = Draft202012Validator(schema['properties']['resources'])
  assert validator.is_valid({current: 4})
  assert not validator.is_valid({retired: 4})


@pytest.mark.parametrize('column', ['annotation', 'fine_label'])
def test_annotations_require_label_column(tmp_path, column):
  annots = tmp_path / 'annotations.csv'
  pl.DataFrame({'cell_id': ['c1'], column: ['Neuron']}).write_csv(annots)
  with pytest.raises(ValueError, match="must have a 'label' column"):
    train_xgboost._load_training_metadata(annots, 'unused.csv', 'sample_id', False, None, set())


def test_only_manifest_models_are_preinstalled():
  spec = importlib.util.spec_from_file_location('setup_without_retired_models', ROOT / 'scripts/setup.py')
  setup = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(setup)
  models = setup.get_xgboost_parameters({})
  assert set(models) == set(train_xgboost.builtin_models())
  assert 'Siletti_hippocampus-2026-06-22' not in models
