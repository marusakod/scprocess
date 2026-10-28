from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
import scprocess_utils


# An explicitly requested annotation step needs at least one classifier.
@pytest.mark.parametrize('config', [{}, {'label_celltypes': []}, {'label_celltypes': None}])
def test_label_target_requires_classifiers(config):
  with pytest.raises(KeyError, match = 'No classifiers configured'):
    scprocess_utils.check_config_ok_for_rule(config = config, rule = 'label_celltypes')


def test_label_target_accepts_configured_model():
  config = {'label_celltypes': [{'labeller': 'scprocess', 'model': 'local'}]}
  scprocess_utils.check_config_ok_for_rule(config = config, rule = 'label_celltypes')


def test_join_without_classifiers_needs_no_registry(tmp_path):
  schema_f = Path(__file__).parents[1] / 'resources/schemas/join.schema.json'
  labeller_entries = scprocess_utils.get_join_labeller_parameters(
    config = {}, schema_f = schema_f, scdata_dir = tmp_path)

  assert labeller_entries == []


def test_join_xgboost_defaults_without_celltypist_registry(tmp_path, monkeypatch):
  # Isolate registry lookup while checking defaults and the returned metadata.
  def model_paths(entry, scdata_dir):
    assert scdata_dir == tmp_path
    entry['label_map_f'] = 'labels.csv'

  monkeypatch.setattr(scprocess_utils, 'get_xgboost_model_paths', model_paths)
  schema_f = Path(__file__).parents[1] / 'resources/schemas/join.schema.json'
  config = {'label_celltypes': [{'labeller': 'scprocess', 'model': 'local'}]}
  labeller_entries = scprocess_utils.get_join_labeller_parameters(
    config = config, schema_f = schema_f, scdata_dir = tmp_path)

  assert labeller_entries[0]['hi_res_cl'] == 'RNA_snn_res.2'
  assert labeller_entries[0]['min_cl_prop'] == 0.5
  assert labeller_entries[0]['save_cluster_names_file'] is False
  assert labeller_entries[0]['label_map_f'] == 'labels.csv'
