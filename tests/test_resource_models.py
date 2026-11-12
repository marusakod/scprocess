import copy
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts.scprocess_utils import (
  RESOURCE_CHEMISTRY_GROUPS, _check_resource_models, _get_chemistry_group,
  _load_resource_models, get_resources, prep_resource_params)


ROOT = Path(__file__).resolve().parents[1]
SCHEMA_F = ROOT / 'resources/schemas/config.schema.json'


def _params(run_chemistries=None):
  params = prep_resource_params({}, ROOT)
  chemistries = run_chemistries or {}
  params['run_chemistries'] = chemistries
  params['RUNS_TO_LIBS'] = {run: run for run in chemistries}
  params['n_batches'] = 10
  return params


def _inputs(tmp_path, size_gb):
  # sparse file whose size is the model input, e.g. 2 -> a 2 GB 'af_h5_f'
  af_h5_f = tmp_path / 'af.h5'
  with af_h5_f.open('wb') as handle:
    handle.truncate(int(size_gb * 1024**3))
  return SimpleNamespace(af_h5_f=str(af_h5_f))


def test_chemistry_groups():
  # 5v3 -> other: chemistries without their own values form the named group 'other'
  assert _get_chemistry_group('3v3') == '3v3'
  assert _get_chemistry_group('3v4') == '3v4'
  assert _get_chemistry_group('flexv1') == 'flex'
  assert _get_chemistry_group('flexv2') == 'flex'
  for chemistry in ['5v3', '5v1', '5v2', '3LT', 'multiome', '3v2']:
    assert _get_chemistry_group(chemistry) == 'other'


def test_shipped_resource_files_are_valid():
  models, floors = _load_resource_models(ROOT)
  assert 'run_cellbender' in models
  assert set(models['run_cellbender']['memory']) == set(RESOURCE_CHEMISTRY_GROUPS)
  assert 'all' in models['run_integration']['memory']


def test_chemistry_specific_model_follows_the_run_chemistry(tmp_path):
  params = _params({'run_a': '3v3', 'run_b': '3v4', 'run_c': '5v3', 'run_d': 'flexv1'})
  rules = SimpleNamespace(run_cellbender=object())
  inputs = _inputs(tmp_path, 2)

  memory = {
    run: get_resources(params, rules, inputs, 'run_cellbender', 'memory', 1, run)
    for run in ['run_a', 'run_b', 'run_c', 'run_d']}

  assert memory['run_a'] != memory['run_b']
  assert memory['run_c'] == memory['run_a']

  # Flex starts at twice the 3v3 model (provisional)
  assert memory['run_d'] == pytest.approx(2 * memory['run_a'], rel=0.01)


def test_chemistry_specific_model_needs_a_run(tmp_path):
  rules = SimpleNamespace(run_cellbender=object())
  with pytest.raises(ValueError, match='pass the run'):
    get_resources(_params(), rules, _inputs(tmp_path, 1), 'run_cellbender', 'memory', 1)


def test_floor_applies_to_a_small_input(tmp_path):
  rules = SimpleNamespace(run_qc_one_run=object())
  params = _params({'run_a': '3v3'})
  memory = get_resources(params, rules, _inputs(tmp_path, 0), 'run_qc_one_run', 'memory', 1, 'run_a')
  assert memory == 2048


def _valid_models():
  model = {'model_var': 'n_smpls_pre_qc', 'rq_intercept': 1.0, 'rq_slope': 2.0}
  models = {
    'rule_a': {'memory': {'all': dict(model)}},
    'rule_b': {'time': {group: dict(model) for group in RESOURCE_CHEMISTRY_GROUPS}}}
  floors = {'rule_a': {'memory': 2}, 'rule_b': {'time': 10}}
  return models, floors


def test_check_accepts_valid_models():
  _check_resource_models(*_valid_models())


def test_check_rejects_a_missing_chemistry_group():
  models, floors = _valid_models()
  del models['rule_b']['time']['flex']
  with pytest.raises(ValueError, match="rule 'rule_b', time"):
    _check_resource_models(models, floors)


def test_check_rejects_all_mixed_with_groups():
  models, floors = _valid_models()
  models['rule_b']['time']['all'] = copy.deepcopy(models['rule_b']['time']['3v3'])
  with pytest.raises(ValueError, match='expected'):
    _check_resource_models(models, floors)


def test_check_rejects_a_missing_floor():
  models, floors = _valid_models()
  del floors['rule_a']['memory']
  with pytest.raises(ValueError, match='no floor'):
    _check_resource_models(models, floors)


def test_check_accepts_a_floor_without_a_model():
  # a rule without a model requests exactly its floor
  models, floors = _valid_models()
  floors['rule_c'] = {'memory': 2, 'time': 10}
  _check_resource_models(models, floors)


@pytest.mark.parametrize('bad_floor', [0, -1, 'eight', None, True])
def test_check_rejects_a_floor_that_is_not_a_positive_number(bad_floor):
  models, floors = _valid_models()
  floors['rule_c'] = {'memory': bad_floor}
  with pytest.raises(ValueError, match='positive number'):
    _check_resource_models(models, floors)


def test_check_rejects_an_unknown_floor_parameter():
  models, floors = _valid_models()
  floors['rule_c'] = {'cpus': 2}
  with pytest.raises(ValueError, match="'memory' or 'time'"):
    _check_resource_models(models, floors)


def test_check_rejects_unexpected_model_fields():
  models, floors = _valid_models()
  models['rule_a']['memory']['all']['floor'] = 1
  with pytest.raises(ValueError, match='expected fields'):
    _check_resource_models(models, floors)


def test_rule_without_a_model_requests_its_floor_and_scales_memory_with_retries():
  params = _params()
  rules = SimpleNamespace(aggregate_tricycle=object())
  # floors: 8 GB and 30 min
  assert get_resources(params, rules, SimpleNamespace(), 'aggregate_tricycle', 'memory', 1) == 8192
  assert get_resources(params, rules, SimpleNamespace(), 'aggregate_tricycle', 'memory', 2) == 12288
  assert get_resources(params, rules, SimpleNamespace(), 'aggregate_tricycle', 'time', 3) == 30


def test_user_value_equal_to_the_floor_is_honoured_over_the_model(tmp_path):
  # run_qc_one_run time has a model with a 10 min floor; a user value of 10 is still an override
  rules = SimpleNamespace(run_qc_one_run=object())
  inputs = _inputs(tmp_path, 5)
  params = _params({'run_a': '3v3'})
  assert get_resources(params, rules, inputs, 'run_qc_one_run', 'time', 1, 'run_a') > 10

  params['user_vals'] = {'mins_run_qc_one_run': 10}
  assert get_resources(params, rules, inputs, 'run_qc_one_run', 'time', 1, 'run_a') == 10


def test_a_rule_with_no_user_value_model_or_floor_is_an_error():
  rules = SimpleNamespace(not_a_real_rule=object())
  with pytest.raises(KeyError, match='No user value, resource model or floor'):
    get_resources(_params(), rules, SimpleNamespace(), 'not_a_real_rule', 'memory', 1)


def test_user_value_does_not_need_a_run_for_a_chemistry_dependent_model(tmp_path):
  rules = SimpleNamespace(run_cellbender=object())
  params = _params()
  params['user_vals'] = {'gb_run_cellbender': 12}
  assert get_resources(params, rules, _inputs(tmp_path, 1), 'run_cellbender', 'memory', 1) == 12288


# zoom rules use the same fitted model as the project rule they follow
ZOOM_MODEL_SOURCES = {
  ('zoom_aggregate_labels', 'memory'): 'merge_labels',
  ('zoom_create_hvg_matrix', 'memory'): 'create_hvg_matrix',
  ('zoom_get_highly_variable_genes', 'memory'): 'get_highly_variable_genes',
  ('zoom_get_stats_for_std_variance_for_sample', 'memory'): 'get_stats_for_std_variance_for_sample',
  ('zoom_make_one_pb_cells', 'memory'): 'make_one_pb_cells',
  ('zoom_make_tmp_csr_matrix', 'memory'): 'make_tmp_csr_matrix',
  ('zoom_make_tmp_csr_matrix', 'time'): 'make_tmp_csr_matrix',
  ('zoom_run_fgsea', 'memory'): 'run_fgsea',
  ('zoom_run_integration', 'memory'): 'run_integration'}


@pytest.mark.parametrize('zoom_rule, param', sorted(ZOOM_MODEL_SOURCES))
def test_zoom_rules_use_the_model_and_floor_of_their_project_rule(zoom_rule, param):
  models, floors = _load_resource_models(ROOT)
  project_rule = ZOOM_MODEL_SOURCES[(zoom_rule, param)]
  assert models[zoom_rule][param] == models[project_rule][param]
  assert floors[zoom_rule][param] == floors[project_rule][param]


def test_flex_models_are_twice_the_3v3_models():
  models, _ = _load_resource_models(ROOT)
  checked = 0
  for rule, params in models.items():
    for param, chem_models in params.items():
      if 'all' in chem_models:
        continue
      for field in ('rq_intercept', 'rq_slope'):
        assert chem_models['flex'][field] == pytest.approx(2 * chem_models['3v3'][field], rel=1e-6)
      assert chem_models['flex']['model_var'] == chem_models['3v3']['model_var']
      checked += 1
  assert checked == 7
