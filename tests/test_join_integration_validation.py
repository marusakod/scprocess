import unittest
import warnings

from scripts.scprocess_utils import _check_join_integration_parameters


class TestJoinIntegrationValidation(unittest.TestCase):
  def test_matching_arrays_are_valid(self):
    config = {'integration': {
      'int_batch_var': ['sample_id', 'project_id'],
      'int_theta': [0.1, 1.0],
    }}
    _check_join_integration_parameters(config)

  def test_mismatched_arrays_are_rejected(self):
    config = {'integration': {
      'int_batch_var': ['sample_id', 'project_id'],
      'int_theta': [0.1],
    }}
    with self.assertRaisesRegex(ValueError, 'must have matching lengths'):
      _check_join_integration_parameters(config)

  def test_scalar_and_array_remain_supported(self):
    for batch_var, theta in [(['sample_id', 'project_id'], 0.1), ('sample_id', [0.1])]:
      with self.subTest(batch_var=batch_var, theta=theta):
        _check_join_integration_parameters({'integration': {
          'int_batch_var': batch_var,
          'int_theta': theta,
        }})

  def test_gpu_multiple_batch_variables_do_not_warn_about_cpu_fallback(self):
    config = {'integration': {
      'int_batch_var': ['project_id', 'sample_id'],
      'int_theta': [1.0, 0.1],
      'int_embedding': 'harmony',
      'int_use_gpu': True,
    }}
    with warnings.catch_warnings(record=True) as caught:
      _check_join_integration_parameters(config)
    self.assertEqual(caught, [])

  def test_cycle_regression_defaults_without_new_score_configuration(self):
    config = {'integration': {'int_batch_var': 'sample_id', 'int_theta': 0.1,
      'int_cell_cycle_correction': 'within_sample'}}
    _check_join_integration_parameters(config)
    self.assertEqual(config['integration']['int_cell_cycle_harmonics'], 2)
    self.assertEqual(config['integration']['int_cell_cycle_ridge_lambda'], 0.1)

  def test_correction_is_opt_in(self):
    config = {'integration': {'int_batch_var': 'sample_id', 'int_theta': 0.1}}
    _check_join_integration_parameters(config)
    self.assertNotIn('int_cell_cycle_correction', config['integration'])
    self.assertNotIn('int_cell_cycle_harmonics', config['integration'])

  def test_cycle_regression_rejects_scanpy(self):
    with self.assertRaisesRegex(ValueError, 'bpcells'):
      _check_join_integration_parameters({'integration': {
        'int_batch_var': 'sample_id', 'int_theta': 0.1,
        'int_cell_cycle_correction': 'global', 'int_pca_method': 'scanpy'}})

  def test_cycle_tuning_requires_correction(self):
    for key in ('int_cell_cycle_harmonics', 'int_cell_cycle_ridge_lambda'):
      with self.subTest(key=key), self.assertRaisesRegex(ValueError, 'require'):
        _check_join_integration_parameters({'integration': {
          'int_batch_var': 'sample_id', 'int_theta': 0.1, key: 1}})

  def test_cycle_schema_constraints(self):
    import json
    from pathlib import Path
    import jsonschema
    schema = json.loads((Path(__file__).resolve().parents[1] /
      'resources/schemas/join.schema.json').read_text())['properties']['integration']
    # Disable PAGA so its conditional required field cannot mask these checks.
    jsonschema.validate({'int_use_paga': False,
      'int_cell_cycle_correction': 'within_sample',
      'int_cell_cycle_harmonics': 2, 'int_cell_cycle_ridge_lambda': 0.1}, schema)
    for params in ({'int_cell_cycle_harmonics': 3},
                   {'int_cell_cycle_ridge_lambda': -0.1},
                   {'int_cell_cycle_correction': 'invalid'}):
      with self.subTest(params=params), self.assertRaises(jsonschema.ValidationError):
        jsonschema.validate({'int_use_paga': False, **params}, schema)


if __name__ == '__main__':
  unittest.main()
