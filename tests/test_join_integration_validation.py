import unittest

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


if __name__ == '__main__':
  unittest.main()
