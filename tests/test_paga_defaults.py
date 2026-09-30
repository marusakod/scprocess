import json
import unittest
from pathlib import Path

import jsonschema

from scripts.scprocess_utils import _get_default_config_from_schema, get_paga_resolution


SCHEMA_DIR = Path(__file__).parents[1] / 'resources' / 'schemas'


class TestPagaDefaults(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    cls.schemas = {}
    for name in ('config', 'zoom', 'join'):
      with open(SCHEMA_DIR / f'{name}.schema.json') as fh:
        cls.schemas[name] = json.load(fh)

  def test_paga_resolution_defaults_to_two_in_all_workflows(self):
    for name, schema in self.schemas.items():
      with self.subTest(schema=name):
        defaults = _get_default_config_from_schema(schema)
        self.assertEqual(defaults['integration']['int_paga_cl_res'], 2)

  def test_marker_resolution_remains_independent(self):
    for name, schema in self.schemas.items():
      with self.subTest(schema=name):
        defaults = _get_default_config_from_schema(schema)
        self.assertEqual(defaults['marker_genes']['mkr_sel_res'], 0.2)

  def test_rule_resolution_uses_schema_or_explicit_value(self):
    for schema in self.schemas.values():
      integration = _get_default_config_from_schema(schema)['integration']
      self.assertEqual(get_paga_resolution(integration), integration['int_paga_cl_res'])
      integration['int_paga_cl_res'] = 0.5
      self.assertEqual(get_paga_resolution(integration), 0.5)

  def test_rule_resolution_is_absent_when_paga_is_disabled(self):
    self.assertIsNone(get_paga_resolution({'int_use_paga': False}))

  def test_rule_resolution_does_not_supply_a_missing_default(self):
    with self.assertRaises(KeyError):
      get_paga_resolution({'int_use_paga': True})

  def test_disabled_paga_does_not_inject_forbidden_resolution(self):
    for name, schema in self.schemas.items():
      with self.subTest(schema=name):
        requested = {'integration': {'int_use_paga': False}}
        defaults = _get_default_config_from_schema(schema, requested)
        integration = defaults['integration'] | requested['integration']
        self.assertNotIn('int_paga_cl_res', integration)
        jsonschema.validate(integration, schema['properties']['integration'])

  def test_explicit_resolution_with_disabled_paga_is_still_invalid(self):
    for name, schema in self.schemas.items():
      with self.subTest(schema=name):
        requested = {'integration': {'int_use_paga': False, 'int_paga_cl_res': 2}}
        defaults = _get_default_config_from_schema(schema, requested)
        integration = defaults['integration'] | requested['integration']
        with self.assertRaises(jsonschema.ValidationError):
          jsonschema.validate(integration, schema['properties']['integration'])

  def test_cell_type_high_resolution_default_remains_two(self):
    for name in ('config', 'join'):
      with self.subTest(schema=name):
        label_props = self.schemas[name]['properties']['label_celltypes']
        default = label_props['items']['properties']['hi_res_cl']['default']
        self.assertEqual(default, 'RNA_snn_res.2')


if __name__ == '__main__':
  unittest.main()
