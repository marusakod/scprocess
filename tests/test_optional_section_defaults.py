import json
import unittest
from pathlib import Path

import jsonschema

from scripts.scprocess_utils import _get_default_config_from_schema


SCHEMA_DIR = Path(__file__).parents[1] / 'resources' / 'schemas'


class TestOptionalSectionDefaults(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    cls.schemas = {}
    for name in ('config', 'zoom', 'join'):
      with open(SCHEMA_DIR / f'{name}.schema.json') as fh:
        cls.schemas[name] = json.load(fh)

  def test_absent_train_xgboost_is_not_materialized(self):
    for name, schema in self.schemas.items():
      with self.subTest(schema=name):
        defaults = _get_default_config_from_schema(schema, {})
        self.assertNotIn('train_xgboost', defaults)

  def test_explicit_train_xgboost_receives_nested_defaults(self):
    for name, schema in self.schemas.items():
      with self.subTest(schema=name):
        config = {'train_xgboost': {'annots_f': 'annotations.csv'}}
        defaults = _get_default_config_from_schema(schema, config)
        self.assertIn('train_xgboost', defaults)
        self.assertTrue(defaults['train_xgboost']['refine_labels'])
        self.assertEqual(defaults['train_xgboost']['seed'], 42)
        self.assertFalse(defaults['train_xgboost']['use_gpu'])

  def test_explicit_empty_train_xgboost_still_requires_annotations(self):
    for name, schema in self.schemas.items():
      with self.subTest(schema=name):
        defaults = _get_default_config_from_schema(
          schema, {'train_xgboost': {}})
        with self.assertRaisesRegex(
            jsonschema.ValidationError, "'annots_f' is a required property"):
          jsonschema.validate(
            defaults['train_xgboost'],
            schema['properties']['train_xgboost'])


if __name__ == '__main__':
  unittest.main()
