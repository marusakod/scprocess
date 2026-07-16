import json
import unittest
from pathlib import Path

from scripts.scprocess_utils import _get_default_config_from_schema


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

  def test_cell_type_high_resolution_default_remains_two(self):
    for name in ('config', 'join'):
      with self.subTest(schema=name):
        label_props = self.schemas[name]['properties']['label_celltypes']
        default = label_props['items']['properties']['hi_res_cl']['default']
        self.assertEqual(default, 'RNA_snn_res.2')


if __name__ == '__main__':
  unittest.main()
