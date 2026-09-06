import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_project_and_zoom_default_to_bpcells_pca():
  for schema_name in ('config.schema.json', 'zoom.schema.json'):
    schema = json.loads((ROOT / 'resources' / 'schemas' / schema_name).read_text())
    pca_method = schema['properties']['integration']['properties']['int_pca_method']
    assert pca_method['default'] == 'bpcells'
    assert pca_method['enum'] == ['bpcells', 'scanpy']


def test_bpcells_pca_uses_lazy_scaled_matrix():
  script = (ROOT / 'scripts' / 'bpcells_pca.R').read_text()
  assert 'write_matrix_dir(mat_norm' in script
  assert 'mat_scaled <- (mat_norm - gene_means) / gene_sds' in script
  assert 'pca_matrix <- t(mat_scaled)' in script
  assert 'irlba(pca_matrix' in script
  assert 'as.matrix(mat_scaled)' not in script


def test_primary_integration_has_persisted_clean_cell_boundary():
  rules = (ROOT / 'rules' / 'integration.smk').read_text()
  assert 'rule classify_additional_doublets:' in rules
  assert 'clean_cells_f = CLEAN_CELLS_F' in rules
  assert 'rule final_pca:' in rules
  assert 'cells_f      = \'{input.clean_cells_f}\'' in rules
