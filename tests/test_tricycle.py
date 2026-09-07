import json
from pathlib import Path

import numpy as np
import polars as pl
import pytest

from scripts.tricycle import (
  aggregate_tricycle,
  _bandwidth_nrd,
  _binned_kde_density,
  _density_equalized_probabilities,
)


ROOT = Path(__file__).resolve().parents[1]


def test_cell_cycle_is_optional_and_not_materialized_by_defaults():
  from scripts.scprocess_utils import _get_default_config_from_schema

  schema = json.loads((ROOT / 'resources/schemas/config.schema.json').read_text())
  defaults = _get_default_config_from_schema(schema, {})
  assert 'cell_cycle' not in defaults
  cell_cycle_schema = schema['properties']['cell_cycle']
  assert 'enabled' not in cell_cycle_schema['properties']
  assert 'enabled' not in cell_cycle_schema['properties']['regression']['properties']


def test_present_cell_cycle_and_regression_blocks_receive_defaults():
  from scripts.scprocess_utils import _get_default_config_from_schema

  schema = json.loads((ROOT / 'resources/schemas/config.schema.json').read_text())
  defaults = _get_default_config_from_schema(
    schema, {'cell_cycle': {'regression': {}}}
  )
  assert defaults['cell_cycle']['origin_method'] == 'kde_density_equalized'


def test_density_equalization_is_capped_and_has_requested_expectation():
  density = np.array([0.1, 0.2, 1.0, 2.0, 10.0])
  probabilities = _density_equalized_probabilities(density, 3)
  assert np.all((probabilities > 0) & (probabilities <= 1))
  assert probabilities.sum() == pytest.approx(3)
  assert probabilities[0] == 1
  assert probabilities[-1] < probabilities[-2]


def test_binned_kde_is_deterministic_and_positive():
  rng = np.random.default_rng(42)
  points = np.vstack([
    rng.normal((-2, 0), (0.2, 0.5), size=(500, 2)),
    rng.normal((2, 0), (0.5, 0.2), size=(500, 2)),
  ])
  density_1, bandwidth_1, bins_1 = _binned_kde_density(points, 1.25)
  density_2, bandwidth_2, bins_2 = _binned_kde_density(points, 1.25)
  np.testing.assert_array_equal(density_1, density_2)
  np.testing.assert_array_equal(bandwidth_1, bandwidth_2)
  assert bins_1 == bins_2
  assert np.all(density_1 > 0)


def test_normal_reference_bandwidth_matches_definition():
  values = np.arange(1, 101, dtype=float)
  robust_scale = min(
    np.std(values, ddof=1),
    (np.quantile(values, 0.75) - np.quantile(values, 0.25)) / 1.34,
  )
  assert _bandwidth_nrd(values) == pytest.approx(
    1.06 * robust_scale * values.size ** (-0.2)
  )


def test_theta_convention_wraps_to_zero_two_pi():
  points = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]], dtype=float)
  theta = np.mod(np.arctan2(points[:, 1], points[:, 0]), 2 * np.pi)
  np.testing.assert_allclose(theta, [0, np.pi / 2, np.pi, 3 * np.pi / 2])


def test_aggregate_assigns_theta_to_all_cells_but_fits_origin_from_kept_cells(tmp_path):
  score_f = tmp_path / 'scores.csv'
  summary_f = tmp_path / 'summary.csv'
  coldata_f = tmp_path / 'coldata.csv'
  out_scores_f = tmp_path / 'all_scores.csv.gz'
  out_origin_f = tmp_path / 'origin.csv'
  out_diagnostics_f = tmp_path / 'diagnostics.csv.gz'
  pl_scores = pl.DataFrame({
    'cell_id': ['a', 'b', 'c', 'known_doublet'],
    'sample_id': ['s1', 's1', 's1', 's1'],
    'tricycle_center_group': ['sample:s1'] * 4,
    'tricycle_pc1': [-1.0, 1.0, 0.0, 100.0],
    'tricycle_pc2': [0.0, 0.0, 1.0, 100.0],
  })
  pl_scores.write_csv(score_f)
  pl.DataFrame({
    'reference': ['tricycle::neuroRef'],
    'tricycle_version': ['1.18.0'],
    'species': ['mouse'],
    'n_reference_genes': [450],
    'n_duplicate_feature_mappings': [12],
  }).write_csv(summary_f)
  pl.DataFrame({
    'cell_id': ['a', 'b', 'c', 'known_doublet'],
    'keep': [True, True, True, False],
  }).write_csv(coldata_f)

  scores, origin, diagnostics = aggregate_tricycle(
    [score_f], [summary_f], coldata_f, [], 1.0, 500, 5000, 7,
    out_scores_f, out_origin_f, out_diagnostics_f,
  )

  assert scores.height == 4
  assert diagnostics.height == 3
  assert origin['n_candidate_cells'][0] == 3
  assert origin['reference_genes_matched_min'][0] == 450
  assert np.all((scores['tricycle_theta'].to_numpy() >= 0) &
                (scores['tricycle_theta'].to_numpy() < 2 * np.pi))


def test_ridge_formula_matches_gene_chunking():
  rng = np.random.default_rng(7)
  y = rng.normal(size=(80, 23))
  x = rng.normal(size=(80, 4))
  gram = x.T @ x
  effective_lambda = 0.1 * np.trace(gram) / x.shape[1]
  inverse = np.linalg.solve(gram + effective_lambda * np.eye(x.shape[1]), np.eye(x.shape[1]))
  corrected = y - x @ inverse @ (x.T @ y)
  chunked = np.column_stack([
    block - x @ inverse @ (x.T @ block)
    for block in np.array_split(y, 5, axis=1)
  ])
  np.testing.assert_allclose(chunked, corrected, rtol=1e-13, atol=1e-13)


def test_bpcells_uses_public_custom_operator_not_deprecated_mult():
  script = (ROOT / 'scripts/bpcells_pca.R').read_text()
  assert 'setClass(\n  "RidgeResidualMatrix"' in script
  assert 'signature(x = "RidgeResidualMatrix", y = "numeric")' in script
  assert 'signature(x = "numeric", y = "RidgeResidualMatrix")' in script
  assert 'mult =' not in script


def test_human_identifier_mapping_uses_declared_rowdata_not_annotation_downloads():
  script = (ROOT / 'scripts' / 'tricycle.R').read_text()
  environment = (ROOT / 'envs' / 'tricycle.yaml').read_text()
  assert 'rowdata_f' in script
  assert 'ensembl_id' in script
  assert 'org.Hs.eg.db' not in script
  assert 'bioconductor-org.hs.eg.db' not in environment
  assert 'bioconductor-org.mm.eg.db' not in environment


def test_cell_cycle_outputs_use_dedicated_directory():
  main_rules = (ROOT / 'rules' / 'scprocess.smk').read_text()
  tricycle_rules = (ROOT / 'rules' / 'tricycle.smk').read_text()
  integration_rules = (ROOT / 'rules' / 'integration.smk').read_text()
  assert 'cell_cycle_dir = f"{PROJ_DIR}/output/{SHORT_TAG}_cell_cycle"' in main_rules
  assert "f'{cell_cycle_dir}/tricycle_scores_" in tricycle_rules
  assert "f'{cell_cycle_dir}/tricycle/sample_" in tricycle_rules
  assert "f'{cell_cycle_dir}/final_cell_cycle_regression_" in integration_rules


def test_tricycle_r_expression_preserves_named_vector_quotes():
  rules = (ROOT / 'rules' / 'tricycle.smk').read_text()
  assert "Rscript -e 'source(\"{scprocess_dir}/scripts/tricycle.R\")" in rules
  assert 'counts_h5_fs = {params.counts_r}' in rules
  assert 'Rscript -e "source(' not in rules


def test_invalid_regression_backend_is_rejected():
  from scripts.scprocess_utils import _check_cell_cycle_parameters

  config = {
    'project': {'ref_txome': 'mouse_2024'},
    'integration': {'int_pca_method': 'scanpy'},
    'cell_cycle': {'regression': {}},
  }
  with pytest.raises(ValueError, match='requires integration.int_pca_method: bpcells'):
    _check_cell_cycle_parameters(config)


def test_species_is_inferred_from_standard_reference():
  from scripts.scprocess_utils import _check_cell_cycle_parameters

  config = {
    'project': {'ref_txome': 'human_2024'},
    'integration': {'int_pca_method': 'bpcells'},
    'cell_cycle': {},
  }
  checked = _check_cell_cycle_parameters(config)
  assert checked['cell_cycle']['species'] == 'human'


def test_present_regression_block_receives_defaults():
  from scripts.scprocess_utils import _check_cell_cycle_parameters

  config = {
    'project': {'ref_txome': 'mouse_2024'},
    'integration': {'int_pca_method': 'bpcells'},
    'cell_cycle': {'regression': {}},
  }
  checked = _check_cell_cycle_parameters(config)
  assert checked['cell_cycle']['regression'] == {
    'harmonics': 2,
    'ridge_lambda': 0.1,
  }
