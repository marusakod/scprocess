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


def test_partitioned_projection_matches_one_pooled_projection():
  rng = np.random.default_rng(11)
  expression = rng.normal(size=(37, 19))
  rotation = rng.normal(size=(19, 2))
  run_projections = [block @ rotation for block in np.array_split(expression, 4)]
  combined = np.vstack(run_projections)
  combined -= combined.mean(axis=0)
  pooled = (expression - expression.mean(axis=0)) @ rotation
  np.testing.assert_allclose(combined, pooled, rtol=1e-13, atol=1e-13)


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
    'tricycle_projection_group': ['run:s1'] * 4,
    'tricycle_raw_pc1': [-1.0, 1.0, 0.0, 100.0],
    'tricycle_raw_pc2': [0.0, 0.0, 1.0, 100.0],
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
  assert origin['n_projection_center_cells'][0] == 4
  assert origin['reference_genes_matched_min'][0] == 450
  assert origin['projection_centring_method'][0] == 'pooled_project_mean'
  np.testing.assert_allclose(
    scores.select('tricycle_pc1', 'tricycle_pc2').mean().to_numpy(), 0,
    atol=1e-14,
  )
  assert np.all((scores['tricycle_theta'].to_numpy() >= 0) &
                (scores['tricycle_theta'].to_numpy() < 2 * np.pi))


def test_hybrid_ridge_formula_matches_gene_chunking_and_preserves_project_mean():
  rng = np.random.default_rng(7)
  y = rng.normal(size=(80, 23))
  theta = rng.uniform(0, 2 * np.pi, size=80)
  raw = np.column_stack((np.sin(theta), np.cos(theta)))
  groups = np.repeat(np.arange(4), 20)
  fit = raw.copy()
  for group in np.unique(groups):
    rows = groups == group
    fit[rows] -= fit[rows].mean(axis=0)
  rms = np.sqrt(np.mean(fit**2, axis=0))
  fit /= rms
  correction = (raw - raw.mean(axis=0)) / rms
  gram = fit.T @ fit
  effective_lambda = 0.1 * np.trace(gram) / fit.shape[1]
  inverse = np.linalg.solve(
    gram + effective_lambda * np.eye(fit.shape[1]), np.eye(fit.shape[1])
  )
  corrected = y - correction @ inverse @ (fit.T @ y)
  chunked = np.column_stack([
    block - correction @ inverse @ (fit.T @ block)
    for block in np.array_split(y, 5, axis=1)
  ])
  np.testing.assert_allclose(chunked, corrected, rtol=1e-13, atol=1e-13)
  np.testing.assert_allclose(corrected.mean(axis=0), y.mean(axis=0), atol=1e-14)
  coefficients = inverse @ (fit.T @ y)
  observed_group_changes = []
  for group in np.unique(groups):
    rows = groups == group
    np.testing.assert_allclose(fit[rows].mean(axis=0), 0, atol=1e-14)
    change = corrected[rows].mean(axis=0) - y[rows].mean(axis=0)
    expected_change = -correction[rows].mean(axis=0) @ coefficients
    np.testing.assert_allclose(change, expected_change, atol=1e-14)
    observed_group_changes.append(change)
  assert np.max(np.abs(observed_group_changes)) > 1e-6


def test_within_sample_fit_does_not_attribute_sample_offset_to_cycle():
  x = np.array([1, 1, 1, -1, -1, -1, -1, 1], dtype=float)
  y = np.array([0, 0, 0, 0, 4, 4, 4, 4], dtype=float)
  groups = np.repeat([0, 1], 4)
  global_slope = np.dot(x - x.mean(), y - y.mean()) / np.sum((x - x.mean())**2)
  fit_x = x.copy()
  for group in np.unique(groups):
    fit_x[groups == group] -= fit_x[groups == group].mean()
  within_sample_slope = np.dot(fit_x, y) / np.dot(fit_x, fit_x)
  assert global_slope == pytest.approx(-1)
  assert within_sample_slope == pytest.approx(0)


def test_bpcells_uses_public_custom_operator_not_deprecated_mult():
  script = (ROOT / 'scripts/bpcells_pca.R').read_text()
  assert 'setClass(\n  "RidgeResidualMatrix"' in script
  assert 'signature(x = "RidgeResidualMatrix", y = "numeric")' in script
  assert 'signature(x = "numeric", y = "RidgeResidualMatrix")' in script
  assert 'fit_design = "matrix"' in script
  assert 'correction_design = "matrix"' in script
  assert 'crossprod(x@fit_design, base_product)' in script
  assert 'mult =' not in script


def test_human_identifier_mapping_uses_declared_rowdata_not_annotation_downloads():
  script = (ROOT / 'scripts' / 'tricycle.R').read_text()
  environment = (ROOT / 'envs' / 'tricycle.yaml').read_text()
  assert 'rowdata_f' in script
  assert 'ensembl_id' in script
  assert 'org.Hs.eg.db' not in script
  assert 'bioconductor-org.hs.eg.db' not in environment
  assert 'bioconductor-org.mm.eg.db' not in environment
  assert 'library(rhdf5)' in script
  assert 'bioconductor-rhdf5' in environment
  assert 'r-r.utils' in environment


def test_cell_cycle_outputs_use_dedicated_directory():
  main_rules = (ROOT / 'rules' / 'scprocess.smk').read_text()
  tricycle_rules = (ROOT / 'rules' / 'tricycle.smk').read_text()
  integration_rules = (ROOT / 'rules' / 'integration.smk').read_text()
  assert 'cell_cycle_dir = f"{PROJ_DIR}/output/{SHORT_TAG}_cell_cycle"' in main_rules
  assert "f'{cell_cycle_dir}/tricycle_scores_" in tricycle_rules
  assert "f'{cell_cycle_dir}/run_{{run}}_" in tricycle_rules
  assert "{cell_cycle_dir}/tricycle/" not in tricycle_rules
  assert "f'{cell_cycle_dir}/final_cell_cycle_regression_" in integration_rules


def test_cell_cycle_has_standalone_diagnostics_report():
  cli = (ROOT / 'scprocess').read_text()
  main_rules = (ROOT / 'rules' / 'scprocess.smk').read_text()
  report_rules = (ROOT / 'rules' / 'render_htmls.smk').read_text()
  renderer = (ROOT / 'scripts' / 'render_htmls.R').read_text()
  template = (ROOT / 'resources/rmd_templates/cell_cycle.Rmd.template').read_text()
  helper = (ROOT / 'scripts/cell_cycle.R').read_text()
  assert 'rule render_html_cell_cycle:' in report_rules
  assert 'rule cell_cycle:' in main_rules
  assert '"hvg", "cell_cycle", "integration"' in cli
  assert "rule_name='cell_cycle'" in report_rules
  assert "f'{docs_dir}/{SHORT_TAG}_cell_cycle.html'" in main_rules
  assert "'cell_cycle', 'integration'" in renderer
  assert 'cell_cycle_link' in renderer
  assert '${cell_cycle_link}' in (ROOT / 'resources/rmd_templates/index.Rmd.template').read_text()
  integration_target = main_rules.split('rule integration:', 1)[1].split('rule marker_genes:', 1)[0]
  assert 'cell_cycle_outs' not in integration_target
  assert 'Pooled projection and estimated origin' in template
  assert 'Observed expression around theta' in template
  assert 'Sample-level diagnostics' in template
  assert 'Approximate phase proportions by sample' in template
  assert 'read_cell_cycle_marker_expression <- function' in helper
  assert 'smooth_periodic_expression <- function' in helper


def test_cell_cycle_target_requires_config_block():
  from scripts.scprocess_utils import check_config_ok_for_rule

  with pytest.raises(KeyError, match="no 'cell_cycle' section"):
    check_config_ok_for_rule({}, 'cell_cycle')
  check_config_ok_for_rule({'cell_cycle': {}}, 'cell_cycle')


def test_bonus_doublet_rule_is_gpu_routed_on_shpc():
  profile = (ROOT / 'profiles/slurm_shpc/config.yaml').read_text()
  assert (
    '  classify_additional_doublets:\n'
    '    slurm_partition: "batch_gpu"\n'
    '    slurm_extra: "\'--gpus=1\'"\n'
    "    qos: '3h'"
  ) in profile


def test_tricycle_rule_runs_once_per_physical_run():
  rules = (ROOT / 'rules' / 'tricycle.smk').read_text()
  assert "Rscript --vanilla -e 'source(\"{scprocess_dir}/scripts/tricycle.R\")" in rules
  assert 'rule estimate_run_tricycle:' in rules
  assert 'run_column = "{params.run_var}"' in rules
  assert 'run_value = "{wildcards.run}"' in rules
  assert 'estimate_sample_tricycle' not in rules
  assert 'estimate_unassigned_doublet_tricycle' not in rules
  assert 'Rscript -e "source(' not in rules


def test_tricycle_reports_incremental_progress():
  script = (ROOT / 'scripts' / 'tricycle.R').read_text()
  assert '.tricycle_log <- function' in script
  assert 'Reading CSC arrays from ' in script
  assert 'Materializing ' in script
  assert 'Calculating uncentred tricycle reference projection' in script
  assert 'Finished ' in script


def test_tricycle_run_jobs_have_dedicated_runtime():
  rules = (ROOT / 'rules' / 'tricycle.smk').read_text()
  assert rules.count('runtime = 60') == 1


def test_projection_is_centered_once_at_project_level():
  r_script = (ROOT / 'scripts' / 'tricycle.R').read_text()
  py_script = (ROOT / 'scripts' / 'tricycle.py').read_text()
  assert 'gene_means <- Matrix::rowMeans(normalized)' not in r_script
  assert 'tricycle_raw_pc1' in r_script
  assert 'projection_center = raw_points.mean(axis=0)' in py_script
  assert 'pooled_project_mean' in py_script


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
