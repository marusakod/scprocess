import unittest
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd
import polars as pl

try:
  import anndata as ad
except ModuleNotFoundError as exc:
  raise unittest.SkipTest("anndata is only available in the integration environment") from exc

import scripts.integration as integration


class TestIntegrationOutput(unittest.TestCase):
  def test_project_id_batch_variable_is_not_selected_twice(self):
    obs = pd.DataFrame({
      'cell_id': ['cell_1', 'cell_2'],
      'project_id': ['project_a', 'project_b'],
      'sample_id': ['sample_1', 'sample_2'],
      'RNA_snn_res.0.2': ['0', '1'],
    })
    adata = ad.AnnData(X=np.zeros((2, 1)), obs=obs)

    result = integration._get_clusts_from_adata(
      adata, 'harmony', ['project_id', 'sample_id']
    )

    self.assertEqual(len(result.columns), len(set(result.columns)))
    self.assertEqual(result.columns.count('project_id'), 1)
    self.assertLess(result.columns.index('project_id'), result.columns.index('sample_id'))
    self.assertEqual(result['project_id'].to_list(), ['project_a', 'project_b'])

  def test_multiple_batch_variables_route_harmony_to_cpu(self):
    gpu_pca = Mock()
    gpu_pca.get.return_value = np.ones((2, 2), dtype=np.float32)
    adata = SimpleNamespace(obsm={'X_pca': gpu_pca})
    cpu_harmony = Mock()
    gpu_harmony = Mock()
    fake_sce = SimpleNamespace(pp=SimpleNamespace(harmony_integrate=cpu_harmony))
    fake_rsc = SimpleNamespace(pp=SimpleNamespace(harmony_integrate=gpu_harmony))

    with patch.object(integration, 'sce', fake_sce, create=True), \
         patch.object(integration, 'sc', fake_rsc, create=True), \
         self.assertWarnsRegex(UserWarning, 'Harmony will run on CPU'):
      integration._run_harmony(
        adata, ['project_id', 'sample_id'], [1.0, 0.1], use_gpu=True
      )

    gpu_pca.get.assert_called_once_with()
    cpu_harmony.assert_called_once_with(
      adata, key=['project_id', 'sample_id'], theta=[1.0, 0.1]
    )
    gpu_harmony.assert_not_called()
    self.assertIsInstance(adata.obsm['X_pca'], np.ndarray)

  def test_single_batch_variable_keeps_gpu_harmony(self):
    adata = SimpleNamespace(obsm={'X_pca': np.ones((2, 2), dtype=np.float32)})
    gpu_harmony = Mock()
    fake_rsc = SimpleNamespace(pp=SimpleNamespace(harmony_integrate=gpu_harmony))
    fake_cp = SimpleNamespace(float32=np.float32)

    with patch.object(integration, 'sc', fake_rsc, create=True), \
         patch.object(integration, 'cp', fake_cp, create=True):
      integration._run_harmony(adata, 'sample_id', 0.1, use_gpu=True)

    gpu_harmony.assert_called_once_with(
      adata, key='sample_id', max_iter_harmony=5,
      dtype=np.float32, theta=0.1
    )

  def test_precomputed_pca_loader_preserves_declared_order(self):
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp_dir:
      pca_f = Path(tmp_dir) / 'pca.csv.gz'
      pl.DataFrame({
        'cell_id': ['cell_b', 'cell_a'],
        'pca_1': [2.0, 1.0],
        'pca_2': [4.0, 3.0],
      }).write_csv(pca_f)
      cells_df = pl.DataFrame({
        'cell_id': ['cell_b', 'cell_a'],
        'sample_id': ['sample_2', 'sample_1'],
      })

      adata = integration._adata_from_precomputed_pca(pca_f, cells_df)

      self.assertEqual(adata.obs['cell_id'].tolist(), ['cell_b', 'cell_a'])
      np.testing.assert_array_equal(
        adata.obsm['X_pca'], np.array([[2.0, 4.0], [1.0, 3.0]], dtype=np.float32)
      )

  def test_precomputed_pca_loader_rejects_reordered_metadata(self):
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp_dir:
      pca_f = Path(tmp_dir) / 'pca.csv.gz'
      pl.DataFrame({'cell_id': ['a', 'b'], 'pca_1': [1.0, 2.0]}).write_csv(pca_f)
      cells_df = pl.DataFrame({
        'cell_id': ['b', 'a'],
        'sample_id': ['sample_1', 'sample_1'],
      })

      with self.assertRaisesRegex(ValueError, 'orders do not match'):
        integration._adata_from_precomputed_pca(pca_f, cells_df)

  def test_tricycle_scores_are_joined_by_cell_id(self):
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp_dir:
      score_f = Path(tmp_dir) / 'tricycle.csv.gz'
      pl.DataFrame({
        'cell_id': ['cell_b', 'cell_a'],
        'tricycle_pc1': [2.0, 1.0],
        'tricycle_pc2': [4.0, 3.0],
        'tricycle_theta': [0.2, 0.1],
      }).write_csv(score_f)
      integrated = pl.DataFrame({
        'cell_id': ['cell_a', 'cell_b'],
        'UMAP1': [1.0, 2.0],
      })

      result = integration._add_tricycle_scores(integrated, score_f)

      self.assertEqual(result['cell_id'].to_list(), ['cell_a', 'cell_b'])
      self.assertEqual(result['tricycle_pc1'].to_list(), [1.0, 2.0])

  def test_tricycle_scores_require_complete_cell_coverage(self):
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp_dir:
      score_f = Path(tmp_dir) / 'tricycle.csv'
      pl.DataFrame({
        'cell_id': ['cell_a'],
        'tricycle_pc1': [1.0],
        'tricycle_pc2': [3.0],
        'tricycle_theta': [0.1],
      }).write_csv(score_f)
      integrated = pl.DataFrame({'cell_id': ['cell_a', 'cell_b']})

      with self.assertRaisesRegex(ValueError, 'lack tricycle scores'):
        integration._add_tricycle_scores(integrated, score_f)

if __name__ == '__main__':
  unittest.main()
