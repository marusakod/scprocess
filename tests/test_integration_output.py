import unittest
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd

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

if __name__ == '__main__':
  unittest.main()
