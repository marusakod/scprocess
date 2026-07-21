import unittest

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

if __name__ == '__main__':
  unittest.main()
