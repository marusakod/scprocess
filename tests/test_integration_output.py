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
  def test_joined_zoom_keeps_parent_harmony_separate_from_new_embedding(self):
    import tempfile
    from pathlib import Path
    from scripts.join_zoom import select_cells

    parent = pl.DataFrame({
      'cell_id': ['a', 'b', 'c'], 'sample_id': ['s1'] * 3,
      'project_id': ['p1'] * 3, 'RNA_snn_res.0.5': ['cl1'] * 3,
      'hmny_pca_01': [10., 11., 12.], 'hmny_pca_02': [20., 21., 22.],
      'hmny_pca_03': [30., 31., 32.], 'UMAP1': [1., 2., 3.],
    })
    selected, _ = select_cells(
      parent, parent, parent.select('sample_id', 'project_id').unique(),
      parent, parent.select('cell_id'),
      {'labels_col': 'RNA_snn_res.0.5', 'sel_labels': ['cl1']},
      {'qc_min_cells': 2}, n_dims=2)

    with tempfile.TemporaryDirectory() as tmp_dir:
      root = Path(tmp_dir)
      pca_f, coldata_f = root / 'pca.csv', root / 'cells.csv'
      pl.DataFrame({'cell_id': ['c', 'a', 'b'],
        'pca_01': [1., 2., 3.], 'pca_02': [4., 5., 6.]}).write_csv(pca_f)
      selected.write_csv(coldata_f)

      for embedding, prefix in [('harmony', 'hmny_pca'), ('pca', 'pca')]:
        with self.subTest(embedding=embedding):
          new_integration = pl.DataFrame({
            'cell_id': ['c', 'a', 'b'], 'sample_id': ['s1'] * 3,
            'project_id': ['p1'] * 3, 'RNA_snn_res.0.5': ['cl2', 'cl1', 'cl2'],
            f'{prefix}_01': [1., 2., 3.], f'{prefix}_02': [4., 5., 6.],
          })
          output_f = root / f'{embedding}.csv.gz'
          with patch.object(integration, '_do_one_integration', return_value=new_integration):
            result = integration.run_zoom_integration(
              hvg_mat_f=None, sample_qc_f=None, coldata_f=coldata_f,
              demux_type='none', exclude_mito=False, embedding=embedding, n_dims=2,
              cl_method='leiden', theta=0.1, res_ls_concat='0.5',
              integration_f=output_f, batch_var='sample_id',
              precomputed_pca_f=pca_f, retain_metadata=True)
          self.assertEqual(result['cell_id'].to_list(), ['c', 'a', 'b'])
          self.assertEqual(result.select(new_integration.columns).to_dicts(),
                           new_integration.to_dicts())
          self.assertNotIn('hmny_pca_03', result.columns)
          for pc in range(1, 4):
            self.assertEqual(result[f'parent_hmny_pca_{pc:02d}'].to_list(),
                             [pc * 10. + 2, pc * 10., pc * 10. + 1])
          self.assertEqual(pl.read_csv(output_f).to_dicts(), result.to_dicts())

  def test_joined_zoom_retains_metadata_in_pca_order(self):
    import tempfile
    from pathlib import Path
    with tempfile.TemporaryDirectory() as tmp_dir:
      root = Path(tmp_dir)
      pca_f, coldata_f, output_f = root / 'pca.csv', root / 'cells.csv', root / 'output.csv.gz'
      pl.DataFrame({'cell_id': ['b', 'a'], 'pca_1': [2., 1.]}).write_csv(pca_f)
      pl.DataFrame({'cell_id': ['a', 'b'], 'sample_id': ['s1', 's2'],
        'project_id': ['p1', 'p2'], 'tissue': ['cortex', 'thalamus'],
        'parent_label': ['cl1', 'cl3']}).write_csv(coldata_f)
      new_clusters = pl.DataFrame({'cell_id': ['b', 'a'], 'sample_id': ['s2', 's1'],
        'project_id': ['p2', 'p1'], 'RNA_snn_res.0.5': ['cl2', 'cl1']})
      with patch.object(integration, '_do_one_integration', return_value=new_clusters):
        result = integration.run_zoom_integration(
          hvg_mat_f=None, sample_qc_f=None, coldata_f=coldata_f,
          demux_type='none', exclude_mito=False, embedding='pca', n_dims=1,
          cl_method='leiden', theta=0.1, res_ls_concat='0.5', integration_f=output_f,
          batch_var='sample_id', precomputed_pca_f=pca_f, retain_metadata=True)
      self.assertEqual(result['cell_id'].to_list(), ['b', 'a'])
      self.assertEqual(result['tissue'].to_list(), ['thalamus', 'cortex'])
      self.assertEqual(result['parent_label'].to_list(), ['cl3', 'cl1'])
      self.assertEqual(result['RNA_snn_res.0.5'].to_list(), ['cl2', 'cl1'])

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

  def test_multiple_batch_variables_keep_gpu_harmony(self):
    gpu_pca = Mock()
    gpu_pca.get.return_value = np.ones((2, 2), dtype=np.float32)
    adata = SimpleNamespace(obsm={'X_pca': gpu_pca})
    cpu_harmony = Mock()
    gpu_harmony = Mock()
    fake_sce = SimpleNamespace(pp=SimpleNamespace(harmony_integrate=cpu_harmony))
    fake_rsc = SimpleNamespace(pp=SimpleNamespace(harmony_integrate=gpu_harmony))

    with patch.object(integration, 'sce', fake_sce, create=True), \
         patch.object(integration, 'sc', fake_rsc, create=True), \
         patch.object(integration, 'cp', SimpleNamespace(float32=np.float32), create=True):
      integration._run_harmony(
        adata, ['project_id', 'sample_id'], [1.0, 0.1], use_gpu=True
      )

    gpu_pca.get.assert_not_called()
    cpu_harmony.assert_not_called()
    gpu_harmony.assert_called_once_with(
      adata, key=['project_id', 'sample_id'], max_iter_harmony=5,
      dtype=np.float32, theta=[1.0, 0.1], flavor='harmony1')

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
      dtype=np.float32, theta=0.1, flavor='harmony1'
    )

  def test_precomputed_pca_loader_preserves_declared_order(self):
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp_dir:
      pca_f = Path(tmp_dir) / 'pca.csv'
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
      pca_f = Path(tmp_dir) / 'pca.csv'
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
      score_f = Path(tmp_dir) / 'tricycle.csv'
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
