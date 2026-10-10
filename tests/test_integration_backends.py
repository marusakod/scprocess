"""Small execution checks of the real CPU and GPU integration backends.

Run unittest discovery in the integration environment with
SCPROCESS_TEST_BACKEND=cpu or gpu to select the allocation's backend.
"""
import os
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd

try:
    import anndata as ad
except ModuleNotFoundError as exc:
    raise unittest.SkipTest('Requires the workflow integration environment') from exc
from scripts import integration


class TestIntegrationBackends(unittest.TestCase):
    def check_backend(self, use_gpu, multiple):
        selected = os.environ.get('SCPROCESS_TEST_BACKEND', 'cpu')
        if selected != ('gpu' if use_gpu else 'cpu'):
            self.skipTest('Backend not selected for this allocation')
        import scanpy.external as sce
        if use_gpu:
            import rapids_singlecell as backend
            import cupy as cp
            # A selected GPU check must fail, rather than skip, without CUDA.
            self.assertGreater(cp.cuda.runtime.getDeviceCount(), 0)
        else:
            import scanpy as backend
            cp = None
        rng = np.random.default_rng(417)
        n_cells, n_dims = 240, 12
        ids = [f'cell{i}' for i in range(n_cells)]
        obs = pd.DataFrame({'cell_id': ids,
                            'sample_id': [f's{i % 4}' for i in range(n_cells)],
                            'project_id': [f'p{(i % 4) // 2}' for i in range(n_cells)]},
                           index=ids)
        data = ad.AnnData(np.zeros((n_cells, 1), dtype=np.float32), obs=obs)
        data.obsm['X_pca'] = rng.normal(size=(n_cells, n_dims)).astype(np.float32)
        batch = ['project_id', 'sample_id'] if multiple else 'sample_id'
        theta = [0.2, 0.1] if multiple else 0.1
        with patch.object(integration, 'sc', backend, create=True), \
             patch.object(integration, 'sce', sce, create=True), \
             patch.object(integration, 'cp', cp, create=True):
            result = integration._do_one_integration(data, batch, 'leiden', n_dims,
                [0.2, 0.5], 'harmony', use_gpu, theta, skip_pca=True, umap_n_iter=100)
        self.assertEqual(result.height, n_cells)
        self.assertEqual(set(result['cell_id']), set(ids))
        self.assertEqual(result['embedding'].unique().to_list(), ['harmony'])
        self.assertTrue(all(c in result.columns for c in ('UMAP1', 'UMAP2',
                           'RNA_snn_res.0.2', 'RNA_snn_res.0.5')))
        pcs = [c for c in result.columns if c.startswith('hmny_pca_')]
        self.assertEqual(len(pcs), n_dims)
        self.assertTrue(np.isfinite(result.select(pcs + ['UMAP1', 'UMAP2']).to_numpy()).all())

    def test_cpu_single_covariate(self):
        self.check_backend(False, False)

    def test_cpu_multiple_covariates(self):
        self.check_backend(False, True)

    def test_gpu_single_covariate(self):
        self.check_backend(True, False)

    def test_gpu_multiple_covariates(self):
        self.check_backend(True, True)
