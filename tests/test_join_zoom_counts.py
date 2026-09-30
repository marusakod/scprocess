"""Small numerical checks; run with unittest in the workflow's HVG environment."""
import tempfile
import unittest
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np
import polars as pl
import yaml

try:
    import anndata as ad
    import h5py
    from scipy import sparse
except ImportError as exc:
    raise unittest.SkipTest('Requires the workflow HVG environment') from exc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from join_zoom import batch_statistics, fit_trend, rank_hvgs, merge_moments, write_csv, prepare
from join import build_joint_matrix, build_joint_coldata, _load_batch_hvg_matrix


class TestJoinedZoomCounts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        rng = np.random.default_rng(47)
        self.x = rng.poisson(np.linspace(0.05, 12, 200), size=(80, 200)).astype(np.float64)
        self.cells = pl.DataFrame({'cell_id': [f'cell{i}' for i in range(80)],
                                  'sample_id': ['a_s1']*40 + ['b_s2']*40})
        self.genes = np.array([f'gene{i}' for i in range(200)])
        self.cells_f = self.root / 'cells.csv.gz'
        write_csv(self.cells, self.cells_f)

    def h5ad(self, name, indices, encoding='csr', gene_order=None):
        gene_order = np.arange(200) if gene_order is None else gene_order
        x = self.x[indices][:, gene_order]
        x = sparse.csc_matrix(x) if encoding == 'csc' else sparse.csr_matrix(x) if encoding == 'csr' else x
        adata = ad.AnnData(x)
        adata.obs_names = self.cells['cell_id'].to_numpy()[indices]
        adata.var_names = self.genes[gene_order]
        path = self.root / name
        adata.write_h5ad(path)
        return path

    def test_chunked_moments_and_clipping_match_direct_sparse_calculation(self):
        from hvgs import calculate_estimated_vars, _calculate_regularized_variance
        # Each biological sample is split across CSR and CSC files.
        paths = [self.h5ad('a.h5ad', np.arange(0, 80, 2)),
                 self.h5ad('b.h5ad', np.arange(1, 80, 2), 'csc')]
        stats = [self.root / f'stats{i}.csv.gz' for i in range(2)]
        clips = [self.root / f'clip{i}.csv.gz' for i in range(2)]
        for path, out in zip(paths, stats):
            batch_statistics(path, self.cells_f, 7, out)
        combined = merge_moments(stats)
        for sample, x in [('a_s1', self.x[:40]), ('b_s2', self.x[40:])]:
            actual = combined.filter(pl.col('sample_id') == sample).sort('gene_id')
            order = np.argsort(self.genes)
            np.testing.assert_allclose(actual['counts_sum'], x.sum(axis=0)[order])
            np.testing.assert_allclose(actual['squared_counts_sum'], (x*x).sum(axis=0)[order])
            self.assertEqual(actual['n_cells'].unique().to_list(), [40])
        trend = self.root / 'trend.csv.gz'
        batch_trends = [self.root / f'trend{i}.csv.gz' for i in range(2)]
        fit_trend(stats, trend, batch_trends)
        for path, out, batch_trend in zip(paths, clips, batch_trends):
            batch_statistics(path, self.cells_f, 11, out, batch_trend)
        actual_clips = merge_moments(clips)
        tr = pl.read_csv(trend)
        for sample, x in [('a_s1', self.x[:40]), ('b_s2', self.x[40:])]:
            means, variances = x.mean(axis=0), x.var(axis=0, ddof=1)
            direct_stats = pl.DataFrame({'sample_id': sample, 'gene_id': self.genes,
                'ensembl_id': self.genes, 'n_cells': 40, 'mean': means, 'variance': variances})
            direct_trend = calculate_estimated_vars(None, 'sample', 'sample_id', mean_var_df=direct_stats)
            direct = _calculate_regularized_variance(sparse.csr_matrix(x.T), self.genes, direct_trend).sort('gene_id')
            actual = actual_clips.filter(pl.col('sample_id') == sample).sort('gene_id')
            np.testing.assert_allclose(actual['counts_sum'], direct['counts_sum'], rtol=1e-10)
            np.testing.assert_allclose(actual['squared_counts_sum'], direct['squared_counts_sum'], rtol=1e-10)
        excluded = self.root / 'exclude.csv'
        pl.DataFrame({'gene_id': ['gene0', 'gene1']}).write_csv(excluded)
        hvgs = self.root / 'hvgs.csv.gz'
        rank_hvgs(clips, trend, hvgs, 20, excluded)
        result = pl.read_csv(hvgs)
        self.assertEqual(result['highly_variable'].sum(), 20)
        self.assertFalse(set(result['gene_id']) & {'gene0', 'gene1'})

    def test_matrix_cell_gene_order_csr_csc_and_dense(self):
        selected = {'cell1', 'cell7', 'cell23'}
        genes = ['gene14', 'gene3', 'gene199']
        expected = self.x[[1, 7, 23]][:, [14, 3, 199]].T
        for encoding in ('csr', 'csc', 'dense'):
            with self.subTest(encoding=encoding):
                path = self.h5ad(encoding+'.h5ad', np.arange(40), encoding,
                                 gene_order=np.arange(200)[::-1])
                actual, ids = _load_batch_hvg_matrix(str(path), selected, genes, 'a', 's1')
                self.assertEqual(ids, ['cell1', 'cell7', 'cell23'])
                np.testing.assert_array_equal(actual.toarray(), expected)

    def test_empty_batch_statistics(self):
        path = self.h5ad('empty-selection.h5ad', np.array([0]))
        empty_cells = self.root / 'other.csv'
        self.cells.tail(1).write_csv(empty_cells)
        out = self.root / 'empty.csv.gz'
        batch_statistics(path, empty_cells, 5, out)
        self.assertEqual(pl.read_csv(out).height, 0)

    def test_csc_statistics_do_not_materialize_the_whole_h5ad(self):
        from anndata._core.sparse_dataset import BaseCompressedSparseDataset
        path = self.h5ad('csc-stream.h5ad', np.arange(80), 'csc')
        out = self.root / 'csc-stats.csv.gz'
        with patch.object(BaseCompressedSparseDataset, 'to_memory', side_effect=AssertionError('full read')):
            batch_statistics(path, self.cells_f, 1, out)
        self.assertEqual(pl.read_csv(out).height, 400)

    def test_missing_hvg_is_an_error(self):
        path = self.h5ad('source.h5ad', np.arange(40))
        with self.assertRaisesRegex(ValueError, 'HVGs missing'):
            _load_batch_hvg_matrix(str(path), {'cell1'}, ['absent'], 'a', 's1')

    def test_prepare_and_build_with_only_one_retained_project(self):
        source = self.h5ad('source.h5ad', np.arange(80))
        parent = self.cells.with_columns(pl.lit('a').alias('project_id'),
            pl.lit('cl1').alias('RNA_snn_res.0.5'), pl.lit(2000).alias('total'))
        # Both samples belong to the one retained project in this test.
        parent = parent.with_columns(pl.col('sample_id').str.replace('b_', 'a_'))
        paths = {k: str(self.root / f'{k}.csv.gz') for k in ('integration', 'coldata', 'metadata', 'labels', 'qc')}
        for k in ('integration', 'coldata', 'labels'):
            write_csv(parent, paths[k])
        meta = parent.select('sample_id', 'project_id').unique()
        write_csv(meta, paths['metadata'])
        write_csv(parent.select('cell_id', 'sample_id', 'project_id').with_columns(pl.lit(3.0).alias('log_counts')), paths['qc'])
        proj = {'project_ids': ['a', 'b'],
                'integrated_fs': [str(self.root / f'{i}.cells.csv.gz') for i in range(2)],
                'sample_meta_fs': [str(self.root / f'{i}.samples.csv') for i in range(2)],
                'qc_all_fs': [str(self.root / f'{i}.qc.csv.gz') for i in range(2)],
                'h5ads_yaml_fs': [str(self.root / f'{i}.yaml') for i in range(2)]}
        ctx = {'paths': paths, 'selection': {'labels_col': 'RNA_snn_res.0.5', 'sel_labels': ['cl1']},
               'qc': {'qc_min_cells': 2}, 'project': proj,
               'batches': {'0': {'path': str(source), 'project_id': 'a', 'key': 'a_batch',
                                  'cells_f': str(self.root / 'batch_cells.csv.gz')}},
               'selected_f': str(self.root / 'selected.csv.gz'), 'summary_f': str(self.root / 'summary.csv'),
               'effective_f': str(self.root / 'effective.yaml')}
        prepare(ctx, {'integration': {'int_n_dims': 2, 'int_batch_var': 'sample_id'}})
        self.assertEqual(pl.read_csv(self.root / 'batch_cells.csv.gz').height, 80)
        self.assertEqual(pl.read_csv(proj['integrated_fs'][1]).height, 0)
        hvgs = self.root / 'hvgs.csv'
        pl.DataFrame({'gene_id': ['gene3', 'gene1']}).write_csv(hvgs)
        counts = self.root / 'matrix.h5'
        build_joint_matrix(hvgs, proj['h5ads_yaml_fs'], ['a', 'b'], proj['integrated_fs'], counts, require_all_cells=True)
        coldata = self.root / 'coldata.csv.gz'
        build_joint_coldata(counts, ['a', 'b'], proj['integrated_fs'], coldata)
        actual = pl.read_csv(coldata)
        self.assertEqual(actual['cell_id'].to_list(), parent['cell_id'].to_list())
        self.assertEqual(actual['sample_id'].to_list(), parent['sample_id'].to_list())
        self.assertEqual(actual.schema['total'], pl.Int64)


if __name__ == '__main__':
    unittest.main()
