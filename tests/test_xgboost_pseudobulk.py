"""The classifier report pseudobulk keeps only plotted genes, with full library sizes."""
from pathlib import Path
import sys

import numpy as np
import polars as pl
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import train_xgboost as training


def test_pseudobulk_has_one_row_per_label_and_batch(tmp_path):
  import anndata as ad

  # 5 cells over 3 genes; only g0 and g2 are plotted, but library sizes count all
  # genes. Label 'Ependymal / ciliated' contains '/', which must not become an HDF5 key.
  counts   = sp.csr_matrix(np.array([[1, 9, 2], [3, 9, 4], [5, 9, 6], [7, 9, 8], [1, 9, 1]], dtype = float))
  cells_df = pl.DataFrame(dict(label = ['Ependymal / ciliated', 'b', 'Ependymal / ciliated', 'b', 'b'],
    batch = ['s1', 's1', 's1', 's2', 's1']))
  lib_sizes = np.asarray(counts.sum(axis = 1)).ravel()
  pb_f = tmp_path / 'pb.h5ad'
  training._compute_pseudobulk(counts, lib_sizes, cells_df, ['g0', 'g1', 'g2'], ['g2', 'g0']).write_h5ad(pb_f)

  pb = ad.read_h5ad(pb_f)
  assert pb.obs['label'].tolist() == ['Ependymal / ciliated', 'b', 'b']
  assert pb.obs['batch'].tolist() == ['s1', 's1', 's2']
  assert pb.obs['n_cells'].tolist() == [2, 2, 1]
  assert pb.obs['lib_size'].tolist() == [32, 27, 24]
  assert pb.var_names.tolist() == ['g2', 'g0']
  assert np.array_equal(pb.X.toarray(), [[8, 6], [5, 4], [8, 7]])


def test_report_genes_are_top_overall_then_top_per_label(monkeypatch):
  monkeypatch.setattr(training, 'REPORT_TOP_GENES', 2)
  monkeypatch.setattr(training, 'REPORT_TOP_GENES_PER_LABEL', 1)
  # 'z' has overall gain but was not selected for the model, so it is skipped.
  gain_df  = pl.DataFrame(dict(gene_id = ['z', 'a', 'b', 'c'], gain = [9., 5., 4., 3.]))
  label_df = pl.DataFrame(dict(label = ['x', 'x', 'y'], gene_id = ['b', 'c', 'd'], gain = [2., 1., 3.]))
  assert training._get_report_genes(gain_df, label_df, ['a', 'b', 'c', 'd']) == ['a', 'b', 'd']
