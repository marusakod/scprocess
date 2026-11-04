"""Model feature matrices for XGBoost prediction."""
from pathlib import Path
import sys

import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from train_xgboost import select_model_features


def _per_column_reference(normalized_counts, gene_ids, selected_genes):
  # The former per-column construction, kept here as the reference result.
  positions = {gene: i for i, gene in enumerate(gene_ids)}
  model_counts = sp.lil_matrix((normalized_counts.shape[0], len(selected_genes)))
  for j, gene in enumerate(selected_genes):
    if gene in positions:
      model_counts[:, j] = normalized_counts[:, positions[gene]]
  return model_counts.tocsr()


def test_features_follow_model_order_with_missing_genes_as_zero():
  normalized_counts = sp.csr_matrix(np.array([[1., 0., 3.], [0., 2., 0.]]))
  model_counts = select_model_features(normalized_counts, ['g1', 'g2', 'g3'], ['g3', 'g9', 'g1'])
  np.testing.assert_array_equal(model_counts.toarray(), [[3., 0., 1.], [0., 0., 0.]])


def test_matches_per_column_construction_and_drops_stored_zeros():
  rng = np.random.default_rng(1)
  normalized_counts = sp.random(50, 40, density=.2, format='csr', random_state=rng)
  normalized_counts.data[::7] = 0.   # stored zeros, as in log1p of stored zero counts
  gene_ids = [f'g{i}' for i in range(40)]
  selected = [f'g{i}' for i in rng.permutation(40)[:25]] + ['absent1', 'absent2']
  model_counts = select_model_features(normalized_counts, gene_ids, selected)
  reference = _per_column_reference(normalized_counts, gene_ids, selected)
  assert (model_counts != reference).nnz == 0
  assert not np.any(model_counts.data == 0)
