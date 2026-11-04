"""Prediction chunks: batches grouped by cell count, predicted with one model load."""
from pathlib import Path
import sys

import numpy as np
import polars as pl
import scipy.sparse as sp
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import label_celltypes


def _write_h5ad(path, n_cells, rng, shift = 0):
  import anndata as ad

  # Two groups of cells, high in genes g0-g1 (label a) or g2-g3 (label b)
  means = np.ones((n_cells, 4))
  means[::2, :2] += 20 + shift
  means[1::2, 2:] += 20 + shift
  data = ad.AnnData(sp.csr_matrix(rng.poisson(means).astype(float)))
  data.obs_names = [f'{path.stem}_c{i}' for i in range(n_cells)]
  data.var_names = [f'g{i}' for i in range(4)]
  data.write_h5ad(path)
  return path


def test_chunks_follow_yaml_order_and_cap(tmp_path, monkeypatch):
  monkeypatch.setattr(label_celltypes, 'PREDICTION_CHUNK_MAX_CELLS', 10)
  rng = np.random.default_rng(1)

  # 6 + 4 fit the cap of 10; 3 starts a new chunk; 12 exceeds the cap alone;
  # the empty placeholder joins the current chunk.
  sizes = {'s1': 6, 's2': 4, 's3': 3, 's4': 12}
  lookup = {batch: str(_write_h5ad(tmp_path / f'{batch}.h5ad', n, rng)) for batch, n in sizes.items()}
  (tmp_path / 'empty.h5ad').touch()
  lookup['empty'] = {'path': str(tmp_path / 'empty.h5ad')}
  h5ads_yaml_f = tmp_path / 'h5ads.yaml'
  h5ads_yaml_f.write_text(yaml.safe_dump(lookup, sort_keys = False))

  chunks_f = tmp_path / 'chunks.csv'
  label_celltypes.assign_prediction_chunks(h5ads_yaml_f, chunks_f)
  chunks_df = pl.read_csv(chunks_f)
  assert chunks_df['batch'].to_list() == ['s1', 's2', 's3', 's4', 'empty']
  assert chunks_df['n_cells'].to_list() == [6, 4, 3, 12, 0]
  assert chunks_df['chunk'].to_list() == ['chunk001', 'chunk001', 'chunk002', 'chunk003', 'chunk003']


def test_chunk_predictions_match_single_batch_predictions(tmp_path):
  import xgboost as xgb
  from train_xgboost import normalize_counts

  # Fit a tiny two-label model on log-normalised counts of one batch.
  rng = np.random.default_rng(2)
  train_f = _write_h5ad(tmp_path / 'train.h5ad', 40, rng)
  import anndata as ad
  train = ad.read_h5ad(train_f)
  labels = np.tile([0, 1], 20)
  booster = xgb.train({'objective': 'multi:softprob', 'num_class': 2, 'max_depth': 2},
    xgb.DMatrix(normalize_counts(sp.csr_matrix(train.X)), label = labels,
      feature_names = [f'g{i}' for i in range(4)]), num_boost_round = 5)
  model_f = tmp_path / 'tiny_xgboost_model.json'
  booster.save_model(model_f)
  cls_f = tmp_path / 'tiny_allowed_cls.csv'
  pl.DataFrame({'class': ['a', 'b']}).write_csv(cls_f)
  genes_f = tmp_path / 'tiny_selected_genes.txt'
  genes_f.write_text('g0\ng1\ng2\ng3\n')

  # Two batches and an empty placeholder in one chunk, against one call per batch
  adata_fs = [_write_h5ad(tmp_path / 'b1.h5ad', 7, rng), _write_h5ad(tmp_path / 'b2.h5ad', 5, rng, shift = 5)]
  (tmp_path / 'b3.h5ad').touch()
  adata_fs.append(tmp_path / 'b3.h5ad')
  batches = ['b1', 'b2', 'b3']
  chunk_df = label_celltypes.run_xgboost_batches(batches, 'sample_id', 'tiny', adata_fs, model_f, cls_f, genes_f)
  single_dfs = [label_celltypes.run_xgboost_batches([b], 'sample_id', 'tiny', [f], model_f, cls_f, genes_f)
    for b, f in zip(batches, adata_fs)]
  assert chunk_df.equals(pl.concat(single_dfs))
  assert chunk_df.group_by('sample_id', maintain_order = True).len()['len'].to_list() == [7, 5]
  assert set(chunk_df['predicted_label_naive']) == {'a', 'b'}
