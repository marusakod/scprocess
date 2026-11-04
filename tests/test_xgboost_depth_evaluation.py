"""Accuracy at fixed lower depths for every trained classifier."""
import json
from pathlib import Path
import sys

import numpy as np
import polars as pl
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT / 'tests'))
import train_xgboost as training
from test_xgboost_downsampling import FIT_SETTINGS


def _write_fixture(tmp_path, base_mean):
  import anndata as ad
  import scipy.sparse as sp

  # 60 cells of two labels: 32 train and 8 validation selected for fitting, 20 unused.
  rng      = np.random.default_rng(5)
  n_cells  = 60
  cell_ids = [f'c{i}' for i in range(n_cells)]
  labels   = ['a' if i % 2 == 0 else 'b' for i in range(n_cells)]
  means    = np.full((n_cells, 12), float(base_mean))
  for i, label in enumerate(labels):
    means[i, :6] += 3 * base_mean if label == 'a' else 0
    means[i, 6:] += 3 * base_mean if label == 'b' else 0
  data = ad.AnnData(sp.csr_matrix(rng.poisson(means).astype(float)))
  data.obs_names = cell_ids
  data.var_names = [f'g{i}' for i in range(12)]
  data.var['gene_type'] = 'protein_coding'
  h5ad_f = tmp_path / 'sample.h5ad'
  data.write_h5ad(h5ad_f)
  h5ads_f = tmp_path / 'h5ads.yaml'
  h5ads_f.write_text(yaml.safe_dump({'sample': str(h5ad_f)}))
  splits   = ['train'] * 40 + ['validation'] * 20
  selected = [True] * 32 + [False] * 8 + [True] * 8 + [False] * 12
  cells_df = pl.DataFrame(dict(cell_id = cell_ids, batch = ['sample'] * n_cells, label = labels,
    split = splits, selected = selected))
  metadata_f = tmp_path / 'training_metadata.csv.gz'
  cells_df.write_csv(metadata_f, compression = 'gzip')
  return metadata_f, h5ads_f, cells_df


def _fit_and_evaluate(tmp_path, base_mean, name = 'eval'):
  metadata_f, h5ads_f, cells_df = _write_fixture(tmp_path, base_mean)
  model_dir = tmp_path / 'model'
  if not (model_dir / 'tiny_xgboost_model.json').exists():
    training.fit_classifier(metadata_f, h5ads_f, model_dir, 1, **FIT_SETTINGS,
      downsample_copies = 0, downsample_min_umis = 300, downsample_max_umis = 3000)
  output_f = tmp_path / f'{name}_depth_evaluation.csv.gz'
  training.evaluate_depths(metadata_f, h5ads_f, model_dir / 'tiny_xgboost_model.json',
    model_dir / 'tiny_allowed_cls.csv', model_dir / 'tiny_selected_genes.txt', output_f, 1, seed = 31)
  return pl.read_csv(output_f, schema_overrides = {'depth': pl.String}), cells_df


def test_unused_deep_cells_are_evaluated_at_every_depth(tmp_path):
  # ~2,400 UMIs per cell: all deeper than DEPTH_EVALUATION_MIN_UMIS
  depth_df, cells_df = _fit_and_evaluate(tmp_path, base_mean = 100)
  unused = set(cells_df.filter(~pl.col('selected'))['cell_id'])
  depths = [str(d) for d in training.DEPTH_EVALUATION_UMIS] + ['full']

  assert set(depth_df['cell_id']) == unused
  assert depth_df.group_by('depth').len().sort('depth')['len'].to_list() == [len(unused)] * len(depths)
  assert set(depth_df['depth']) == set(depths)
  assert set(depth_df['cell_group']) == {'remaining'}

  # thinned totals sit near their targets and never exceed the original
  full = dict(depth_df.filter(pl.col('depth') == 'full').select('cell_id', 'n_umis').iter_rows())
  for depth in training.DEPTH_EVALUATION_UMIS:
    thinned = depth_df.filter(pl.col('depth') == str(depth))
    assert thinned['n_umis'].mean() == pytest.approx(depth, rel = .15)
    assert all(row['n_umis'] <= full[row['cell_id']] for row in thinned.iter_rows(named = True))


def test_depth_evaluation_is_reproducible(tmp_path):
  first, _  = _fit_and_evaluate(tmp_path, base_mean = 100, name = 'first')
  second, _ = _fit_and_evaluate(tmp_path, base_mean = 100, name = 'second')
  assert first.equals(second)


def test_shallow_datasets_give_an_empty_table(tmp_path):
  # ~480 UMIs per cell: no cell is deep enough, so the table has only its header
  depth_df, _ = _fit_and_evaluate(tmp_path, base_mean = 20)
  assert depth_df.is_empty()
  assert depth_df.columns == ['cell_id', 'label', 'cell_group', 'depth', 'n_umis',
    'predicted_label', 'probability']


@pytest.mark.parametrize('rule_f,rule', [('train_xgboost.smk', 'train_xgboost_depth_evaluation'),
  ('join.smk', 'join_train_xgboost_depth_evaluation'), ('zoom.smk', 'zoom_train_xgboost_depth_evaluation')])
def test_rules_evaluate_depth_and_feed_the_report(rule_f, rule):
  text = ' '.join((ROOT / 'rules' / rule_f).read_text().split())
  assert f'rule {rule}:' in text and 'train_xgboost.py depth' in text
  report_rules = {'train_xgboost.smk': 'render_htmls.smk', 'join.smk': 'join.smk', 'zoom.smk': 'zoom.smk'}
  report_text = ' '.join((ROOT / 'rules' / report_rules[rule_f]).read_text().split())
  assert "depth_evaluation_f = '{input.depth_evaluation_f}'" in report_text
