"""Depth augmentation (binomial thinning) for XGBoost training."""
import json
from pathlib import Path
import sys

import numpy as np
import polars as pl
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import train_xgboost as training


FIT_SETTINGS = dict(xgboost_tag = 'tiny', use_gpu = False, gene_exclude_re = None,
  min_cells_expressed = 1, seed = 31, pass1_subsample = 1., pass1_colsample_bytree = .5,
  pass1_learning_rate = .3, pass1_nrounds = 8, pass1_early_stopping = 3,
  gain_threshold = .9, min_genes = 2, max_genes = 12, pass2_colsample_bytree = .5,
  pass2_learning_rate = .3, pass2_nrounds = 8, pass2_early_stopping = 3)


def _counts(totals, n_genes = 20, seed = 3):
  import scipy.sparse as sp

  # one row per cell, with roughly the requested total spread over the genes
  rng = np.random.default_rng(seed)
  rows = [rng.multinomial(total, np.full(n_genes, 1 / n_genes)) for total in totals]
  return sp.csr_matrix(np.array(rows, dtype = float))


def test_thinned_totals_follow_targets_and_never_grow():
  counts_csr  = _counts([20000, 20000, 5000, 800])
  target_umis = np.array([1000., 300., 2500., 400.])
  thinned     = training._thin_counts(counts_csr, target_umis, np.random.default_rng(1))

  totals_before = np.asarray(counts_csr.sum(axis = 1)).ravel()
  totals_after  = np.asarray(thinned.sum(axis = 1)).ravel()
  assert np.all(thinned.toarray() <= counts_csr.toarray())
  assert np.allclose(thinned.data, np.round(thinned.data))
  # binomial totals: within 5 standard deviations of each target
  sd = np.sqrt(target_umis * (1 - target_umis / totals_before))
  assert np.all(np.abs(totals_after - target_umis) < 5 * sd)


def test_augmented_copies_are_reproducible_and_skip_shallow_cells():
  # totals 10,000 / 2,000 / 150 / 0 UMIs; only the first two exceed min_umis = 200
  counts_csr = _counts([10000, 2000, 150, 0])
  settings   = dict(copies = 2, min_umis = 200, max_umis = 3000, seed = 7)
  first, first_rows   = training._make_depth_augmented_counts(counts_csr, **settings)
  second, second_rows = training._make_depth_augmented_counts(counts_csr, **settings)

  assert first_rows.tolist() == [0, 1, 0, 1]
  np.testing.assert_array_equal(first_rows, second_rows)
  np.testing.assert_array_equal(first.toarray(), second.toarray())
  assert first.shape == (4, counts_csr.shape[1])

  # disabled, or no cell deep enough: no copies
  assert training._make_depth_augmented_counts(counts_csr, copies = 0, min_umis = 200,
    max_umis = 3000, seed = 7)[0] is None
  assert training._make_depth_augmented_counts(counts_csr, copies = 1, min_umis = 20000,
    max_umis = 30000, seed = 7)[0] is None


def test_non_integer_counts_are_rejected():
  counts_csr = _counts([1000, 1000])
  counts_csr.data = counts_csr.data + .25
  with pytest.raises(ValueError, match = 'integer raw counts'):
    training._make_depth_augmented_counts(counts_csr, copies = 1, min_umis = 10,
      max_umis = 100, seed = 1)


def _write_fixture(tmp_path):
  import anndata as ad
  import scipy.sparse as sp

  # Two labels with distinct high-expression genes; cells hold ~1,000 UMIs.
  rng      = np.random.default_rng(31)
  cell_ids = [f'c{i}' for i in range(48)]
  labels   = ['a' if i % 2 == 0 else 'b' for i in range(48)]
  means    = np.full((48, 12), 40.)
  for i, label in enumerate(labels):
    means[i, :6] += 100 if label == 'a' else 0
    means[i, 6:] += 100 if label == 'b' else 0
  data = ad.AnnData(sp.csr_matrix(rng.poisson(means).astype(float)))
  data.obs_names = cell_ids
  data.var_names = [f'g{i}' for i in range(12)]
  data.var['gene_type'] = 'protein_coding'
  h5ad_f = tmp_path / 'sample.h5ad'
  data.write_h5ad(h5ad_f)
  h5ads_f = tmp_path / 'h5ads.yaml'
  h5ads_f.write_text(yaml.safe_dump({'sample': str(h5ad_f)}))
  cells_df = pl.DataFrame(dict(cell_id = cell_ids, batch = ['sample'] * 48, label = labels,
    split = ['train'] * 32 + ['validation'] * 16, selected = [True] * 48))
  metadata_f = tmp_path / 'training_metadata.csv.gz'
  cells_df.write_csv(metadata_f, compression = 'gzip')
  return metadata_f, h5ads_f, cells_df


def _fit(tmp_path, name, metadata_f, h5ads_f, **augmentation):
  model_dir = tmp_path / name
  training.fit_classifier(metadata_f, h5ads_f, model_dir, 1, **FIT_SETTINGS, **augmentation)
  return model_dir


def test_disabled_augmentation_matches_no_eligible_cells(tmp_path):
  # Copies requested but no cell deeper than min_umis must give the unaugmented fit.
  metadata_f, h5ads_f, _ = _write_fixture(tmp_path)
  off  = _fit(tmp_path, 'off', metadata_f, h5ads_f, downsample_copies = 0,
    downsample_min_umis = 300, downsample_max_umis = 3000)
  none = _fit(tmp_path, 'none', metadata_f, h5ads_f, downsample_copies = 1,
    downsample_min_umis = 100000, downsample_max_umis = 200000)
  assert (off / 'tiny_xgboost_model.json').read_bytes() == (none / 'tiny_xgboost_model.json').read_bytes()
  card = json.loads((off / 'tiny_model_card.json').read_text())
  assert card['n_train_augmented'] == card['n_validation_augmented'] == 0
  assert card['depth_augmentation'] == 'none'


def test_enabled_outputs_keep_one_row_per_original_cell(tmp_path):
  metadata_f, h5ads_f, cells_df = _write_fixture(tmp_path)
  off = _fit(tmp_path, 'off', metadata_f, h5ads_f, downsample_copies = 0,
    downsample_min_umis = 50, downsample_max_umis = 500)
  on  = _fit(tmp_path, 'on', metadata_f, h5ads_f, downsample_copies = 1,
    downsample_min_umis = 50, downsample_max_umis = 500)

  # predictions: same original cells, same order, no augmented rows
  calls_off = pl.read_csv(off / 'tiny_predictions.csv.gz')
  calls_on  = pl.read_csv(on / 'tiny_predictions.csv.gz')
  assert calls_on['cell_id'].to_list() == calls_off['cell_id'].to_list()
  assert calls_on.height == cells_df.height

  # every cell is deep enough, so each gets one copy in its own split
  card_on  = json.loads((on / 'tiny_model_card.json').read_text())
  card_off = json.loads((off / 'tiny_model_card.json').read_text())
  assert card_on['n_train'] == 32 and card_on['n_validation'] == 16
  assert card_on['n_train_augmented'] == 32 and card_on['n_validation_augmented'] == 16
  assert card_on['configuration']['downsample_copies'] == 1
  assert card_on['depth_augmentation'] != 'none'

  # the detection filter uses original training cells, so the tested genes match
  assert card_on['n_genes_tested'] == card_off['n_genes_tested']


@pytest.mark.parametrize('schema_name', ['config', 'join', 'zoom'])
def test_schemas_define_augmentation_defaults(schema_name):
  schema = json.loads((ROOT / f'resources/schemas/{schema_name}.schema.json').read_text())
  properties = schema['properties']['train_xgboost']['properties']
  assert properties['downsample_copies']['default'] == 0
  assert properties['downsample_min_umis']['default'] == 300
  assert properties['downsample_max_umis']['default'] == 3000


@pytest.mark.parametrize('rule_f', ['train_xgboost.smk', 'join.smk', 'zoom.smk'])
def test_rules_pass_augmentation_settings(rule_f):
  text = (ROOT / 'rules' / rule_f).read_text()
  for name in ['downsample_copies', 'downsample_min_umis', 'downsample_max_umis']:
    assert f'--{name} {{params.{name}:q}}' in ' '.join(text.split())


def test_paired_depth_bounds_are_validated():
  import scprocess_utils

  ok = dict(downsample_copies = 1, downsample_min_umis = 300, downsample_max_umis = 3000)
  scprocess_utils._check_xgboost_depth_augmentation(ok)
  with pytest.raises(ValueError, match = 'below downsample_max_umis'):
    scprocess_utils._check_xgboost_depth_augmentation(dict(ok, downsample_min_umis = 3000))
  # bounds are not checked while augmentation is off
  scprocess_utils._check_xgboost_depth_augmentation(dict(ok, downsample_copies = 0,
    downsample_min_umis = 3000))
