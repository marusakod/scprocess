"""Grouped splits, seeded sampling and read-only count loading."""
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
from train_xgboost import _balanced_sample, _read_selected_h5ad, _best_iteration_model
from train_xgboost import combine_predictions


def frame():
  rows = []
  for donor in ['d1', 'd2', 'd3']:
    for tube in range(5):
      for label in ['a', 'b']:
        for cell in range(10):
          rows.append(dict(cell_id = f'{donor}_{tube}_{label}_{cell}', label = label,
            donor = donor, tube = f'{donor}_t{tube}', project = 'p', batch = f'{donor}_{tube}'))
  return pl.DataFrame(rows)


def cfg():
  return dict(test_group_var = 'donor', test_groups = ['d3'], validation_group_var = 'tube',
    validation_fraction = .4, validation_variable_combinations = ['donor', 'project'],
    min_cells_per_label = 5, excluded_labels = [], train_cells_per_label = 15,
    validation_cells_per_label = 8, sampling_variable_combinations = ['donor', 'tube'], seed = 7,
    cluster_majority_vote_resolution = None, cluster_majority_vote_min_purity = .9)


def assign(data, config):
  keys = ['test_group_var', 'test_groups', 'validation_group_var', 'validation_fraction',
    'validation_variable_combinations', 'excluded_labels', 'min_cells_per_label']
  return training._assign_splits(data, np.random.default_rng(config['seed']), **{key: config[key] for key in keys})


def training_metadata_fixture(tmp_path, data, config):
  annotations = tmp_path / 'annotations.csv'
  data.write_csv(annotations)
  integration = tmp_path / 'integration.csv'
  data.select('cell_id', pl.col('batch').alias('sample_id'), pl.lit('cl1').alias('RNA_snn_res.2')).write_csv(integration)
  manifest = tmp_path / 'manifest.yaml'
  manifest.write_text(yaml.safe_dump({batch: 'unused.h5ad' for batch in data['batch'].unique()}))
  return training.prepare_training_metadata(annotations, integration, manifest, 'sample_id',
    tmp_path / 'training_metadata.csv.gz', tmp_path / 'coverage.csv', tmp_path / 'counts.csv', is_join = False, **config)


def test_donor_and_tube_disjoint():
  result, _, _ = assign(frame(), cfg())
  test = result.filter(pl.col('split') == 'test')
  train = result.filter(pl.col('split') == 'train')
  val = result.filter(pl.col('split') == 'validation')
  assert set(test['donor']) == {'d3'}
  assert 'd3' not in set(train['donor']) | set(val['donor'])
  assert not set(train['tube']) & set(val['tube'])
  assert set(result['split']) == {'train', 'validation', 'test'}


def test_group_spanning_combinations_fails():
  data = frame().with_columns(
    pl.when(pl.col('donor').is_in(['d1', 'd2'])).then(pl.lit('shared')).otherwise(pl.col('tube')).alias('tube'))
  with pytest.raises(ValueError, match = 'spans'):
    assign(data, cfg())


def test_unknown_test_group_fails():
  config = cfg()
  config['test_groups'] = ['absent']
  with pytest.raises(ValueError, match = 'absent'):
    assign(frame(), config)


def test_joint_combination_sampling():
  data = pl.DataFrame({'cell': range(105), 'group': ['large'] * 100 + ['rare1'] * 3 + ['rare2'] * 2})
  sample = _balanced_sample(data, 9, ['group'], np.random.default_rng(8))
  assert len(sample) == 9 and sample['cell'].n_unique() == 9
  assert dict(sample.group_by('group').len().iter_rows()) == {'large': 4, 'rare1': 3, 'rare2': 2}


def test_prepare_training_metadata_and_no_default_smoothing(tmp_path):
  config = cfg()
  config['excluded_labels'] = ['b']
  result = training_metadata_fixture(tmp_path, frame(), config)
  assert set(result.filter(pl.col('selected'))['label']) == {'a'}
  assert set(result.filter(pl.col('label') == 'b')['split']) == {'excluded'}
  assert result['label'].equals(result['original_label'])
  assert result.filter(pl.col('selected') & (pl.col('split') == 'train')).height == 15
  assert result.filter(pl.col('selected') & (pl.col('split') == 'validation')).height == 8


def test_schema_renames_and_test_group_pair():
  from jsonschema import Draft202012Validator
  for context in ['config', 'join', 'zoom']:
    schema = json.loads((ROOT / f'resources/schemas/{context}.schema.json').read_text())['properties']['train_xgboost']
    settings = cfg()
    settings.pop('cluster_majority_vote_resolution')
    good = dict(annots_f = 'annotations.csv', xgboost_tag = 'model', **settings)
    Draft202012Validator(schema).validate(good)
    assert schema['properties']['validation_group_var']['default'] == 'batch'
    for bad in [dict(annots_f = 'a', ref_tag = 'old'), dict(annots_f = 'a', refine_labels = True),
      dict(annots_f = 'a', test_group_var = 'donor')]:
      assert list(Draft202012Validator(schema).iter_errors(bad))


def test_backed_rows_align_and_preserve_bytes(tmp_path):
  import anndata as ad
  import scipy.sparse as sp
  paths = []
  for name, genes, raw in [('a', ['g1', 'g2'], [[1, 2], [3, 4]]), ('b', ['g2', 'g1'], [[6, 5], [8, 7]])]:
    path = tmp_path / f'{name}.h5ad'
    data = ad.AnnData(sp.csr_matrix(raw))
    data.obs_names = ['c1', 'c2']
    data.var_names = genes
    data.write_h5ad(path)
    paths.append(path)
  before = [(p.read_bytes(), p.stat().st_mtime_ns) for p in paths]
  raw, genes, _ = _read_selected_h5ad(paths[0], ['c2', 'c1'])
  assert raw.toarray().tolist() == [[3, 4], [1, 2]]
  with pytest.warns(UserWarning, match = 'different order'):
    raw, _, _ = _read_selected_h5ad(paths[1], ['c1'], genes)
  assert raw.toarray().tolist() == [[5, 6]]
  assert before == [(p.read_bytes(), p.stat().st_mtime_ns) for p in paths]
  with pytest.raises(ValueError, match = 'Missing'):
    _read_selected_h5ad(paths[0], ['absent'])
  with pytest.raises(ValueError, match = 'Different'):
    _read_selected_h5ad(paths[1], ['c1'], ['g1', 'absent'])


def test_best_iteration_export():
  import xgboost as xgb
  x = np.arange(80, dtype = float).reshape(40, 2)
  y = (x[:, 0] > 30).astype(int)
  model = xgb.train({'objective': 'binary:logistic', 'nthread': 1, 'seed': 7},
    xgb.DMatrix(x, label = y), num_boost_round = 6)
  model.set_attr(best_iteration = '2')
  trimmed = _best_iteration_model(model)
  assert trimmed.num_boosted_rounds() == 3
  np.testing.assert_array_equal(trimmed.predict(xgb.DMatrix(x)), model.predict(xgb.DMatrix(x), iteration_range = (0, 3)))


def test_aggregation_preserves_test_and_counts_unsupported(tmp_path):
  plan = tmp_path / 'training_metadata.csv.gz'
  pl.DataFrame(dict(cell_id = ['a', 'b', 'c'], label = ['x', 'unseen', 'x'], original_label = ['x', 'unseen', 'x'],
    split = ['train', 'test', 'validation'], selected = [True, False, False],
    label_supported_in_training = [True, False, True])).write_csv(plan, compression = 'gzip')
  pred = tmp_path / 'pred.csv.gz'
  pl.DataFrame(dict(cell_id = ['a', 'b', 'c', 'outside'], predicted_label_naive = ['x'] * 4,
    probability_naive = [.9, .8, .9, .7])).write_csv(pred, compression = 'gzip')
  integration = tmp_path / 'integration.csv'
  pl.DataFrame({'cell_id': ['a', 'b', 'c']}).write_csv(integration)
  output = tmp_path / 'demo_fulldata_predictions.csv.gz'
  combine_predictions([pred], plan, output, integration)
  result = pl.read_csv(output)
  assert set(result['cell_id']) == {'a', 'b', 'c'}
  assert result.filter(pl.col('cell_id') == 'b')['split'][0] == 'test'
  metrics = pl.read_csv(tmp_path / 'demo_evaluation.csv')
  test = metrics.filter(pl.col('split') == 'test')
  assert test['accuracy'][0] == 0 and test['macro_f1'][0] == 0
  assert 'development_holdout' in metrics['split']
  assert 'coarse_predicted' not in result.columns



def test_validation_only_genes_do_not_pass_detection_filter():
  import scipy.sparse as sp
  raw = sp.csr_matrix([[1, 0, 1], [2, 0, 0], [3, 0, 0], [0, 4, 0], [0, 5, 0]])
  normalized_counts = raw.astype(float) * .5
  filtered_norm, filtered_raw, gene_names = training._filter_training_genes(
    normalized_counts, raw,
    pl.Series('gene_id', ['g1', 'g2', 'g3']),
    pl.Series('gene_type', ['protein_coding'] * 3),
    np.array([True, True, True, False, False]),
    gene_exclude_re = None, min_cells_expressed = 2)
  assert gene_names == ['g1']
  np.testing.assert_array_equal(filtered_raw.toarray(), raw[:, [0]].toarray())
  np.testing.assert_array_equal(filtered_norm.toarray(), normalized_counts[:, [0]].toarray())


def test_majority_votes_use_training_labels_only(tmp_path):
  data = frame().with_columns(
    pl.when((pl.col('donor') == 'd3') | pl.col('cell_id').str.ends_with('_b_9'))
      .then(pl.lit('b')).otherwise(pl.lit('a')).alias('label'))
  config = cfg()
  config.update(cluster_majority_vote_resolution = 2, cluster_majority_vote_min_purity = .9)
  result = training_metadata_fixture(tmp_path, data, config)
  assert set(result.filter(pl.col('split') == 'train')['label']) == {'a'}
  assert set(result.filter(pl.col('split') == 'test')['label']) == {'b'}
  assert set(result.filter(pl.col('split') == 'test')['original_label']) == {'b'}


def test_gpu_request_rejects_cpu_build(monkeypatch):
  import xgboost as xgb
  monkeypatch.setattr(xgb, 'build_info', lambda: {'USE_CUDA': False})
  with pytest.raises(RuntimeError, match = 'without CUDA'):
    training._verify_gpu(2)


@pytest.mark.parametrize('batch_var', ['sample_id', 'pool_id'])
@pytest.mark.parametrize('is_join,qualified', [(False, False), (True, False), (True, True)])
def test_planning_context_is_independent_of_batch_column(tmp_path, batch_var, is_join, qualified):
  data = frame()
  if is_join:
    parts = []
    for project_id in ['p1', 'p2']:
      parts.append(data.with_columns(
        pl.lit(project_id).alias('project_id'),
        (pl.lit(project_id + '_') + pl.col('cell_id')).alias('cell_id'),
        (pl.lit(project_id + '_') + pl.col('tube')).alias('tube')))
    data = pl.concat(parts)

  annotations = tmp_path / 'annotations.csv'
  data.select('cell_id', 'label', 'donor', 'tube', 'project').write_csv(annotations)
  batch = pl.col('batch')
  if qualified:
    batch = pl.col('project_id') + '_' + batch
  columns = ['cell_id', batch.alias(batch_var)]
  if is_join:
    columns.append('project_id')
  integration = tmp_path / 'integration.csv'
  data.select(columns).write_csv(integration)

  keys = data['batch']
  if is_join:
    keys = data['project_id'] + '_' + keys
  manifest = tmp_path / 'manifest.yaml'
  manifest.write_text(yaml.safe_dump({key: 'unused.h5ad' for key in keys.unique()}))
  training_metadata_f = tmp_path / 'training_metadata.csv.gz'
  result = training.prepare_training_metadata(annotations, integration, manifest, batch_var,
    training_metadata_f, tmp_path / 'coverage.csv', tmp_path / 'counts.csv', is_join = is_join, **cfg())

  assert result['cell_id'].to_list() == data['cell_id'].to_list()
  assert result['batch'].to_list() == keys.to_list()
  assert json.loads(Path(str(training_metadata_f) + '.json').read_text())['configuration']['is_join'] == is_join

  if is_join and not qualified:
    # The same raw IDs must fail when incorrectly declared as a single project.
    with pytest.raises(ValueError, match = 'batch identifiers'):
      training.prepare_training_metadata(annotations, integration, manifest, batch_var,
        training_metadata_f, tmp_path / 'coverage.csv', tmp_path / 'counts.csv', is_join = False, **cfg())


def test_feature_selection_rejects_no_positive_gains():
  from train_xgboost import _select_features_by_gain

  class NoInformativeFeatures:
    def get_score(self, importance_type):
      return {}

  with pytest.raises(ValueError, match = 'no genes with positive feature gain'):
    _select_features_by_gain(NoInformativeFeatures(), ['g1', 'g2'], .9, 1, 2)


@pytest.mark.parametrize('a_cells,b_cells,expected_unit', [(10, 10, 'a'), (8, 10, 'b')])
def test_training_support_prefers_largest_unit_then_unit_id(a_cells, b_cells, expected_unit):
  # Choose the two large units for validation so training needs one moved back.
  class FixedValidationSelection:
    def choice(self, values, n, replace):
      assert n == 2 and not replace
      return np.array(['b', 'a'])

  cells_df = pl.DataFrame({
    'tube': ['a'] * a_cells + ['b'] * b_cells + ['c', 'd'],
    'label': ['rare'] * (a_cells + b_cells + 2)})
  result_df, _, moved_groups = training._assign_splits(cells_df, FixedValidationSelection(),
    test_group_var = None, test_groups = [], validation_group_var = 'tube',
    validation_fraction = .5, validation_variable_combinations = [],
    excluded_labels = [], min_cells_per_label = 10)
  assert moved_groups == [{'label': 'rare', 'group': expected_unit}]
  assert result_df.filter(pl.col('tube') == expected_unit)['split'].unique().to_list() == ['train']


@pytest.mark.parametrize('pattern,expected_indices', [
  (None, [0, 1, 2, 3]),
  ('', [0, 1, 2, 3]),
  ('lncRNA|pseudogene', [0, 3]),
  ('^protein_coding$', [1, 2, 3]),
  ('absent', [0, 1, 2, 3]),
  ('.*', [])])
def test_gene_type_filter_preserves_matrix_and_metadata_alignment(pattern, expected_indices):
  import scipy.sparse as sp

  gene_names = pl.Series('gene_id', ['g3', 'g1', 'g4', 'g2'])
  gene_types = pl.Series('gene_type', ['protein_coding', 'lncRNA', 'processed_pseudogene', 'unknown'])
  raw_counts = sp.csr_matrix(np.arange(8).reshape(2, 4))
  normalized_counts = raw_counts.astype(float) * .5
  if not expected_indices:
    with pytest.raises(ValueError, match = 'No genes pass training-only filtering'):
      training._filter_training_genes(normalized_counts, raw_counts, gene_names, gene_types,
        np.array([True, True]), gene_exclude_re = pattern, min_cells_expressed = 1)
    return

  filtered_norm, filtered_raw, filtered_names = training._filter_training_genes(
    normalized_counts, raw_counts, gene_names, gene_types, np.array([True, True]),
    gene_exclude_re = pattern, min_cells_expressed = 1)

  np.testing.assert_array_equal(filtered_raw.toarray(), raw_counts[:, expected_indices].toarray())
  np.testing.assert_array_equal(filtered_norm.toarray(), normalized_counts[:, expected_indices].toarray())
  assert filtered_names == [gene_names[index] for index in expected_indices]


def test_two_pass_fit_and_separate_coarse_specifications(tmp_path):
  import anndata as ad
  import scipy.sparse as sp

  # A small synthetic fixture checks the complete fitting/export interface.
  rng = np.random.default_rng(31)
  cell_ids = [f'c{i}' for i in range(48)]
  labels = ['a' if i % 2 == 0 else 'b' for i in range(48)]
  means = np.full((48, 12), 2.)
  for i, label in enumerate(labels):
    means[i, :6] += 20 if label == 'a' else 0
    means[i, 6:] += 20 if label == 'b' else 0
  data = ad.AnnData(sp.csr_matrix(rng.poisson(means)))
  data.obs_names = cell_ids
  data.var_names = [f'g{i}' for i in range(12)]
  data.var['gene_type'] = 'protein_coding'
  h5ad_f = tmp_path / 'sample.h5ad'
  data.write_h5ad(h5ad_f)
  h5ads_f = tmp_path / 'h5ads.yaml'
  h5ads_f.write_text(yaml.safe_dump({'sample': str(h5ad_f)}))
  cells_df = pl.DataFrame(dict(cell_id=cell_ids, batch=['sample'] * 48, label=labels,
    split=['train'] * 32 + ['validation'] * 16, selected=[True] * 48))
  metadata_f = tmp_path / 'training_metadata.csv.gz'
  cells_df.write_csv(metadata_f, compression='gzip')
  model_dir = tmp_path / 'model'
  training.fit_classifier(metadata_f, h5ads_f, model_dir, 1,
    xgboost_tag='tiny', use_gpu=False, gene_exclude_re=None,
    min_cells_expressed=1, seed=31, pass1_subsample=1., pass1_colsample_bytree=.5,
    pass1_learning_rate=.3, pass1_nrounds=8, pass1_early_stopping=3,
    gain_threshold=.9, min_genes=2, max_genes=12, pass2_colsample_bytree=.5,
    pass2_learning_rate=.3, pass2_nrounds=8, pass2_early_stopping=3)
  model_f = model_dir / 'tiny_xgboost_model.json'
  model_before = model_f.read_bytes()
  calls_df = pl.read_csv(model_dir / 'tiny_predictions.csv.gz')
  assert calls_df.height == 48 and 'coarse_predicted' not in calls_df.columns

  # Coarse evaluation uses the specification stage without changing fitted weights.
  calls_df = calls_df.with_columns(pl.lit(True).alias('selected'))
  predictions_f = tmp_path / 'full_predictions.csv.gz'
  calls_df.write_csv(predictions_f, compression='gzip')
  counts_df = cells_df.group_by('label').len()
  counts_f = tmp_path / 'counts.csv'
  counts_df.write_csv(counts_f)
  specs_f = tmp_path / 'label_specs.csv'
  specs_f.write_text('label,coarse_label\na,Neuron\nb,Neuron\n')
  training.write_label_specifications(model_dir / 'tiny_allowed_cls.csv', counts_f,
    predictions_f, model_dir, xgboost_tag='tiny', label_specs_f=specs_f,
    coarse_label_specs_f=None)
  assert model_f.read_bytes() == model_before
  metrics_df = pl.read_csv(model_dir / 'tiny_coarse_evaluation.csv')
  assert set(metrics_df['split']) == {'train', 'validation'}
  assert metrics_df['accuracy'].to_list() == [1., 1.]
