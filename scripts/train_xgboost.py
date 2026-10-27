# Train an XGBoost classifier to predict cell type labels from scRNA-seq data.

import argparse
import gc
import hashlib
import json
from pathlib import Path
import re
import shutil
import tarfile
import tempfile
import urllib.request
import numpy as np
import polars as pl
import yaml


# Training plans

def prepare_plan(annots_f, integration_f, h5ads_yaml, batch_var, plan_f, coverage_f, counts_f,
  *, is_join, excluded_labels, label_map_f, cluster_majority_vote_resolution,
  cluster_majority_vote_min_purity, test_group_var, test_groups, validation_group_var,
  validation_fraction, validation_variable_combinations, sampling_variable_combinations,
  train_cells_per_type, validation_cells_per_type, min_cells_per_type, seed):

  # Read cell/batch IDs, join project IDs and optional voting clusters from integration_f.
  # Requested grouping/sampling metadata comes from annotations when available,
  # otherwise from integration_f.
  vote_res = cluster_majority_vote_resolution
  vote_col = f'RNA_snn_res.{float(vote_res):g}' if vote_res is not None else None
  needed   = set(validation_variable_combinations + sampling_variable_combinations)
  needed |= {v for v in [test_group_var, validation_group_var] if v}

  # Join labels and integration metadata, preserving integration cell order.
  plan_df = _load_plan_metadata(annots_f, integration_f, batch_var, is_join, vote_col, needed)
  plan_df = _resolve_plan_batches(plan_df, h5ads_yaml, batch_var, is_join, needed)

  # Split whole groups before exclusions, optional voting and cell selection.
  rng = np.random.default_rng(seed)
  plan_df, labels_with_minimum_development_cells, validation_groups_moved_to_training = _assign_splits(plan_df, rng,
    test_group_var                   = test_group_var,
    test_groups                      = test_groups,
    validation_group_var             = validation_group_var,
    validation_fraction              = validation_fraction,
    validation_variable_combinations = validation_variable_combinations,
    excluded_labels                  = excluded_labels,
    min_cells_per_type               = min_cells_per_type)

  plan_df = _prepare_plan_labels(plan_df, labels_with_minimum_development_cells, excluded_labels, vote_col,
    cluster_majority_vote_min_purity)
  plan_df = _select_plan_cells(plan_df, rng, sampling_variable_combinations,
    train_cells_per_type, validation_cells_per_type)

  # Record the explicit planning settings alongside the selected-cell tables.
  configuration = dict(
    is_join                          = bool(is_join),
    excluded_labels                  = excluded_labels,
    label_map_f                      = label_map_f,
    cluster_majority_vote_resolution = cluster_majority_vote_resolution,
    cluster_majority_vote_min_purity = cluster_majority_vote_min_purity,
    test_group_var                   = test_group_var,
    test_groups                      = test_groups,
    validation_group_var             = validation_group_var,
    validation_fraction              = validation_fraction,
    validation_variable_combinations = validation_variable_combinations,
    sampling_variable_combinations   = sampling_variable_combinations,
    train_cells_per_type             = train_cells_per_type,
    validation_cells_per_type        = validation_cells_per_type,
    min_cells_per_type               = min_cells_per_type,
    seed                             = seed)
  _write_training_plan(plan_df, plan_f, coverage_f, counts_f, configuration, validation_groups_moved_to_training)
  return plan_df


def _load_plan_metadata(annots_f, integration_f, batch_var, is_join, vote_col, needed):
  # Read labels and only the integration columns required for planning.
  integration_columns = ['cell_id', batch_var]
  if is_join:
    integration_columns.append('project_id')
  if vote_col:
    integration_columns.append(vote_col)
  annotations_df = pl.read_csv(annots_f)
  # Accept one standard or legacy label column, then use the canonical name.
  label_columns = [name for name in ['label', 'annotation', 'fine_label'] if name in annotations_df.columns]
  if len(label_columns) != 1:
    raise ValueError('annots_f needs one label column (label; legacy annotation/fine_label accepted)')
  annotations_df = annotations_df.rename({label_columns[0]: 'label'})
  if _missing_values(annotations_df, ['cell_id']) or annotations_df['cell_id'].n_unique() != len(annotations_df):
    raise ValueError('Annotations must have unique nonmissing cell_id values')
  # Grouping metadata may come from either input table.
  header = pl.read_csv(integration_f, n_rows = 0).columns
  integration_columns += [
    v for v in sorted(needed)
    if v in header and v not in annotations_df.columns and v not in integration_columns]
  int_df = pl.read_csv(integration_f, columns = list(dict.fromkeys(integration_columns)))
  if _missing_values(int_df, ['cell_id']) or int_df['cell_id'].n_unique() != len(int_df):
    raise ValueError('Integration must have unique nonmissing cell_id values')

  # Reject conflicting shared metadata before dropping duplicate columns.
  overlap = [v for v in annotations_df.columns if v != 'cell_id' and v in int_df.columns]
  if overlap:
    comparison_df = int_df.join(
      annotations_df.select(['cell_id'] + overlap),
      on     = 'cell_id',
      suffix = '_annotation')
    for var in overlap:
      conflict = (pl.col(var).is_not_null() & pl.col(var + '_annotation').is_not_null()
        & (pl.col(var).cast(pl.String) != pl.col(var + '_annotation').cast(pl.String)))
      if comparison_df.filter(conflict).height:
        raise ValueError(f'Conflicting integration and annotation values: {var}')
    annotations_df = annotations_df.drop(overlap)

  # Keep integration order; the matrix loader later groups selected rows by batch.
  plan_df = int_df.join(
    annotations_df,
    on             = 'cell_id',
    how            = 'inner',
    maintain_order = 'left')
  if plan_df.is_empty():
    raise ValueError('No overlapping original labels and integration cells')

  return plan_df


def _resolve_plan_batches(plan_df, h5ads_yaml, batch_var, is_join, needed):
  # Joined integrations may already have project-qualified batch identifiers.
  h5ad_lookup     = yaml.safe_load(Path(h5ads_yaml).read_text())
  batch_col       = plan_df[batch_var].cast(pl.String)
  missing_batches = set(batch_col) - set(h5ad_lookup)
  if is_join and missing_batches:
    batch_col       = plan_df['project_id'].cast(pl.String) + '_' + batch_col
    missing_batches = set(batch_col) - set(h5ad_lookup)
  if missing_batches:
    raise ValueError('Input batch identifiers do not match H5AD manifest')
  plan_df    = plan_df.with_columns(batch_col.alias('batch'))
  # Check grouping metadata before any split or sampling decisions.
  required = needed - {'batch'}
  missing  = required - set(plan_df.columns)
  if missing:
    raise ValueError(f'Missing grouping/sampling columns: {sorted(missing)}')
  if _missing_values(plan_df, required):
    raise ValueError('Missing values in grouping/sampling columns')

  return plan_df


def _assign_splits(plan_df, rng, *, test_group_var, test_groups, validation_group_var,
  validation_fraction, validation_variable_combinations, excluded_labels,
  min_cells_per_type):


  # Reserve explicitly held-out groups for the independent test set.
  plan_df = plan_df.with_columns(pl.lit('development').alias('split'))
  if bool(test_group_var) != bool(test_groups):
    raise ValueError('test_group_var and nonempty test_groups must be supplied together')
  if test_group_var:
    if _missing_values(plan_df, [test_group_var]):
      raise ValueError(f'Missing test grouping values: {test_group_var}')
    groups  = plan_df[test_group_var].cast(pl.String)
    missing = set(test_groups) - set(groups)
    if missing:
      raise ValueError(f'Test groups absent from annotations/integration: {sorted(missing)}')
    plan_df = plan_df.with_columns(
      pl.when(groups.is_in(test_groups)).then(pl.lit('test'))
        .otherwise(pl.col('split')).alias('split'))

  # Validation units must be disjoint from test groups and have consistent metadata.
  development_df = plan_df.filter(pl.col('split') != 'test')
  group_var      = validation_group_var
  test_df        = plan_df.filter(pl.col('split') == 'test')
  if set(development_df[group_var]) & set(test_df[group_var]):
    raise ValueError('Validation groups overlap the held-out test groups')
  variables = validation_variable_combinations
  if _missing_values(development_df, [group_var] + variables):
    raise ValueError('Missing validation grouping/combination values')
  group_units_df = development_df.select(list(dict.fromkeys([group_var] + variables))).unique(maintain_order = True)
  if group_units_df[group_var].n_unique() != len(group_units_df):
    raise ValueError('A validation group spans several variable combinations; choose compatible variables')

  # Select whole validation units within each variable combination.
  validation = set()
  groups     = _sorted_groups(group_units_df, variables) if variables else [group_units_df]
  for group_df in groups:
    n = min(len(group_df) - 1, max(1, int(len(group_df) * validation_fraction)))
    if n > 0:
      validation.update(rng.choice(group_df[group_var].to_numpy(), n, replace = False))
  if not validation:
    raise ValueError('No whole groups available for validation; supply more groups or different validation combinations')

  # Move whole validation groups to training until rare labels have enough cells.
  eligible_df                           = development_df.filter(pl.col('label').is_not_null() & ~pl.col('label').is_in(excluded_labels))
  label_counts_df                       = eligible_df.group_by('label').len()
  labels_with_minimum_development_cells = set(label_counts_df.filter(pl.col('len') >= min_cells_per_type)['label'])
  validation_groups_moved_to_training   = []
  for label in sorted(labels_with_minimum_development_cells):
    group_df = development_df.filter(pl.col('label') == label)
    while sum(value not in validation for value in group_df[group_var]) < min_cells_per_type:
      group_sizes_df = (group_df.filter(pl.col(group_var).is_in(list(validation)))
        .group_by(group_var).len().sort(group_var))
      if group_sizes_df.is_empty():
        raise ValueError(f'Unable to retain training support for {label}')
      # Preserve the previous descending quicksort tie order and seeded selections.
      indices = np.arange(len(group_sizes_df))[::-1]
      order   = indices[np.argsort(group_sizes_df['len'].to_numpy()[::-1], kind = 'quicksort')][::-1]
      unit    = group_sizes_df[group_var][int(order[0])]
      validation.remove(unit)
      validation_groups_moved_to_training.append({'label': label, 'group': str(unit)})

  # Apply the final group assignments after moving any groups back to training.
  plan_df = plan_df.with_columns(
    pl.when(pl.col('split') == 'test').then(pl.lit('test'))
      .when(pl.col(group_var).is_in(list(validation))).then(pl.lit('validation'))
      .otherwise(pl.lit('train')).alias('split'))
  if not plan_df.filter(pl.col('split') == 'validation').height:
    raise ValueError('No validation groups remain after rescuing rare labels')
  return plan_df, labels_with_minimum_development_cells, validation_groups_moved_to_training


def _prepare_plan_labels(plan_df, labels_with_minimum_development_cells, excluded_labels, vote_col,
  cluster_majority_vote_min_purity):

  # Preserve reference labels and flag excluded, unlabelled or unsupported cells.
  plan_df = plan_df.with_columns(
    pl.col('label').alias('original_label'),
    pl.when(pl.col('label').is_null()).then(pl.lit('unlabeled'))
      .when(pl.col('label').is_in(excluded_labels)).then(pl.lit('excluded'))
      .when(~pl.col('label').is_in(list(labels_with_minimum_development_cells)) & pl.col('split').is_in(['train', 'validation']))
      .then(pl.lit('unsupported')).otherwise(pl.col('split')).alias('split'))

  # Optional voting uses training labels only, leaving the test labels untouched.
  if vote_col:
    cluster_counts_df = plan_df.filter(pl.col('split') == 'train').group_by([vote_col, 'label']).len()
    cluster_labels_df = (cluster_counts_df.with_columns(pl.col('len').sum().over(vote_col).alias('total'))
      .sort([vote_col, 'len', 'label'], descending = [False, True, False])
      .unique(vote_col, keep = 'first', maintain_order = True)
      .filter((pl.col('len') / pl.col('total') >= cluster_majority_vote_min_purity) & (pl.col('total') >= 10)))
    mapping = dict(cluster_labels_df.select(vote_col, 'label').iter_rows())
    plan_df = plan_df.with_columns(
      pl.when(pl.col('split').is_in(['train', 'validation']))
        .then(pl.col(vote_col).replace_strict(mapping, default = pl.col('label')))
        .otherwise(pl.col('label')).alias('label'))

  return plan_df


def _select_plan_cells(plan_df, rng, sampling_variable_combinations,
  train_cells_per_type, validation_cells_per_type):

  # Sample independently within each label and split, retaining plan row order.
  plan_df       = plan_df.with_row_index('_plan_row')
  selected_rows = []
  eligible_df   = plan_df.filter(pl.col('split').is_in(['train', 'validation']))
  for group_df in _sorted_groups(eligible_df, ['label', 'split']):
    cap       = train_cells_per_type if group_df['split'][0] == 'train' else validation_cells_per_type
    sample_df = _balanced_sample(group_df, cap, sampling_variable_combinations, rng)
    selected_rows.extend(sample_df['_plan_row'].to_list())
  plan_df    = plan_df.with_columns(pl.col('_plan_row').is_in(selected_rows).alias('selected')).drop('_plan_row')
  # Require both splits and training support for every validation label.
  selected_df = plan_df.filter(pl.col('selected'))
  if not {'train', 'validation'} <= set(selected_df['split']):
    raise ValueError('Both training and validation cells are required')
  train_labels = set(selected_df.filter(pl.col('split') == 'train')['label'])
  if not set(selected_df.filter(pl.col('split') == 'validation')['label']) <= train_labels:
    raise ValueError('Validation label absent from training')
  plan_df = plan_df.with_columns(
    pl.col('label').is_in(list(train_labels)).fill_null(False).alias('label_supported_in_training'))

  return plan_df


def _write_training_plan(plan_df, plan_f, coverage_f, counts_f, configuration, validation_groups_moved_to_training):
  selected_df  = plan_df.filter(pl.col('selected'))
  train_labels = set(selected_df.filter(pl.col('split') == 'train')['label'])

  # Save the full plan, selection coverage and decisions needed for provenance.
  plan_df.write_csv(plan_f, compression = 'gzip')
  coverage_df = (plan_df.group_by(['label', 'split']).agg(
    pl.len().alias('n_cells'), pl.col('selected').sum().alias('n_selected'))
    .sort(['label', 'split'], nulls_last = True))
  coverage_df.write_csv(coverage_f)
  (selected_df.filter(pl.col('split') == 'train')
    .group_by('label').len().sort('label')
    .write_csv(counts_f))
  # Record planning settings and any validation groups moved back to training.
  metadata = dict(
    configuration               = configuration,
    rescued_validation_groups   = validation_groups_moved_to_training,
    n_annotated                 = len(plan_df),
    n_train                     = selected_df.filter(pl.col('split') == 'train').height,
    n_validation                = selected_df.filter(pl.col('split') == 'validation').height,
    n_test                      = plan_df.filter(pl.col('split') == 'test').height,
    n_classes                   = len(train_labels),
    majority_vote_training_only = True)
  Path(str(plan_f) + '.json').write_text(json.dumps(metadata, indent = 2) + '\n')
  print(coverage_df, flush = True)


def _balanced_sample(input_df, n, variables, rng):

  # Keep all cells below the cap; sample uniformly when no balancing is requested.
  if len(input_df) <= n:
    return input_df.clone()
  if not variables:
    return input_df[rng.choice(len(input_df), n, replace = False)]

  # Balance joint combinations, redistributing unused quotas to larger groups.
  indexed_df = input_df.with_row_index('_sample_row')
  groups     = [group_df['_sample_row'].to_numpy() for group_df in _sorted_groups(indexed_df, variables)]
  order      = rng.permutation(len(groups))
  quota      = np.zeros(len(groups), dtype = int)
  remaining  = n
  while remaining:
    active = [i for i in order if quota[i] < len(groups[i])]
    share  = max(1, remaining // len(active))
    for i in active:
      take = min(share, len(groups[i]) - quota[i], remaining)
      quota[i] += take
      remaining -= take
      if not remaining:
        break
  idx = np.concatenate([rng.choice(indices, q, replace = False) for indices, q in zip(groups, quota) if q])
  return input_df[idx]


def _sorted_groups(input_df, variables):
  # Keep original row order within lexicographically ordered combinations.
  if not variables:
    return [input_df]
  return (
    input_df.sort(variables, nulls_last = True, maintain_order = True)
      .partition_by(variables, maintain_order = True))


def _missing_values(input_df, columns):
  # Treat floating-point NaN as missing alongside Polars null values.
  for name in columns:
    values = input_df[name]
    if values.null_count() or (values.dtype.is_float() and values.is_nan().any()):
      return True
  return False


# Two-pass training and feature selection

def fit_plan(plan_f, h5ads_yaml, output_dir, n_cores, *, xgboost_tag, use_gpu,
  gene_exclude_re, min_cells_expressed, seed, pass1_subsample, pass1_colsample_bytree,
  pass1_learning_rate, pass1_nrounds, pass1_early_stopping, gain_threshold, min_genes,
  max_genes, pass2_colsample_bytree, pass2_learning_rate, pass2_nrounds,
  pass2_early_stopping, label_map_f):

  import xgboost as xgb

  out = Path(output_dir)
  out.mkdir(parents = True, exist_ok = True)
  gpu_device = _verify_gpu(n_cores) if use_gpu else 'cpu'

  # Load only the cells selected by the saved training plan.
  plan_df     = pl.read_csv(plan_f, schema_overrides = {'label': pl.String})
  plan_df     = plan_df.filter(pl.col('selected') & pl.col('split').is_in(['train', 'validation']))
  h5ad_lookup = yaml.safe_load(Path(h5ads_yaml).read_text())
  raw, norm, genes, types, plan_df = _load_plan_expression(plan_df, h5ad_lookup, normalize_counts)

  # Encode labels consistently across the training and validation matrices.
  train   = plan_df['split'].to_numpy() == 'train'
  val     = ~train
  classes = sorted(plan_df.filter(train)['label'].unique())
  if not set(plan_df.filter(val)['label']) <= set(classes):
    raise ValueError('Validation label absent from training')
  class_indices = {name: i for i, name in enumerate(classes)}
  y             = np.array([class_indices[name] for name in plan_df['label']])

  # Detection filtering uses training cells only; held-out cells cannot select genes.
  norm, raw, genes, types = _filter_genes_by_type(norm, raw, genes, types, gene_exclude_re)
  keep  = _training_gene_mask(raw, train, min_cells_expressed)
  norm  = norm[:, keep].tocsr()
  raw   = raw[:, keep].tocsr()
  genes = [g for g, k in zip(genes, keep) if k]
  if not genes:
    raise ValueError('No genes pass training-only filtering')
  xtrain = norm[train]
  xval   = norm[val]

  # Pass 1 considers all retained genes with the lower column-sampling fraction.
  model1 = _run_xgboost_pass1(
    X_train          = xtrain,
    y_train          = y[train],
    X_val            = xval,
    y_val            = y[val],
    n_cores          = n_cores,
    seed             = seed,
    use_gpu          = use_gpu,
    subsample        = pass1_subsample,
    colsample_bytree = pass1_colsample_bytree,
    learning_rate    = pass1_learning_rate,
    nrounds          = pass1_nrounds,
    early_stopping   = pass1_early_stopping)
  best1  = int(model1.best_iteration)
  model1 = _best_iteration_model(model1)

  # Retain genes contributing the requested cumulative gain, within size bounds.
  selected, names, gene_importance_df = _select_features_by_gain(
    model          = model1,
    gene_names     = genes,
    gain_threshold = gain_threshold,
    min_genes      = min_genes,
    max_genes      = max_genes)
  del model1, norm
  gc.collect()
  xtrain = xtrain[:, selected].tocsr()
  xval   = xval[:, selected].tocsr()

  # Pass 2 fits the final classifier on those genes, stopping on validation loss.
  model2 = _run_xgboost_pass2(
    X_train          = xtrain,
    y_train          = y[train],
    X_val            = xval,
    y_val            = y[val],
    n_cores          = n_cores,
    seed             = seed,
    use_gpu          = use_gpu,
    subsample        = pass1_subsample,
    colsample_bytree = pass2_colsample_bytree,
    learning_rate    = pass2_learning_rate,
    nrounds          = pass2_nrounds,
    early_stopping   = pass2_early_stopping)
  best2  = int(model2.best_iteration)
  model2 = _best_iteration_model(model2)

  # Export the best iteration, selected-cell predictions and raw-count summaries.
  pb = _compute_pseudobulk(raw, plan_df, genes, selected)
  _save_outputs(
    model              = model2,
    class_names        = classes,
    selected_genes     = names,
    gene_importance_df = gene_importance_df,
    output_dir         = out,
    ref_tag            = xgboost_tag,
    cells_df           = plan_df,
    label_map          = _load_label_mapping(label_map_f),
    X_train            = xtrain,
    y_train            = y[train],
    X_val              = xval,
    y_val              = y[val],
    pseudobulk_adata   = pb)

  # Record fitting settings separately from display specifications.
  configuration = dict(
    xgboost_tag            = xgboost_tag,
    use_gpu                = bool(use_gpu),
    gene_exclude_re        = gene_exclude_re,
    min_cells_expressed    = min_cells_expressed,
    seed                   = seed,
    pass1_subsample        = pass1_subsample,
    pass1_colsample_bytree = pass1_colsample_bytree,
    pass1_learning_rate    = pass1_learning_rate,
    pass1_nrounds          = pass1_nrounds,
    pass1_early_stopping   = pass1_early_stopping,
    gain_threshold         = gain_threshold,
    min_genes              = min_genes,
    max_genes              = max_genes,
    pass2_colsample_bytree = pass2_colsample_bytree,
    pass2_learning_rate    = pass2_learning_rate,
    pass2_nrounds          = pass2_nrounds,
    pass2_early_stopping   = pass2_early_stopping,
    label_map_f            = label_map_f)
  model_f = out / f'{xgboost_tag}_xgboost_model.json'
  card    = dict(
    configuration                  = configuration,
    n_train                        = int(train.sum()),
    n_validation                   = int(val.sum()),
    n_classes                      = len(classes),
    n_genes                        = len(names),
    n_genes_tested                 = len(genes),
    pass1_best_iteration           = best1,
    pass2_best_iteration           = best2,
    detection_filter_training_only = True,
    best_iteration_model_sliced    = True,
    model_sha256                   = hashlib.sha256(model_f.read_bytes()).hexdigest(),
    plan_sha256                    = hashlib.sha256(Path(plan_f).read_bytes()).hexdigest(),
    xgboost_version                = xgb.__version__,
    xgboost_build                  = xgb.build_info(),
    training_device                = gpu_device,
    normalization                  = 'log1p(10000 * counts / all-input-gene library size)',
    status                         = 'Unreviewed candidate; not registered or published')
  (out / f'{xgboost_tag}_model_card.json').write_text(json.dumps(card, indent = 2) + '\n')


def _verify_gpu(n_cores):
  import xgboost as xgb

  # Run a tiny CUDA fit and reject silent fallback to CPU.
  if not xgb.build_info().get('USE_CUDA'):
    raise RuntimeError('GPU training requested but XGBoost was built without CUDA')
  probe = xgb.DMatrix(np.array([[0., 1.], [1., 0.], [1., 1.], [0., 0.]], dtype = np.float32),
    label = np.array([0, 1, 1, 0]), nthread = min(n_cores, 2))
  model = xgb.train({'objective': 'binary:logistic', 'device': 'cuda',
    'tree_method': 'hist', 'nthread': min(n_cores, 2)}, probe, num_boost_round = 1)
  device = json.loads(model.save_config())['learner']['generic_param']['device']
  if not device.startswith('cuda'):
    raise RuntimeError('GPU training requested but XGBoost fell back to CPU')
  print('Verified XGBoost GPU execution on ' + device, flush = True)
  return device


def _load_plan_expression(plan_df, h5ad_lookup, normalize):
  import scipy.sparse as sp

  # Read batches in deterministic order and normalise before gene filtering.
  raw_parts, norm_parts, ordered = [], [], []
  genes = types = None
  for batch_df in _sorted_groups(plan_df, ['batch']):
    batch = batch_df['batch'][0]
    entry = h5ad_lookup[batch]
    path  = entry['path'] if isinstance(entry, dict) else entry
    raw, current, current_types = _read_selected_h5ad(path, batch_df['cell_id'].to_list(), genes)
    if genes is None:
      genes, types = current, current_types
    raw_parts.append(raw)
    norm_parts.append(normalize(raw).astype(np.float32))
    ordered.append(batch_df)
    print(f'Loading {batch}: {len(batch_df)} selected cells', flush = True)
  # Stack count matrices and metadata in the same cell order.
  raw = sp.vstack(raw_parts, format = 'csr')
  raw_parts.clear()
  norm = sp.vstack(norm_parts, format = 'csr')
  norm_parts.clear()
  return raw, norm, genes, types, pl.concat(ordered)


def _read_selected_h5ad(path, cells, reference = None):
  import scipy.sparse as sp

  import anndata as ad

  # Backed read-only access keeps count inputs unchanged.
  before = (Path(path).stat().st_size, Path(path).stat().st_mtime_ns)
  data   = ad.read_h5ad(path, backed = 'r')
  try:
    ids = (data.obs['cell_id'] if 'cell_id' in data.obs else data.obs_names).astype(str).tolist()
    if len(set(ids)) != len(ids):
      raise ValueError(f'Duplicate cell IDs: {path}')
    positions = {cell: i for i, cell in enumerate(ids)}
    idx       = np.array([positions.get(cell, -1) for cell in cells])
    if np.any(idx < 0) or len(set(cells)) != len(cells):
      raise ValueError(f'Missing/duplicated selected cells: {path}')
    genes = data.var_names.astype(str).tolist()
    if len(set(genes)) != len(genes) or any(not g.strip() for g in genes):
      raise ValueError(f'Invalid gene IDs: {path}')
    types = (data.var['gene_type'].astype(object).fillna('unknown').astype(str).tolist()
      if 'gene_type' in data.var else ['unknown'] * len(genes))

    # HDF5 needs sorted row indices; restore the requested order after reading.
    order  = np.argsort(idx)
    matrix = sp.csr_matrix(data.X[idx[order], :])[np.argsort(order)].tocsr()

    # Align identical gene sets; mismatched contents remain an error.
    if reference is not None and genes != reference:
      if set(genes) != set(reference):
        raise ValueError(f'Different gene sets: {path}')
      import warnings
      warnings.warn(f'Same genes, different order: aligning {path}', stacklevel = 2)
      positions = {gene: i for i, gene in enumerate(genes)}
      indices   = [positions[gene] for gene in reference]
      matrix    = matrix[:, indices].tocsr()
      types     = [types[i] for i in indices]
      genes     = list(reference)
  finally:
    data.file.close()
  if before != (Path(path).stat().st_size, Path(path).stat().st_mtime_ns):
    raise RuntimeError(f'Reading changed input signature: {path}')
  return matrix, genes, types


def normalize_counts(X, scale_factor = 10000):
  import scipy.sparse as sp

  # Scale by each cell's full-input library size, then apply log1p.
  X_csr     = X.tocsr().astype(np.float64)
  lib_sizes = np.array(X_csr.sum(axis = 1)).flatten()
  lib_sizes[lib_sizes == 0] = 1.0
  inv_lib = sp.diags(1.0 / lib_sizes)
  X_norm  = inv_lib @ X_csr
  X_norm *= scale_factor
  X_norm.data = np.log1p(X_norm.data)
  return X_norm


def _filter_genes_by_type(X_norm, X_raw, gene_names, gene_types, gene_exclude_re):
  # remove genes whose gene_type matches the exclusion regex.

  if gene_exclude_re is None or gene_exclude_re == "":
    return X_norm, X_raw, gene_names, gene_types

  pattern   = re.compile(gene_exclude_re)
  keep_mask = np.array([not pattern.search(gt) for gt in gene_types])
  n_before  = len(gene_names)
  n_after   = int(keep_mask.sum())

  X_norm_filtered     = X_norm[:, keep_mask]
  X_raw_filtered      = X_raw[:, keep_mask]
  gene_names_filtered = [g for g, keep in zip(gene_names, keep_mask) if keep]
  gene_types_filtered = [gt for gt, keep in zip(gene_types, keep_mask) if keep]

  print(f"  Gene type filter (exclude '{gene_exclude_re}'): {n_before} → {n_after} genes "
        f"({n_before - n_after} removed)")

  return X_norm_filtered, X_raw_filtered, gene_names_filtered, gene_types_filtered


def _training_gene_mask(raw, train, minimum):

  # Count detection in training cells only, avoiding held-out information.
  detected = np.asarray(raw[train].getnnz(axis = 0)).ravel()
  return detected >= minimum


def _run_xgboost_pass1(X_train, y_train, X_val, y_val,
  n_cores, seed, use_gpu, subsample, colsample_bytree, learning_rate,
  nrounds, early_stopping):
  import xgboost as xgb

  num_class = int(y_train.max()) + 1

  # Fit the broad-feature model, with validation log loss controlling early stopping.
  params = _build_xgb_params(n_cores, seed, use_gpu, num_class,
    subsample        = subsample,
    colsample_bytree = colsample_bytree,
    learning_rate    = learning_rate,
  )

  dtrain = xgb.DMatrix(X_train, label = y_train)
  dval   = xgb.DMatrix(X_val, label = y_val)

  model = xgb.train(
    params                = params, dtrain = dtrain,
    num_boost_round       = nrounds,
    evals                 = [(dtrain, "train"), (dval, "val")],
    early_stopping_rounds = early_stopping,
    verbose_eval          = 50,
  )

  print(f"  Best iteration: {model.best_iteration}")
  return model


def _run_xgboost_pass2(X_train, y_train, X_val, y_val, n_cores, seed, use_gpu,
  subsample, colsample_bytree, learning_rate, nrounds, early_stopping):

  import xgboost as xgb

  num_class = int(y_train.max()) + 1

  params = _build_xgb_params(n_cores, seed, use_gpu, num_class,
    subsample        = subsample,
    colsample_bytree = colsample_bytree,
    learning_rate    = learning_rate,
  )

  # Use stable names for the selected features in the saved model.
  feat_names = [f"g{i}" for i in range(X_train.shape[1])]
  dtrain     = xgb.DMatrix(X_train, label = y_train, feature_names = feat_names)
  dval       = xgb.DMatrix(X_val, label = y_val, feature_names = feat_names)

  model = xgb.train(
    params                = params, dtrain = dtrain,
    num_boost_round       = nrounds,
    evals                 = [(dtrain, "train"), (dval, "val")],
    early_stopping_rounds = early_stopping,
    verbose_eval          = 50,
  )

  print(f"  Best iteration: {model.best_iteration}")
  return model


def _build_xgb_params(n_cores, seed, use_gpu, num_class, **overrides):
  # Use multiclass probabilities and log loss in both training passes.
  params = {
    "objective": "multi:softprob",
    "num_class": num_class,
    "eval_metric": "mlogloss",
    "nthread": n_cores,
    "seed": seed,
  }
  if use_gpu:
    params["device"] = "cuda"
    params["tree_method"] = "hist"
  params.update(overrides)
  return params


def _select_features_by_gain(model, gene_names,
  gain_threshold, min_genes, max_genes):


  # Map XGBoost feature indices back to the input genes.
  importance = model.get_score(importance_type = "gain")

  gene_gains = np.zeros(len(gene_names))
  for feat, gain in importance.items():
    gene_gains[int(feat[1:])] = gain

  sorted_idx   = np.argsort(-gene_gains)
  sorted_gains = gene_gains[sorted_idx]

  # Accumulate gain over features actually used by the first model.
  nonzero_mask = sorted_gains > 0
  sorted_idx   = sorted_idx[nonzero_mask]
  sorted_gains = sorted_gains[nonzero_mask]

  total_gain = sorted_gains.sum()
  cum_gain   = np.cumsum(sorted_gains) / total_gain

  # Find the cumulative-gain cutoff and apply the configured size bounds.
  n_sel = int(np.searchsorted(cum_gain, gain_threshold)) + 1
  n_sel = max(n_sel, min_genes)
  n_sel = min(n_sel, max_genes)
  n_sel = min(n_sel, len(sorted_idx))

  sel_indices    = sorted_idx[:n_sel]
  sel_gene_names = [gene_names[i] for i in sel_indices]

  gene_importance_df = pl.DataFrame({
    "gene_id": [gene_names[i] for i in sorted_idx],
    "gain": sorted_gains.tolist(),
    "cumulative_gain_frac": cum_gain.tolist(),
  })

  print(f"  Total genes with non-zero gain: {len(sorted_idx)}")
  print(f"  Genes at {gain_threshold*100:.0f}% cumulative gain: {n_sel}")
  print(f"  Top 10 genes: {sel_gene_names[:10]}")

  return sel_indices, sel_gene_names, gene_importance_df


def _best_iteration_model(model):
  # Drop rounds after the best validation iteration before scoring or export.
  return model[:int(model.best_iteration) + 1]


# Save model outputs and selected-cell summaries

def _save_outputs(model, class_names, selected_genes, gene_importance_df,
  output_dir, ref_tag, cells_df, label_map = None, X_train = None, y_train = None,
  X_val = None, y_val = None, pseudobulk_adata = None):

  # save model, gene list, classes, importance, predictions, and pseudobulk
  out_dir = Path(output_dir)

  model_path = out_dir / f"{ref_tag}_xgboost_model.json"
  model.save_model(str(model_path))
  print(f"  Model: {model_path}")

  cls_path = out_dir / f"{ref_tag}_allowed_cls.csv"
  pl.DataFrame({"class": class_names}).write_csv(str(cls_path))
  print(f"  Classes: {cls_path}")

  genes_path = out_dir / f"{ref_tag}_selected_genes.txt"
  with open(genes_path, "w") as f:
    f.write("\n".join(selected_genes) + "\n")
  print(f"  Genes: {genes_path}")

  imp_path = out_dir / f"{ref_tag}_gene_importance.csv"
  gene_importance_df.write_csv(str(imp_path))
  print(f"  Importance: {imp_path}")

  if label_map is not None:
    map_path = out_dir / f"{ref_tag}_label_map.csv"
    pl.DataFrame({
      "label": list(label_map.keys()),
      "coarse_label": list(label_map.values()),
    }).write_csv(str(map_path))
    print(f"  Label map: {map_path}")

  # save predictions for all cells (train + validation)
  if X_val is not None and X_train is not None:
    preds_path = out_dir / f"{ref_tag}_predictions.csv.gz"
    preds_df   = _make_predictions_df(
      model, X_train, y_train, X_val, y_val, cells_df, class_names, label_map
    )
    preds_df.write_csv(str(preds_path), compression = "gzip")
    print(f"  Predictions: {preds_path}")

  if pseudobulk_adata is not None:
    pb_path = out_dir / f"{ref_tag}_pseudobulk.h5ad"
    pseudobulk_adata.write_h5ad(str(pb_path))
    print(f"  Pseudobulk: {pb_path}")


def _make_predictions_df(model, X_train, y_train, X_val, y_val, cells_df, class_names, label_map = None):
  # Generate predictions in the selected cell order.
  import xgboost as xgb

  feat_names = [f"g{i}" for i in range(X_val.shape[1])]

  # Predict on validation set
  dval      = xgb.DMatrix(X_val, feature_names = feat_names)
  val_probs = model.predict(dval)
  val_pred  = val_probs.argmax(axis = 1)
  val_pmax  = val_probs.max(axis = 1)

  # Predict on training set
  dtrain      = xgb.DMatrix(X_train, feature_names = feat_names)
  train_probs = model.predict(dtrain)
  train_pred  = train_probs.argmax(axis = 1)
  train_pmax  = train_probs.max(axis = 1)

  # Fill predictions in cells_df order using the split mask.
  train_mask = cells_df["split"].to_numpy() == "train"
  all_pred   = np.empty(len(cells_df), dtype = int)
  all_pmax   = np.empty(len(cells_df))
  all_pred[train_mask] = train_pred
  all_pred[~train_mask] = val_pred
  all_pmax[train_mask] = train_pmax
  all_pmax[~train_mask] = val_pmax

  pred_labels = [class_names[i] for i in all_pred]

  predictions_df = cells_df.select(["cell_id", "batch", "label", "split"]).with_columns([
    pl.Series("predicted_label", pred_labels),
    pl.Series("probability", all_pmax),
  ])

  # Add coarse labels if mapping is available
  if label_map is not None:
    predictions_df = predictions_df.with_columns([
      pl.col("label").replace(label_map).alias("coarse_true"),
      pl.col("predicted_label").replace(label_map).alias("coarse_predicted"),
    ])

  return predictions_df


def _compute_pseudobulk(X_raw, cells_df, gene_names, sel_indices):
  # AnnData requires pandas for obs, var and labelled obsm tables.

  import anndata as ad

  X_sel          = X_raw[:, sel_indices]
  sel_gene_names = [gene_names[i] for i in sel_indices]

  labels  = cells_df["label"].to_list()
  batches = cells_df["batch"].to_list()

  unique_labels  = sorted(set(labels))
  unique_batches = sorted(set(batches))
  n_batches      = len(unique_batches)
  n_genes        = len(sel_gene_names)
  batch_to_idx   = {s: i for i, s in enumerate(unique_batches)}

  # Index cells by label and batch for raw-count aggregation.
  groups: dict[tuple[str, str], list[int]] = {}
  for i, (lbl, bid) in enumerate(zip(labels, batches)):
    key = (lbl, bid)
    if key not in groups:
      groups[key] = []
    groups[key].append(i)

  # Store one batch-by-gene count matrix per label, including empty combinations.
  layers = {}
  n_cells_dict: dict[str, list[int]] = {}
  for lbl in unique_labels:
    count_mat         = np.zeros((n_batches, n_genes), dtype = np.float64)
    n_cells_per_batch = [0] * n_batches
    for bid in unique_batches:
      key = (lbl, bid)
      if key in groups:
        indices = groups[key]
        count_mat[batch_to_idx[bid], :] = np.array(
          X_sel[indices, :].sum(axis = 0)).flatten()
        n_cells_per_batch[batch_to_idx[bid]] = len(indices)

    layers[lbl] = count_mat
    n_cells_dict[lbl] = n_cells_per_batch

  obs_df = pl.DataFrame({"batch": unique_batches}).to_pandas()
  var_df = pl.DataFrame({"gene_id": sel_gene_names}).to_pandas()
  obs_df.index = unique_batches
  var_df.index = sel_gene_names

  pb_adata = ad.AnnData(
    X      = np.zeros((n_batches, n_genes)),
    obs    = obs_df,
    var    = var_df,
    layers = layers,
  )

  # Keep contributing cell counts alongside the pseudobulk sums.
  n_cells_mat = np.column_stack([n_cells_dict[lbl] for lbl in unique_labels])
  n_cells_df  = pl.DataFrame(n_cells_mat, schema = unique_labels, orient = 'row').to_pandas()
  n_cells_df.index = unique_batches
  pb_adata.obsm["n_cells"] = n_cells_df

  print(f"  Pseudobulk raw counts: {n_batches} batches, {n_genes} genes, "
        f"{len(unique_labels)} labels (stored as layers)")

  return pb_adata


def _load_label_mapping(label_map_f):

  # Leave labels unchanged when no legacy mapping file is supplied.
  if label_map_f is None:
    return None

  records, _ = read_label_specs(label_map_f = label_map_f)
  return label_mapping(records)


# Full-data evaluation and label specifications

def aggregate_plan_predictions(pred_fs, plan_f, output_f, integration_f):
  # Combine per-batch calls and require one prediction per integrated cell.
  predictions_df = pl.concat([pl.read_csv(path, columns = ['cell_id', 'predicted_label_naive', 'probability_naive'],
    schema_overrides = {'cell_id': pl.String, 'predicted_label_naive': pl.String,
      'probability_naive': pl.Float64}) for path in pred_fs])
  if predictions_df['cell_id'].n_unique() != len(predictions_df):
    raise ValueError('Duplicate full-data predictions')
  predictions_df = predictions_df.select('cell_id', pl.col('predicted_label_naive').alias('predicted_label'),
    pl.col('probability_naive').alias('probability'))
  integration_cells_df = pl.read_csv(integration_f, columns = ['cell_id'])
  predictions_df       = predictions_df.join(integration_cells_df, on = 'cell_id', how = 'inner', maintain_order = 'left')
  if len(predictions_df) != len(integration_cells_df):
    raise ValueError('Integrated cells missing from full-data predictions')

  # Preserve all integrated cells, including those without a reference label.
  plan_df   = pl.read_csv(plan_f, schema_overrides = {'label': pl.String, 'original_label': pl.String})
  fields    = ['cell_id', 'label', 'original_label', 'split', 'selected', 'label_supported_in_training']
  result_df = predictions_df.join(plan_df.select(fields), on = 'cell_id', how = 'left', maintain_order = 'left').with_columns(
    pl.col('split').fill_null('unlabeled'), pl.col('selected').fill_null(False),
    pl.col('label_supported_in_training').fill_null(False))
  if not set(plan_df['cell_id']) <= set(result_df['cell_id']):
    raise ValueError('Planned cells missing from full-data predictions')
  metadata_f    = Path(str(plan_f) + '.json')
  configuration = json.loads(metadata_f.read_text())['configuration'] if metadata_f.exists() else {}
  if configuration.get('label_map_f'):
    mapping   = _load_label_mapping(configuration['label_map_f'])
    result_df = result_df.with_columns(pl.col('label').replace(mapping).alias('coarse_true'),
      pl.col('predicted_label').replace(mapping).alias('coarse_predicted'))
  result_df.write_csv(output_f, compression = 'gzip')

  # Score unselected development cells separately from fitting/validation cells.
  prefix = str(output_f).removesuffix('_fulldata_predictions.csv.gz')
  _write_evaluation(_evaluation_counts(result_df), prefix, 'label')


def write_label_specifications(cls_f, counts_f, predictions_f, output_dir, *,
  xgboost_tag, label_specs_f, coarse_label_specs_f, label_map_f):

  # Use training abundance for the default display order unless a table supplies one.
  classes       = pl.read_csv(cls_f)['class'].to_list()
  label_counts  = dict(pl.read_csv(counts_f).select('label', 'len').iter_rows())
  default_order = sorted(classes, key = lambda name: (-label_counts.get(name, 0), name))
  specs, coarse = read_label_specs(label_specs_f, coarse_label_specs_f, default_order, label_map_f)
  out = Path(output_dir)
  out.mkdir(parents = True, exist_ok = True)
  prefix = out / xgboost_tag
  pl.DataFrame(specs, schema = {'label': pl.String, 'colour': pl.String, 'coarse_label': pl.String}).write_csv(
    str(prefix) + '_label_specs.csv')
  pl.DataFrame(coarse, schema = {'coarse_label': pl.String, 'coarse_colour': pl.String}).write_csv(
    str(prefix) + '_coarse_label_specs.csv')

  # Coarse evaluation is a direct lookup of fine labels, without another vote.
  mapping = label_mapping(specs)
  if mapping:
    evaluation_df = _evaluation_counts(pl.scan_csv(predictions_f,
      schema_overrides = {'label': pl.String, 'predicted_label': pl.String}), mapping)
  else:
    evaluation_df = pl.DataFrame(schema = {'split': pl.String, 'label': pl.String,
      'predicted_label': pl.String, 'len': pl.UInt32})
  _write_evaluation(evaluation_df, str(prefix) + '_coarse', 'coarse_label')


def _evaluation_counts(predictions_df, mapping = None):
  # Summarise lazily so full prediction tables need not be converted or collected.
  if isinstance(predictions_df, pl.DataFrame):
    predictions_df = predictions_df.lazy()
  if 'selected' not in predictions_df.collect_schema():
    predictions_df = predictions_df.with_columns(pl.lit(True).alias('selected'))
  predictions_df = predictions_df.filter(
    pl.col('label').is_not_null() & ~pl.col('split').is_in(['excluded', 'unsupported', 'unlabeled']))
  predictions_df = predictions_df.with_columns(
    pl.when(~pl.col('selected') & pl.col('split').is_in(['train', 'validation']))
      .then(pl.lit('development_holdout')).otherwise(pl.col('split')).alias('split'))
  if mapping:
    predictions_df = predictions_df.with_columns(pl.col('label').replace(mapping), pl.col('predicted_label').replace(mapping))
  return predictions_df.group_by(['split', 'label', 'predicted_label']).len().collect(engine = 'streaming')


def _write_evaluation(counts_df, prefix, label_column_name):

  # Evaluate each split using both observed and predicted labels.
  metrics   = []
  per_label = []
  for split_df in _sorted_groups(counts_df, ['split']):
    split     = split_df['split'][0]
    labels    = sorted(set(split_df['label']) | set(split_df['predicted_label']))
    confusion = {(a, b): n for _, a, b, n in split_df.select('split', 'label', 'predicted_label', 'len').iter_rows()}
    scores    = []
    for name in labels:
      support   = sum(n for (a, b), n in confusion.items() if a == name)
      predicted = sum(n for (a, b), n in confusion.items() if b == name)
      tp        = confusion.get((name, name), 0)
      precision = tp / max(predicted, 1)
      recall    = tp / max(support, 1)
      f1        = 2 * precision * recall / (precision + recall) if precision + recall else 0.
      per_label.append(dict(split = split, **{label_column_name: name},
        precision = precision, recall = recall, f1 = f1, n_test = support))
      # Average over labels with reference cells in this split.
      if support:
        scores.append((recall, f1))
    n_cells = sum(confusion.values())
    metrics.append(dict(split = split, n_cells = n_cells,
      accuracy          = sum(n for (a, b), n in confusion.items() if a == b) / n_cells,
      macro_f1          = float(np.mean([s[1] for s in scores])),
      balanced_accuracy = float(np.mean([s[0] for s in scores]))))

  # Write stable table columns even when there are no evaluable cells.
  metric_schema = {'split': pl.String, 'n_cells': pl.Int64, 'accuracy': pl.Float64,
    'macro_f1': pl.Float64, 'balanced_accuracy': pl.Float64}
  label_schema = {'split': pl.String, label_column_name: pl.String, 'precision': pl.Float64,
    'recall': pl.Float64, 'f1': pl.Float64, 'n_test': pl.Int64}
  summary_df = pl.DataFrame(metrics, schema = metric_schema)
  summary_df.write_csv(str(prefix) + '_evaluation.csv')
  pl.DataFrame(per_label, schema = label_schema).write_csv(str(prefix) + '_per_label_evaluation.csv')
  print(summary_df, flush = True)


# Model installation, registry identities and label specifications

def resolve_xgboost_entry(entry, scdata_dir):
  # Resolve a single installed classifier and check its required files.
  registry_df = pl.read_csv(Path(scdata_dir) / 'xgboost/xgboost_models.csv')
  matches_df  = registry_df.filter(pl.col('model') == entry['model'])
  if len(matches_df) != 1:
    raise ValueError(f"Expected one registered XGBoost model named {entry['model']}; found {len(matches_df)}")
  row = matches_df.row(0, named = True)
  for key in ['model_f', 'cls_f', 'genes_f']:
    entry[key] = row[key]
    if not Path(entry[key]).is_file():
      raise FileNotFoundError(f"XGBoost file missing: {entry[key]}; rerun scprocess setup")
  # Older registry entries may omit display and coarse-label specifications.
  entry['label_specs_f'] = row.get('label_specs_f') or ''
  entry['coarse_label_specs_f'] = row.get('coarse_label_specs_f') or ''
  entry['label_map_f'] = row.get('label_map_f') or ''
  if entry['label_specs_f']:
    classes = pl.read_csv(entry['cls_f'])['class'].to_list()
    specs, coarse = read_label_specs(entry['label_specs_f'], entry['coarse_label_specs_f'], classes)
    entry['label_map_f'] = entry['label_specs_f'] if label_mapping(specs) else ''
  # Attach the installed model identity after verifying its checksums.
  entry.update(model_identity(entry['model_f'], entry['model']))
  return entry


def install_builtin(name, dest_dir):
  manifest = builtin_models()[name]
  dest_dir = Path(dest_dir)
  dest_dir.parent.mkdir(parents = True, exist_ok = True)
  # Validate the complete download before replacing any installed files.
  with tempfile.TemporaryDirectory(dir = dest_dir.parent) as directory:
    stage   = Path(directory)
    archive = stage / 'download.tar.gz'
    urllib.request.urlretrieve(manifest['url'], archive)
    if sha256(archive) != manifest['archive_sha256']:
      raise ValueError(f'Archive checksum mismatch for {name}')
    # Accept only the expected regular files within this model archive.
    found = set()
    with tarfile.open(archive, 'r:gz') as tar:
      for member in tar.getmembers():
        path = Path(member.name)
        if member.isdir() and path.parts == (manifest['model_id'],):
          continue
        if (not member.isfile() or len(path.parts) != 2 or
            path.parts[0] != manifest['model_id'] or
            path.name not in manifest['files'] or path.name in found):
          raise ValueError(f'Unexpected model archive member: {member.name}')
        with tar.extractfile(member) as source, (stage / path.name).open('wb') as target:
          shutil.copyfileobj(source, target)
        if sha256(stage / path.name) != manifest['files'][path.name]:
          raise ValueError(f'Model member checksum mismatch: {path.name}')
        found.add(path.name)
    if found != set(manifest['files']):
      raise ValueError(f'Incomplete model archive for {name}')
    # Install files only after the entire staged archive passes validation.
    dest_dir.mkdir(exist_ok = True)
    for filename in found:
      shutil.copy2(stage / filename, dest_dir / filename)
    (dest_dir / f'{name}_manifest.json').write_text(json.dumps(manifest, indent = 2) + '\n')


def model_identity(model_f, name):
  model_f    = Path(model_f)
  checksum   = sha256(model_f)
  manifest_f = model_f.with_name(model_f.name.removesuffix('_xgboost_model.json') + '_manifest.json')
  model_id   = name
  # Installed manifests identify the model and protect all companion files.
  if manifest_f.exists():
    manifest = json.loads(manifest_f.read_text())
    expected = manifest['files'][model_f.name]
    if checksum != expected:
      raise ValueError(f'XGBoost model checksum mismatch: {model_f}')
    for filename, expected_hash in manifest['files'].items():
      if filename != model_f.name and sha256(model_f.with_name(filename)) != expected_hash:
        raise ValueError(f'XGBoost file checksum mismatch: {filename}')
    model_id = manifest['model_id']
  return {'model_id': model_id, 'model_sha256': checksum}


def read_label_specs(label_specs_f = None, coarse_label_specs_f = None, labels = None, label_map_f = None):
  import csv

  if label_specs_f and label_map_f:
    raise ValueError('Use label_specs_f or legacy label_map_f, not both')

  # Read standard labels in the supplied display order.
  source  = label_specs_f or label_map_f
  records = []
  if source:
    with Path(source).open(newline = '') as stream:
      reader        = csv.DictReader(stream)
      columns       = reader.fieldnames or []
      label_columns = [name for name in ['label', 'annotation', 'fine_label'] if name in columns]
      if len(label_columns) != 1:
        raise ValueError(f'{source} needs one label column (label; legacy annotation/fine_label accepted)')
      column = label_columns[0]
      for row in reader:
        records.append({'label': (row.get(column) or '').strip(),
          'colour': (row.get('colour') or '').strip(),
          'coarse_label': (row.get('coarse_label') or '').strip()})
    if not records:
      raise ValueError(f'Empty label specifications: {source}')
  elif labels is not None:
    records = [{'label': str(label), 'colour': '', 'coarse_label': ''} for label in labels]

  # Validate label coverage and complete optional colour/group specifications.
  names = [row['label'] for row in records]
  if any(not name for name in names) or len(set(names)) != len(names):
    raise ValueError('Label specifications contain missing or duplicate labels')
  if labels is not None and not set(labels) <= set(names):
    raise ValueError(f"Label specifications miss model labels: {sorted(set(labels)-set(names))}")
  colours = [row['colour'] for row in records]
  if any(colours) and not all(colours):
    raise ValueError('Supply colours for every label, or omit the colour column')
  mapping = [row['coarse_label'] for row in records]
  if any(mapping) and not all(mapping):
    raise ValueError('Supply coarse_label for every label, or omit the coarse_label column')

  # A separate table can give coarse labels their own order and colours.
  groups = list(dict.fromkeys(value for value in mapping if value))
  coarse = []
  if coarse_label_specs_f:
    with Path(coarse_label_specs_f).open(newline = '') as stream:
      reader = csv.DictReader(stream)
      if 'coarse_label' not in (reader.fieldnames or []):
        raise ValueError('Coarse specifications need a coarse_label column')
      coarse = [{'coarse_label': (row.get('coarse_label') or '').strip(),
        'coarse_colour': (row.get('coarse_colour') or '').strip()} for row in reader]
    if coarse and not groups:
      raise ValueError('Coarse specifications require a coarse_label mapping in label_specs_f')
    names_coarse = [row['coarse_label'] for row in coarse]
    if any(not name for name in names_coarse) or len(names_coarse) != len(set(names_coarse)):
      raise ValueError('Coarse specifications contain missing or duplicate labels')
    if not set(groups) <= set(names_coarse):
      raise ValueError(f"Coarse specifications miss groups: {sorted(set(groups)-set(names_coarse))}")
  else:
    coarse = [{'coarse_label': name, 'coarse_colour': ''} for name in groups]
  coarse_colours = [row['coarse_colour'] for row in coarse]
  if any(coarse_colours) and not all(coarse_colours):
    raise ValueError('Supply colours for every coarse label, or omit coarse_colour')

  # Explicit hex colours are portable between Python, R and installed models.
  for colour in colours + coarse_colours:
    if colour and not re.fullmatch(r'#[0-9a-fA-F]{6}(?:[0-9a-fA-F]{2})?', colour):
      raise ValueError(f'Invalid colour {colour!r}; use #RRGGBB or #RRGGBBAA')
  return records, coarse


def builtin_models():
  # Load the centrally maintained model manifests shipped with this checkout.
  root = Path(__file__).resolve().parents[1] / 'resources/xgboost'
  return {record['name']: record for path in sorted(root.glob('*.json'))
    for record in [json.loads(path.read_text())]}


def sha256(path):
  # Hash large model files without reading them fully into memory.
  with Path(path).open('rb') as stream:
    return hashlib.file_digest(stream, 'sha256').hexdigest()


def label_mapping(records):
  # Include only explicitly assigned coarse labels.
  return {row['label']: row['coarse_label'] for row in records if row['coarse_label']}


if __name__ == '__main__':
  parser = argparse.ArgumentParser(
    description = 'Train an XGBoost classifier from scprocess outputs')
  commands = parser.add_subparsers(dest = 'command', required = True)

  # The workflow supplies validated settings; defaults belong in the schemas.
  # Planning uses metadata only, before expression matrices are loaded.
  plan = commands.add_parser('plan')
  for name in ('annots_f', 'integration_f', 'h5ads_yaml', 'batch_var',
    'plan_f', 'coverage_f', 'counts_f', 'validation_group_var'):
    plan.add_argument('--' + name, required = True)

  for name in ('label_map_f', 'test_group_var'):
    plan.add_argument('--' + name)

  for name in ('train_cells_per_type', 'validation_cells_per_type', 'min_cells_per_type', 'seed'):
    plan.add_argument('--' + name, required = True, type = int)

  for name in ('cluster_majority_vote_min_purity', 'validation_fraction'):
    plan.add_argument('--' + name, required = True, type = float)

  for name in ('excluded_labels', 'test_groups', 'validation_variable_combinations',
    'sampling_variable_combinations'):
    plan.add_argument('--' + name, nargs = '*', default = [])

  plan.add_argument('--cluster_majority_vote_resolution',
    type = lambda value: float(value) if value else None)
  plan.add_argument('--is_join', type = int, choices = [0, 1], required = True)

  # Both training passes use the saved cell selection and these fitting settings.
  fit = commands.add_parser('train')
  for name in ('plan_f', 'h5ads_yaml', 'output_dir', 'xgboost_tag', 'gene_exclude_re'):
    fit.add_argument('--' + name, required = True)

  fit.add_argument('--label_map_f')

  for name in ('n_cores', 'seed', 'min_cells_expressed', 'min_genes', 'max_genes',
    'pass1_nrounds', 'pass1_early_stopping', 'pass2_nrounds', 'pass2_early_stopping'):
    fit.add_argument('--' + name, required = True, type = int)

  for name in ('pass1_subsample', 'pass1_colsample_bytree', 'pass1_learning_rate',
    'pass2_colsample_bytree', 'pass2_learning_rate', 'gain_threshold'):
    fit.add_argument('--' + name, required = True, type = float)

  fit.add_argument('--use_gpu', type = int, choices = [0, 1], required = True)

  # Full-data evaluation combines per-batch predictions with the original plan.
  aggregate = commands.add_parser('aggregate_fulldata')
  for name in ('plan_f', 'integration_f', 'output_f'):
    aggregate.add_argument('--' + name, required = True)

  aggregate.add_argument('pred_fs', nargs = '+')

  # Display specifications and coarse summaries can change without refitting.
  specs = commands.add_parser('specs')
  for name in ('cls_f', 'counts_f', 'predictions_f', 'output_dir', 'xgboost_tag'):
    specs.add_argument('--' + name, required = True)

  for name in ('label_specs_f', 'coarse_label_specs_f', 'label_map_f'):
    specs.add_argument('--' + name)

  args    = vars(parser.parse_args())
  command = args.pop('command')

  # Empty strings represent omitted optional features, not alternative defaults.
  for name in ('test_group_var', 'label_map_f', 'label_specs_f',
    'coarse_label_specs_f', 'gene_exclude_re'):
    if name in args and args[name] == '':
      args[name] = None

  if command == 'plan':
    prepare_plan(**args)
  elif command == 'train':
    fit_plan(**args)
  elif command == 'specs':
    write_label_specifications(**args)
  else:
    aggregate_plan_predictions(**args)
