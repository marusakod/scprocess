# Train an XGBoost classifier to predict cell type labels from scRNA-seq data.

import argparse
import gc
import hashlib
import json
from pathlib import Path
import re
import warnings
import numpy as np
import polars as pl
import yaml


# XGBoost, SciPy and AnnData imports stay local: setup also imports this module.

# Minimum training cells supporting a cluster's majority-vote label.
MIN_CLUSTER_MAJORITY_VOTE_CELLS = 10


# Training metadata

def prepare_training_metadata(annots_f, integration_f, h5ads_yaml, batch_var, training_metadata_f, coverage_f, counts_f,
  *, is_join, excluded_labels, label_map_f, cluster_majority_vote_resolution,
  cluster_majority_vote_min_purity, test_group_var, test_groups, validation_group_var,
  validation_fraction, validation_variable_combinations, sampling_variable_combinations,
  train_cells_per_label, validation_cells_per_label, min_cells_per_label, seed):

  # Read cell/batch IDs, join project IDs and optional voting clusters from integration_f.
  # Requested grouping/sampling metadata comes from annotations when available,
  # otherwise from integration_f.
  vote_res    = cluster_majority_vote_resolution
  vote_col    = f'RNA_snn_res.{float(vote_res):g}' if vote_res is not None else None
  cols_needed = set(validation_variable_combinations + sampling_variable_combinations)
  cols_needed |= {v for v in [test_group_var, validation_group_var] if v}

  # training_metadata_df has one row per cell present in both input tables, in integration order.
  training_metadata_df = _load_training_metadata(annots_f, integration_f, batch_var, is_join, vote_col, cols_needed)
  training_metadata_df = _check_training_batches(training_metadata_df, h5ads_yaml, batch_var, is_join, cols_needed)

  # Split whole groups before exclusions, optional voting and cell selection.
  rng = np.random.default_rng(seed)
  training_metadata_df, labels_ok, validation_groups_moved_to_training = _assign_splits(training_metadata_df, rng,
    test_group_var                   = test_group_var,
    test_groups                      = test_groups,
    validation_group_var             = validation_group_var,
    validation_fraction              = validation_fraction,
    validation_variable_combinations = validation_variable_combinations,
    excluded_labels                  = excluded_labels,
    min_cells_per_label              = min_cells_per_label)

  training_metadata_df = _prepare_training_labels(training_metadata_df, labels_ok, excluded_labels, vote_col,
    cluster_majority_vote_min_purity)
  training_metadata_df = _select_training_cells(training_metadata_df, rng, sampling_variable_combinations,
    train_cells_per_label, validation_cells_per_label)

  # Record the explicit preparation settings alongside the selected-cell tables.
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
    train_cells_per_label            = train_cells_per_label,
    validation_cells_per_label       = validation_cells_per_label,
    min_cells_per_label              = min_cells_per_label,
    seed                             = seed)
  _write_training_metadata(training_metadata_df, training_metadata_f, coverage_f, counts_f, configuration, validation_groups_moved_to_training)

  return training_metadata_df


def _load_training_metadata(annots_f, integration_f, batch_var, is_join, vote_col, cols_needed):
  # Read labels and only the integration columns required for preparation.
  integration_columns = ['cell_id', batch_var]
  if is_join:
    integration_columns.append('project_id')
  if vote_col:
    integration_columns.append(vote_col)

  # labels_df: one row per annotated cell, including any supplied cell metadata.
  labels_df = pl.read_csv(annots_f)

  # Accept one standard or legacy label column, then use the canonical name.
  label_columns = [name for name in ['label', 'annotation', 'fine_label'] if name in labels_df.columns]
  if len(label_columns) != 1:
    raise ValueError('annots_f needs one label column (label; legacy annotation/fine_label accepted)')
  labels_df = labels_df.rename({label_columns[0]: 'label'})
  if _missing_values(labels_df, ['cell_id']) or labels_df['cell_id'].n_unique() != len(labels_df):
    raise ValueError('Annotations must have unique nonmissing cell_id values')

  # Grouping metadata may come from either input table.
  header = pl.read_csv(integration_f, n_rows = 0).columns
  integration_columns += [
    v for v in sorted(cols_needed)
    if v in header and v not in labels_df.columns and v not in integration_columns]

  # int_df: one row per integrated cell; selecting columns does not aggregate cells.
  int_df = pl.read_csv(integration_f, columns = list(dict.fromkeys(integration_columns)))
  if _missing_values(int_df, ['cell_id']) or int_df['cell_id'].n_unique() != len(int_df):
    raise ValueError('Integration must have unique nonmissing cell_id values')

  # comparison_df has one row per shared cell, comparing metadata from both inputs.
  overlap = [v for v in labels_df.columns if v != 'cell_id' and v in int_df.columns]
  if overlap:
    comparison_df = int_df.join(
      labels_df.select(['cell_id'] + overlap),
      on     = 'cell_id',
      suffix = '_annotation')
    for var in overlap:
      conflict = (pl.col(var).is_not_null() & pl.col(var + '_annotation').is_not_null()
        & (pl.col(var).cast(pl.String) != pl.col(var + '_annotation').cast(pl.String)))
      if comparison_df.filter(conflict).height:
        raise ValueError(f'Conflicting integration and annotation values: {var}')
    labels_df = labels_df.drop(overlap)

  # Keep integration order; the matrix loader later groups selected rows by batch.
  training_metadata_df = int_df.join(
    labels_df,
    on             = 'cell_id',
    how            = 'inner',
    maintain_order = 'left')
  if training_metadata_df.is_empty():
    raise ValueError('No overlapping original labels and integration cells')

  return training_metadata_df


def _check_training_batches(training_metadata_df, h5ads_yaml, batch_var, is_join, cols_needed):
  # Joined integrations may already have project-qualified batch identifiers.
  h5ad_lookup     = yaml.safe_load(Path(h5ads_yaml).read_text())
  batch_col       = training_metadata_df[batch_var].cast(pl.String)
  missing_batches = set(batch_col) - set(h5ad_lookup)
  if is_join and missing_batches:
    batch_col       = training_metadata_df['project_id'].cast(pl.String) + '_' + batch_col
    missing_batches = set(batch_col) - set(h5ad_lookup)
  if missing_batches:
    raise ValueError('Input batch identifiers do not match H5AD manifest')
  training_metadata_df = training_metadata_df.with_columns(batch_col.alias('batch'))

  # Check grouping metadata before any split or sampling decisions.
  required     = cols_needed - {'batch'}
  cols_missing = required - set(training_metadata_df.columns)
  if cols_missing:
    raise ValueError(f'Missing grouping/sampling columns: {sorted(cols_missing)}')
  if _missing_values(training_metadata_df, required):
    raise ValueError('Missing values in grouping/sampling columns')

  return training_metadata_df


def _assign_splits(training_metadata_df, rng, *, test_group_var, test_groups, validation_group_var,
  validation_fraction, validation_variable_combinations, excluded_labels,
  min_cells_per_label):

  # Reserve explicitly held-out groups for the independent test set.
  training_metadata_df = training_metadata_df.with_columns(pl.lit('unassigned').alias('split'))
  if bool(test_group_var) != bool(test_groups):
    raise ValueError('test_group_var and nonempty test_groups must be supplied together')
  if test_group_var:
    if _missing_values(training_metadata_df, [test_group_var]):
      raise ValueError(f'Missing test grouping values: {test_group_var}')
    groups  = training_metadata_df[test_group_var].cast(pl.String)
    missing = set(test_groups) - set(groups)
    if missing:
      raise ValueError(f'Test groups absent from annotations/integration: {sorted(missing)}')
    training_metadata_df = training_metadata_df.with_columns(
      pl.when(groups.is_in(test_groups)).then(pl.lit('test'))
        .otherwise(pl.col('split')).alias('split'))

  # Split whole validation_group_var units (e.g. tissue tubes), keeping their cells together.
  # No unit may have cells in both the held-out test set and the train/validation pool.
  # Each remaining unit must also have one consistent combination of the requested
  # validation_variable_combinations values, so it belongs to exactly one sampling group.
  # Both filtered tables still have one row per cell; no aggregation occurs here.
  train_validation_df = training_metadata_df.filter(pl.col('split') != 'test')
  group_var           = validation_group_var
  test_df             = training_metadata_df.filter(pl.col('split') == 'test')
  if set(train_validation_df[group_var]) & set(test_df[group_var]):
    raise ValueError('Validation groups overlap the held-out test groups')
  variables = validation_variable_combinations
  if _missing_values(train_validation_df, [group_var] + variables):
    raise ValueError('Missing validation grouping/combination values')

  # group_units_df has one row per distinct (validation unit, metadata combination),
  # rather than per cell. The check below requires exactly one combination per unit.
  # For example,
  # validation_group_var = 'tube_id' and validation_variable_combinations = ['donor', 'region']:
  #   tube_id  donor  region
  #   tube1    d1     cortex
  #   tube2    d1     cortex
  #   tube3    d2     cerebellum
  #   tube4    d2     cerebellum
  # Validation selection below samples whole tubes within each (donor, region) combination.

  # dict.fromkeys removes duplicate column names in order: ['tube_id', 'donor', 'tube_id'] -> ['tube_id', 'donor'].
  group_units_df = train_validation_df.select(list(dict.fromkeys([group_var] + variables))).unique(maintain_order = True)
  if group_units_df[group_var].n_unique() != len(group_units_df):
    raise ValueError('A validation group spans several variable combinations; choose compatible variables')

  # Select whole validation units within each variable combination.
  validation = set()
  groups     = _sorted_groups(group_units_df, variables)
  for validation_units_df in groups:
    # validation_units_df contains validation-unit rows sharing the same metadata combination.
    # Apply the fraction to validation units (e.g. tubes), not cells, rounding down.
    # Select at least one unit where possible, leaving at least one for training;
    # a combination with only one unit contributes none to validation.
    n = min(len(validation_units_df) - 1, max(1, int(len(validation_units_df) * validation_fraction)))
    if n > 0:
      validation.update(rng.choice(validation_units_df[group_var].to_numpy(), n, replace = False))
  if not validation:
    raise ValueError('No whole groups available for validation; supply more groups or different validation combinations')

  # Move whole validation groups to training until rare labels have enough cells.
  eligible_df = train_validation_df.filter(pl.col('label').is_not_null() & ~pl.col('label').is_in(excluded_labels))

  # eligible_df is cell-level; label_counts_df has one row per label and its cell count.
  label_counts_df                     = eligible_df.group_by('label').len()
  labels_ok                           = set(label_counts_df.filter(pl.col('len') >= min_cells_per_label)['label'])
  validation_groups_moved_to_training = []
  for label in sorted(labels_ok):
    # Here label_cells_df is cell-level: all train/validation-pool cells with this label.
    label_cells_df = train_validation_df.filter(pl.col('label') == label)

    # Count this label's cells outside validation; these are the available training cells.
    while label_cells_df[group_var].is_in(list(validation)).not_().sum() < min_cells_per_label:
      # group_sizes_df: one row per validation unit, counting cells of this label.
      group_sizes_df = (label_cells_df.filter(pl.col(group_var).is_in(list(validation)))
        .group_by(group_var).len()
        .sort(['len', group_var], descending = [True, False]))
      if group_sizes_df.is_empty():
        raise ValueError(f'Unable to retain training support for {label}')

      # Move the unit with the most cells; break ties by validation-unit ID.
      unit = group_sizes_df[group_var][0]
      validation.remove(unit)
      validation_groups_moved_to_training.append({'label': label, 'group': str(unit)})

  # Apply the final group assignments after moving any groups back to training.
  training_metadata_df = training_metadata_df.with_columns(
    pl.when(pl.col('split') == 'test').then(pl.lit('test'))
      .when(pl.col(group_var).is_in(list(validation))).then(pl.lit('validation'))
      .otherwise(pl.lit('train')).alias('split'))
  if not training_metadata_df.filter(pl.col('split') == 'validation').height:
    raise ValueError('No validation groups remain after rescuing rare labels')
  return training_metadata_df, labels_ok, validation_groups_moved_to_training


def _prepare_training_labels(training_metadata_df, labels_ok, excluded_labels, vote_col,
  cluster_majority_vote_min_purity):

  # Preserve reference labels and flag excluded, unlabelled or unsupported cells.
  training_metadata_df = training_metadata_df.with_columns(
    pl.col('label').alias('original_label'),
    pl.when(pl.col('label').is_null()).then(pl.lit('unlabeled'))
      .when(pl.col('label').is_in(excluded_labels)).then(pl.lit('excluded'))
      .when(~pl.col('label').is_in(list(labels_ok)) & pl.col('split').is_in(['train', 'validation']))
      .then(pl.lit('unsupported')).otherwise(pl.col('split')).alias('split'))

  # Optional voting uses training labels only, leaving the test labels untouched.
  if vote_col:
    # cluster_counts_df: one row per (cluster, label), with its training-cell count.
    cluster_counts_df = training_metadata_df.filter(pl.col('split') == 'train').group_by([vote_col, 'label']).len()

    # cluster_labels_df: one row per qualifying cluster, with its majority label.
    cluster_labels_df = (cluster_counts_df.with_columns(pl.col('len').sum().over(vote_col).alias('total'))
      .sort([vote_col, 'len', 'label'], descending = [False, True, False])
      .unique(vote_col, keep = 'first', maintain_order = True)
      .filter((pl.col('len') / pl.col('total') >= cluster_majority_vote_min_purity) & (pl.col('total') >= MIN_CLUSTER_MAJORITY_VOTE_CELLS)))
    mapping = dict(cluster_labels_df.select(vote_col, 'label').iter_rows())
    training_metadata_df = training_metadata_df.with_columns(
      pl.when(pl.col('split').is_in(['train', 'validation']))
        .then(pl.col(vote_col).replace_strict(mapping, default = pl.col('label')))
        .otherwise(pl.col('label')).alias('label'))

  return training_metadata_df


def _select_training_cells(training_metadata_df, rng, sampling_variable_combinations,
  train_cells_per_label, validation_cells_per_label):

  # training_metadata_df, eligible_df, label_split_cells_df and sample_df all have one row per cell here.
  # Sample independently within each label and split, retaining cell row order.
  training_metadata_df       = training_metadata_df.with_row_index('_tmp_row_idx')
  selected_rows = []
  eligible_df   = training_metadata_df.filter(pl.col('split').is_in(['train', 'validation']))
  for label_split_cells_df in _sorted_groups(eligible_df, ['label', 'split']):
    cap       = train_cells_per_label if label_split_cells_df['split'][0] == 'train' else validation_cells_per_label
    sample_df = _balanced_sample(label_split_cells_df, cap, sampling_variable_combinations, rng)
    selected_rows.extend(sample_df['_tmp_row_idx'].to_list())
  training_metadata_df = training_metadata_df.with_columns(pl.col('_tmp_row_idx').is_in(selected_rows).alias('selected')).drop('_tmp_row_idx')

  # Require both splits and training support for every validation label.
  selected_df = training_metadata_df.filter(pl.col('selected'))
  if not {'train', 'validation'} <= set(selected_df['split']):
    raise ValueError('Both training and validation cells are required')
  train_labels = set(selected_df.filter(pl.col('split') == 'train')['label'])
  if not set(selected_df.filter(pl.col('split') == 'validation')['label']) <= train_labels:
    raise ValueError('Validation label absent from training')
  training_metadata_df = training_metadata_df.with_columns(
    pl.col('label').is_in(list(train_labels)).fill_null(False).alias('label_supported_in_training'))

  return training_metadata_df


def _write_training_metadata(training_metadata_df, training_metadata_f, coverage_f, counts_f, configuration, validation_groups_moved_to_training):
  selected_df  = training_metadata_df.filter(pl.col('selected'))
  train_labels = set(selected_df.filter(pl.col('split') == 'train')['label'])

  # Save the training metadata, selection coverage and decisions needed for provenance.
  training_metadata_df.write_csv(training_metadata_f, compression = 'gzip')

  # coverage_df: one row per (label, split), with total and selected cell counts.
  coverage_df = (training_metadata_df.group_by(['label', 'split']).agg(
    pl.len().alias('n_cells'), pl.col('selected').sum().alias('n_selected'))
    .sort(['label', 'split'], nulls_last = True))
  coverage_df.write_csv(coverage_f)

  # training_counts_df: one row per label, counting selected training cells.
  training_counts_df = (selected_df.filter(pl.col('split') == 'train')
    .group_by('label').len().sort('label'))
  training_counts_df.write_csv(counts_f)

  # Record preparation settings and any validation groups moved back to training.
  metadata = dict(
    configuration               = configuration,
    rescued_validation_groups   = validation_groups_moved_to_training,
    n_annotated                 = len(training_metadata_df),
    n_train                     = selected_df.filter(pl.col('split') == 'train').height,
    n_validation                = selected_df.filter(pl.col('split') == 'validation').height,
    n_test                      = training_metadata_df.filter(pl.col('split') == 'test').height,
    n_classes                   = len(train_labels),
    majority_vote_training_only = True)
  Path(str(training_metadata_f) + '.json').write_text(json.dumps(metadata, indent = 2) + '\n')
  print(coverage_df, flush = True)


def _balanced_sample(input_df, n, variables, rng):

  # Keep all cells below the cap; sample uniformly when no balancing is requested.
  if len(input_df) <= n:
    return input_df.clone()
  if not variables:
    return input_df[rng.choice(len(input_df), n, replace = False)]

  # Balance joint combinations, redistributing unused quotas to larger groups.
  # For groups with 2, 10 and 10 cells and a cap of 8, quotas become 2, 3 and 3.
  # Randomise group priority so any remainder is not always given to the first group.
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
  # Partition rows into lexicographically ordered combinations, preserving row order
  # within each table. This does not aggregate: cell-level input stays cell-level,
  # while validation-unit input stays one row per unit.
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

def fit_classifier(training_metadata_f, h5ads_yaml, output_dir, n_cores, *, xgboost_tag, use_gpu,
  gene_exclude_re, min_cells_expressed, seed, pass1_subsample, pass1_colsample_bytree,
  pass1_learning_rate, pass1_nrounds, pass1_early_stopping, gain_threshold, min_genes,
  max_genes, pass2_colsample_bytree, pass2_learning_rate, pass2_nrounds,
  pass2_early_stopping, label_map_f):

  import xgboost as xgb

  output_path = Path(output_dir)
  output_path.mkdir(parents = True, exist_ok = True)
  gpu_device = _verify_gpu(n_cores) if use_gpu else 'cpu'

  # training_metadata_df remains one row per cell, filtered to selected training/validation cells.
  training_metadata_df     = pl.read_csv(training_metadata_f, schema_overrides = {'label': pl.String})
  training_metadata_df     = training_metadata_df.filter(pl.col('selected') & pl.col('split').is_in(['train', 'validation']))
  h5ad_lookup = yaml.safe_load(Path(h5ads_yaml).read_text())
  raw_counts, normalized_counts, gene_ids, gene_types, training_metadata_df = _load_training_expression(training_metadata_df, h5ad_lookup)

  # Encode labels and filter genes before the two fitting passes.
  class_names, label_indices, train_mask, validation_mask = _encode_training_labels(training_metadata_df)
  normalized_counts, raw_counts, gene_names = _filter_training_genes(
    normalized_counts, raw_counts, gene_ids, gene_types, train_mask,
    gene_exclude_re = gene_exclude_re,
    min_cells_expressed = min_cells_expressed)
  X_train      = normalized_counts[train_mask]
  X_validation = normalized_counts[validation_mask]

  # Pass 1 considers all retained genes with the lower column-sampling fraction.
  pass1_model = _fit_xgboost(
    X_train          = X_train,
    y_train          = label_indices[train_mask],
    X_val            = X_validation,
    y_val            = label_indices[validation_mask],
    n_cores          = n_cores,
    seed             = seed,
    use_gpu          = use_gpu,
    subsample        = pass1_subsample,
    colsample_bytree = pass1_colsample_bytree,
    learning_rate    = pass1_learning_rate,
    nrounds          = pass1_nrounds,
    early_stopping   = pass1_early_stopping)
  pass1_best_iteration = int(pass1_model.best_iteration)
  pass1_model          = _best_iteration_model(pass1_model)

  # Retain genes contributing the requested cumulative gain, within size bounds.
  selected_gene_indices, selected_gene_names, gene_importance_df = _select_features_by_gain(
    model          = pass1_model,
    gene_names     = gene_names,
    gain_threshold = gain_threshold,
    min_genes      = min_genes,
    max_genes      = max_genes)
  del pass1_model, normalized_counts
  gc.collect()
  X_train      = X_train[:, selected_gene_indices].tocsr()
  X_validation = X_validation[:, selected_gene_indices].tocsr()

  # Pass 2 fits the final classifier on those genes, stopping on validation loss.
  pass2_model = _fit_xgboost(
    X_train          = X_train,
    y_train          = label_indices[train_mask],
    X_val            = X_validation,
    y_val            = label_indices[validation_mask],
    n_cores          = n_cores,
    seed             = seed,
    use_gpu          = use_gpu,
    subsample        = pass1_subsample,
    colsample_bytree = pass2_colsample_bytree,
    learning_rate    = pass2_learning_rate,
    nrounds          = pass2_nrounds,
    early_stopping   = pass2_early_stopping,
    feature_names    = [f"g{i}" for i in range(X_train.shape[1])])
  pass2_best_iteration = int(pass2_model.best_iteration)
  pass2_model          = _best_iteration_model(pass2_model)

  # Export the best iteration, selected-cell predictions and raw-count summaries.
  pseudobulk_adata = _compute_pseudobulk(raw_counts, training_metadata_df, gene_names, selected_gene_indices)
  _save_outputs(
    model              = pass2_model,
    class_names        = class_names,
    selected_genes     = selected_gene_names,
    gene_importance_df = gene_importance_df,
    output_dir         = output_path,
    xgboost_tag        = xgboost_tag,
    cells_df           = training_metadata_df,
    label_map          = _load_label_mapping(label_map_f),
    X_train            = X_train,
    X_val              = X_validation,
    pseudobulk_adata   = pseudobulk_adata)

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
  model_f    = output_path / f'{xgboost_tag}_xgboost_model.json'
  model_card = dict(
    configuration                  = configuration,
    n_train                        = int(train_mask.sum()),
    n_validation                   = int(validation_mask.sum()),
    n_classes                      = len(class_names),
    n_genes                        = len(selected_gene_names),
    n_genes_tested                 = len(gene_names),
    pass1_best_iteration           = pass1_best_iteration,
    pass2_best_iteration           = pass2_best_iteration,
    detection_filter_training_only = True,
    best_iteration_model_sliced    = True,
    model_sha256                   = sha256(model_f),
    training_metadata_sha256       = sha256(training_metadata_f),
    xgboost_version                = xgb.__version__,
    xgboost_build                  = xgb.build_info(),
    training_device                = gpu_device,
    normalization                  = 'log1p(10000 * counts / all-input-gene library size)',
    status                         = 'Unreviewed candidate; not registered or published')
  (output_path / f'{xgboost_tag}_model_card.json').write_text(json.dumps(model_card, indent = 2) + '\n')


def _encode_training_labels(training_metadata_df):

  # training_metadata_df is cell-level; use the same class indices in both fitting splits.
  train_mask      = training_metadata_df['split'].to_numpy() == 'train'
  validation_mask = ~train_mask
  class_names     = sorted(training_metadata_df.filter(train_mask)['label'].unique())
  if not set(training_metadata_df.filter(validation_mask)['label']) <= set(class_names):
    raise ValueError('Validation label absent from training')
  class_indices = {name: i for i, name in enumerate(class_names)}
  label_indices = np.array([class_indices[name] for name in training_metadata_df['label']])
  return class_names, label_indices, train_mask, validation_mask


def _filter_training_genes(normalized_counts, raw_counts, gene_ids, gene_types,
  train_mask, *, gene_exclude_re, min_cells_expressed):

  # Remove excluded gene types, keeping gene IDs aligned with matrix columns.
  if gene_exclude_re:
    keep      = gene_types.str.contains(gene_exclude_re).not_()
    keep_mask = keep.to_numpy()
    n_before  = len(gene_ids)
    n_after   = int(keep.sum())
    normalized_counts = normalized_counts[:, keep_mask]
    raw_counts        = raw_counts[:, keep_mask]
    gene_ids          = gene_ids.filter(keep)
    print(f"  Gene type filter (exclude '{gene_exclude_re}'): {n_before} → {n_after} genes "
          f"({n_before - n_after} removed)")

  # Detection uses training cells only; held-out cells cannot select genes.
  detected            = np.asarray(raw_counts[train_mask].getnnz(axis = 0)).ravel()
  expressed_gene_mask = detected >= min_cells_expressed
  normalized_counts   = normalized_counts[:, expressed_gene_mask].tocsr()
  raw_counts          = raw_counts[:, expressed_gene_mask].tocsr()
  gene_ids            = gene_ids.filter(expressed_gene_mask)
  if gene_ids.is_empty():
    raise ValueError('No genes pass training-only filtering')

  # Model feature selection and pseudobulk output helpers consume gene-name lists.
  gene_names = gene_ids.to_list()
  return normalized_counts, raw_counts, gene_names


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


def _load_training_expression(training_metadata_df, h5ad_lookup):
  import scipy.sparse as sp

  # Read batches in deterministic order and normalise before gene filtering.
  raw_count_parts, normalized_count_parts, batch_cells_dfs = [], [], []
  gene_ids = gene_types = None
  for batch_cells_df in _sorted_groups(training_metadata_df, ['batch']):
    batch_id = batch_cells_df['batch'][0]
    entry    = h5ad_lookup[batch_id]
    h5ad_f   = entry['path'] if isinstance(entry, dict) else entry
    raw_counts, batch_gene_ids, batch_gene_types = _read_selected_h5ad(
      h5ad_f, batch_cells_df['cell_id'].to_list(), gene_ids)
    if gene_ids is None:
      gene_ids, gene_types = batch_gene_ids, batch_gene_types
    raw_count_parts.append(raw_counts)
    normalized_count_parts.append(normalize_counts(raw_counts).astype(np.float32))
    batch_cells_dfs.append(batch_cells_df)
    print(f'Loading {batch_id}: {len(batch_cells_df)} selected cells', flush = True)

  # Stack count matrices and cell-level metadata in the same order.
  raw_counts = sp.vstack(raw_count_parts, format = 'csr')
  raw_count_parts.clear()
  normalized_counts = sp.vstack(normalized_count_parts, format = 'csr')
  normalized_count_parts.clear()
  cells_df = pl.concat(batch_cells_dfs)

  # Gene metadata stays as Series through filtering, aligned to matrix columns.
  gene_ids   = pl.Series('gene_id', gene_ids, dtype = pl.String)
  gene_types = pl.Series('gene_type', gene_types, dtype = pl.String)
  return raw_counts, normalized_counts, gene_ids, gene_types, cells_df


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


def _fit_xgboost(X_train, y_train, X_val, y_val, *, n_cores, seed, use_gpu,
  subsample, colsample_bytree, learning_rate, nrounds, early_stopping,
  feature_names = None):

  import xgboost as xgb

  params = _build_xgb_params(n_cores, seed, use_gpu, int(y_train.max()) + 1,
    subsample        = subsample,
    colsample_bytree = colsample_bytree,
    learning_rate    = learning_rate)

  # Both passes stop on validation log loss; only pass two needs saved feature names.
  dtrain = xgb.DMatrix(X_train, label = y_train, feature_names = feature_names)
  dval   = xgb.DMatrix(X_val, label = y_val, feature_names = feature_names)
  model  = xgb.train(
    params                = params,
    dtrain                = dtrain,
    num_boost_round       = nrounds,
    evals                 = [(dtrain, "train"), (dval, "val")],
    early_stopping_rounds = early_stopping,
    verbose_eval          = 50)

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

  if not len(sorted_idx):
    raise ValueError('First-pass XGBoost model has no genes with positive feature gain')

  total_gain = sorted_gains.sum()
  cum_gain   = np.cumsum(sorted_gains) / total_gain

  # Find the cumulative-gain cutoff and apply the configured size bounds.
  n_sel = int(np.searchsorted(cum_gain, gain_threshold)) + 1
  n_sel = max(n_sel, min_genes)
  n_sel = min(n_sel, max_genes)
  n_sel = min(n_sel, len(sorted_idx))

  sel_indices    = sorted_idx[:n_sel]
  sel_gene_names = [gene_names[i] for i in sel_indices]

  # gene_importance_df: one row per gene with positive gain, in descending gain order.
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
  output_dir, xgboost_tag, cells_df, X_train, X_val, pseudobulk_adata, label_map = None):

  # Save model files and the metadata required to use the classifier.
  out_dir = Path(output_dir)

  model_path = out_dir / f"{xgboost_tag}_xgboost_model.json"
  model.save_model(str(model_path))
  print(f"  Model: {model_path}")

  cls_path   = out_dir / f"{xgboost_tag}_allowed_cls.csv"
  classes_df = pl.DataFrame({"class": class_names})
  classes_df.write_csv(cls_path)
  print(f"  Classes: {cls_path}")

  genes_path = out_dir / f"{xgboost_tag}_selected_genes.txt"
  with open(genes_path, "w") as f:
    f.write("\n".join(selected_genes) + "\n")
  print(f"  Genes: {genes_path}")

  imp_path = out_dir / f"{xgboost_tag}_gene_importance.csv"
  gene_importance_df.write_csv(str(imp_path))
  print(f"  Importance: {imp_path}")

  if label_map is not None:
    map_path     = out_dir / f"{xgboost_tag}_label_map.csv"
    label_map_df = pl.DataFrame({
      "label": list(label_map.keys()),
      "coarse_label": list(label_map.values()),
    })
    label_map_df.write_csv(map_path)
    print(f"  Label map: {map_path}")

  # The fitting workflow always supplies selected-cell predictions and pseudobulks.
  predictions_f  = out_dir / f"{xgboost_tag}_predictions.csv.gz"
  predictions_df = _make_predictions_df(
    model       = model,
    X_train     = X_train,
    X_val       = X_val,
    cells_df    = cells_df,
    class_names = class_names,
    label_map   = label_map)
  predictions_df.write_csv(predictions_f, compression = "gzip")
  print(f"  Predictions: {predictions_f}")

  pseudobulk_f = out_dir / f"{xgboost_tag}_pseudobulk.h5ad"
  pseudobulk_adata.write_h5ad(pseudobulk_f)
  print(f"  Pseudobulk: {pseudobulk_f}")


def _make_predictions_df(model, X_train, X_val, cells_df, class_names, label_map = None):
  import xgboost as xgb

  feat_names = [f"g{i}" for i in range(X_val.shape[1])]

  # Score each split separately, then restore the interleaved selected-cell order.
  dval      = xgb.DMatrix(X_val, feature_names = feat_names)
  val_probs = model.predict(dval)
  val_pred  = val_probs.argmax(axis = 1)
  val_pmax  = val_probs.max(axis = 1)

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

  # AnnData metadata has one row per batch (obs_df) or selected gene (var_df).
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

def combine_predictions(pred_fs, training_metadata_f, output_f, integration_f):

  # Prediction tables have one row per cell; concatenation preserves that row unit.
  # Read consistent types, including empty per-batch prediction files.
  prediction_columns = ['cell_id', 'predicted_label_naive', 'probability_naive']
  prediction_schema  = {'cell_id': pl.String, 'predicted_label_naive': pl.String,
    'probability_naive': pl.Float64}
  prediction_tables = []
  for path in pred_fs:
    batch_predictions_df = pl.read_csv(path,
      columns          = prediction_columns,
      schema_overrides = prediction_schema)
    prediction_tables.append(batch_predictions_df)
  predictions_df = pl.concat(prediction_tables)
  if predictions_df['cell_id'].n_unique() != len(predictions_df):
    raise ValueError('Duplicate full-data predictions')
  predictions_df = predictions_df.rename({
    'predicted_label_naive': 'predicted_label',
    'probability_naive': 'probability'})

  # Require one prediction for every integrated cell, retaining prediction order.
  integration_cells_df = pl.read_csv(integration_f, columns = ['cell_id'])
  predictions_df       = predictions_df.join(integration_cells_df,
    on             = 'cell_id',
    how            = 'inner',
    maintain_order = 'left')
  if len(predictions_df) != len(integration_cells_df):
    raise ValueError('Integrated cells missing from full-data predictions')

  # training_metadata_df is cell-level; the left join keeps result_df one row per integrated
  # cell, including cells absent from the training metadata because they have no reference label.
  training_metadata_df = pl.read_csv(training_metadata_f,
    schema_overrides = {'label': pl.String, 'original_label': pl.String})
  training_metadata_columns = ['cell_id', 'label', 'original_label', 'split', 'selected',
    'label_supported_in_training']
  result_df = predictions_df.join(training_metadata_df.select(training_metadata_columns),
    on             = 'cell_id',
    how            = 'left',
    maintain_order = 'left')
  result_df = result_df.with_columns(
    pl.col('split').fill_null('unlabeled'),
    pl.col('selected').fill_null(False),
    pl.col('label_supported_in_training').fill_null(False))
  if not set(training_metadata_df['cell_id']) <= set(result_df['cell_id']):
    raise ValueError('Training-metadata cells missing from full-data predictions')

  # Training metadata distinguishes an omitted mapping from a lost mapping setting.
  metadata_f = Path(str(training_metadata_f) + '.json')
  if not metadata_f.is_file():
    raise FileNotFoundError(
      f'Training metadata JSON missing: {metadata_f}; cannot determine whether coarse-label mapping was requested')
  configuration = json.loads(metadata_f.read_text())['configuration']
  if configuration.get('label_map_f'):
    label_map = _load_label_mapping(configuration['label_map_f'])
    result_df = result_df.with_columns(
      pl.col('label').replace(label_map).alias('coarse_true'),
      pl.col('predicted_label').replace(label_map).alias('coarse_predicted'))
  result_df.write_csv(output_f, compression = 'gzip')

  # Score unselected development cells separately from fitting/validation cells.
  prefix = str(output_f).removesuffix('_fulldata_predictions.csv.gz')
  _write_evaluation(_evaluation_counts(result_df), prefix, 'label')


def write_label_specifications(cls_f, counts_f, predictions_f, output_dir, *,
  xgboost_tag, label_specs_f, coarse_label_specs_f, label_map_f):

  # One row per model class or training label; abundance sets the default order.
  classes_df      = pl.read_csv(cls_f)
  label_counts_df = pl.read_csv(counts_f)
  class_names     = classes_df['class'].to_list()
  label_counts    = dict(label_counts_df.select('label', 'len').iter_rows())
  default_label_order = sorted(class_names,
    key = lambda label: (-label_counts.get(label, 0), label))

  # Supplied specification tables override the default display order.
  fine_specs, coarse_specs = read_label_specs(
    label_specs_f        = label_specs_f,
    coarse_label_specs_f = coarse_label_specs_f,
    labels               = default_label_order,
    label_map_f          = label_map_f)
  output_path = Path(output_dir)
  output_path.mkdir(parents = True, exist_ok = True)
  output_prefix = output_path / xgboost_tag

  # Write one row per fine or coarse label, in its specified display order.
  fine_specs_df = pl.DataFrame(fine_specs,
    schema = {'label': pl.String, 'colour': pl.String, 'coarse_label': pl.String})
  coarse_specs_df = pl.DataFrame(coarse_specs,
    schema = {'coarse_label': pl.String, 'coarse_colour': pl.String})
  fine_specs_df.write_csv(str(output_prefix) + '_label_specs.csv')
  coarse_specs_df.write_csv(str(output_prefix) + '_coarse_label_specs.csv')

  # Read cell-level predictions lazily; only the grouped counts are collected.
  # Coarse labels are a direct lookup of fine labels, without another vote.
  label_map = label_mapping(fine_specs)
  if label_map:
    predictions_df = pl.scan_csv(predictions_f,
      schema_overrides = {'label': pl.String, 'predicted_label': pl.String})
    evaluation_df = _evaluation_counts(predictions_df, label_map = label_map)
  else:
    # No mapping is a valid choice; write empty coarse-evaluation tables.
    evaluation_df = pl.DataFrame(schema = {
      'split': pl.String, 'label': pl.String,
      'predicted_label': pl.String, 'len': pl.UInt32})
  _write_evaluation(
    counts_df         = evaluation_df,
    prefix            = str(output_prefix) + '_coarse',
    label_column_name = 'coarse_label')


def _evaluation_counts(predictions_df, label_map = None):

  # Input has one row per cell. Keep processing lazy until the final count table.
  if isinstance(predictions_df, pl.DataFrame):
    predictions_df = predictions_df.lazy()
  if 'selected' not in predictions_df.collect_schema():
    raise ValueError('Evaluation predictions need a selected column to distinguish fitted and unselected cells')

  # Score only cells with reference labels that are eligible for evaluation.
  has_reference_label = pl.col('label').is_not_null()
  is_included_split   = ~pl.col('split').is_in(['excluded', 'unsupported', 'unlabeled'])
  predictions_df = predictions_df.filter(has_reference_label & is_included_split)

  # Separate unselected train/validation cells from cells used during fitting.
  # Independent test cells keep their test assignment regardless of selection.
  is_holdout = ~pl.col('selected') & pl.col('split').is_in(['train', 'validation'])
  predictions_df = predictions_df.with_columns(
    pl.when(is_holdout)
      .then(pl.lit('development_holdout'))
      .otherwise(pl.col('split'))
      .alias('split'))

  # Optional coarse evaluation maps both reference and predicted labels directly.
  if label_map:
    predictions_df = predictions_df.with_columns(
      pl.col('label').replace(label_map),
      pl.col('predicted_label').replace(label_map))

  # Output has one row per (split, reference label, predicted label), with a cell count.
  counts_df = predictions_df.group_by(['split', 'label', 'predicted_label']).len()
  return counts_df.collect(engine = 'streaming')


def _write_evaluation(counts_df, prefix, label_column_name):

  # counts_df and split_df contain counted label pairs, not individual cells.
  metrics   = []
  per_label = []
  for split_df in _sorted_groups(counts_df, ['split']):
    split = split_df['split'][0]

    # Each count table has one row per label: reference, predicted or correct cells.
    reference_counts_df = split_df.group_by('label').agg(pl.col('len').sum().alias('n_reference'))
    predicted_counts_df = (split_df.group_by('predicted_label')
      .agg(pl.col('len').sum().alias('n_predicted'))
      .rename({'predicted_label': 'label'}))
    correct_counts_df = (split_df.filter(pl.col('label') == pl.col('predicted_label'))
      .group_by('label').agg(pl.col('len').sum().alias('n_correct')))

    # Include labels seen only in predictions, with zero reference support.
    label_totals_df = pl.concat([
      split_df.select('label'),
      split_df.select(pl.col('predicted_label').alias('label'))]).unique().sort('label')
    label_totals_df = (label_totals_df
      .join(reference_counts_df, on = 'label', how = 'left', maintain_order = 'left')
      .join(predicted_counts_df, on = 'label', how = 'left', maintain_order = 'left')
      .join(correct_counts_df, on = 'label', how = 'left', maintain_order = 'left')
      .fill_null(0))

    # Compute per-label metrics from the grouped counts, keeping formulas explicit.
    scores = []
    for label, n_reference, n_predicted, n_correct in label_totals_df.iter_rows():
      precision = n_correct / max(n_predicted, 1)
      recall    = n_correct / max(n_reference, 1)
      f1        = 2 * precision * recall / (precision + recall) if precision + recall else 0.
      per_label.append(dict(
        split     = split, **{label_column_name: label},
        precision = precision, recall = recall, f1 = f1, n_test = n_reference))

      # Average over labels with reference cells in this split.
      if n_reference:
        scores.append((recall, f1))
    n_cells = int(split_df['len'].sum())
    metrics.append(dict(
      split             = split,
      n_cells           = n_cells,
      accuracy          = int(label_totals_df['n_correct'].sum()) / n_cells,
      macro_f1          = float(np.mean([score[1] for score in scores])),
      balanced_accuracy = float(np.mean([score[0] for score in scores]))))

  # Write stable table columns even when there are no evaluable cells.
  metric_schema = {'split': pl.String, 'n_cells': pl.Int64, 'accuracy': pl.Float64,
    'macro_f1': pl.Float64, 'balanced_accuracy': pl.Float64}
  label_schema = {'split': pl.String, label_column_name: pl.String, 'precision': pl.Float64,
    'recall': pl.Float64, 'f1': pl.Float64, 'n_test': pl.Int64}
  summary_df           = pl.DataFrame(metrics, schema = metric_schema)
  per_label_metrics_df = pl.DataFrame(per_label, schema = label_schema)
  summary_df.write_csv(str(prefix) + '_evaluation.csv')
  per_label_metrics_df.write_csv(str(prefix) + '_per_label_evaluation.csv')
  print(summary_df, flush = True)


# Model identity and label specifications

def get_xgboost_model_info(model_f, name):
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

  if label_specs_f and label_map_f:
    raise ValueError('Use label_specs_f or legacy label_map_f, not both')

  # Keep fine and coarse display order independent, with a direct label mapping.
  fine_specs_f = label_specs_f or label_map_f
  fine_specs   = _read_fine_label_specs(fine_specs_f, labels)
  coarse_names = list(dict.fromkeys(row['coarse_label'] for row in fine_specs if row['coarse_label']))
  coarse_specs = _read_coarse_label_specs(coarse_label_specs_f, coarse_names)

  # Explicit hex colours are portable between Python, R and installed models.
  colours = [row['colour'] for row in fine_specs] + [row['coarse_colour'] for row in coarse_specs]
  for colour in colours:
    if colour and not re.fullmatch(r'#[0-9a-fA-F]{6}(?:[0-9a-fA-F]{2})?', colour):
      raise ValueError(f'Invalid colour {colour!r}; use #RRGGBB or #RRGGBBAA')
  return fine_specs, coarse_specs


def _read_fine_label_specs(label_specs_f, labels):

  # Read all fields as strings, preserving label names and the supplied row order.
  fine_specs = []
  if label_specs_f:
    specs_df      = pl.read_csv(label_specs_f, infer_schema = False, raise_if_empty = False)
    label_columns = [name for name in ['label', 'annotation', 'fine_label'] if name in specs_df.columns]
    if len(label_columns) != 1:
      raise ValueError(f'{label_specs_f} needs one label column (label; legacy annotation/fine_label accepted)')
    label_col = label_columns[0]
    for row in specs_df.iter_rows(named = True):
      fine_specs.append({
        'label': (row.get(label_col) or '').strip(),
        'colour': (row.get('colour') or '').strip(),
        'coarse_label': (row.get('coarse_label') or '').strip()})
    if not fine_specs:
      raise ValueError(f'Empty label specifications: {label_specs_f}')
  elif labels is not None:
    fine_specs = [{'label': str(label), 'colour': '', 'coarse_label': ''} for label in labels]

  # Coverage is required; colours and coarse labels must be complete when present.
  label_names = [row['label'] for row in fine_specs]
  if any(not name for name in label_names) or len(set(label_names)) != len(label_names):
    raise ValueError('Label specifications contain missing or duplicate labels')
  if labels is not None and not set(labels) <= set(label_names):
    raise ValueError(f"Label specifications miss model labels: {sorted(set(labels) - set(label_names))}")
  colours = [row['colour'] for row in fine_specs]
  if any(colours) and not all(colours):
    raise ValueError('Supply colours for every label, or omit the colour column')
  coarse_labels = [row['coarse_label'] for row in fine_specs]
  if any(coarse_labels) and not all(coarse_labels):
    raise ValueError('Supply coarse_label for every label, or omit the coarse_label column')
  return fine_specs


def _read_coarse_label_specs(coarse_label_specs_f, coarse_names):

  # Without a separate table, retain first occurrence in the fine-label mapping.
  if not coarse_label_specs_f:
    return [{'coarse_label': name, 'coarse_colour': ''} for name in coarse_names]

  specs_df = pl.read_csv(coarse_label_specs_f, infer_schema = False, raise_if_empty = False)
  if 'coarse_label' not in specs_df.columns:
    raise ValueError('Coarse specifications need a coarse_label column')
  coarse_specs = []
  for row in specs_df.iter_rows(named = True):
    coarse_specs.append({
      'coarse_label': (row.get('coarse_label') or '').strip(),
      'coarse_colour': (row.get('coarse_colour') or '').strip()})

  # Every mapped group must be represented, with unique names and optional colours.
  if coarse_specs and not coarse_names:
    raise ValueError('Coarse specifications require a coarse_label mapping in label_specs_f')
  label_names = [row['coarse_label'] for row in coarse_specs]
  if any(not name for name in label_names) or len(label_names) != len(set(label_names)):
    raise ValueError('Coarse specifications contain missing or duplicate labels')
  if not set(coarse_names) <= set(label_names):
    raise ValueError(f"Coarse specifications miss groups: {sorted(set(coarse_names) - set(label_names))}")
  colours = [row['coarse_colour'] for row in coarse_specs]
  if any(colours) and not all(colours):
    raise ValueError('Supply colours for every coarse label, or omit coarse_colour')
  return coarse_specs


def builtin_models():
  # Load the centrally maintained model manifests shipped with this checkout.
  root   = Path(__file__).resolve().parents[1] / 'resources/xgboost'
  models = {}
  for path in sorted(root.glob('*.json')):
    model_manifest = json.loads(path.read_text())
    models[model_manifest['name']] = model_manifest
  return models


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
  # Preparation uses metadata only, before expression matrices are loaded.
  prepare = commands.add_parser('prepare')
  for name in ('annots_f', 'integration_f', 'h5ads_yaml', 'batch_var',
    'training_metadata_f', 'coverage_f', 'counts_f', 'validation_group_var'):
    prepare.add_argument('--' + name, required = True)

  for name in ('label_map_f', 'test_group_var'):
    prepare.add_argument('--' + name)

  for name in ('train_cells_per_label', 'validation_cells_per_label', 'min_cells_per_label', 'seed'):
    prepare.add_argument('--' + name, required = True, type = int)

  for name in ('cluster_majority_vote_min_purity', 'validation_fraction'):
    prepare.add_argument('--' + name, required = True, type = float)

  for name in ('excluded_labels', 'test_groups', 'validation_variable_combinations',
    'sampling_variable_combinations'):
    prepare.add_argument('--' + name, nargs = '*', default = [])

  prepare.add_argument('--cluster_majority_vote_resolution',
    type = lambda value: float(value) if value else None)
  prepare.add_argument('--is_join', type = int, choices = [0, 1], required = True)

  # Both training passes use the saved cell selection and these fitting settings.
  fit = commands.add_parser('train')
  for name in ('training_metadata_f', 'h5ads_yaml', 'output_dir', 'xgboost_tag', 'gene_exclude_re'):
    fit.add_argument('--' + name, required = True)

  fit.add_argument('--label_map_f')

  for name in ('n_cores', 'seed', 'min_cells_expressed', 'min_genes', 'max_genes',
    'pass1_nrounds', 'pass1_early_stopping', 'pass2_nrounds', 'pass2_early_stopping'):
    fit.add_argument('--' + name, required = True, type = int)

  for name in ('pass1_subsample', 'pass1_colsample_bytree', 'pass1_learning_rate',
    'pass2_colsample_bytree', 'pass2_learning_rate', 'gain_threshold'):
    fit.add_argument('--' + name, required = True, type = float)

  fit.add_argument('--use_gpu', type = int, choices = [0, 1], required = True)

  # Full-data evaluation combines per-batch predictions with the training metadata.
  combine = commands.add_parser('combine')
  for name in ('training_metadata_f', 'integration_f', 'output_f'):
    combine.add_argument('--' + name, required = True)

  combine.add_argument('pred_fs', nargs = '+')

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

  if command == 'prepare':
    prepare_training_metadata(**args)
  elif command == 'train':
    fit_classifier(**args)
  elif command == 'combine':
    combine_predictions(**args)
  else:
    write_label_specifications(**args)
