"""Plan group-disjoint XGBoost training before reading expression matrices."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import polars as pl
import yaml


def balanced_sample(frame, n, variables, rng):
    if len(frame) <= n:
        return frame.copy()
    if not variables:
        return frame.iloc[rng.choice(len(frame), n, replace=False)].copy()
    # Balance joint combinations; combinations exceeding the budget are sampled.
    groups = list(frame.groupby(variables, sort=True, dropna=False, observed=True).indices.values())
    order = rng.permutation(len(groups))
    quota = np.zeros(len(groups), dtype=int)
    remaining = n
    while remaining:
        active = [i for i in order if quota[i] < len(groups[i])]
        share = max(1, remaining // len(active))
        for i in active:
            take = min(share, len(groups[i]) - quota[i], remaining)
            quota[i] += take
            remaining -= take
            if not remaining:
                break
    idx = np.concatenate([rng.choice(g, q, replace=False) for g, q in zip(groups, quota) if q])
    return frame.iloc[idx].copy()


def assign_splits(frame, config, rng):
    result = frame.copy()
    result['split'] = 'development'
    test_var = config.get('test_group_var')
    test_groups = config.get('test_groups', [])
    if bool(test_var) != bool(test_groups):
        raise ValueError('test_group_var and nonempty test_groups must be supplied together')
    if test_var:
        if result[test_var].isna().any():
            raise ValueError(f'Missing test grouping values: {test_var}')
        missing = set(map(str, test_groups)) - set(result[test_var].astype(str))
        if missing:
            raise ValueError(f'Test groups absent from annotations/integration: {sorted(missing)}')
        result.loc[result[test_var].astype(str).isin(map(str, test_groups)), 'split'] = 'test'
    pool = result[result.split != 'test']
    group_var = config.get('validation_group_var') or 'batch'
    if set(pool[group_var]) & set(result.loc[result.split == 'test', group_var]):
        raise ValueError('Validation groups overlap the held-out test groups')
    variables = config.get('validation_variable_combinations', [])
    if pool[group_var].isna().any() or pool[variables].isna().any().any():
        raise ValueError('Missing validation grouping/combination values')
    units = pool[[group_var] + [v for v in variables if v != group_var]].drop_duplicates()
    if units[group_var].duplicated().any():
        raise ValueError('A validation group spans several variable combinations; choose compatible variables')
    validation = set()
    strata = [g for _, g in units.groupby(variables, sort=True, observed=True, dropna=False)] if variables else [units]
    for group in strata:
        n = min(len(group) - 1, max(1, int(len(group) * config.get('validation_fraction', .2))))
        if n > 0:
            validation.update(rng.choice(group[group_var].to_numpy(), n, replace=False))
    if not validation:
        raise ValueError('No whole groups available for validation; supply more groups or different validation combinations')
    result.loc[result.split != 'test', 'split'] = np.where(pool[group_var].isin(validation), 'validation', 'train')
    # Preserve rare-label training support by rescuing whole validation groups.
    excluded = config.get('excluded_labels', [])
    eligible = pool[pool.label.notna() & ~pool.label.isin(excluded)]
    support = eligible.groupby('label').size()
    supported = set(support[support >= config.get('min_cells_per_type', 20)].index)
    rescued = []
    for label in sorted(supported):
        g = result[(result.label == label) & (result.split != 'test')]
        target = config.get('min_cells_per_type', 20)
        while sum(g[group_var].isin(validation) == False) < target:
            sizes = g[g[group_var].isin(validation)].groupby(group_var).size().sort_values(ascending=False)
            if sizes.empty:
                raise ValueError(f'Unable to retain training support for {label}')
            unit = sizes.index[0]
            validation.remove(unit)
            rescued.append({'label': label, 'group': str(unit)})
    mask = result.split != 'test'
    result.loc[mask, 'split'] = np.where(result.loc[mask, group_var].isin(validation), 'validation', 'train')
    if not (result.split == 'validation').any():
        raise ValueError('No validation groups remain after rescuing rare labels')
    return result, supported, rescued


def prepare_plan(annots_f, cluster_csv, h5ads_yaml, config, batch_var, plan_f, coverage_f, counts_f):
    columns = ['cell_id', batch_var]
    if batch_var != 'sample_id':
        columns.append('project_id')
    vote_res = config.get('cluster_majority_vote_resolution')
    vote_col = f'RNA_snn_res.{float(vote_res):g}' if vote_res is not None else None
    if vote_col:
        columns.append(vote_col)
    annotations = pl.read_csv(annots_f).to_pandas()
    if annotations.cell_id.isna().any() or annotations.cell_id.duplicated().any():
        raise ValueError('Annotations must have unique nonmissing cell_id values')
    header = pl.read_csv(cluster_csv, n_rows=0).columns
    needed = set(config.get('validation_variable_combinations', []) + config.get('sampling_variable_combinations', []))
    needed |= {v for v in [config.get('test_group_var'), config.get('validation_group_var')] if v}
    columns += [v for v in needed if v in header and v not in annotations.columns and v not in columns]
    # Metadata in annotations is read independently of the batch identifier.
    integration = pl.read_csv(cluster_csv, columns=list(dict.fromkeys(columns))).to_pandas()
    if integration.cell_id.isna().any() or integration.cell_id.duplicated().any():
        raise ValueError('Integration must have unique nonmissing cell_id values')
    overlap = [v for v in annotations if v != 'cell_id' and v in integration]
    if overlap:
        check = integration.merge(annotations[['cell_id'] + overlap], on='cell_id', suffixes=('', '_annotation'))
        for var in overlap:
            both = check[var].notna() & check[var + '_annotation'].notna()
            if (check.loc[both, var].astype(str) != check.loc[both, var + '_annotation'].astype(str)).any():
                raise ValueError(f'Conflicting integration and annotation values: {var}')
        annotations = annotations.drop(columns=overlap)
    # Only annotated rows are needed for training/evaluation; full predictions
    # later retain unannotated cells independently.
    frame = integration.merge(annotations, on='cell_id', how='inner').rename(columns={'annotation': 'label'})
    if frame.empty:
        raise ValueError('No overlapping original labels and integration cells')
    manifest = yaml.safe_load(Path(h5ads_yaml).read_text())
    batch = frame[batch_var].astype(str)
    if not set(batch) <= set(manifest):
        if 'project_id' not in frame:
            raise ValueError('Input batch identifiers do not match H5AD manifest')
        batch = frame.project_id.astype(str) + '_' + batch
    if not set(batch) <= set(manifest):
        raise ValueError('Input batch identifiers do not match H5AD manifest')
    frame['batch'] = batch
    required = needed - {'batch'}
    missing = required - set(frame)
    if missing:
        raise ValueError(f'Missing grouping/sampling columns: {sorted(missing)}')
    if required and frame[list(required)].isna().any().any():
        raise ValueError('Missing values in grouping/sampling columns')
    rng = np.random.default_rng(config.get('seed', 42))
    frame, supported, rescued = assign_splits(frame, config, rng)
    frame['original_label'] = frame.label
    excluded = frame.label.isin(config.get('excluded_labels', []))
    frame.loc[excluded, 'split'] = 'excluded'
    frame.loc[frame.label.isna(), 'split'] = 'unlabeled'
    frame.loc[(~frame.label.isin(supported)) & frame.split.isin(['train', 'validation']), 'split'] = 'unsupported'
    if vote_col:
        training = frame[frame.split == 'train']
        tab = training.groupby([vote_col, 'label']).size().rename('n').reset_index()
        tab['total'] = tab.groupby(vote_col).n.transform('sum')
        top = tab.sort_values([vote_col, 'n', 'label'], ascending=[True, False, True]).drop_duplicates(vote_col)
        top = top[(top.n / top.total >= config.get('cluster_majority_vote_min_purity', .9)) & (top.total >= 10)]
        mapping = top.set_index(vote_col).label
        dev = frame.split.isin(['train', 'validation'])
        frame.loc[dev, 'label'] = frame.loc[dev, vote_col].map(mapping).fillna(frame.loc[dev, 'label'])
    frame['selected'] = False
    variables = config.get('sampling_variable_combinations', [])
    for (label, split), group in frame[frame.split.isin(['train', 'validation'])].groupby(['label', 'split'], sort=True):
        cap = config['train_cells_per_type'] if split == 'train' else config['validation_cells_per_type']
        sample = balanced_sample(group, cap, variables, rng)
        frame.loc[sample.index, 'selected'] = True
    selected = frame[frame.selected]
    if not {'train', 'validation'} <= set(selected.split):
        raise ValueError('Both training and validation cells are required')
    train_labels = set(selected[selected.split == 'train'].label)
    if not set(selected[selected.split == 'validation'].label) <= train_labels:
        raise ValueError('Validation label absent from training')
    frame['label_supported_in_training'] = frame.label.isin(train_labels)
    frame.to_csv(plan_f, index=False, compression='gzip')
    coverage = frame.groupby(['label', 'split'], dropna=False).agg(n_cells=('cell_id', 'size'), n_selected=('selected', 'sum')).reset_index()
    coverage.to_csv(coverage_f, index=False)
    selected[selected.split == 'train'].groupby('label').size().rename('len').reset_index().to_csv(counts_f, index=False)
    Path(str(plan_f) + '.json').write_text(json.dumps(dict(configuration=config, rescued_validation_groups=rescued,
        n_annotated=len(frame), n_train=int(sum(selected.split == 'train')),
        n_validation=int(sum(selected.split == 'validation')), n_test=int(sum(frame.split == 'test')),
        n_classes=len(train_labels), majority_vote_training_only=True), indent=2) + '\n')
    print(coverage.to_string(index=False), flush=True)
    return frame
