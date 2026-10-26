"""Fit and evaluate classifiers from an explicit scprocess training plan."""
from pathlib import Path
import gc
import hashlib
import json
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp
import xgboost as xgb
import yaml


def read_selected_h5ad(path, cells, reference=None):
    import anndata as ad
    before = (Path(path).stat().st_size, Path(path).stat().st_mtime_ns)
    data = ad.read_h5ad(path, backed='r')
    try:
        ids = pd.Index(data.obs['cell_id'].astype(str) if 'cell_id' in data.obs else data.obs_names.astype(str))
        if ids.has_duplicates:
            raise ValueError(f'Duplicate cell IDs: {path}')
        idx = ids.get_indexer(cells)
        if np.any(idx < 0) or len(set(cells)) != len(cells):
            raise ValueError(f'Missing/duplicated selected cells: {path}')
        genes = data.var_names.astype(str).tolist()
        if len(set(genes)) != len(genes) or any(not g.strip() for g in genes):
            raise ValueError(f'Invalid gene IDs: {path}')
        types = data.var['gene_type'].astype(object).fillna('unknown').astype(str).tolist() if 'gene_type' in data.var else ['unknown'] * len(genes)
        order = np.argsort(idx)
        matrix = sp.csr_matrix(data.X[idx[order], :])[np.argsort(order)].tocsr()
        if reference is not None and genes != reference:
            if set(genes) != set(reference):
                raise ValueError(f'Different gene sets: {path}')
            import warnings
            warnings.warn(f'Same genes, different order: aligning {path}', stacklevel=2)
            indices = pd.Index(genes).get_indexer(reference)
            matrix = matrix[:, indices].tocsr()
            types = [types[i] for i in indices]
            genes = list(reference)
    finally:
        data.file.close()
    if before != (Path(path).stat().st_size, Path(path).stat().st_mtime_ns):
        raise RuntimeError(f'Reading changed input signature: {path}')
    return matrix, genes, types


def load_plan_expression(plan, manifest, normalize):
    raw_parts, norm_parts, ordered = [], [], []
    genes = types = None
    for batch, group in plan.groupby('batch', sort=True):
        entry = manifest[batch]
        path = entry['path'] if isinstance(entry, dict) else entry
        raw, current, current_types = read_selected_h5ad(path, group.cell_id.tolist(), genes)
        if genes is None:
            genes, types = current, current_types
        raw_parts.append(raw)
        norm_parts.append(normalize(raw).astype(np.float32))
        ordered.append(group)
        print(f'Loading {batch}: {len(group)} selected cells', flush=True)
    raw = sp.vstack(raw_parts, format='csr'); raw_parts.clear()
    norm = sp.vstack(norm_parts, format='csr'); norm_parts.clear()
    return raw, norm, genes, types, pd.concat(ordered, ignore_index=True)


def best_iteration_model(model):
    return model[:int(model.best_iteration) + 1]


def training_gene_mask(raw, train, minimum):
    detected = np.asarray(raw[train].getnnz(axis=0)).ravel()
    return detected >= minimum


def verify_gpu(n_cores):
    if not xgb.build_info().get('USE_CUDA'):
        raise RuntimeError('GPU training requested but XGBoost was built without CUDA')
    probe = xgb.DMatrix(np.array([[0., 1.], [1., 0.], [1., 1.], [0., 0.]], dtype=np.float32),
        label=np.array([0, 1, 1, 0]), nthread=min(n_cores, 2))
    model = xgb.train({'objective': 'binary:logistic', 'device': 'cuda',
        'tree_method': 'hist', 'nthread': min(n_cores, 2)}, probe, num_boost_round=1)
    device = json.loads(model.save_config())['learner']['generic_param']['device']
    if not device.startswith('cuda'):
        raise RuntimeError('GPU training requested but XGBoost fell back to CPU')
    print('Verified XGBoost GPU execution on ' + device, flush=True)
    return device


def fit_plan(plan_f, h5ads_yaml, config, output_dir, n_cores):
    import train_xgboost as native
    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    gpu_device = verify_gpu(n_cores) if config.get('use_gpu', False) else 'cpu'
    tag = config['xgboost_tag']
    plan = pd.read_csv(plan_f)
    plan = plan[plan.selected & plan.split.isin(['train', 'validation'])].copy()
    manifest = yaml.safe_load(Path(h5ads_yaml).read_text())
    raw, norm, genes, types, plan = load_plan_expression(plan, manifest, native.normalize_counts)
    train = plan.split.to_numpy() == 'train'; val = ~train
    classes = sorted(plan.loc[train, 'label'].unique())
    if not set(plan.loc[val, 'label']) <= set(classes):
        raise ValueError('Validation label absent from training')
    y = pd.Index(classes).get_indexer(plan.label)
    norm, raw, genes, types = native.filter_genes_by_type(norm, raw, genes, types, config.get('gene_exclude_re'))
    keep = training_gene_mask(raw, train, config.get('min_cells_expressed', 10))
    norm = norm[:, keep].tocsr(); raw = raw[:, keep].tocsr()
    genes = [g for g, k in zip(genes, keep) if k]
    if not genes:
        raise ValueError('No genes pass training-only filtering')
    xtrain, xval = norm[train], norm[val]
    model1 = native.run_xgboost_pass1(xtrain, y[train], xval, y[val], n_cores,
        config.get('seed', 42), config.get('use_gpu', False), config.get('pass1_subsample', .632),
        config.get('pass1_colsample_bytree', .1), config.get('pass1_learning_rate', .1),
        config.get('pass1_nrounds', 300), config.get('pass1_early_stopping', 10))
    best1 = int(model1.best_iteration)
    model1 = best_iteration_model(model1)
    selected, names, importance = native.select_features_by_gain(model1, genes,
        config.get('gain_threshold', .9), config.get('min_genes', 100), config.get('max_genes', 3000))
    del model1, norm; gc.collect()
    xtrain, xval = xtrain[:, selected].tocsr(), xval[:, selected].tocsr()
    model2 = native.run_xgboost_pass2(xtrain, y[train], xval, y[val], n_cores,
        config.get('seed', 42), config.get('use_gpu', False), config.get('pass1_subsample', .632),
        config.get('pass2_colsample_bytree', .5), config.get('pass2_learning_rate', .05),
        config.get('pass2_nrounds', 500), config.get('pass2_early_stopping', 10))
    best2 = int(model2.best_iteration)
    model2 = best_iteration_model(model2)
    plan_pl = pl.from_pandas(plan)
    pb = native.compute_pseudobulk(raw, plan_pl, genes, selected)
    native.save_outputs(model2, classes, names, importance, out, tag, plan_pl,
        native.load_label_mapping(config.get('label_map_f')), xtrain, y[train], xval, y[val], pb)
    path = out / f'{tag}_xgboost_model.json'
    (out / f'{tag}_model_card.json').write_text(json.dumps(dict(configuration=config,
        n_train=int(train.sum()), n_validation=int(val.sum()), n_classes=len(classes),
        n_genes=len(names), n_genes_tested=len(genes), pass1_best_iteration=best1, pass2_best_iteration=best2,
        detection_filter_training_only=True, best_iteration_model_sliced=True,
        model_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        plan_sha256=hashlib.sha256(Path(plan_f).read_bytes()).hexdigest(),
        xgboost_version=xgb.__version__, xgboost_build=xgb.build_info(), training_device=gpu_device, normalization='log1p(10000 * counts / all-input-gene library size)',
        status='Unreviewed candidate; not registered or published'), indent=2) + '\n')


def aggregate_plan_predictions(pred_fs, plan_f, output_f, cluster_csv):
    preds = pl.concat([pl.read_csv(p, columns=['cell_id', 'predicted_label_naive',
        'probability_naive'], schema_overrides={'cell_id': pl.String,
        'predicted_label_naive': pl.String, 'probability_naive': pl.Float64}) for p in pred_fs])
    if preds['cell_id'].n_unique() != len(preds):
        raise ValueError('Duplicate full-data predictions')
    preds = preds.select('cell_id', pl.col('predicted_label_naive').alias('predicted_label'),
        pl.col('probability_naive').alias('probability'))
    join_cells = pl.read_csv(cluster_csv, columns=['cell_id'])
    preds = preds.join(join_cells, on='cell_id', how='inner')
    if len(preds) != len(join_cells):
        raise ValueError('Integrated cells missing from full-data predictions')
    plan = pl.read_csv(plan_f)
    fields = ['cell_id', 'label', 'original_label', 'split', 'selected', 'label_supported_in_training']
    result = preds.join(plan.select(fields), on='cell_id', how='left').with_columns(
        pl.col('split').fill_null('unlabeled'), pl.col('selected').fill_null(False),
        pl.col('label_supported_in_training').fill_null(False))
    if not set(plan['cell_id']) <= set(result['cell_id']):
        raise ValueError('Planned cells missing from full-data predictions')
    configuration = json.loads(Path(str(plan_f) + '.json').read_text())['configuration'] if Path(str(plan_f) + '.json').exists() else {}
    if configuration.get('label_map_f'):
        mapping = pd.read_csv(configuration['label_map_f']).set_index('annotation').coarse_label.to_dict()
        result = result.with_columns(pl.col('label').replace(mapping).alias('coarse_true'),
            pl.col('predicted_label').replace(mapping).alias('coarse_predicted'))
    result.write_csv(output_f, compression='gzip')
    from sklearn.metrics import precision_recall_fscore_support, accuracy_score
    metrics, per_label = [], []
    frame = result.to_pandas()
    frame.loc[(~frame.selected) & frame.split.isin(['train', 'validation']), 'split'] = 'development_holdout'
    for split, group in frame[frame.label.notna() & ~frame.split.isin(['excluded', 'unsupported', 'unlabeled'])].groupby('split'):
        labels = sorted(set(group.label) | set(group.predicted_label))
        precision, recall, f1, support = precision_recall_fscore_support(group.label, group.predicted_label, labels=labels, zero_division=0)
        observed = support > 0
        metrics.append(dict(split=split, n_cells=len(group), accuracy=accuracy_score(group.label, group.predicted_label),
            macro_f1=float(f1[observed].mean()), balanced_accuracy=float(recall[observed].mean())))
        per_label.append(pd.DataFrame(dict(split=split, label=labels, precision=precision, recall=recall, f1=f1, n_test=support)))
    prefix = str(output_f).removesuffix('_fulldata_predictions.csv.gz')
    pd.DataFrame(metrics).to_csv(prefix + '_evaluation.csv', index=False)
    pd.concat(per_label, ignore_index=True).to_csv(prefix + '_per_label_evaluation.csv', index=False)
    print(pd.DataFrame(metrics).to_string(index=False), flush=True)
