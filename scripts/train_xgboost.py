# Train an XGBoost classifier to predict cell type labels from scRNA-seq data.

import argparse
import pathlib
import re
import sys
import pandas as pd
import anndata as ad
import numpy as np
import polars as pl
import scipy.sparse as sp
import xgboost as xgb
import yaml


def load_label_mapping(label_map_f):

  if label_map_f is None:
    return None

  map_df = pl.read_csv(label_map_f)
  assert "annotation" in map_df.columns, "label_map_f must have 'annotation' column"
  assert "coarse_label" in map_df.columns, "label_map_f must have 'coarse_label' column"

  return dict(zip(
    map_df["annotation"].to_list(),
    map_df["coarse_label"].to_list(),
  ))


def normalize_counts(X, scale_factor=10000):
  X_csr = X.tocsr().astype(np.float64)
  lib_sizes = np.array(X_csr.sum(axis=1)).flatten()
  lib_sizes[lib_sizes == 0] = 1.0
  inv_lib = sp.diags(1.0 / lib_sizes)
  X_norm = inv_lib @ X_csr
  X_norm *= scale_factor
  X_norm.data = np.log1p(X_norm.data)
  return X_norm


def filter_genes_by_type(X_norm, X_raw, gene_names, gene_types,gene_exclude_re):
  # remove genes whose gene_type matches the exclusion regex.
  if gene_exclude_re is None or gene_exclude_re == "":
    return X_norm, X_raw, gene_names, gene_types

  pattern = re.compile(gene_exclude_re)
  keep_mask = np.array([not pattern.search(gt) for gt in gene_types])
  n_before = len(gene_names)
  n_after = int(keep_mask.sum())

  X_norm_filtered = X_norm[:, keep_mask]
  X_raw_filtered = X_raw[:, keep_mask]
  gene_names_filtered = [g for g, keep in zip(gene_names, keep_mask) if keep]
  gene_types_filtered = [gt for gt, keep in zip(gene_types, keep_mask) if keep]

  print(f"  Gene type filter (exclude '{gene_exclude_re}'): {n_before} → {n_after} genes "
        f"({n_before - n_after} removed)")

  return X_norm_filtered, X_raw_filtered, gene_names_filtered, gene_types_filtered


def filter_uninformative_genes(X_norm, X_raw, gene_names, min_cells=10):
  # remove genes expressed in fewer than min_cells cells.
  cells_per_gene = np.array(X_raw.getnnz(axis=0)).flatten()
  keep_mask = cells_per_gene >= min_cells
  n_before = len(gene_names)
  n_after = int(keep_mask.sum())

  X_norm_filtered = X_norm[:, keep_mask]
  X_raw_filtered = X_raw[:, keep_mask]
  gene_names_filtered = [g for g, keep in zip(gene_names, keep_mask) if keep]

  print(f"  Gene filter (>= {min_cells} cells): {n_before} → {n_after} genes "
        f"({n_before - n_after} removed)")

  return X_norm_filtered, X_raw_filtered, gene_names_filtered


def compute_pseudobulk(X_raw, cells_df, gene_names, sel_indices):

  X_sel = X_raw[:, sel_indices]
  sel_gene_names = [gene_names[i] for i in sel_indices]

  labels = cells_df["label"].to_list()
  batches = cells_df["batch"].to_list()

  unique_labels = sorted(set(labels))
  unique_batches = sorted(set(batches))
  n_batches = len(unique_batches)
  n_genes = len(sel_gene_names)
  batch_to_idx = {s: i for i, s in enumerate(unique_batches)}

  groups: dict[tuple[str, str], list[int]] = {}
  for i, (lbl, bid) in enumerate(zip(labels, batches)):
    key = (lbl, bid)
    if key not in groups:
      groups[key] = []
    groups[key].append(i)

  layers = {}
  n_cells_dict: dict[str, list[int]] = {}
  for lbl in unique_labels:
    count_mat = np.zeros((n_batches, n_genes), dtype=np.float64)
    n_cells_per_batch = [0] * n_batches
    for bid in unique_batches:
      key = (lbl, bid)
      if key in groups:
        indices = groups[key]
        count_mat[batch_to_idx[bid], :] = np.array(
          X_sel[indices, :].sum(axis=0)).flatten()
        n_cells_per_batch[batch_to_idx[bid]] = len(indices)

    layers[lbl] = count_mat
    n_cells_dict[lbl] = n_cells_per_batch

  obs_df = pd.DataFrame({"batch": unique_batches}, index=unique_batches)
  var_df = pd.DataFrame({"gene_id": sel_gene_names}, index=sel_gene_names)

  pb_adata = ad.AnnData(
    X=np.zeros((n_batches, n_genes)),
    obs=obs_df,
    var=var_df,
    layers=layers,
  )

  n_cells_mat = np.column_stack([n_cells_dict[lbl] for lbl in unique_labels])
  pb_adata.obsm["n_cells"] = pd.DataFrame(
    n_cells_mat, index=unique_batches, columns=unique_labels)

  print(f"  Pseudobulk raw counts: {n_batches} batches, {n_genes} genes, "
        f"{len(unique_labels)} labels (stored as layers)")

  return pb_adata


# xgboost training and feature selection
def _build_xgb_params(n_cores, seed, use_gpu, num_class, **overrides):
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


def run_xgboost_pass1(X_train, y_train, X_val, y_val,
  n_cores, seed, use_gpu, subsample, colsample_bytree, learning_rate,
  nrounds, early_stopping):
  num_class = int(y_train.max()) + 1

  params = _build_xgb_params(n_cores, seed, use_gpu, num_class,
    subsample=subsample,
    colsample_bytree=colsample_bytree,
    learning_rate=learning_rate,
  )

  dtrain = xgb.DMatrix(X_train, label=y_train)
  dval = xgb.DMatrix(X_val, label=y_val)

  model = xgb.train(
    params=params, dtrain=dtrain,
    num_boost_round=nrounds,
    evals=[(dtrain, "train"), (dval, "val")],
    early_stopping_rounds=early_stopping,
    verbose_eval=50,
  )

  print(f"  Best iteration: {model.best_iteration}")
  return model


def select_features_by_gain(model, gene_names,
  gain_threshold, min_genes, max_genes):

  importance = model.get_score(importance_type="gain")

  gene_gains = np.zeros(len(gene_names))
  for feat, gain in importance.items():
    gene_gains[int(feat[1:])] = gain

  sorted_idx = np.argsort(-gene_gains)
  sorted_gains = gene_gains[sorted_idx]

  nonzero_mask = sorted_gains > 0
  sorted_idx = sorted_idx[nonzero_mask]
  sorted_gains = sorted_gains[nonzero_mask]

  total_gain = sorted_gains.sum()
  cum_gain = np.cumsum(sorted_gains) / total_gain

  n_sel = int(np.searchsorted(cum_gain, gain_threshold)) + 1
  n_sel = max(n_sel, min_genes)
  n_sel = min(n_sel, max_genes)
  n_sel = min(n_sel, len(sorted_idx))

  sel_indices = sorted_idx[:n_sel]
  sel_gene_names = [gene_names[i] for i in sel_indices]

  gene_importance = pl.DataFrame({
    "gene_id": [gene_names[i] for i in sorted_idx],
    "gain": sorted_gains.tolist(),
    "cumulative_gain_frac": cum_gain.tolist(),
  })

  print(f"  Total genes with non-zero gain: {len(sorted_idx)}")
  print(f"  Genes at {gain_threshold*100:.0f}% cumulative gain: {n_sel}")
  print(f"  Top 10 genes: {sel_gene_names[:10]}")

  return sel_indices, sel_gene_names, gene_importance


def run_xgboost_pass2(X_train, y_train, X_val, y_val, n_cores, seed, use_gpu,
  subsample, colsample_bytree, learning_rate, nrounds, early_stopping):

  num_class = int(y_train.max()) + 1

  params = _build_xgb_params(n_cores, seed, use_gpu, num_class,
    subsample=subsample,
    colsample_bytree=colsample_bytree,
    learning_rate=learning_rate,
  )

  feat_names = [f"g{i}" for i in range(X_train.shape[1])]
  dtrain = xgb.DMatrix(X_train, label=y_train, feature_names=feat_names)
  dval = xgb.DMatrix(X_val, label=y_val, feature_names=feat_names)

  model = xgb.train(
    params=params, dtrain=dtrain,
    num_boost_round=nrounds,
    evals=[(dtrain, "train"), (dval, "val")],
    early_stopping_rounds=early_stopping,
    verbose_eval=50,
  )

  print(f"  Best iteration: {model.best_iteration}")
  return model


# model evaluation
def evaluate_model(model, X_val, y_val, class_names):
  
  feat_names = [f"g{i}" for i in range(X_val.shape[1])]
  dval = xgb.DMatrix(X_val, feature_names=feat_names)
  probs = model.predict(dval)
  y_pred = probs.argmax(axis=1)
  p_max = probs.max(axis=1)

  overall_acc = (y_pred == y_val).mean()
  print(f"  Overall accuracy: {overall_acc:.3f}")

  hc_mask = p_max > 0.5
  if hc_mask.sum() > 0:
    hc_acc = (y_pred[hc_mask] == y_val[hc_mask]).mean()
    print(f"  High-confidence accuracy (p>0.5): {hc_acc:.3f} ({hc_mask.sum()} cells)")

  print(f"\n  {'Class':<30} {'N':>6} {'Acc':>6} {'F1':>6}")
  print("  " + "-" * 52)

  f1_scores = []
  for i, name in enumerate(class_names):
    mask_true = y_val == i
    n_true = mask_true.sum()
    if n_true == 0:
      continue

    acc = (y_pred[mask_true] == i).mean()
    pred_pos = y_pred == i
    tp = ((y_pred == i) & (y_val == i)).sum()
    precision = tp / pred_pos.sum() if pred_pos.sum() > 0 else 0
    recall = tp / n_true if n_true > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    f1_scores.append(f1)

    print(f"  {name:<30} {n_true:>6} {acc:>6.3f} {f1:>6.3f}")

  macro_f1 = np.mean(f1_scores) if f1_scores else 0
  print(f"\n  Macro F1: {macro_f1:.3f}")

  print("\n  Top misclassifications:")
  print(f"  {'True':<20} {'Predicted':<20} {'N':>5} {'%':>6}")
  print("  " + "-" * 55)

  confusion_pairs = {}
  for i in range(len(y_val)):
    if y_pred[i] != y_val[i]:
      key = (class_names[y_val[i]], class_names[y_pred[i]])
      confusion_pairs[key] = confusion_pairs.get(key, 0) + 1

  sorted_pairs = sorted(confusion_pairs.items(), key=lambda x: -x[1])[:15]
  for (true_name, pred_name), count in sorted_pairs:
    n_true_class = (y_val == class_names.index(true_name)).sum()
    pct = count / n_true_class * 100 if n_true_class > 0 else 0
    print(f"  {true_name:<20} {pred_name:<20} {count:>5} {pct:>5.1f}%")


def evaluate_model_coarse(model, X_val, y_val, class_names, label_map):

  feat_names = [f"g{i}" for i in range(X_val.shape[1])]
  dval = xgb.DMatrix(X_val, feature_names=feat_names)
  probs = model.predict(dval)
  y_pred_fine = probs.argmax(axis=1)

  def to_coarse(fine_idx):
    return label_map.get(class_names[fine_idx], class_names[fine_idx])

  y_true_coarse = np.array([to_coarse(i) for i in y_val])
  y_pred_coarse = np.array([to_coarse(i) for i in y_pred_fine])

  coarse_classes = sorted(set(y_true_coarse))
  overall_acc = (y_true_coarse == y_pred_coarse).mean()
  print(f"  Overall coarse accuracy: {overall_acc:.3f}")

  print(f"\n  {'Coarse class':<30} {'N':>6} {'Acc':>6} {'F1':>6}")
  print("  " + "-" * 52)

  f1_scores = []
  for cls in coarse_classes:
    mask_true = y_true_coarse == cls
    n_true = mask_true.sum()
    if n_true == 0:
      continue

    acc = (y_pred_coarse[mask_true] == cls).mean()
    pred_pos = y_pred_coarse == cls
    tp = ((y_pred_coarse == cls) & (y_true_coarse == cls)).sum()
    precision = tp / pred_pos.sum() if pred_pos.sum() > 0 else 0
    recall = tp / n_true if n_true > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    f1_scores.append(f1)

    print(f"  {cls:<30} {n_true:>6} {acc:>6.3f} {f1:>6.3f}")

  macro_f1 = np.mean(f1_scores) if f1_scores else 0
  print(f"\n  Coarse Macro F1: {macro_f1:.3f}")


# save outputs
def save_outputs(model, class_names, selected_genes, gene_importance,
  output_dir, ref_tag, cells_df, label_map=None, X_train=None, y_train=None,
  X_val=None, y_val=None, pseudobulk_adata=None):

  # save model, gene list, classes, importance, predictions, and pseudobulk
  out_dir = pathlib.Path(output_dir)

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
  gene_importance.write_csv(str(imp_path))
  print(f"  Importance: {imp_path}")

  if label_map is not None:
    map_path = out_dir / f"{ref_tag}_label_map.csv"
    pl.DataFrame({
      "fine_label": list(label_map.keys()),
      "coarse_label": list(label_map.values()),
    }).write_csv(str(map_path))
    print(f"  Label map: {map_path}")

  # save predictions for all cells (train + validation)
  if X_val is not None and X_train is not None:
    preds_path = out_dir / f"{ref_tag}_predictions.csv.gz"
    preds_df = _make_predictions_df(
      model, X_train, y_train, X_val, y_val, cells_df, class_names, label_map
    )
    preds_df.write_csv(str(preds_path), compression="gzip")
    print(f"  Predictions: {preds_path}")

  if pseudobulk_adata is not None:
    pb_path = out_dir / f"{ref_tag}_pseudobulk.h5ad"
    pseudobulk_adata.write_h5ad(str(pb_path))
    print(f"  Pseudobulk: {pb_path}")


def _make_predictions_df(model, X_train, y_train, X_val, y_val, cells_df, class_names, label_map=None):
  # genertae predictions for all cells
  feat_names = [f"g{i}" for i in range(X_val.shape[1])]

  # Predict on validation set
  dval = xgb.DMatrix(X_val, feature_names=feat_names)
  val_probs = model.predict(dval)
  val_pred = val_probs.argmax(axis=1)
  val_pmax = val_probs.max(axis=1)

  # Predict on training set
  dtrain = xgb.DMatrix(X_train, feature_names=feat_names)
  train_probs = model.predict(dtrain)
  train_pred = train_probs.argmax(axis=1)
  train_pmax = train_probs.max(axis=1)

  # Build arrays in cells_df order (train first, then validation)
  train_mask = cells_df["split"].to_numpy() == "train"
  all_pred = np.empty(len(cells_df), dtype=int)
  all_pmax = np.empty(len(cells_df))
  all_pred[train_mask] = train_pred
  all_pred[~train_mask] = val_pred
  all_pmax[train_mask] = train_pmax
  all_pmax[~train_mask] = val_pmax

  pred_labels = [class_names[i] for i in all_pred]

  result = cells_df.select(["cell_id", "batch", "label", "split"]).with_columns([
    pl.Series("predicted_label", pred_labels),
    pl.Series("probability", all_pmax),
  ])

  # Add coarse labels if mapping is available
  if label_map is not None:
    result = result.with_columns([
      pl.col("label").replace(label_map).alias("coarse_true"),
      pl.col("predicted_label").replace(label_map).alias("coarse_predicted"),
    ])

  return result


def aggregate_fulldata_predictions(pred_fs, annots_f, subsample_preds_f,
  label_map_f, output_f):

  print("=" * 60)
  print("Aggregating full-data predictions")
  print("=" * 60)

  all_preds = pl.concat([pl.read_csv(f) for f in pred_fs])
  print(f"  Total cells predicted: {all_preds.shape[0]}")

  all_preds = all_preds.select([
    "cell_id",
    pl.col("predicted_label_naive").alias("predicted_label"),
    pl.col("probability_naive").alias("probability"),
  ])

  annots_df = pl.read_csv(annots_f)
  assert "cell_id" in annots_df.columns, "annots_f must have 'cell_id' column"
  assert "annotation" in annots_df.columns, "annots_f must have 'annotation' column"
  annots_df = annots_df.select(["cell_id", pl.col("annotation").alias("label")])

  subsample_df = pl.read_csv(subsample_preds_f).select(["cell_id", "split"])

  result = all_preds.join(annots_df, on="cell_id", how="left")
  result = result.join(subsample_df, on="cell_id", how="left")

  result = result.with_columns(
    pl.when(pl.col("split").is_not_null())
      .then(pl.col("split"))
      .when(pl.col("label").is_not_null())
      .then(pl.lit("holdout"))
      .otherwise(pl.lit("unlabeled"))
      .alias("split")
  )

  label_map = load_label_mapping(label_map_f)
  if label_map is not None:
    result = result.with_columns([
      pl.col("label").replace(label_map).alias("coarse_true"),
      pl.col("predicted_label").replace(label_map).alias("coarse_predicted"),
    ])

  split_summary = result.group_by("split").len().sort("split")
  for row in split_summary.iter_rows(named=True):
    print(f"  {row['split']}: {row['len']} cells")

  result.write_csv(output_f, compression="gzip")
  print(f"  Saved to {output_f}")


if __name__ == "__main__":
  import json
  from train_xgboost_plan import prepare_plan
  from train_xgboost_workflow import fit_plan, aggregate_plan_predictions
  parser = argparse.ArgumentParser(description="Train an XGBoost classifier from scprocess outputs")
  commands = parser.add_subparsers(dest="command", required=True)
  plan = commands.add_parser("plan")
  for name in ("annots_f", "cluster_csv", "h5ads_yaml", "batch_var", "config_json", "plan_f", "coverage_f", "counts_f"):
    plan.add_argument("--" + name, required=True)
  fit = commands.add_parser("train")
  for name in ("plan_f", "h5ads_yaml", "config_json", "output_dir"):
    fit.add_argument("--" + name, required=True)
  fit.add_argument("--n_cores", type=int, required=True)
  aggregate = commands.add_parser("aggregate_fulldata")
  aggregate.add_argument("pred_fs", nargs="+")
  aggregate.add_argument("--plan_f", required=True)
  aggregate.add_argument("--cluster_csv", required=True)
  aggregate.add_argument("--output_f", required=True)
  args = parser.parse_args()
  if args.command == "plan":
    prepare_plan(args.annots_f, args.cluster_csv, args.h5ads_yaml, json.loads(args.config_json),
      args.batch_var, args.plan_f, args.coverage_f, args.counts_f)
  elif args.command == "train":
    fit_plan(args.plan_f, args.h5ads_yaml, json.loads(args.config_json), args.output_dir, args.n_cores)
  else:
    aggregate_plan_predictions(args.pred_fs, args.plan_f, args.output_f, args.cluster_csv)
