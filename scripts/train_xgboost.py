"""Train an XGBoost classifier to predict cell type labels from scRNA-seq data.

Workflow:
  1. Plan training cells from cluster CSV + annotations (no H5AD loaded)
  2. Load expression data one H5AD at a time, subset to selected cells
  3. XGBoost pass 1: broad exploration with all genes
  4. Feature selection via cumulative gain
  5. XGBoost pass 2: final model on selected genes
  6. Evaluate and save model
"""

import argparse
import pathlib
import sys
from typing import Optional
import anndata as ad
import numpy as np
import polars as pl
import scipy.sparse as sp
import xgboost as xgb
import yaml


def make_classifier(
  annots_f: str, cluster_csv: str, h5ads_yaml: str,
  ref_tag: str, output_dir: str, batch_var: str, int_res_ls: str,
  label_map_f: Optional[str],
  refine_labels: bool, purity_threshold: float,
  n_cells_per_type: int, min_cells_per_type: int,
  min_cells_expressed: int, seed: int, n_cores: int, use_gpu: bool,
  pass1_subsample: float, pass1_colsample_bytree: float,
  pass1_learning_rate: float, pass1_nrounds: int, pass1_early_stopping: int,
  pass2_colsample_bytree: float, pass2_learning_rate: float,
  pass2_nrounds: int, pass2_early_stopping: int,
  gain_threshold: float, min_genes: int, max_genes: int,
) -> None:
  
  print("=" * 60)
  print("XGBoost cell type classifier training")
  print("=" * 60)

  out_dir = pathlib.Path(output_dir)
  out_dir.mkdir(parents=True, exist_ok=True)
  (out_dir / "plots").mkdir(exist_ok=True)

  with open(h5ads_yaml) as f:
    h5ad_dict = yaml.safe_load(f)

  res_ls = [float(x) for x in int_res_ls.split()]
  max_res = max(res_ls)
  max_res_str = str(int(max_res)) if max_res == int(max_res) else str(max_res)
  hi_res_col = f"RNA_snn_res.{max_res_str}"

  label_map = load_label_mapping(label_map_f)
  if label_map is not None:
    print(f"  Label mapping loaded: {len(label_map)} fine→coarse entries")
    print(f"  Coarse categories: {sorted(set(label_map.values()))}")

  # Phase 1: plan which cells to use (CSV only, no H5AD)
  print("\n--- Planning training cells ---")
  cells_df = plan_training_cells(
    annots_f, cluster_csv, batch_var, hi_res_col,
    refine_labels, purity_threshold,
    min_cells_per_type, n_cells_per_type, seed,
  )

  # Phase 2: load expression data one H5AD at a time
  print("\n--- Loading expression data ---")
  X, gene_names, cell_ids = load_expression_matrix(cells_df, h5ad_dict)
  X, gene_names = filter_uninformative_genes(X, gene_names, min_cells_expressed)

  # Align cells_df row order with the matrix
  cells_df = cells_df.filter(pl.col("cell_id").is_in(cell_ids))
  id_order = pl.DataFrame({"cell_id": cell_ids, "_row_idx": range(len(cell_ids))})
  cells_df = cells_df.join(id_order, on="cell_id").sort("_row_idx").drop("_row_idx")

  # Encode labels as integers (sorted alphabetically)
  class_names = sorted(cells_df["label"].unique().to_list())
  label_to_int = {name: i for i, name in enumerate(class_names)}
  y = np.array([label_to_int[lbl] for lbl in cells_df["label"].to_list()])

  train_mask = cells_df["split"].to_numpy() == "train"
  val_mask = ~train_mask
  X_train, y_train = X[train_mask], y[train_mask]
  X_val, y_val = X[val_mask], y[val_mask]

  print(f"  Train: {X_train.shape[0]} cells, Validation: {X_val.shape[0]} cells")
  print(f"  Features (genes after filtering): {X_train.shape[1]}")
  print(f"  Classes: {len(class_names)}")

  # Phase 3: XGBoost pass 1 — broad gene exploration
  print("\n--- XGBoost Pass 1 (all genes) ---")
  model_pass1 = run_xgboost_pass1(
    X_train, y_train, X_val, y_val,
    n_cores, seed, use_gpu,
    pass1_subsample, pass1_colsample_bytree, pass1_learning_rate,
    pass1_nrounds, pass1_early_stopping,
  )

  # Phase 4: feature selection via gain scores
  print("\n--- Feature selection ---")
  sel_indices, sel_gene_names, gene_importance = select_features_by_gain(
    model_pass1, gene_names, gain_threshold, min_genes, max_genes,
  )
  print(f"  Selected {len(sel_indices)} genes (of {len(gene_names)})")

  X_train_sub = X_train[:, sel_indices]
  X_val_sub = X_val[:, sel_indices]

  # Phase 5: XGBoost pass 2 — final training on curated gene set
  print("\n--- XGBoost Pass 2 (selected genes) ---")
  model_pass2 = run_xgboost_pass2(
    X_train_sub, y_train, X_val_sub, y_val,
    n_cores, seed, use_gpu,
    pass1_subsample, pass2_colsample_bytree, pass2_learning_rate,
    pass2_nrounds, pass2_early_stopping,
  )

  # Phase 6: evaluate and save
  print("\n--- Evaluation (fine labels) ---")
  evaluate_model(model_pass2, X_val_sub, y_val, class_names)

  if label_map is not None:
    print("\n--- Evaluation (coarse labels) ---")
    evaluate_model_coarse(model_pass2, X_val_sub, y_val, class_names, label_map)

  print("\n--- Saving outputs ---")
  save_outputs(
    model_pass2, class_names, sel_gene_names, gene_importance,
    output_dir, ref_tag, cells_df, label_map,
    X_train_sub, y_train, X_val_sub, y_val,
  )

  print("\nDone.")


# planning phase (CSV only)

def plan_training_cells(
  annots_f: str, cluster_csv: str, batch_var: str, hi_res_col: str,
  refine_labels: bool, purity_threshold: float,
  min_cells_per_type: int, n_cells_per_type: int, seed: int,
) -> pl.DataFrame:

  cols_to_read = ["cell_id", batch_var, hi_res_col]
  cluster_df = pl.read_csv(cluster_csv, columns=cols_to_read)
  print(f"  Loaded cluster CSV: {cluster_df.shape[0]} cells")

  annots_df = pl.read_csv(annots_f)
  assert "cell_id" in annots_df.columns, "annots_f must have 'cell_id' column"
  assert "annotation" in annots_df.columns, "annots_f must have 'annotation' column"
  annots_df = annots_df.select(["cell_id", "annotation"])

  df = cluster_df.join(annots_df, on="cell_id", how="left")
  n_annotated = df.filter(pl.col("annotation").is_not_null()).shape[0]
  print(f"  Cells with annotations: {n_annotated} / {df.shape[0]}")

  if refine_labels:
    print("  Refining labels by cluster majority voting...")
    df = refine_labels_by_cluster(df, hi_res_col, purity_threshold)
  else:
    df = df.with_columns(pl.col("annotation").alias("label"))

  df = df.filter(pl.col("label").is_not_null())
  print(f"  Cells with labels after processing: {df.shape[0]}")

  # Exclude rare cell types
  type_counts = df.group_by("label").len()
  keep_types = (
    type_counts
    .filter(pl.col("len") >= min_cells_per_type)
    ["label"].to_list()
  )
  excluded = type_counts.filter(pl.col("len") < min_cells_per_type)
  if excluded.shape[0] > 0:
    print(f"  Excluding {excluded.shape[0]} cell types with < {min_cells_per_type} cells:")
    for row in excluded.iter_rows(named=True):
      print(f"    {row['label']}: {row['len']} cells")
  df = df.filter(pl.col("label").is_in(keep_types))

  print("  Downsampling...")
  df = downsample_per_type(df, n_cells_per_type, seed)

  print("  Assigning train/val split...")
  df = assign_train_val_split(df, seed)

  # Batch column determines which H5AD file contains each cell
  if batch_var == "sample_id":
    df = df.with_columns(pl.col("sample_id").alias("batch"))
  else:
    df = df.with_columns(pl.col(batch_var).alias("batch"))

  split_summary = df.group_by("split").len()
  for row in split_summary.iter_rows(named=True):
    print(f"    {row['split']}: {row['len']} cells")

  return df.select(["cell_id", "sample_id", "label", "split", "batch"])


def refine_labels_by_cluster(
  df: pl.DataFrame, hi_res_col: str, purity_threshold: float,
) -> pl.DataFrame:

  min_cluster_size = 10

  annotated = df.filter(pl.col("annotation").is_not_null())

  cluster_annot_counts = (
    annotated.group_by([hi_res_col, "annotation"])
    .len().rename({"len": "n"})
  )
  cluster_totals = (
    annotated.group_by(hi_res_col)
    .len().rename({"len": "n_total"})
  )

  cluster_stats = (
    cluster_annot_counts
    .join(cluster_totals, on=hi_res_col)
    .with_columns((pl.col("n") / pl.col("n_total")).alias("purity"))
    .sort([hi_res_col, "purity"], descending=[False, True])
    .group_by(hi_res_col)
    .first()
    .filter(
      (pl.col("purity") >= purity_threshold) & (pl.col("n_total") >= min_cluster_size)
    )
    .select([hi_res_col, pl.col("annotation").alias("majority_label")])
  )

  n_refined = cluster_stats.shape[0]
  n_total_clusters = df[hi_res_col].n_unique()
  print(f"    {n_refined}/{n_total_clusters} clusters pass purity threshold")

  df = df.join(cluster_stats, on=hi_res_col, how="left")
  df = df.with_columns(
    pl.when(pl.col("majority_label").is_not_null())
    .then(pl.col("majority_label")).otherwise(pl.col("annotation"))
    .alias("label")
  ).drop("majority_label")

  return df


def load_label_mapping(label_map_f: Optional[str]) -> Optional[dict[str, str]]:

  if label_map_f is None:
    return None

  map_df = pl.read_csv(label_map_f)
  assert "annotation" in map_df.columns, "label_map_f must have 'annotation' column"
  assert "coarse_label" in map_df.columns, "label_map_f must have 'coarse_label' column"

  return dict(zip(
    map_df["annotation"].to_list(),
    map_df["coarse_label"].to_list(),
  ))


def downsample_per_type(
  df: pl.DataFrame, n_cells_per_type: int, seed: int,
) -> pl.DataFrame:
  sampled = (
    df
    .with_columns(pl.lit(1).alias("_dummy"))
    .group_by("label")
    .map_groups(
      lambda group: group.sample(n=min(n_cells_per_type, group.shape[0]), seed=seed)
    )
    .drop("_dummy")
  )

  type_summary = sampled.group_by("label").len().sort("len", descending=True)
  print(f"    {sampled.shape[0]} cells after downsampling ({type_summary.shape[0]} types)")

  return sampled


def assign_train_val_split(df: pl.DataFrame, seed: int) -> pl.DataFrame:

  rng = np.random.default_rng(seed)
  samples = df["sample_id"].drop_nulls().unique().to_list()
  n_samples = len(samples)

  if n_samples > 4:
    samples = sorted(samples)
    n_val = max(1, len(samples) // 5)
    val_samples = rng.choice(samples, size=n_val, replace=False).tolist()
    print(f"    Sample-level holdout: {n_val} validation samples of {len(samples)}")
    print(f"    Validation samples: {val_samples}")

    df = df.with_columns(
      pl.when(pl.col("sample_id").is_in(val_samples))
      .then(pl.lit("validation")).otherwise(pl.lit("train"))
      .alias("split")
    )
  else:
    print(f"    Cell-level stratified split (<=4 samples)")

    labels = df["label"].to_list()
    assignments = ["train"] * len(labels)

    for label in df["label"].unique().to_list():
      label_indices = [i for i, lbl in enumerate(labels) if lbl == label]
      n_val = max(1, int(len(label_indices) * 0.2))
      val_indices = rng.choice(label_indices, size=n_val, replace=False).tolist()
      for idx in val_indices:
        assignments[idx] = "validation"

    df = df.with_columns(pl.Series("split", assignments))

  return df


# data loading

def normalize_counts(X: sp.spmatrix, scale_factor: int = 10000) -> sp.csr_matrix:
  """Total-count normalize to scale_factor, then log1p. Matches label_celltypes.R."""
  X_csr = X.tocsr().astype(np.float64)
  lib_sizes = np.array(X_csr.sum(axis=1)).flatten()
  lib_sizes[lib_sizes == 0] = 1.0
  inv_lib = sp.diags(1.0 / lib_sizes)
  X_norm = inv_lib @ X_csr
  X_norm *= scale_factor
  X_norm.data = np.log1p(X_norm.data)
  return X_norm


def filter_uninformative_genes(
  X: sp.csr_matrix, gene_names: list[str], min_cells: int = 10
) -> tuple[sp.csr_matrix, list[str]]:
  """Remove genes expressed in fewer than min_cells cells."""
  cells_per_gene = np.array(X.getnnz(axis=0)).flatten()
  keep_mask = cells_per_gene >= min_cells
  n_before = len(gene_names)
  n_after = keep_mask.sum()

  X_filtered = X[:, keep_mask]
  gene_names_filtered = [g for g, keep in zip(gene_names, keep_mask) if keep]

  print(f"  Gene filter (>= {min_cells} cells): {n_before} → {n_after} genes "
        f"({n_before - n_after} removed)")

  return X_filtered, gene_names_filtered


def load_expression_matrix(
  cells_df: pl.DataFrame, h5ad_dict: dict,
) -> tuple[sp.csr_matrix, list[str], list[str]]:
  
  batches_needed = cells_df["batch"].drop_nulls().unique().to_list()

  matrices = []
  cell_ids_ordered = []
  gene_names: Optional[list[str]] = None

  for batch in sorted(batches_needed):
    if batch not in h5ad_dict:
      print(f"  WARNING: batch '{batch}' not found in h5ad_dict, skipping")
      continue

    h5ad_path = h5ad_dict[batch]
    batch_cell_ids = (
      cells_df
      .filter(pl.col("batch") == batch)
      ["cell_id"].to_list()
    )

    if len(batch_cell_ids) == 0:
      continue

    print(f"  Loading batch '{batch}': {len(batch_cell_ids)} cells")
    adata = ad.read_h5ad(h5ad_path)

    cell_mask = adata.obs["cell_id"].isin(batch_cell_ids)
    adata = adata[cell_mask].copy()

    current_genes = adata.var_names.tolist()
    if gene_names is None:
      gene_names = current_genes
    else:
      assert current_genes == gene_names, (
        f"Gene mismatch in batch '{batch}': expected {len(gene_names)} genes, "
        f"got {len(current_genes)}")

    X_norm = normalize_counts(adata.X, scale_factor=10000)
    matrices.append(X_norm)
    cell_ids_ordered.extend(adata.obs["cell_id"].tolist())

    del adata

  X = sp.vstack(matrices, format="csr")
  print(f"  Final matrix: {X.shape[0]} cells x {X.shape[1]} genes")

  return X, gene_names, cell_ids_ordered

# xgboost training and feature selection

def _build_xgb_params(
  n_cores: int, seed: int, use_gpu: bool, num_class: int, **overrides,
) -> dict:
  """Build XGBoost parameter dict. Adds GPU params if use_gpu=True."""
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


def run_xgboost_pass1(
  X_train: sp.csr_matrix, y_train: np.ndarray,
  X_val: sp.csr_matrix, y_val: np.ndarray,
  n_cores: int, seed: int, use_gpu: bool,
  subsample: float, colsample_bytree: float, learning_rate: float,
  nrounds: int, early_stopping: int,
) -> xgb.Booster:
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


def select_features_by_gain(
  model: xgb.Booster, gene_names: list[str],
  gain_threshold: float, min_genes: int, max_genes: int,
) -> tuple[np.ndarray, list[str], pl.DataFrame]:

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
    "gene": [gene_names[i] for i in sorted_idx],
    "gain": sorted_gains.tolist(),
    "cumulative_gain_frac": cum_gain.tolist(),
  })

  print(f"  Total genes with non-zero gain: {len(sorted_idx)}")
  print(f"  Genes at {gain_threshold*100:.0f}% cumulative gain: {n_sel}")
  print(f"  Top 10 genes: {sel_gene_names[:10]}")

  return sel_indices, sel_gene_names, gene_importance


def run_xgboost_pass2(
  X_train: sp.csr_matrix, y_train: np.ndarray,
  X_val: sp.csr_matrix, y_val: np.ndarray,
  n_cores: int, seed: int, use_gpu: bool,
  subsample: float, colsample_bytree: float, learning_rate: float,
  nrounds: int, early_stopping: int,
) -> xgb.Booster:
  
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

def evaluate_model(
  model: xgb.Booster, X_val: sp.csr_matrix,
  y_val: np.ndarray, class_names: list[str],
) -> None:
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


def evaluate_model_coarse(
  model: xgb.Booster, X_val: sp.csr_matrix, y_val: np.ndarray,
  class_names: list[str], label_map: dict[str, str],
) -> None:
  
  feat_names = [f"g{i}" for i in range(X_val.shape[1])]
  dval = xgb.DMatrix(X_val, feature_names=feat_names)
  probs = model.predict(dval)
  y_pred_fine = probs.argmax(axis=1)

  def to_coarse(fine_idx: int) -> str:
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

def save_outputs(
  model: xgb.Booster, class_names: list[str], selected_genes: list[str],
  gene_importance: pl.DataFrame,
  output_dir: str, ref_tag: str, cells_df: pl.DataFrame,
  label_map: Optional[dict[str, str]] = None,
  X_train: sp.csr_matrix = None, y_train: np.ndarray = None,
  X_val: sp.csr_matrix = None, y_val: np.ndarray = None,
) -> None:
  """Save model, gene list, classes, importance, predictions, and label map."""
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

  # Save predictions for all cells (train + validation) — consumed by Rmd report
  if X_val is not None and X_train is not None:
    preds_path = out_dir / f"{ref_tag}_predictions.csv.gz"
    preds_df = _make_predictions_df(
      model, X_train, y_train, X_val, y_val, cells_df, class_names, label_map
    )
    preds_df.write_csv(str(preds_path), compression="gzip")
    print(f"  Predictions: {preds_path}")


def _make_predictions_df(
  model: xgb.Booster,
  X_train: sp.csr_matrix, y_train: np.ndarray,
  X_val: sp.csr_matrix, y_val: np.ndarray,
  cells_df: pl.DataFrame, class_names: list[str],
  label_map: Optional[dict[str, str]] = None,
) -> pl.DataFrame:
  """Generate predictions for all cells and return as a DataFrame."""
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

  result = cells_df.select(["cell_id", "sample_id", "label", "split"]).with_columns([
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


if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    description="Train XGBoost cell type classifier from scprocess output"
  )
  parser.add_argument("--annots_f",          type=str, required=True)
  parser.add_argument("--cluster_csv",       type=str, required=True)
  parser.add_argument("--h5ads_yaml",        type=str, required=True)
  parser.add_argument("--ref_tag",           type=str, required=True)
  parser.add_argument("--output_dir",        type=str, required=True)
  parser.add_argument("--batch_var",         type=str, required=True)
  parser.add_argument("--int_res_ls",        type=str, required=True)
  parser.add_argument("--label_map_f",       type=str, default=None)
  parser.add_argument("--refine_labels",     type=str, required=True)
  parser.add_argument("--purity_threshold",  type=float, required=True)
  parser.add_argument("--n_cells_per_type",  type=int, required=True)
  parser.add_argument("--min_cells_per_type", type=int, required=True)
  parser.add_argument("--min_cells_expressed", type=int, required=True)
  parser.add_argument("--seed",             type=int, required=True)
  parser.add_argument("--n_cores",          type=int, required=True)
  parser.add_argument("--use_gpu",          action="store_true")
  parser.add_argument("--pass1_subsample",      type=float, required=True)
  parser.add_argument("--pass1_colsample_bytree", type=float, required=True)
  parser.add_argument("--pass1_learning_rate",  type=float, required=True)
  parser.add_argument("--pass1_nrounds",        type=int, required=True)
  parser.add_argument("--pass1_early_stopping", type=int, required=True)
  parser.add_argument("--pass2_colsample_bytree", type=float, required=True)
  parser.add_argument("--pass2_learning_rate",  type=float, required=True)
  parser.add_argument("--pass2_nrounds",        type=int, required=True)
  parser.add_argument("--pass2_early_stopping", type=int, required=True)
  parser.add_argument("--gain_threshold",   type=float, required=True)
  parser.add_argument("--min_genes",        type=int, required=True)
  parser.add_argument("--max_genes",        type=int, required=True)

  args = parser.parse_args()

  make_classifier(
    annots_f=args.annots_f,
    cluster_csv=args.cluster_csv,
    h5ads_yaml=args.h5ads_yaml,
    ref_tag=args.ref_tag,
    output_dir=args.output_dir,
    batch_var=args.batch_var,
    int_res_ls=args.int_res_ls,
    label_map_f=args.label_map_f,
    refine_labels=args.refine_labels.lower() == "true",
    purity_threshold=args.purity_threshold,
    n_cells_per_type=args.n_cells_per_type,
    min_cells_per_type=args.min_cells_per_type,
    min_cells_expressed=args.min_cells_expressed,
    seed=args.seed,
    n_cores=args.n_cores,
    use_gpu=args.use_gpu,
    pass1_subsample=args.pass1_subsample,
    pass1_colsample_bytree=args.pass1_colsample_bytree,
    pass1_learning_rate=args.pass1_learning_rate,
    pass1_nrounds=args.pass1_nrounds,
    pass1_early_stopping=args.pass1_early_stopping,
    pass2_colsample_bytree=args.pass2_colsample_bytree,
    pass2_learning_rate=args.pass2_learning_rate,
    pass2_nrounds=args.pass2_nrounds,
    pass2_early_stopping=args.pass2_early_stopping,
    gain_threshold=args.gain_threshold,
    min_genes=args.min_genes,
    max_genes=args.max_genes,
  )
