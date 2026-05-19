"""Train an XGBoost classifier to predict cell type labels from scRNA-seq data.

Workflow:
  1. Plan training cells from cluster CSV + annotations (no H5AD loaded)
  2. Load expression data one H5AD at a time, subset to selected cells
  3. XGBoost pass 1: broad exploration with all genes
  4. Feature selection via cumulative gain
  5. XGBoost pass 2: final model on selected genes
  6. Evaluate and save model + diagnostic plots

Usage:
  python scripts/train_xgboost.py path/to/train_config.yaml
"""

import argparse
import pathlib
import sys
from typing import Optional

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import scipy.sparse as sp
import xgboost as xgb
import yaml


# test files
annotation_f = '/projects/site/pred/neurogenomics/users/kodermam/scprocess/test_xgboost_annotation.csv.gz'
label_map_f  = '/projects/site/pred/neurogenomics/users/kodermam/scprocess/test_xgboost_label_map.csv.gz'
scprocess_config_f = '/pmount/projects/site/pred/brain-sc-analysis/configs/config-siletti_2023_hippocampus.yaml'

DEFAULTS = {
  "label_map_f": None,
  "refine_labels": True,
  "purity_threshold": 0.65,
  "n_cells_per_type": 1000,
  "min_cells_per_type": 20,
  "seed": 42,
  "n_cores": 16,
  "min_cells_expressed": 10,
  "use_gpu": False,
  # XGBoost Pass 1 — low colsample_bytree forces broad gene exploration
  "pass1_subsample": 0.632,
  "pass1_colsample_bytree": 0.1,
  "pass1_learning_rate": 0.1,
  "pass1_nrounds": 300,
  "pass1_early_stopping": 10,
  # XGBoost Pass 2 — focused training on selected genes
  "pass2_colsample_bytree": 0.5,
  "pass2_learning_rate": 0.05,
  "pass2_nrounds": 500,
  "pass2_early_stopping": 10,
  # Feature selection — genes contributing to top fraction of cumulative gain
  "gain_threshold": 0.9,
  "min_genes": 100,
  "max_genes": 3000,
}


def make_classifier(yaml_f: str) -> None:
  """Main entry point. Orchestrates the full training pipeline."""
  print("=" * 60)
  print("XGBoost cell type classifier training")
  print("=" * 60)

  config = load_config(yaml_f)
  paths = resolve_scprocess_paths(config)

  # Label mapping is not used for training — saved alongside model for prediction time
  label_map = load_label_mapping(config)
  if label_map is not None:
    print(f"  Label mapping loaded: {len(label_map)} fine→coarse entries")
    print(f"  Coarse categories: {sorted(set(label_map.values()))}")

  # Phase 1: plan which cells to use (CSV only, no H5AD)
  print("\n--- Planning training cells ---")
  cells_df = plan_training_cells(config, paths)

  # Phase 2: load expression data one H5AD at a time
  print("\n--- Loading expression data ---")
  X, gene_names, cell_ids = load_expression_matrix(cells_df, paths)

  min_cells_expressed = config.get("min_cells_expressed", 10)
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
  model_pass1 = run_xgboost_pass1(X_train, y_train, X_val, y_val, config)

  # Phase 4: feature selection via gain scores
  print("\n--- Feature selection ---")
  sel_indices, sel_gene_names, gene_importance = select_features_by_gain(
    model_pass1, gene_names, config
  )
  print(f"  Selected {len(sel_indices)} genes (of {len(gene_names)})")

  X_train_sub = X_train[:, sel_indices]
  X_val_sub = X_val[:, sel_indices]

  # Phase 5: XGBoost pass 2 — final training on curated gene set
  print("\n--- XGBoost Pass 2 (selected genes) ---")
  model_pass2 = run_xgboost_pass2(X_train_sub, y_train, X_val_sub, y_val, config)

  # Phase 6: evaluate and save
  print("\n--- Evaluation (fine labels) ---")
  evaluate_model(model_pass2, X_val_sub, y_val, class_names)

  if label_map is not None:
    print("\n--- Evaluation (coarse labels) ---")
    evaluate_model_coarse(model_pass2, X_val_sub, y_val, class_names, label_map)

  print("\n--- Saving outputs ---")
  save_outputs(
    model_pass2, class_names, sel_gene_names, gene_importance,
    config, cells_df, paths, label_map
  )

  print("\nDone.")


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------


def load_config(yaml_f: str) -> dict:
  """Load training config YAML + referenced scprocess config. Fill defaults."""
  yaml_path = pathlib.Path(yaml_f)
  assert yaml_path.is_file(), f"Config file not found: {yaml_f}"

  with open(yaml_path) as f:
    config = yaml.safe_load(f)

  required = ["scprocess_config_f", "annots_f", "output_dir", "ref_tag"]
  for key in required:
    assert key in config, f"Required key '{key}' missing from config"

  for key, default in DEFAULTS.items():
    if key not in config:
      config[key] = default

  assert pathlib.Path(config["annots_f"]).is_file(), (
    f"Annotations file not found: {config['annots_f']}")
  assert pathlib.Path(config["scprocess_config_f"]).is_file(), (
    f"scprocess config not found: {config['scprocess_config_f']}")
  if config["label_map_f"] is not None:
    assert pathlib.Path(config["label_map_f"]).is_file(), (
      f"Label map file not found: {config['label_map_f']}")

  with open(config["scprocess_config_f"]) as f:
    config["scprocess"] = yaml.safe_load(f)

  out_dir = pathlib.Path(config["output_dir"])
  out_dir.mkdir(parents=True, exist_ok=True)
  (out_dir / "plots").mkdir(exist_ok=True)

  return config


def resolve_scprocess_paths(config: dict) -> dict:
  """Build paths to scprocess integration outputs."""
  scp = config["scprocess"]
  proj_dir = scp["project"]["proj_dir"]
  short_tag = scp["project"]["short_tag"]
  full_tag = scp["project"]["full_tag"]
  date_stamp = scp["project"]["date_stamp"]
  int_dir = f"{proj_dir}/output/{short_tag}_integration"

  paths = {
    "cluster_csv": f"{int_dir}/integrated_dt_{full_tag}_{date_stamp}.csv.gz",
    "h5ads_yaml": f"{int_dir}/h5ads_clean_paths_{full_tag}_{date_stamp}.yaml",
    "int_dir": int_dir,
    "batch_var": scp.get("integration", {}).get("int_batch_var", "sample_id"),
  }

  assert pathlib.Path(paths["cluster_csv"]).is_file(), (
    f"Cluster CSV not found: {paths['cluster_csv']}")
  assert pathlib.Path(paths["h5ads_yaml"]).is_file(), (
    f"H5AD paths YAML not found: {paths['h5ads_yaml']}")

  with open(paths["h5ads_yaml"]) as f:
    paths["h5ad_dict"] = yaml.safe_load(f)

  return paths


# ---------------------------------------------------------------------------
# Planning phase (CSV only)
# ---------------------------------------------------------------------------


def plan_training_cells(config: dict, paths: dict) -> pl.DataFrame:
  """Decide which cells to use for training/validation using only CSV files.

  Label mapping (fine→coarse) is NOT applied here — the model trains on
  fine-grained labels. The mapping is saved alongside the model for prediction.
  """
  batch_var = paths["batch_var"]

  # Highest clustering resolution for label refinement (typically RNA_snn_res.2)
  res_ls = config["scprocess"].get("integration", {}).get("int_res_ls", [0.1, 0.2, 0.5, 1, 2])
  hi_res_col = f"RNA_snn_res.{max(res_ls)}"

  cols_to_read = ["cell_id", batch_var, hi_res_col]
  cluster_df = pl.read_csv(paths["cluster_csv"], columns=cols_to_read)
  print(f"  Loaded cluster CSV: {cluster_df.shape[0]} cells")

  annots_df = pl.read_csv(config["annots_f"])
  assert "cell_id" in annots_df.columns, "annots_f must have 'cell_id' column"
  assert "annotation" in annots_df.columns, "annots_f must have 'annotation' column"
  annots_df = annots_df.select(["cell_id", "annotation"])

  df = cluster_df.join(annots_df, on="cell_id", how="left")
  n_annotated = df.filter(pl.col("annotation").is_not_null()).shape[0]
  print(f"  Cells with annotations: {n_annotated} / {df.shape[0]}")

  if config["refine_labels"]:
    print("  Refining labels by cluster majority voting...")
    df = refine_labels_by_cluster(df, hi_res_col, config)
  else:
    df = df.with_columns(pl.col("annotation").alias("label"))

  df = df.filter(pl.col("label").is_not_null())
  print(f"  Cells with labels after processing: {df.shape[0]}")

  # Exclude rare cell types
  type_counts = df.group_by("label").len()
  keep_types = (
    type_counts
    .filter(pl.col("len") >= config["min_cells_per_type"])
    ["label"].to_list()
  )
  excluded = type_counts.filter(pl.col("len") < config["min_cells_per_type"])
  if excluded.shape[0] > 0:
    print(f"  Excluding {excluded.shape[0]} cell types with < {config['min_cells_per_type']} cells:")
    for row in excluded.iter_rows(named=True):
      print(f"    {row['label']}: {row['len']} cells")
  df = df.filter(pl.col("label").is_in(keep_types))

  print("  Downsampling...")
  df = downsample_per_type(df, config)

  print("  Assigning train/val split...")
  df = assign_train_val_split(df, config)

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
  df: pl.DataFrame, hi_res_col: str, config: dict
) -> pl.DataFrame:
  """Smooth annotations via cluster majority voting.

  For each high-res cluster, if majority annotation >= purity_threshold and
  cluster has >= 10 annotated cells, relabel all cells. Otherwise keep originals.
  """
  purity_thr = config["purity_threshold"]
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
      (pl.col("purity") >= purity_thr) & (pl.col("n_total") >= min_cluster_size)
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


def load_label_mapping(config: dict) -> Optional[dict[str, str]]:
  """Load fine→coarse label mapping CSV. Not applied during training —
  saved alongside model and applied at prediction time.

  Expected CSV columns: annotation, coarse_label
  """
  if config["label_map_f"] is None:
    return None

  map_df = pl.read_csv(config["label_map_f"])
  assert "annotation" in map_df.columns, "label_map_f must have 'annotation' column"
  assert "coarse_label" in map_df.columns, "label_map_f must have 'coarse_label' column"

  return dict(zip(
    map_df["annotation"].to_list(),
    map_df["coarse_label"].to_list(),
  ))


def downsample_per_type(df: pl.DataFrame, config: dict) -> pl.DataFrame:
  """Downsample to at most n_cells_per_type cells per label."""
  n_max = config["n_cells_per_type"]
  seed = config["seed"]

  sampled = (
    df
    .with_columns(pl.lit(1).alias("_dummy"))
    .group_by("label")
    .map_groups(
      lambda group: group.sample(n=min(n_max, group.shape[0]), seed=seed)
    )
    .drop("_dummy")
  )

  type_summary = sampled.group_by("label").len().sort("len", descending=True)
  print(f"    {sampled.shape[0]} cells after downsampling ({type_summary.shape[0]} types)")

  return sampled


def assign_train_val_split(df: pl.DataFrame, config: dict) -> pl.DataFrame:
  """Assign cells to train/validation. Sample-level holdout if >4 samples,
  otherwise stratified random split at cell level.
  """
  seed = config["seed"]
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


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


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
  cells_df: pl.DataFrame, paths: dict
) -> tuple[sp.csr_matrix, list[str], list[str]]:
  """Load and normalize expression data, one H5AD at a time.

  Peak memory is proportional to one batch — each H5AD is loaded, subsetted
  to selected cells, normalized, and freed before loading the next.
  """
  h5ad_dict = paths["h5ad_dict"]
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


# ---------------------------------------------------------------------------
# XGBoost training
# ---------------------------------------------------------------------------


def _build_xgb_params(config: dict, num_class: int, **overrides) -> dict:
  """Build XGBoost parameter dict. Adds GPU params if use_gpu=True."""
  params = {
    "objective": "multi:softprob",
    "num_class": num_class,
    "eval_metric": "mlogloss",
    "nthread": config["n_cores"],
    "seed": config["seed"],
  }
  if config.get("use_gpu", False):
    params["device"] = "cuda"
    params["tree_method"] = "hist"
  params.update(overrides)
  return params


def run_xgboost_pass1(
  X_train: sp.csr_matrix, y_train: np.ndarray,
  X_val: sp.csr_matrix, y_val: np.ndarray, config: dict,
) -> xgb.Booster:
  """Pass 1: broad exploration with all genes (low colsample_bytree)."""
  num_class = int(y_train.max()) + 1

  params = _build_xgb_params(config, num_class,
    subsample=config["pass1_subsample"],
    colsample_bytree=config["pass1_colsample_bytree"],
    learning_rate=config["pass1_learning_rate"],
  )

  dtrain = xgb.DMatrix(X_train, label=y_train)
  dval = xgb.DMatrix(X_val, label=y_val)

  model = xgb.train(
    params=params, dtrain=dtrain,
    num_boost_round=config["pass1_nrounds"],
    evals=[(dtrain, "train"), (dval, "val")],
    early_stopping_rounds=config["pass1_early_stopping"],
    verbose_eval=50,
  )

  print(f"  Best iteration: {model.best_iteration}")
  return model


def select_features_by_gain(
  model: xgb.Booster, gene_names: list[str], config: dict
) -> tuple[np.ndarray, list[str], pl.DataFrame]:
  """Select genes contributing to top gain_threshold fraction of cumulative gain.
  Clamped to [min_genes, max_genes].
  """
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

  n_sel = int(np.searchsorted(cum_gain, config["gain_threshold"])) + 1
  n_sel = max(n_sel, config["min_genes"])
  n_sel = min(n_sel, config["max_genes"])
  n_sel = min(n_sel, len(sorted_idx))

  sel_indices = sorted_idx[:n_sel]
  sel_gene_names = [gene_names[i] for i in sel_indices]

  gene_importance = pl.DataFrame({
    "gene": [gene_names[i] for i in sorted_idx],
    "gain": sorted_gains.tolist(),
    "cumulative_gain_frac": cum_gain.tolist(),
  })

  print(f"  Total genes with non-zero gain: {len(sorted_idx)}")
  print(f"  Genes at {config['gain_threshold']*100:.0f}% cumulative gain: {n_sel}")
  print(f"  Top 10 genes: {sel_gene_names[:10]}")

  return sel_indices, sel_gene_names, gene_importance


def run_xgboost_pass2(
  X_train: sp.csr_matrix, y_train: np.ndarray,
  X_val: sp.csr_matrix, y_val: np.ndarray, config: dict,
) -> xgb.Booster:
  """Pass 2: fresh retrain on selected genes with relaxed parameters."""
  num_class = int(y_train.max()) + 1

  params = _build_xgb_params(config, num_class,
    subsample=config["pass1_subsample"],
    colsample_bytree=config["pass2_colsample_bytree"],
    learning_rate=config["pass2_learning_rate"],
  )

  feat_names = [f"g{i}" for i in range(X_train.shape[1])]
  dtrain = xgb.DMatrix(X_train, label=y_train, feature_names=feat_names)
  dval = xgb.DMatrix(X_val, label=y_val, feature_names=feat_names)

  model = xgb.train(
    params=params, dtrain=dtrain,
    num_boost_round=config["pass2_nrounds"],
    evals=[(dtrain, "train"), (dval, "val")],
    early_stopping_rounds=config["pass2_early_stopping"],
    verbose_eval=50,
  )

  print(f"  Best iteration: {model.best_iteration}")
  return model


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------


def evaluate_model(
  model: xgb.Booster, X_val: sp.csr_matrix,
  y_val: np.ndarray, class_names: list[str],
) -> None:
  """Print per-class accuracy, F1, and top misclassifications on validation set."""
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
  """Evaluate after collapsing fine predictions to coarse labels via mapping."""
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


# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------


def save_outputs(
  model: xgb.Booster, class_names: list[str], selected_genes: list[str],
  gene_importance: pl.DataFrame, config: dict, cells_df: pl.DataFrame,
  paths: dict, label_map: Optional[dict[str, str]] = None,
) -> None:
  """Save model (.json), allowed classes, selected genes, importance, and plots."""
  out_dir = pathlib.Path(config["output_dir"])
  ref_tag = config["ref_tag"]

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

  print("  Generating plots...")
  make_diagnostic_plots(cells_df, model, class_names, paths, config)


def make_diagnostic_plots(
  cells_df: pl.DataFrame, model: xgb.Booster,
  class_names: list[str], paths: dict, config: dict,
) -> None:
  """UMAP plots using pre-computed coordinates from scprocess."""
  out_dir = pathlib.Path(config["output_dir"]) / "plots"
  ref_tag = config["ref_tag"]

  umap_df = pl.read_csv(paths["cluster_csv"], columns=["cell_id", "UMAP1", "UMAP2"])

  plot_df = cells_df.join(umap_df, on="cell_id", how="left")
  plot_df = plot_df.filter(
    pl.col("UMAP1").is_not_null() & pl.col("UMAP2").is_not_null()
  )

  labels = sorted(plot_df["label"].unique().to_list())
  cmap = plt.cm.get_cmap("tab20", len(labels))
  label_colors = {lbl: cmap(i) for i, lbl in enumerate(labels)}

  # Plot 1: UMAP colored by label
  fig, ax = plt.subplots(1, 1, figsize=(10, 8))
  for lbl in labels:
    subset = plot_df.filter(pl.col("label") == lbl)
    ax.scatter(
      subset["UMAP1"].to_numpy(), subset["UMAP2"].to_numpy(),
      c=[label_colors[lbl]], s=1, alpha=0.5, label=lbl, rasterized=True,
    )
  ax.set_xlabel("UMAP1")
  ax.set_ylabel("UMAP2")
  ax.set_title("True labels")
  ax.legend(markerscale=5, fontsize=7, loc="center left", bbox_to_anchor=(1, 0.5))
  ax.set_aspect("equal")
  ax.set_xticks([])
  ax.set_yticks([])
  plt.tight_layout()
  plt.savefig(out_dir / f"{ref_tag}_umap_true_labels.png", dpi=150, bbox_inches="tight")
  plt.close()

  # Plot 2: UMAP colored by train/validation split
  fig, ax = plt.subplots(1, 1, figsize=(8, 7))
  split_colors = {"train": "#1f77b4", "validation": "#ff7f0e"}
  for split in ["train", "validation"]:
    subset = plot_df.filter(pl.col("split") == split)
    ax.scatter(
      subset["UMAP1"].to_numpy(), subset["UMAP2"].to_numpy(),
      c=split_colors[split], s=1, alpha=0.4, label=split, rasterized=True,
    )
  ax.set_xlabel("UMAP1")
  ax.set_ylabel("UMAP2")
  ax.set_title("Train / Validation split")
  ax.legend(markerscale=5)
  ax.set_aspect("equal")
  ax.set_xticks([])
  ax.set_yticks([])
  plt.tight_layout()
  plt.savefig(out_dir / f"{ref_tag}_umap_train_val_split.png", dpi=150, bbox_inches="tight")
  plt.close()

  print(f"    Saved plots to {out_dir}/")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


if __name__ == "__main__" and len(sys.argv) > 1:
  parser = argparse.ArgumentParser(
    description="Train XGBoost cell type classifier from scprocess output"
  )
  parser.add_argument("config_yaml", type=str, help="Path to training config YAML")
  args = parser.parse_args()
  make_classifier(args.config_yaml)
