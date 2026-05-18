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


# Default values for optional config parameters. 
DEFAULTS = {
    # Label refinement
    "label_map_f": None,            # optional CSV mapping fine→coarse labels
    "refine_labels": True,          # whether to smooth labels via cluster majority voting
    "purity_threshold": 0.65,       # min fraction of majority annotation to relabel a cluster

    # Downsampling
    "n_cells_per_type": 1000,       # max cells per cell type in training set
    "min_cells_per_type": 20,       # cell types below this are excluded entirely

    # General
    "seed": 42,
    "n_cores": 16,

    # XGBoost Pass 1 — broad exploration across all genes
    # Low colsample_bytree forces the model to explore many genes per tree,
    # effectively discovering which genes are informative for classification.
    "pass1_subsample": 0.632,       # bootstrap fraction of cells per tree
    "pass1_colsample_bytree": 0.1,  # fraction of genes sampled per tree (low = broad)
    "pass1_learning_rate": 0.1,
    "pass1_nrounds": 300,
    "pass1_early_stopping": 10,

    # XGBoost Pass 2 — focused training on selected genes
    # Higher colsample_bytree since the gene set is already curated.
    # Lower learning rate for a tighter fit.
    "pass2_colsample_bytree": 0.5,
    "pass2_learning_rate": 0.05,
    "pass2_nrounds": 500,
    "pass2_early_stopping": 10,

    # Feature selection between passes
    # Selects genes that contribute to the top fraction of total gain from pass 1.
    # This is adaptive: complex datasets select more genes, simple ones fewer.
    "gain_threshold": 0.9,          # cumulative gain fraction cutoff
    "min_genes": 100,               # floor on number of selected genes
    "max_genes": 3000,              # ceiling on number of selected genes
}


def make_classifier(yaml_f: str) -> None:
    """Main entry point. Orchestrates the full training pipeline.

    The pipeline has two major phases designed for memory efficiency:
      - Planning phase: uses only lightweight CSV files to decide which cells
        go into training vs. validation, apply label refinement, and downsample.
        No H5AD files are loaded during this phase.
      - Data phase: loads H5AD files one batch at a time, subsets to only the
        cells selected in the planning phase, normalizes, and appends to a
        sparse matrix. This keeps peak memory proportional to one batch.

    After data loading, two passes of XGBoost are run:
      - Pass 1 uses ALL genes with aggressive column subsampling (colsample=0.1)
        to discover which genes are informative across the full transcriptome.
      - Feature selection picks genes contributing to 90% of cumulative gain.
      - Pass 2 retrains from scratch on the selected gene subset with relaxed
        parameters for a tighter final model.
    """
    print("=" * 60)
    print("XGBoost cell type classifier training")
    print("=" * 60)

    config = load_config(yaml_f)
    paths  = resolve_scprocess_paths(config)

    # -------------------------------------------------------------------------
    # Phase 1: plan which cells to use (CSV only, no H5AD)
    # This decides label refinement, downsampling, and train/val split using
    # only the cluster CSV and annotations file.
    # -------------------------------------------------------------------------
    print("\n--- Planning training cells ---")
    cells_df = plan_training_cells(config, paths)

    # -------------------------------------------------------------------------
    # Phase 2: load expression data one H5AD at a time
    # Each batch's H5AD is loaded, subsetted to selected cells, normalized
    # (total-count to 10k + log1p), and freed before loading the next.
    # -------------------------------------------------------------------------
    print("\n--- Loading expression data ---")
    X, gene_names, cell_ids = load_expression_matrix(cells_df, paths)

    # Align cells_df row order with the matrix row order
    cells_df = cells_df.filter(pl.col("cell_id").is_in(cell_ids))
    id_order = pl.DataFrame({"cell_id": cell_ids, "_row_idx": range(len(cell_ids))})
    cells_df = cells_df.join(id_order, on="cell_id").sort("_row_idx").drop("_row_idx")

    # Encode string labels as integers for XGBoost (sorted alphabetically)
    class_names = sorted(cells_df["label"].unique().to_list())
    label_to_int = {name: i for i, name in enumerate(class_names)}
    y = np.array([label_to_int[lbl] for lbl in cells_df["label"].to_list()])

    # Split into train/val matrices using the pre-assigned split column
    train_mask = cells_df["split"].to_numpy() == "train"
    val_mask = ~train_mask

    X_train, y_train = X[train_mask], y[train_mask]
    X_val, y_val = X[val_mask], y[val_mask]

    print(f"  Train: {X_train.shape[0]} cells, Val: {X_val.shape[0]} cells")
    print(f"  Features (all genes): {X_train.shape[1]}")
    print(f"  Classes: {len(class_names)}")

    # -------------------------------------------------------------------------
    # Phase 3: XGBoost pass 1 — broad gene exploration
    # With colsample_bytree=0.1, each tree only sees 10% of genes, forcing the
    # ensemble to explore the full transcriptome across many trees. This lets us
    # identify informative genes without prior HVG selection.
    # -------------------------------------------------------------------------
    print("\n--- XGBoost Pass 1 (all genes) ---")
    model_pass1 = run_xgboost_pass1(X_train, y_train, X_val, y_val, config)

    # -------------------------------------------------------------------------
    # Phase 4: feature selection via gain scores
    # Genes are ranked by their total gain (reduction in loss) across all trees.
    # We keep genes contributing to the top 90% of cumulative gain. This is
    # adaptive: complex datasets with many informative genes keep more; simple
    # datasets keep fewer.
    # -------------------------------------------------------------------------
    print("\n--- Feature selection ---")
    sel_indices, sel_gene_names, gene_importance = select_features_by_gain(
        model_pass1, gene_names, config
    )
    print(f"  Selected {len(sel_indices)} genes (of {len(gene_names)})")

    # Subset expression matrices to selected genes only
    X_train_sub = X_train[:, sel_indices]
    X_val_sub = X_val[:, sel_indices]

    # -------------------------------------------------------------------------
    # Phase 5: XGBoost pass 2 — final training on curated gene set
    # Now that we have a refined feature set, we retrain from scratch with:
    #   - Higher colsample (0.5) since genes are already informative
    #   - Lower learning rate (0.05) for finer convergence
    #   - More rounds (500) to compensate for the lower learning rate
    # -------------------------------------------------------------------------
    print("\n--- XGBoost Pass 2 (selected genes) ---")
    model_pass2 = run_xgboost_pass2(X_train_sub, y_train, X_val_sub, y_val, config)

    # -------------------------------------------------------------------------
    # Phase 6: evaluate and save outputs
    # -------------------------------------------------------------------------
    print("\n--- Evaluation ---")
    evaluate_model(model_pass2, X_val_sub, y_val, class_names)

    print("\n--- Saving outputs ---")
    save_outputs(
        model_pass2, class_names, sel_gene_names, gene_importance,
        config, cells_df, paths
    )

    print("\nDone.")


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------


def load_config(yaml_f: str) -> dict:
    """Load and validate training config + referenced scprocess config.

    The training config YAML must contain:
      - scprocess_config_f: path to the scprocess pipeline config YAML
        (this gives us proj_dir, short_tag, full_tag, date_stamp to locate outputs)
      - annots_f: path to annotations CSV (columns: cell_id, annotation)
      - output_dir: directory where model and plots will be saved
      - ref_tag: short name for this model (used in output filenames)

    Optional keys (see DEFAULTS dict for values):
      - label_map_f: CSV mapping fine annotations to coarse labels
      - refine_labels, purity_threshold: label smoothing parameters
      - n_cells_per_type, min_cells_per_type: downsampling parameters
      - pass1_*/pass2_*/gain_*: XGBoost and feature selection parameters
    """
    yaml_path = pathlib.Path(yaml_f)
    assert yaml_path.is_file(), f"Config file not found: {yaml_f}"

    with open(yaml_path) as f:
        config = yaml.safe_load(f)

    # check required keys
    required = ["scprocess_config_f", "annots_f", "output_dir", "ref_tag"]
    for key in required:
        assert key in config, f"Required key '{key}' missing from config"

    # fill defaults for any unspecified optional parameters
    for key, default in DEFAULTS.items():
        if key not in config:
            config[key] = default

    # Validate that input files exist
    assert pathlib.Path(config["annots_f"]).is_file(), (
        f"Annotations file not found: {config['annots_f']}"
    )
    assert pathlib.Path(config["scprocess_config_f"]).is_file(), (
        f"scprocess config not found: {config['scprocess_config_f']}"
    )
    if config["label_map_f"] is not None:
        assert pathlib.Path(config["label_map_f"]).is_file(), (
            f"Label map file not found: {config['label_map_f']}"
        )

    # load the scprocess config
    with open(config["scprocess_config_f"]) as f:
        config["scprocess"] = yaml.safe_load(f)

    # create output directory structure
    out_dir = pathlib.Path(config["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "plots").mkdir(exist_ok=True)

    return config


def resolve_scprocess_paths(config: dict) -> dict:
    """Build paths to scprocess integration outputs from the pipeline config.

    scprocess outputs are located at:
      {proj_dir}/output/{short_tag}_integration/

    Key files:
      - integrated_dt_{full_tag}_{date_stamp}.csv.gz
          Cluster assignments, UMAP coords, and sample IDs for all cells.
      - h5ads_clean_paths_{full_tag}_{date_stamp}.yaml
          YAML dict mapping batch_name → H5AD file path.

    The 'batch_var' (int_batch_var in scprocess config) determines how cells
    are split across H5AD files: if "sample_id", one H5AD per sample; if
    "pool_id", one H5AD per pool (which may contain multiple samples).
    """
    scp        = config["scprocess"]
    proj_dir   = scp["project"]["proj_dir"]
    short_tag  = scp["project"]["short_tag"]
    full_tag   = scp["project"]["full_tag"]
    date_stamp = scp["project"]["date_stamp"]
    int_dir    = f"{proj_dir}/output/{short_tag}_integration"

    paths = {
      "cluster_csv": f"{int_dir}/integrated_dt_{full_tag}_{date_stamp}.csv.gz",
      "h5ads_yaml": f"{int_dir}/h5ads_clean_paths_{full_tag}_{date_stamp}.yaml",
      "int_dir": int_dir,
      "batch_var": scp.get("integration", {}).get("int_batch_var", "sample_id"),
    }

    # validate that scprocess has been run and outputs exist
    assert pathlib.Path(paths["cluster_csv"]).is_file(), (
        f"Cluster CSV not found: {paths['cluster_csv']}"
    )
    assert pathlib.Path(paths["h5ads_yaml"]).is_file(), (
        f"H5AD paths YAML not found: {paths['h5ads_yaml']}"
    )

    # load the batch→h5ad mapping (e.g. {"sample_A": "/path/to/sample_A.h5ad", ...})
    with open(paths["h5ads_yaml"]) as f:
        paths["h5ad_dict"] = yaml.safe_load(f)

    return paths


# ---------------------------------------------------------------------------
# Planning phase (CSV only)
# ---------------------------------------------------------------------------


def plan_training_cells(config: dict, paths: dict) -> pl.DataFrame:
    """Decide which cells to use for training/validation using only CSV files.

    This is the "planning phase" — it operates entirely on lightweight tabular
    data (the cluster CSV and annotations CSV). No H5AD files are loaded here.

    Steps:
      1. Load the cluster CSV (has cell_id, sample_id, cluster assignments)
      2. Load annotations CSV and join to clusters
      3. (Optional) Refine labels via cluster majority voting
      4. (Optional) Apply a fine→coarse label mapping
      5. Exclude cell types with too few cells
      6. Downsample to n_cells_per_type per cell type
      7. Assign train/val split (sample-level holdout or stratified random)
      8. Map each cell to its batch (determines which H5AD to load later)

    Returns:
      DataFrame with columns: cell_id, sample_id, label, split, batch
    """
    batch_var = paths["batch_var"]

    # Use the highest clustering resolution from scprocess for label refinement.
    # scprocess defaults to [0.1, 0.2, 0.5, 1, 2], so this is typically res=2.
    res_ls = config["scprocess"].get("integration", {}).get("int_res_ls", [0.1, 0.2, 0.5, 1, 2])
    hi_res = max(res_ls)
    hi_res_col = f"RNA_snn_res.{hi_res}"

    # Load only the columns we need from the cluster CSV (can be large)
    cols_to_read = ["cell_id", "sample_id", hi_res_col]
    if batch_var != "sample_id" and batch_var not in cols_to_read:
        cols_to_read.append(batch_var)

    cluster_df = pl.read_csv(paths["cluster_csv"], columns=cols_to_read)
    print(f"  Loaded cluster CSV: {cluster_df.shape[0]} cells")

    # Load annotations — must have cell_id and annotation columns
    annots_df = pl.read_csv(config["annots_f"])
    assert "cell_id" in annots_df.columns, "annots_f must have 'cell_id' column"
    assert "annotation" in annots_df.columns, "annots_f must have 'annotation' column"
    annots_df = annots_df.select(["cell_id", "annotation"])

    # Left join: keep all cells from cluster CSV, add annotations where available.
    # Cells without annotations will have annotation=null.
    df = cluster_df.join(annots_df, on="cell_id", how="left")
    n_annotated = df.filter(pl.col("annotation").is_not_null()).shape[0]
    print(f"  Cells with annotations: {n_annotated} / {df.shape[0]}")

    # Step 3: Label refinement — smooth per-cell annotations using cluster consensus.
    # For each high-resolution cluster, if the majority annotation makes up >=65%
    # of annotated cells, all cells in that cluster get that label. Mixed clusters
    # keep their original per-cell annotations.
    if config["refine_labels"]:
        print("  Refining labels by cluster majority voting...")
        df = refine_labels_by_cluster(df, hi_res_col, config)
    else:
        df = df.with_columns(pl.col("annotation").alias("label"))

    # Step 4: Apply label mapping (e.g. "excitatory_L5" → "neuron")
    if config["label_map_f"] is not None:
        print("  Applying label mapping...")
        df = apply_label_mapping(df, config)

    # Drop cells that still have no label (were never annotated and not in a
    # pure-enough cluster to receive a refined label)
    df = df.filter(pl.col("label").is_not_null())
    print(f"  Cells with labels after processing: {df.shape[0]}")

    # Exclude rare cell types that don't have enough cells for meaningful training
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

    # Downsample each cell type to at most n_cells_per_type
    print("  Downsampling...")
    df = downsample_per_type(df, config)

    # Assign cells to train or validation split
    print("  Assigning train/val split...")
    df = assign_train_val_split(df, config)

    # Add batch column — this determines which H5AD file contains each cell.
    # When batch_var="sample_id", each sample has its own H5AD.
    # When batch_var="pool_id", each pool (potentially containing multiple samples)
    # has one H5AD.
    if batch_var == "sample_id":
        df = df.with_columns(pl.col("sample_id").alias("batch"))
    else:
        df = df.with_columns(pl.col(batch_var).alias("batch"))

    # Print summary
    split_summary = df.group_by("split").len()
    for row in split_summary.iter_rows(named=True):
        print(f"    {row['split']}: {row['len']} cells")

    return df.select(["cell_id", "sample_id", "label", "split", "batch"])


def refine_labels_by_cluster(
    df: pl.DataFrame, hi_res_col: str, config: dict
) -> pl.DataFrame:
    """Use highest-resolution clustering to smooth annotations via majority voting.

    Rationale: raw per-cell annotations can be noisy (misassigned cells, ambiguous
    boundaries). By looking at Harmony clusters at high resolution, we can identify
    groups of cells that are transcriptionally similar. If a cluster is dominated
    (>= purity_threshold) by one annotation, we trust that consensus and apply it
    to all cells in the cluster — including any that were mislabeled individually.

    Clusters that are mixed (below purity threshold) or too small (< 10 annotated
    cells) are left alone — their cells keep their original per-cell annotations.
    This avoids imposing false consensus on genuinely heterogeneous regions.
    """
    purity_thr = config["purity_threshold"]
    min_cluster_size = 10

    # Only use annotated cells to compute cluster purity
    annotated = df.filter(pl.col("annotation").is_not_null())

    # Count (cluster, annotation) pairs
    cluster_annot_counts = (
        annotated
        .group_by([hi_res_col, "annotation"])
        .len()
        .rename({"len": "n"})
    )

    # Total annotated cells per cluster
    cluster_totals = (
        annotated
        .group_by(hi_res_col)
        .len()
        .rename({"len": "n_total"})
    )

    # For each cluster, find the majority annotation and its purity (proportion).
    # Keep only clusters that pass both the purity threshold and minimum size.
    cluster_stats = (
        cluster_annot_counts
        .join(cluster_totals, on=hi_res_col)
        .with_columns((pl.col("n") / pl.col("n_total")).alias("purity"))
        .sort([hi_res_col, "purity"], descending=[False, True])
        .group_by(hi_res_col)
        .first()  # takes the top-purity annotation per cluster
        .filter(
            (pl.col("purity") >= purity_thr) & (pl.col("n_total") >= min_cluster_size)
        )
        .select([hi_res_col, pl.col("annotation").alias("majority_label")])
    )

    n_refined = cluster_stats.shape[0]
    n_total_clusters = df[hi_res_col].n_unique()
    print(f"    {n_refined}/{n_total_clusters} clusters pass purity threshold")

    # Apply: cells in pure clusters get the majority label; others keep original
    df = df.join(cluster_stats, on=hi_res_col, how="left")
    df = df.with_columns(
        pl.when(pl.col("majority_label").is_not_null())
        .then(pl.col("majority_label"))
        .otherwise(pl.col("annotation"))
        .alias("label")
    ).drop("majority_label")

    return df


def apply_label_mapping(df: pl.DataFrame, config: dict) -> pl.DataFrame:
    """Map fine-grained annotations to coarse labels.

    Reads a CSV with columns 'annotation' and 'coarse_label'. Any label that
    appears in the mapping gets replaced; labels not in the mapping are kept as-is.

    Example mapping CSV:
      annotation,coarse_label
      excitatory_L2/3,Neuron
      excitatory_L5,Neuron
      inhibitory_PV,Neuron
      Astrocyte,Astrocyte

    This allows training on diverse subtypes while predicting coarse categories.
    """
    map_df = pl.read_csv(config["label_map_f"])
    assert "annotation" in map_df.columns, "label_map_f must have 'annotation' column"
    assert "coarse_label" in map_df.columns, "label_map_f must have 'coarse_label' column"

    map_df = map_df.select(["annotation", "coarse_label"]).rename(
        {"annotation": "label", "coarse_label": "mapped_label"}
    )

    # Left join: labels in the map get a mapped_label; others get null
    df = df.join(map_df, on="label", how="left")
    df = df.with_columns(
        pl.when(pl.col("mapped_label").is_not_null())
        .then(pl.col("mapped_label"))
        .otherwise(pl.col("label"))
        .alias("label")
    ).drop("mapped_label")

    return df


def downsample_per_type(df: pl.DataFrame, config: dict) -> pl.DataFrame:
    """Downsample to at most n_cells_per_type cells per label.

    This ensures the training set is balanced and not too large. Cell types with
    fewer cells than the cap are kept entirely (no upsampling).
    """
    n_max = config["n_cells_per_type"]
    seed = config["seed"]

    sampled = (
        df
        .with_columns(pl.lit(1).alias("_dummy"))
        .group_by("label")
        .map_groups(
            lambda group: group.sample(
                n=min(n_max, group.shape[0]),
                seed=seed,
            )
        )
        .drop("_dummy")
    )

    type_summary = sampled.group_by("label").len().sort("len", descending=True)
    print(f"    {sampled.shape[0]} cells after downsampling ({type_summary.shape[0]} types)")

    return sampled


def assign_train_val_split(df: pl.DataFrame, config: dict) -> pl.DataFrame:
    """Assign cells to train or validation split.

    Strategy depends on the number of samples in the dataset:
      - >4 samples: hold out ~20% of samples entirely. This gives the strongest
        test of generalization since cells from the same sample share technical
        artifacts (batch effects, library prep noise).
      - <=4 samples: fall back to stratified random split at the cell level
        (80% train / 20% val), preserving label proportions. This is weaker
        (overestimates performance) but necessary when too few samples exist
        to hold any out entirely.
    """
    seed = config["seed"]
    rng = np.random.default_rng(seed)
    n_samples = df["sample_id"].n_unique()

    if n_samples > 4:
        # Sample-level holdout: entire samples go to validation
        samples = sorted(df["sample_id"].unique().to_list())
        n_val = max(1, len(samples) // 5)
        val_samples = rng.choice(samples, size=n_val, replace=False).tolist()
        print(f"    Sample-level holdout: {n_val} val samples of {len(samples)}")
        print(f"    Val samples: {val_samples}")

        df = df.with_columns(
            pl.when(pl.col("sample_id").is_in(val_samples))
            .then(pl.lit("val"))
            .otherwise(pl.lit("train"))
            .alias("split")
        )
    else:
        # Cell-level stratified split: 20% of cells per label go to validation
        print(f"    Cell-level stratified split (<=4 samples)")

        labels = df["label"].to_list()
        assignments = ["train"] * len(labels)

        # For each cell type, randomly assign 20% of its cells to validation
        for label in df["label"].unique().to_list():
            label_indices = [i for i, lbl in enumerate(labels) if lbl == label]
            n_val = max(1, int(len(label_indices) * 0.2))
            val_indices = rng.choice(label_indices, size=n_val, replace=False).tolist()
            for idx in val_indices:
                assignments[idx] = "val"

        df = df.with_columns(pl.Series("split", assignments))

    return df


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def normalize_counts(X: sp.spmatrix, scale_factor: int = 10000) -> sp.csr_matrix:
    """Total-count normalize, scale by factor, and log1p-transform.

    This matches the normalization in scprocess's label_celltypes.R:
      1. Compute library size per cell (sum of ALL genes, not just a subset)
      2. Divide each cell's counts by its library size
      3. Multiply by scale_factor (default 10k)
      4. Apply log1p

    Importantly, log1p(0) = 0, so the sparse structure is preserved — zero
    entries stay zero and don't get stored.
    """
    X_csr = X.tocsr().astype(np.float64)
    lib_sizes = np.array(X_csr.sum(axis=1)).flatten()
    lib_sizes[lib_sizes == 0] = 1.0  # avoid division by zero for empty cells
    # Multiply each row by (1 / lib_size) via diagonal matrix
    inv_lib = sp.diags(1.0 / lib_sizes)
    X_norm = inv_lib @ X_csr
    X_norm *= scale_factor
    # log1p on the .data array only touches non-zero entries (preserves sparsity)
    X_norm.data = np.log1p(X_norm.data)
    return X_norm


def load_expression_matrix(
    cells_df: pl.DataFrame, paths: dict
) -> tuple[sp.csr_matrix, list[str], list[str]]:
    """Load and normalize expression data, one H5AD at a time.

    Memory efficiency strategy:
      - Only batches (H5AD files) that contain selected cells are loaded.
      - Within each batch, only the selected cells are kept.
      - After normalization, the AnnData object is freed before loading the next.
      - Peak memory is thus proportional to one batch, not the full dataset.

    The H5AD files contain raw counts in .X (sparse CSC). Gene names come from
    .var_names (gene symbols, set in make_clean_h5ad.py). All batches from the
    same scprocess run share identical gene sets.

    Returns:
      X: sparse CSR matrix (n_cells x n_genes), log-normalized
      gene_names: list of gene symbols (column labels)
      cell_ids_ordered: list of cell_id strings (row labels, matches X row order)
    """
    h5ad_dict = paths["h5ad_dict"]
    batches_needed = cells_df["batch"].unique().to_list()

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

        # Subset to only the cells we selected in the planning phase
        cell_mask = adata.obs["cell_id"].isin(batch_cell_ids)
        adata = adata[cell_mask].copy()

        # All H5AD files from one scprocess run must have identical gene sets
        current_genes = adata.var_names.tolist()
        if gene_names is None:
            gene_names = current_genes
        else:
            assert current_genes == gene_names, (
                f"Gene mismatch in batch '{batch}': expected {len(gene_names)} genes, "
                f"got {len(current_genes)}"
            )

        # Normalize raw counts: total-count → scale to 10k → log1p
        X_norm = normalize_counts(adata.X, scale_factor=10000)
        matrices.append(X_norm)
        cell_ids_ordered.extend(adata.obs["cell_id"].tolist())

        # Free memory before loading next batch
        del adata

    # Vertically stack all batch matrices into one (cells x genes) matrix
    X = sp.vstack(matrices, format="csr")
    print(f"  Final matrix: {X.shape[0]} cells x {X.shape[1]} genes")

    return X, gene_names, cell_ids_ordered


# ---------------------------------------------------------------------------
# XGBoost training
# ---------------------------------------------------------------------------


def run_xgboost_pass1(
    X_train: sp.csr_matrix,
    y_train: np.ndarray,
    X_val: sp.csr_matrix,
    y_val: np.ndarray,
    config: dict,
) -> xgb.Booster:
    """XGBoost pass 1: broad exploration with all genes.

    Purpose: discover which genes (out of ~20-30k) are informative for cell type
    classification, WITHOUT prior HVG selection.

    Key parameter choices:
      - colsample_bytree=0.1: each tree only sees 10% of genes. Over 300 trees,
        the ensemble explores the full transcriptome. Genes that appear in many
        trees and produce high gain are the informative ones.
      - subsample=0.632: bootstrap fraction — standard "bagging" proportion.
      - learning_rate=0.1: moderate step size since individual trees are weak
        (only 10% of features).

    The validation set is used for early stopping (no improvement in mlogloss
    for pass1_early_stopping rounds triggers stopping).
    """
    num_class = int(y_train.max()) + 1

    params = {
        "objective": "multi:softprob",
        "num_class": num_class,
        "eval_metric": "mlogloss",
        "subsample": config["pass1_subsample"],
        "colsample_bytree": config["pass1_colsample_bytree"],
        "learning_rate": config["pass1_learning_rate"],
        "nthread": config["n_cores"],
        "seed": config["seed"],
    }

    dtrain = xgb.DMatrix(X_train, label=y_train)
    dval = xgb.DMatrix(X_val, label=y_val)
    evals = [(dtrain, "train"), (dval, "val")]

    model = xgb.train(
        params=params,
        dtrain=dtrain,
        num_boost_round=config["pass1_nrounds"],
        evals=evals,
        early_stopping_rounds=config["pass1_early_stopping"],
        verbose_eval=50,
    )

    print(f"  Best iteration: {model.best_iteration}")
    return model


def select_features_by_gain(
    model: xgb.Booster, gene_names: list[str], config: dict
) -> tuple[np.ndarray, list[str], pl.DataFrame]:
    """Select genes contributing to the top fraction of cumulative gain.

    "Gain" is XGBoost's measure of how much a feature reduces the loss function
    across all trees and splits where it was used. Higher gain = more informative.

    Selection logic:
      1. Rank all genes by their total gain (descending).
      2. Compute cumulative gain as a fraction of total.
      3. Select genes until the cumulative fraction reaches gain_threshold (0.9).
      4. Clamp the count to [min_genes, max_genes].

    This is adaptive: datasets with many informative genes (complex) will select
    more; datasets dominated by a few markers (simple) will select fewer.

    Returns:
      sel_indices: numpy array of column indices into the original gene matrix
      sel_gene_names: corresponding gene symbol names
      gene_importance: full DataFrame of all genes with non-zero gain (for output)
    """
    # get_score returns {feature_name: gain} where names are "f0", "f1", etc.
    importance = model.get_score(importance_type="gain")

    # Map feature index strings back to numeric indices
    gene_gains = np.zeros(len(gene_names))
    for feat, gain in importance.items():
        idx = int(feat[1:])  # "f123" → 123
        gene_gains[idx] = gain

    # Sort genes by gain, highest first
    sorted_idx = np.argsort(-gene_gains)
    sorted_gains = gene_gains[sorted_idx]

    # Only consider genes that were actually used by the model
    nonzero_mask = sorted_gains > 0
    sorted_idx = sorted_idx[nonzero_mask]
    sorted_gains = sorted_gains[nonzero_mask]

    # Cumulative gain fraction (what percentage of total gain do the top N genes explain?)
    total_gain = sorted_gains.sum()
    cum_gain = np.cumsum(sorted_gains) / total_gain

    # Find how many genes are needed to reach the threshold
    n_sel = int(np.searchsorted(cum_gain, config["gain_threshold"])) + 1
    n_sel = max(n_sel, config["min_genes"])   # never select fewer than 100
    n_sel = min(n_sel, config["max_genes"])   # never select more than 3000
    n_sel = min(n_sel, len(sorted_idx))       # can't exceed available genes

    sel_indices = sorted_idx[:n_sel]
    sel_gene_names = [gene_names[i] for i in sel_indices]

    # Build importance table (saved as output for biological interpretation)
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
    X_train: sp.csr_matrix,
    y_train: np.ndarray,
    X_val: sp.csr_matrix,
    y_val: np.ndarray,
    config: dict,
) -> xgb.Booster:
    """XGBoost pass 2: final model trained on the selected gene subset.

    This is a fresh retrain (not continuing from pass 1) on the curated feature
    set identified by select_features_by_gain. Since the gene set is already
    informative, we use:
      - Higher colsample_bytree (0.5): more features per tree since they're all good
      - Lower learning_rate (0.05): smaller steps for finer convergence
      - More rounds (500): compensates for the lower learning rate

    This model is the one that gets saved and used for downstream prediction.
    """
    num_class = int(y_train.max()) + 1

    params = {
        "objective": "multi:softprob",
        "num_class": num_class,
        "eval_metric": "mlogloss",
        "subsample": config["pass1_subsample"],
        "colsample_bytree": config["pass2_colsample_bytree"],
        "learning_rate": config["pass2_learning_rate"],
        "nthread": config["n_cores"],
        "seed": config["seed"],
    }

    # Feature names are generic ("g0", "g1", ...) — the actual gene symbols are
    # saved separately in the selected_genes.txt output file.
    dtrain = xgb.DMatrix(X_train, label=y_train, feature_names=[f"g{i}" for i in range(X_train.shape[1])])
    dval = xgb.DMatrix(X_val, label=y_val, feature_names=[f"g{i}" for i in range(X_val.shape[1])])
    evals = [(dtrain, "train"), (dval, "val")]

    model = xgb.train(
        params=params,
        dtrain=dtrain,
        num_boost_round=config["pass2_nrounds"],
        evals=evals,
        early_stopping_rounds=config["pass2_early_stopping"],
        verbose_eval=50,
    )

    print(f"  Best iteration: {model.best_iteration}")
    return model


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------


def evaluate_model(
    model: xgb.Booster,
    X_val: sp.csr_matrix,
    y_val: np.ndarray,
    class_names: list[str],
) -> None:
    """Print confusion metrics on the held-out validation set.

    Reports:
      - Overall accuracy
      - High-confidence accuracy (cells where max probability > 0.5)
      - Per-class accuracy and F1 score
      - Macro-averaged F1
      - Top misclassification pairs (which true classes get confused with which)
    """
    dval = xgb.DMatrix(X_val, feature_names=[f"g{i}" for i in range(X_val.shape[1])])
    probs = model.predict(dval)  # shape: (n_cells, n_classes)
    y_pred = probs.argmax(axis=1)
    p_max = probs.max(axis=1)

    n_classes = len(class_names)
    overall_acc = (y_pred == y_val).mean()
    print(f"  Overall accuracy: {overall_acc:.3f}")

    # High-confidence accuracy (p > 0.5)
    hc_mask = p_max > 0.5
    if hc_mask.sum() > 0:
        hc_acc = (y_pred[hc_mask] == y_val[hc_mask]).mean()
        print(f"  High-confidence accuracy (p>0.5): {hc_acc:.3f} ({hc_mask.sum()} cells)")

    # Per-class metrics
    print(f"\n  {'Class':<30} {'N':>6} {'Acc':>6} {'F1':>6}")
    print("  " + "-" * 52)

    f1_scores = []
    for i, name in enumerate(class_names):
        mask_true = y_val == i
        n_true = mask_true.sum()
        if n_true == 0:
            continue

        # Accuracy for this class
        acc = (y_pred[mask_true] == i).mean()

        # F1: precision and recall
        pred_pos = y_pred == i
        tp = ((y_pred == i) & (y_val == i)).sum()
        precision = tp / pred_pos.sum() if pred_pos.sum() > 0 else 0
        recall = tp / n_true if n_true > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        f1_scores.append(f1)

        print(f"  {name:<30} {n_true:>6} {acc:>6.3f} {f1:>6.3f}")

    macro_f1 = np.mean(f1_scores) if f1_scores else 0
    print(f"\n  Macro F1: {macro_f1:.3f}")

    # Top misclassifications
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


# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------


def save_outputs(
    model: xgb.Booster,
    class_names: list[str],
    selected_genes: list[str],
    gene_importance: pl.DataFrame,
    config: dict,
    cells_df: pl.DataFrame,
    paths: dict,
) -> None:
    """Save all model artifacts to the output directory.

    Output files:
      - {ref_tag}_xgboost_model.json: XGBoost model in native JSON format.
        Portable — can be loaded in both Python (xgb.Booster().load_model)
        and R (xgb.load.model).
      - {ref_tag}_allowed_cls.csv: one column 'class' with the class names in
        the order matching the model's integer encoding (0 = first row, etc.).
      - {ref_tag}_selected_genes.txt: gene symbols, one per line, in the order
        matching the model's feature columns. At prediction time, the new data
        must be subset to these genes in this order.
      - {ref_tag}_gene_importance.csv: full gene ranking by gain from pass 1
        (useful for biological interpretation).
      - plots/: UMAP diagnostic plots.
    """
    out_dir = pathlib.Path(config["output_dir"])
    ref_tag = config["ref_tag"]

    # Save model in native JSON format (portable across Python/R)
    model_path = out_dir / f"{ref_tag}_xgboost_model.json"
    model.save_model(str(model_path))
    print(f"  Model: {model_path}")

    # Save class names — row order matches integer label encoding
    cls_path = out_dir / f"{ref_tag}_allowed_cls.csv"
    pl.DataFrame({"class": class_names}).write_csv(str(cls_path))
    print(f"  Classes: {cls_path}")

    # Save selected gene list — order matches model feature columns
    genes_path = out_dir / f"{ref_tag}_selected_genes.txt"
    with open(genes_path, "w") as f:
        f.write("\n".join(selected_genes) + "\n")
    print(f"  Genes: {genes_path}")

    # Save full gene importance table from pass 1
    imp_path = out_dir / f"{ref_tag}_gene_importance.csv"
    gene_importance.write_csv(str(imp_path))
    print(f"  Importance: {imp_path}")

    # Diagnostic plots
    print("  Generating plots...")
    make_diagnostic_plots(cells_df, model, class_names, paths, config)


def make_diagnostic_plots(
    cells_df: pl.DataFrame,
    model: xgb.Booster,
    class_names: list[str],
    paths: dict,
    config: dict,
) -> None:
    """Generate UMAP diagnostic plots using pre-computed coordinates from scprocess.

    Uses UMAP1/UMAP2 from the integration cluster CSV (no re-computation needed).
    Produces:
      1. UMAP colored by assigned training label — sanity check that labels
         correspond to distinct transcriptomic clusters.
      2. UMAP colored by train/val split — verify that spatial coverage is
         reasonable in both sets (especially for sample-level holdout).
    """
    out_dir = pathlib.Path(config["output_dir"]) / "plots"
    ref_tag = config["ref_tag"]

    # Load UMAP coordinates from the scprocess cluster CSV
    umap_df = pl.read_csv(paths["cluster_csv"], columns=["cell_id", "UMAP1", "UMAP2"])

    # Join with our training cells to get labels and split assignments
    plot_df = cells_df.join(umap_df, on="cell_id", how="left")
    plot_df = plot_df.filter(
        pl.col("UMAP1").is_not_null() & pl.col("UMAP2").is_not_null()
    )

    # Assign colors to labels
    labels = sorted(plot_df["label"].unique().to_list())
    cmap = plt.cm.get_cmap("tab20", len(labels))
    label_colors = {lbl: cmap(i) for i, lbl in enumerate(labels)}

    # Plot 1: UMAP colored by true (refined) label
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    for lbl in labels:
        subset = plot_df.filter(pl.col("label") == lbl)
        ax.scatter(
            subset["UMAP1"].to_numpy(),
            subset["UMAP2"].to_numpy(),
            c=[label_colors[lbl]],
            s=1,
            alpha=0.5,
            label=lbl,
            rasterized=True,
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

    # Plot 2: UMAP showing train/val split
    fig, ax = plt.subplots(1, 1, figsize=(8, 7))
    split_colors = {"train": "#1f77b4", "val": "#ff7f0e"}
    for split in ["train", "val"]:
        subset = plot_df.filter(pl.col("split") == split)
        ax.scatter(
            subset["UMAP1"].to_numpy(),
            subset["UMAP2"].to_numpy(),
            c=split_colors[split],
            s=1,
            alpha=0.4,
            label=split,
            rasterized=True,
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Train XGBoost cell type classifier from scprocess output"
    )
    parser.add_argument("config_yaml", type=str, help="Path to training config YAML")
    args = parser.parse_args()
    make_classifier(args.config_yaml)
