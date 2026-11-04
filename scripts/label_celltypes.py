import os
import argparse
import pathlib
import polars as pl
import gzip
import requests

# Prediction jobs process consecutive batches until adding the next batch would
# exceed this many cells; a larger single batch forms its own chunk. Each job
# loads the model once, so fewer, larger jobs avoid per-job start-up costs.
PREDICTION_CHUNK_MAX_CELLS = 150_000


def download_celltypist_models(models_f):
  # lazy import: only available in the celltypist conda env
  import celltypist
  # download
  celltypist.models.download_models()
  models_dir = pathlib.Path(celltypist.models.models_path)

  # define Allen Brain Immunology models
  allen_base_url = "https://allenimmunology.org/public/publication/download/84792154-cdfb-42d0-8e42-39e210e980b4/filesets/c5300f8b-f5ff-4010-9371-edc33d489143"
  allen_urls = {
    "AIFI_L1": f"{allen_base_url}/ref_pbmc_clean_celltypist_model_AIFI_L1_2024-04-18.pkl",
    "AIFI_L2": f"{allen_base_url}/ref_pbmc_clean_celltypist_model_AIFI_L2_2024-04-19.pkl",
    "AIFI_L3": f"{allen_base_url}/ref_pbmc_clean_celltypist_model_AIFI_L3_2024-04-19.pkl"
  }

  # download them
  for model, url in allen_urls.items():
    # Extract the filename from the URL (or use the model)
    file_path = models_dir / (model + ".pkl")
    
    print(f"Downloading {model} to {file_path}...")
    try:
      response = requests.get(url, stream = True)
      response.raise_for_status() # Check for HTTP errors (404, 500, etc.)

      with open(file_path, "wb") as f:
        for chunk in response.iter_content(chunk_size = 8192):
          f.write(chunk)
      print(f"Successfully downloaded {model}.")

    except requests.exceptions.RequestException as e:
      print(f"Failed to download {model}: {e}")

  # record their names
  models_ls = [ f.replace(".pkl", "") for f in os.listdir(models_dir) if f.endswith(".pkl") ]

  # save per-model allowed_cls.csv
  cls_dir = pathlib.Path(models_f).parent / "allowed_cls"
  os.makedirs(cls_dir, exist_ok = True)
  
  for model_name in models_ls:
    try:
      m = celltypist.models.Model.load(model = model_name + ".pkl")
      cls_path = cls_dir / f"{model_name}_allowed_cls.csv"
      classes_df = pl.DataFrame({"class": sorted(m.cell_types.tolist())})
      classes_df.write_csv(cls_path)
      print(f"  Saved classes: {cls_path}")
    except Exception as e:
      print(f"  Failed to extract classes for {model_name}: {e}")

  # make dataframe, save
  models_df = pl.DataFrame({ "model": models_ls }).sort("model")
  models_df.write_csv(models_f)


  return


def assign_prediction_chunks(h5ads_yaml_f, chunks_f):
  """Group batches into prediction chunks of up to PREDICTION_CHUNK_MAX_CELLS cells.

  Batches keep their order in the H5AD YAML; zero-byte placeholders for failed
  batches count as 0 cells. E.g. 100k, 40k, 30k and 200k cells give chunks
  [chunk001: 100k + 40k], [chunk002: 30k], [chunk003: 200k]. Writes one row per batch.
  """
  import h5py
  import yaml

  with open(h5ads_yaml_f) as f:
    h5ad_lookup = yaml.safe_load(f)

  # Count cells from each H5AD's obs index without reading the matrix.
  batches, n_cells = [], []
  for batch, entry in h5ad_lookup.items():
    h5ad_f = entry['path'] if isinstance(entry, dict) else entry
    if pathlib.Path(h5ad_f).stat().st_size == 0:
      batch_cells = 0
    else:
      with h5py.File(h5ad_f, 'r') as h5:
        batch_cells = h5['obs'][h5['obs'].attrs['_index']].shape[0]
    batches.append(str(batch))
    n_cells.append(batch_cells)

  # Start a new chunk when the next batch would take the current one over the
  # cap; empty placeholders join the current chunk.
  chunk_ids, chunk_idx, chunk_cells = [], 1, 0
  for batch_cells in n_cells:
    if chunk_cells > 0 and batch_cells > 0 and chunk_cells + batch_cells > PREDICTION_CHUNK_MAX_CELLS:
      chunk_idx  += 1
      chunk_cells = 0
    chunk_cells += batch_cells
    chunk_ids.append(f'chunk{chunk_idx:03d}')

  chunks_df = pl.DataFrame({'batch': batches, 'n_cells': n_cells, 'chunk': chunk_ids})
  chunks_df.write_csv(chunks_f)
  print(f'Assigned {len(batches)} batches ({sum(n_cells)} cells) to {chunk_idx} prediction chunks', flush = True)


def run_celltypist_batches(batches, batch_var, model_name, adata_fs):
  """CellTypist calls for several batches, loading the model once."""
  # lazy imports: only available in the celltypist conda env, not scprocess_local
  import celltypist

  model = celltypist.models.Model.load(model_name + ".pkl")
  pred_dfs = [_run_celltypist_one_batch(batch, batch_var, model_name, model, adata_f)
    for batch, adata_f in zip(batches, adata_fs)]
  return pl.concat(pred_dfs, how = 'vertical')


def _run_celltypist_one_batch(sel_batch, batch_var, model_name, model, adata_f):
  import scanpy as sc
  import celltypist

  # read anndata
  print(f'Predicting {sel_batch}', flush = True)
  adata = sc.read_h5ad(adata_f)
  adata.obs_names = adata.obs_names.to_list()
  adata.var_names = adata.var['symbol'].to_list()

  # normalize
  sc.pp.normalize_total(adata)
  sc.pp.log1p(adata)

  # make predictions
  predictions = celltypist.annotate(adata, model = model, majority_voting = False)

  # turn into nice output
  pred_df = pl.from_pandas( predictions.predicted_labels.reset_index(names = "cell_id") ).rename({"predicted_labels": "predicted_label_naive"})
  all_probs_df = pl.from_pandas( predictions.probability_matrix.reset_index(names = "cell_id") )
  probs_df = all_probs_df.select([
    pl.col('cell_id'),
    pl.max_horizontal(
        pl.exclude('cell_id')
    ).alias('probability_naive')
  ])

  # join together
  pred_df = pred_df.join( probs_df, on = "cell_id" )
  pred_df = pred_df.with_columns(
    pl.lit("celltypist").alias("labeller"),
    pl.lit(sel_batch).alias(batch_var),
    pl.lit(model_name).alias("model")
  )

  return pred_df


def run_xgboost_batches(batches, batch_var, model_name, adata_fs, model_f, cls_f, genes_f):
  """XGBoost calls for several batches, loading the model once."""
  # Model and expression libraries are needed only for inference, not downloads.
  import xgboost as xgb
  from train_xgboost import get_xgboost_model_info

  model_info = get_xgboost_model_info(model_f, model_name)
  model = xgb.Booster()
  model.load_model(model_f)

  # Class and gene order are part of the saved model, not alphabetical order.
  class_names = pl.read_csv(cls_f)['class'].to_list()
  with open(genes_f, 'r') as f:
    selected_genes = [line.strip() for line in f if line.strip()]

  pred_dfs = [_run_xgboost_one_batch(batch, batch_var, model_name, adata_f, model = model,
    model_info = model_info, class_names = class_names, selected_genes = selected_genes)
    for batch, adata_f in zip(batches, adata_fs)]
  return pl.concat(pred_dfs, how = 'vertical')


def _run_xgboost_one_batch(sel_batch, batch_var, model_name, adata_f, *, model, model_info,
  class_names, selected_genes):
  # Failed/excluded batches can be represented by a zero-byte H5AD placeholder;
  # they contribute no rows.
  pred_schema = {
    'cell_id': pl.String, 'predicted_label_naive': pl.String,
    'probability_naive': pl.Float64, 'labeller': pl.String,
    batch_var: pl.String, 'model': pl.String,
    'model_id': pl.String, 'model_sha256': pl.String,
  }
  if pathlib.Path(adata_f).stat().st_size == 0:
    print(f"Skipping empty H5AD placeholder: {adata_f}", flush = True)
    return pl.DataFrame(schema = pred_schema)

  import anndata
  import scipy.sparse as sp
  import xgboost as xgb
  from train_xgboost import normalize_counts, select_model_features

  print(f'Predicting {sel_batch}', flush = True)
  adata = anndata.read_h5ad(adata_f)
  cell_ids = adata.obs_names.tolist()
  gene_ids = adata.var_names.tolist()

  raw_counts = adata.X
  if not sp.issparse(raw_counts):
    raw_counts = sp.csr_matrix(raw_counts)

  # Normalise with all input genes, then put the model genes in its saved order;
  # missing model genes remain zero.
  normalized_counts = normalize_counts(raw_counts)
  model_counts      = select_model_features(normalized_counts, gene_ids, selected_genes)

  feat_names = [f"g{i}" for i in range(len(selected_genes))]
  dmat = xgb.DMatrix(model_counts, feature_names = feat_names)
  # Report the most probable label and its score, without confidence filtering.
  probs = model.predict(dmat)
  all_pred_labels = [class_names[idx] for idx in probs.argmax(axis = 1)]
  all_pred_probs = probs.max(axis = 1).tolist()

  pred_df = pl.DataFrame({
    'cell_id': cell_ids,
    'predicted_label_naive': all_pred_labels,
    'probability_naive': all_pred_probs,
  }).with_columns(
    pl.lit("scprocess").alias("labeller"),
    pl.lit(sel_batch).alias(batch_var),
    pl.lit(model_name).alias("model"),
    pl.lit(model_info['model_id']).alias('model_id'),
    pl.lit(model_info['model_sha256']).alias('model_sha256'),
  )

  return pred_df.cast(pred_schema)


def aggregate_predictions(pred_fs, int_f, hi_res_cl, min_cl_prop, batch_var, label_map_f = None):
  # load integration, check cluster column is there
  int_df = pl.read_csv(int_f)
  if not hi_res_cl in int_df:
    raise KeyError(f"specified high resolution cluster column {hi_res_cl} is not the integration file:\n  {int_f}")
  
  # restrict to just cells that are not doublets
  int_df = int_df.filter( pl.col(hi_res_cl).is_not_null() )
  int_cols = ["cell_id", hi_res_cl]
  if batch_var in int_df.columns:
    int_cols.append(batch_var)
  int_df = int_df.select(int_cols)
  int_df = int_df.rename({hi_res_cl: "hi_res_cl"})


  # get all prediction files
  preds_df = pl.concat([ pl.read_csv(f, schema_overrides = {
    'probability_naive': pl.Float64}) for f in pred_fs ], how = "diagonal_relaxed")

  # A shared model name is insufficient when installed weights have changed.
  if 'model_sha256' in preds_df.columns and len(preds_df):
    checksums = preds_df['model_sha256']
    if checksums.null_count() or checksums.n_unique() != 1:
      raise ValueError('Prediction files contain missing or mixed model checksums; rerun annotation')

  # Keep only cell-level prediction columns. This also allows an already
  # aggregated labels file to be re-aggregated against a different integration
  # (for example, the high-resolution clusters from a zoom analysis) without
  # retaining its previous hi_res_cl/predicted_label_agg columns.
  pred_cols = ['model', 'model_id', 'model_sha256', batch_var, 'cell_id', 'predicted_label_naive',
    'probability_naive']
  pred_cols = [col for col in pred_cols if col in preds_df.columns]
  preds_df = preds_df.select(pred_cols)
  data_df = preds_df.join(int_df, on = "cell_id")

  # join to int_df
  counts_df = data_df.group_by("hi_res_cl", "predicted_label_naive").agg(
    pl.len().alias("N")
  )
  counts_df = counts_df.with_columns(
    (pl.col("N") / pl.col("N").sum().over("hi_res_cl")).alias("prop")
  )
  counts_df = counts_df.sort("hi_res_cl", "prop", descending = [False, True])

  # take top prediction for each cluster
  hi_res_lu = counts_df.group_by("hi_res_cl").first().select(
    "hi_res_cl",
    predicted_label_agg = pl.when(pl.col("prop") < min_cl_prop)
      .then(pl.lit("ambiguous"))
      .otherwise(pl.col("predicted_label_naive")),
    prop_hi_res_cl = pl.col("prop")
  )

  # Return one row per retained cell, with its cluster-majority and naive calls.
  out_cols = ['model', 'cell_id', 'hi_res_cl', 'predicted_label_agg',
    'prop_hi_res_cl', 'predicted_label_naive', 'probability_naive']
  out_cols += [name for name in ['model_id', 'model_sha256'] if name in data_df.columns]
  if batch_var in data_df.columns:
    out_cols.insert(1, batch_var)
  agg_df = data_df.join(hi_res_lu, on = "hi_res_cl", how = "left").select(out_cols)

  if label_map_f and pathlib.Path(label_map_f).is_file():
    from train_xgboost import read_label_specs, label_mapping
    records, _ = read_label_specs(label_specs_f = label_map_f)
    label_map = label_mapping(records)
    if not label_map:
      return agg_df
    agg_df = agg_df.with_columns(
      pl.col('predicted_label_naive').replace(label_map).alias('coarse_predicted_label_naive'),
      pl.col('predicted_label_agg').replace(label_map).alias('coarse_predicted_label_agg'),
    )

  return agg_df


def save_cluster_names(labels_f, integration_f, mkr_sel_res, output_f):
  """Generate a cluster_name CSV from aggregated labels at the marker gene resolution."""
  cluster_col = f"RNA_snn_res.{mkr_sel_res}"

  labels_df = pl.read_csv(labels_f)
  int_df = pl.read_csv(integration_f, columns = ["cell_id", cluster_col])
  int_df = int_df.filter(pl.col(cluster_col).is_not_null())

  merged_df = labels_df.join(int_df, on = "cell_id")

  # majority vote per cluster
  counts_df = (
    merged_df
    .group_by(cluster_col, "predicted_label_naive")
    .agg(pl.len().alias("N"))
    .with_columns(
      (pl.col("N") / pl.col("N").sum().over(cluster_col)).alias("prop")
    )
    .sort(cluster_col, "prop", descending = [False, True])
  )
  top_df = (
    counts_df
    .group_by(cluster_col).first()
    .select(
      pl.col(cluster_col).alias("cluster"),
      pl.col("predicted_label_naive").alias("cluster_name"),
    )
    .sort("cluster")
  )

  # make cluster_name unique by appending numbers for duplicates
  names = top_df["cluster_name"].to_list()
  name_counts = {}
  for name in names:
    name_counts[name] = name_counts.get(name, 0) + 1
  seen = {}
  unique_names = []
  for name in names:
    if name_counts[name] > 1:
      seen[name] = seen.get(name, 0) + 1
      unique_names.append(f"{name} {seen[name]}")
    else:
      unique_names.append(name)

  cluster_names_df = top_df.with_columns(pl.Series("cluster_name", unique_names))
  cluster_names_df.write_csv(output_f)


def extract_naive_predictions(source_labels_fs, int_f, model_name, pred_f, model_sha256 = None):
  """Extract naive predictions from source project labels for cells in the join integration."""
  int_df = pl.read_csv(int_f, columns = ["cell_id"])
  valid_cells = set(int_df["cell_id"].to_list())

  prediction_dfs = []
  for path in source_labels_fs:
    source_df = pl.read_csv(path)
    source_df = source_df.filter(pl.col("cell_id").is_in(valid_cells))
    if model_sha256:
      if ('model_sha256' not in source_df.columns
          or source_df['model_sha256'].null_count()
          or not (source_df['model_sha256'] == model_sha256).all()):
        raise ValueError('Source predictions use different model weights; rerun annotation')
    identity_columns = [name for name in ('model_id', 'model_sha256') if name in source_df.columns]
    source_df = source_df.select(
      pl.lit(model_name).alias('model'), 'cell_id', 'predicted_label_naive',
      'probability_naive', *identity_columns)
    prediction_dfs.append(source_df)

  predictions_df = pl.concat(prediction_dfs)

  # Write the header once and serialize bounded chunks into the gzip stream.
  with gzip.open(pred_f, 'wt') as stream:
    stream.write(predictions_df.head(0).write_csv())
    for batch_df in predictions_df.iter_slices(n_rows = 100_000):
      stream.write(batch_df.write_csv(include_header = False))


def main():
  # Define the available annotation and aggregation commands.
  parser = argparse.ArgumentParser()

  # define subparsers
  subparsers      = parser.add_subparsers(dest = 'subcommand', required = True)
  downloader_prsr = subparsers.add_parser('download_models')
  chunks_prsr     = subparsers.add_parser('assign_prediction_chunks')
  typist_prsr     = subparsers.add_parser('celltypist_batches')
  xgb_prsr        = subparsers.add_parser('xgboost_batches')
  agg_prsr        = subparsers.add_parser('aggregate_predictions')
  names_prsr      = subparsers.add_parser('save_cluster_names')
  extract_prsr    = subparsers.add_parser('extract_naive_predictions')

  # get arguments
  downloader_prsr.add_argument("models_f",type = str)  
  
  # prediction chunks: batches grouped by cell count
  chunks_prsr.add_argument("--h5ads_yaml_f", type = str, required = True)
  chunks_prsr.add_argument("--chunks_f",     type = str, required = True)

  # celltypist: one model, several batches with matching H5ADs
  typist_prsr.add_argument(  "batch_var", type = str)
  typist_prsr.add_argument(  "model",     type = str)
  typist_prsr.add_argument("--pred_f",    type = str, required = True)
  typist_prsr.add_argument("--batches",   type = str, nargs = "+", required = True)
  typist_prsr.add_argument("--adata_fs",  type = str, nargs = "+", required = True)

  # xgboost: one model, several batches with matching H5ADs
  xgb_prsr.add_argument(  "batch_var", type = str)
  xgb_prsr.add_argument(  "model",     type = str)
  xgb_prsr.add_argument("--model_f",   type = str)
  xgb_prsr.add_argument("--cls_f",     type = str)
  xgb_prsr.add_argument("--genes_f",   type = str)
  xgb_prsr.add_argument("--pred_f",    type = str, required = True)
  xgb_prsr.add_argument("--batches",   type = str, nargs = "+", required = True)
  xgb_prsr.add_argument("--adata_fs",  type = str, nargs = "+", required = True)

  # get arguments
  agg_prsr.add_argument(  "pred_fs",      type = str, nargs = "+")
  agg_prsr.add_argument("--int_f",        type = str)
  agg_prsr.add_argument("--hi_res_cl",    type = str)
  agg_prsr.add_argument("--min_cl_prop",  type = float)
  agg_prsr.add_argument("--batch_var",    type = str)
  agg_prsr.add_argument("--label_map_f",  type = str, default = None)
  agg_prsr.add_argument("--agg_f",        type = str)

  # save_cluster_names arguments
  names_prsr.add_argument("--labels_f",      type = str, required = True)
  names_prsr.add_argument("--integration_f", type = str, required = True)
  names_prsr.add_argument("--mkr_sel_res",   type = str, required = True)
  names_prsr.add_argument("--output_f",      type = str, required = True)

  # extract_naive_predictions arguments
  extract_prsr.add_argument("--model_sha256")
  extract_prsr.add_argument(  "source_labels_fs", type = str, nargs = "+")
  extract_prsr.add_argument("--int_f",        type = str, required = True)
  extract_prsr.add_argument("--model",        type = str, required = True)
  extract_prsr.add_argument("--pred_f",       type = str, required = True)

  # Parse the arguments
  args = parser.parse_args()

  # Dispatch the requested command.
  if args.subcommand == "assign_prediction_chunks":
    assign_prediction_chunks(args.h5ads_yaml_f, args.chunks_f)

  elif args.subcommand in ("celltypist_batches", "xgboost_batches"):
    # Batches and H5ADs are paired by position.
    if len(args.batches) != len(args.adata_fs):
      raise ValueError(f"Got {len(args.batches)} batches but {len(args.adata_fs)} H5AD files")
    for adata_f in args.adata_fs:
      if not pathlib.Path(adata_f).is_file():
        raise FileNotFoundError(f"adata_f is not a valid file:\n  {adata_f}")

    if args.subcommand == "celltypist_batches":
      pred_df = run_celltypist_batches(args.batches, args.batch_var, args.model, args.adata_fs)
    else:
      pred_df = run_xgboost_batches(args.batches, args.batch_var, args.model, args.adata_fs,
        args.model_f, args.cls_f, args.genes_f)

    with gzip.open(args.pred_f, 'wb') as f:
      pred_df.write_csv(f)

  elif args.subcommand == "aggregate_predictions":
    # set up some locations
    pred_fs = [ pathlib.Path(f) for f in args.pred_fs ]
    for f in pred_fs:
      if not f.is_file():
        raise FileNotFoundError(f"this file does not exist:\n{str(f)}")
    int_f = pathlib.Path(args.int_f)
    if not int_f.is_file():
      raise FileNotFoundError(f"this file does not exist:\n{str(int_f)}")

    # do some checks
    if (args.min_cl_prop < 0) or (args.min_cl_prop >= 1):
      raise ValueError("min_cl_prop must be greater than or equal to 0 and strictly less than 1")

    # run
    agg_df = aggregate_predictions(pred_fs, int_f, args.hi_res_cl, args.min_cl_prop, args.batch_var, label_map_f = args.label_map_f)

    # save
    with gzip.open(args.agg_f, 'wb') as f:
      agg_df.write_csv(f)
  
  elif args.subcommand == "save_cluster_names":
    save_cluster_names(args.labels_f, args.integration_f, args.mkr_sel_res, args.output_f)

  elif args.subcommand == "extract_naive_predictions":
    extract_naive_predictions(
      args.source_labels_fs, args.int_f, args.model, args.pred_f, model_sha256 = args.model_sha256)

  elif args.subcommand == "download_models":
    download_celltypist_models(args.models_f)


if __name__ == '__main__':
  main()
