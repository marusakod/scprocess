import os
import argparse
import pathlib
import polars as pl
import gzip
import requests


def download_celltypist_models(models_f):
  import celltypist
  # download
  celltypist.models.download_models()
  models_dir  = pathlib.Path(celltypist.models.models_path)

  # define Allen Brain Immunology models
  allen_base_url  = "https://allenimmunology.org/public/publication/download/84792154-cdfb-42d0-8e42-39e210e980b4/filesets/c5300f8b-f5ff-4010-9371-edc33d489143"
  allen_urls  = {
    "AIFI_L1": f"{allen_base_url}/ref_pbmc_clean_celltypist_model_AIFI_L1_2024-04-18.pkl",
    "AIFI_L2": f"{allen_base_url}/ref_pbmc_clean_celltypist_model_AIFI_L2_2024-04-19.pkl",
    "AIFI_L3": f"{allen_base_url}/ref_pbmc_clean_celltypist_model_AIFI_L3_2024-04-19.pkl"
  }

  # download them
  for model, url in allen_urls.items():
    # Extract the filename from the URL (or use the model)
    file_path   = models_dir / (model + ".pkl")
    
    print(f"Downloading {model} to {file_path}...")
    try:
      response = requests.get(url, stream=True)
      response.raise_for_status() # Check for HTTP errors (404, 500, etc.)

      with open(file_path, "wb") as f:
        for chunk in response.iter_content(chunk_size=8192):
          f.write(chunk)
      print(f"Successfully downloaded {model}.")

    except requests.exceptions.RequestException as e:
      print(f"Failed to download {model}: {e}")

  # record their names
  models_ls   = [ f.replace(".pkl", "") for f in os.listdir(models_dir) if f.endswith(".pkl") ]

  # make dataframe, save
  models_df   = pl.DataFrame({ "model": models_ls }).sort("model")
  models_df.write_csv(models_f)

  return


def run_celltypist(sel_batch, batch_var, model_name, adata_f):
  import scanpy as sc
  import celltypist

  # read anndata
  adata = sc.read_h5ad(adata_f)
  adata.obs_names = adata.obs_names.to_list()
  adata.var_names = adata.var['symbol'].to_list()
  
  # normalize
  sc.pp.normalize_total(adata)
  sc.pp.log1p(adata)

  # make predictions
  predictions = celltypist.annotate(adata, model = model_name + ".pkl", majority_voting = False)

  # turn into nice output
  pred_df     = pl.from_pandas( predictions.predicted_labels.reset_index(names = "cell_id") ).rename({"predicted_labels": "predicted_label_naive"})
  all_probs   = pl.from_pandas( predictions.probability_matrix.reset_index(names = "cell_id") )
  probs_df    = all_probs.select([
    pl.col('cell_id'),
    pl.max_horizontal(
        pl.exclude('cell_id')
    ).alias('probability_naive')
  ])
  
  # join together
  pred_df     = pred_df.join( probs_df, on = "cell_id" )
  pred_df     = pred_df.with_columns(
    pl.lit("celltypist").alias("labeller"),
    pl.lit(sel_batch).alias(batch_var),
    pl.lit(model_name).alias("model")
  )

  return pred_df


def aggregate_predictions(pred_fs, int_f, hi_res_cl, min_cl_prop, batch_var):
  # load integration, check cluster column is there
  int_df      = pl.read_csv(int_f)
  if not hi_res_cl in int_df:
    raise KeyError(f"specified high resolution cluster column {hi_res_cl} is not the integration file:\n  {int_f}")
  
  # restrict to just cells that are not doublets
  int_df      = int_df.filter( pl.col(hi_res_cl).is_not_null() )
  int_cols    = ["cell_id", hi_res_cl]
  if batch_var in int_df.columns:
    int_cols.append(batch_var)
  int_df      = int_df.select(int_cols)
  int_df      = int_df.rename({hi_res_cl: "hi_res_cl"})


  # get all prediction files
  preds_df    = pl.concat([ pl.read_csv(f) for f in pred_fs ], how = "vertical")
  data_df     = preds_df.join(int_df, on = "cell_id")

  # join to int_df
  counts_df   = data_df.group_by("hi_res_cl", "predicted_label_naive").agg(
    pl.len().alias("N")
  )
  counts_df   = counts_df.with_columns(
    (pl.col("N") / pl.col("N").sum().over("hi_res_cl")).alias("prop")
  )
  counts_df   = counts_df.sort("hi_res_cl", "prop", descending = [False, True])

  # take top prediction for each cluster
  hi_res_lu   = counts_df.group_by("hi_res_cl").first().select(
    "hi_res_cl",
    predicted_label_agg = pl.when(pl.col("prop") < min_cl_prop)
      .then(pl.lit("ambiguous"))
      .otherwise(pl.col("predicted_label_naive")),
    prop_hi_res_cl      = pl.col("prop")
  )

  # join to
  out_cols = ['model', 'cell_id', 'hi_res_cl', 'predicted_label_agg',
    'prop_hi_res_cl', 'predicted_label_naive', 'probability_naive']
  if batch_var in data_df.columns:
    out_cols.insert(1, batch_var)
  agg_df      = data_df.join(hi_res_lu, on = "hi_res_cl", how = "left").select(out_cols)

  return agg_df


def save_cluster_names(labels_f, integration_f, mkr_sel_res, output_f):
  """Generate a cluster_name CSV from aggregated labels at the marker gene resolution."""
  cluster_col = f"RNA_snn_res.{mkr_sel_res}"

  labels_df = pl.read_csv(labels_f)
  int_df = pl.read_csv(integration_f, columns=["cell_id", cluster_col])
  int_df = int_df.filter(pl.col(cluster_col).is_not_null())

  merged_df = labels_df.join(int_df, on="cell_id")

  # majority vote per cluster
  counts_df = (
    merged_df
    .group_by(cluster_col, "predicted_label_naive")
    .agg(pl.len().alias("N"))
    .with_columns(
      (pl.col("N") / pl.col("N").sum().over(cluster_col)).alias("prop")
    )
    .sort(cluster_col, "prop", descending=[False, True])
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

  annotation_df = top_df.with_columns(pl.Series("cluster_name", unique_names))
  annotation_df.write_csv(output_f)


def extract_naive_predictions(source_labels_fs, int_f, model_name, pred_f):
  """Extract naive predictions from source project labels for cells in the join integration."""
  int_df = pl.read_csv(int_f, columns=["cell_id"])
  valid_cells = set(int_df["cell_id"].to_list())

  dfs = []
  for f in source_labels_fs:
    df = pl.read_csv(f)
    df = df.filter(pl.col("cell_id").is_in(valid_cells)).select(
      pl.lit(model_name).alias("model"),
      "cell_id",
      "predicted_label_naive",
      "probability_naive",
    )
    dfs.append(df)

  result = pl.concat(dfs)
  with gzip.open(pred_f, 'wb') as fh:
    result.write_csv(fh)


if __name__ == "__main__":
  # define arguments
  parser      = argparse.ArgumentParser()

  # define subparsers
  subparsers  = parser.add_subparsers(dest='subcommand', required=True)
  downloader_prsr = subparsers.add_parser('download_models')
  typist_prsr     = subparsers.add_parser('celltypist_one_batch')
  agg_prsr        = subparsers.add_parser('aggregate_predictions')
  names_prsr      = subparsers.add_parser('save_cluster_names')
  extract_prsr    = subparsers.add_parser('extract_naive_predictions')

  # get arguments
  downloader_prsr.add_argument("models_f",type = str)  
  
  # get arguments
  typist_prsr.add_argument(  "batch",     type=str)
  typist_prsr.add_argument(  "batch_var", type=str)
  typist_prsr.add_argument(  "model",     type=str)
  typist_prsr.add_argument("--adata_f",  type=str)
  typist_prsr.add_argument("--pred_f",    type=str)

  # get arguments
  agg_prsr.add_argument(  "pred_fs",      type=str, nargs="+")
  agg_prsr.add_argument("--int_f",        type=str)
  agg_prsr.add_argument("--hi_res_cl",    type=str)
  agg_prsr.add_argument("--min_cl_prop",  type=float)
  agg_prsr.add_argument("--batch_var",    type=str)
  agg_prsr.add_argument("--agg_f",        type=str)

  # save_cluster_names arguments
  names_prsr.add_argument("--labels_f",      type=str, required=True)
  names_prsr.add_argument("--integration_f", type=str, required=True)
  names_prsr.add_argument("--mkr_sel_res",   type=str, required=True)
  names_prsr.add_argument("--output_f",      type=str, required=True)

  # extract_naive_predictions arguments
  extract_prsr.add_argument(  "source_labels_fs", type=str, nargs="+")
  extract_prsr.add_argument("--int_f",        type=str, required=True)
  extract_prsr.add_argument("--model",        type=str, required=True)
  extract_prsr.add_argument("--pred_f",       type=str, required=True)

  # Parse the arguments
  args      = parser.parse_args()

  # create new project folder
  if args.subcommand == "celltypist_one_batch":
    # set up some locations
    adata_f   = pathlib.Path(args.adata_f)

    # check that they exist
    if not adata_f.is_file():
      raise FileNotFoundError(f"adata_f is not a valid file:\n  {adata_f}")
    
    # run
    pred_df   = run_celltypist(args.batch, args.batch_var, args.model, args.adata_f)

    # save
    with gzip.open(args.pred_f, 'wb') as f: 
      pred_df.write_csv(f)

  elif args.subcommand == "aggregate_predictions":
    # set up some locations
    pred_fs   = [ pathlib.Path(f) for f in args.pred_fs ]
    for f in pred_fs:
      if not f.is_file():
        raise FileNotFoundError(f"this file does not exist:\n{str(f)}")
    int_f     = pathlib.Path(args.int_f)
    if not int_f.is_file():
      raise FileNotFoundError(f"this file does not exist:\n{str(int_f)}")

    # do some checks
    if (args.min_cl_prop < 0) or (args.min_cl_prop >= 1):
      raise ValueError("min_cl_prop must be greater than or equal to 0 and strictly less than 1")

    # run
    agg_df    = aggregate_predictions(pred_fs, int_f, args.hi_res_cl, args.min_cl_prop, args.batch_var)

    # save
    with gzip.open(args.agg_f, 'wb') as f:
      agg_df.write_csv(f)
  
  elif args.subcommand == "save_cluster_names":
    save_cluster_names(args.labels_f, args.integration_f, args.mkr_sel_res, args.output_f)

  elif args.subcommand == "extract_naive_predictions":
    extract_naive_predictions(
      args.source_labels_fs, args.int_f, args.model, args.pred_f)

  elif args.subcommand == "download_models":
    download_celltypist_models(args.models_f)

