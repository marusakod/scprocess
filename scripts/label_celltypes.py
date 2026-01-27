import os
import argparse
import pathlib
import polars as pl
import scanpy as sc
import celltypist
import gzip


def download_celltypist_models(models_f):
  # download
  celltypist.models.download_models()
    
  # record their names
  models_dir  = celltypist.models.models_path
  models_ls   = [ f.replace(".pkl", "") for f in os.listdir(models_dir) if f.endswith(".pkl") ]

  # make dataframe, save
  models_df   = pl.DataFrame({ "model": models_ls })
  models_df.write_csv(models_f)

  return


def run_celltypist(sel_batch, batch_var, model_name, adata_f):
 
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
  pred_df     = pl.from_pandas( predictions.predicted_labels.reset_index(names = "cell_id") ).rename({"predicted_labels": "predicted_label"})
  all_probs   = pl.from_pandas( predictions.probability_matrix.reset_index(names = "cell_id") )
  probs_df    = all_probs.select([
    pl.col('cell_id'),
    pl.max_horizontal(
        pl.exclude('cell_id')
    ).alias('probability')
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
  int_df      = int_df.select( pl.col("cell_id"), pl.col(hi_res_cl).alias("hi_res_cl") )


  # get all prediction files
  preds_df    = pl.concat([ pl.read_csv(f) for f in pred_fs ], how = "vertical")
  data_df     = preds_df.join(int_df, on = "cell_id")

  # join to int_df
  counts_df   = data_df.group_by("hi_res_cl", "predicted_label").agg(
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
      .otherwise(pl.col("predicted_label")),
    prop_hi_res_cl      = pl.col("prop")
  )

  # join to 
  agg_df      = data_df.join(hi_res_lu, on = "hi_res_cl", how = "left").select(
    'model', batch_var, 'cell_id', 'hi_res_cl', 'predicted_label_agg', 
    'prop_hi_res_cl',
    predicted_label_naive = pl.col('predicted_label'),
    probability_naive = pl.col('probability'),
  )

  return agg_df


if __name__ == "__main__":
  # define arguments
  parser      = argparse.ArgumentParser()

  # define subparsers
  subparsers  = parser.add_subparsers(dest='subcommand', required=True)
  downloader_prsr = subparsers.add_parser('download_models')
  typist_prsr     = subparsers.add_parser('celltypist_one_batch')
  agg_prsr        = subparsers.add_parser('aggregate_predictions')

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
  
  elif args.subcommand == "download_models":
    download_celltypist_models(args.models_f)

