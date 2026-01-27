import polars as pl
import anndata as ad
import scanpy as sc
from scipy.sparse import csc_matrix
import os
import h5py
import argparse
import re
import numpy as np


def make_clean_h5ad(sel_b, sel_run, integration_f, h5_paths_f, 
  coldata_f, rowdata_f, run_var, batch_var, clean_h5ad_f):
  
  # load, exclude doublets
  int_dt = pl.read_csv(integration_f)
  int_dt = int_dt.filter((pl.col("is_dbl") == False) & (pl.col("in_dbl_cl") == False))
    
  # check if ok
  ok_batches = int_dt[batch_var].unique().to_list()
  if sel_b not in ok_batches:
    print(f'excluded {sel_b}; creating empty file.')
    open(clean_h5ad_f, 'a').close()
    return None
    
  print(f'creating clean h5ad file for {sel_b}')

  # get some subsets
  print(f'  getting values specific to {sel_b}')
  batch_int  = int_dt.filter(pl.col(batch_var) == sel_b)
  batch_ids  = batch_int["cell_id"].to_list()

  # get adata objects
  print('  loading counts into AnnData')
  h5_paths   = pl.read_csv(h5_paths_f)
  filtered_f = h5_paths.filter(pl.col(run_var) == sel_run).select("amb_filt_f").unique().item()
  assert os.path.exists(filtered_f), f"File {filtered_f} not found"
  
  adata = _get_clean_adata(filtered_f, sel_run, rowdata_f, coldata_f, batch_ids)
  
  # add integration outputs
  print('  adding integration variables')
  adata = _add_int_variables(adata, batch_int)
  
  # save
  print('  saving h5ad')
  
  # check for strange types that can't be saved
  obj_cols = adata.obs.select_dtypes(include=['object']).columns
  for col in obj_cols:
    adata.obs[col] = adata.obs[col].fillna("").astype(str)

  cat_cols = adata.obs.select_dtypes(include=['category']).columns
  for col in cat_cols:
    if adata.obs[col].isnull().any():
      if "" not in adata.obs[col].cat.categories:
        adata.obs[col] = adata.obs[col].cat.add_categories("")
      adata.obs[col] = adata.obs[col].fillna("")

  adata.write_h5ad(clean_h5ad_f)
  print('done!')



def _get_clean_adata(filtered_f, sel_run, rowdata_f, coldata_f, sub_cols ):
  
  with h5py.File(filtered_f, 'r') as f:
    indptr      = f['matrix/indptr'][:]
    indices     = f['matrix/indices'][:]
    data        = f['matrix/data'][:]
    features    = f['matrix/features/name'][:]
    barcodes    = f['matrix/barcodes'][:].astype('U16')
    num_rows    = f['matrix/shape'][0]
    num_cols    = f['matrix/shape'][1]

  # make a csc sparse matrix
  sua_csc_mat = csc_matrix((data, indices, indptr), shape=(num_rows, num_cols))
  csc_mat, uniq_features = sum_SUA(sua_csc_mat, features)
  
  barcodes    = barcodes.astype('<U21')
  barcodes    = [f"{sel_run}:{bc}" for bc in barcodes]  
  barcodes    = np.array(barcodes)

  # create anndata object
  adata = ad.AnnData(X=csc_mat.T)
  adata.obs_names = barcodes
  adata.var_names = uniq_features
  
  # filter features and add feature metadata
  feats_dt = pl.read_csv(rowdata_f)
  keep_ids = feats_dt["ensembl_id"].to_list()
  adata    = adata[:, adata.var_names.isin(keep_ids)].copy()

  var_df = (pl.DataFrame({"ensembl_id": adata.var_names})
    .join(feats_dt, on="ensembl_id", how="left")
  )
  adata.var = var_df.to_pandas().set_index("gene_id")
  
  # filter cells and add cell metadata
  batch_cols = pl.read_csv(coldata_f).filter(pl.col("cell_id").is_in(sub_cols))
  adata      = adata[adata.obs_names.isin(sub_cols), :].copy()

  obs_df = (pl.DataFrame({"cell_id": adata.obs_names})
    .join(batch_cols, on = "cell_id", how = 'left')
  )
  adata.obs = obs_df.to_pandas().set_index('cell_id')
  
  return adata


def _add_int_variables(adata, batch_int):
  # make sure things are aligned
  batch_int_pd = batch_int.to_pandas().set_index("cell_id")
  batch_int_pd = batch_int_pd.reindex(adata.obs_names)

  # find umap and clustering cols
  int_vs = ["UMAP1", "UMAP2"] + [col for col in batch_int_pd.columns if "RNA_snn_res" in col]
    
  for v in int_vs:
    if "RNA_snn_res" in v:
      # convert clusters to categorical
      adata.obs[v] = batch_int_pd[v].astype('category')
    else:
      adata.obs[v] = batch_int_pd[v]

  return adata


# same function as in hvgs.py
def sum_SUA(sua_mat, row_names):
  # define some objects
  row_names = row_names.astype(str)
  types = ['_S$', '_U$', '_A$']
  mats = []

  # sum up each type of gene mapping
  for t in types:
    indices = [i for i, name in enumerate(row_names) if re.search(t,name)]
    mat = sua_mat[indices, :]
    mats.append(mat)

  # give them names
  genes = [re.sub(r'_[SUA]$', '', name) for name in row_names]
  uniq_genes = list(dict.fromkeys(genes))
  uniq_genes = np.array(uniq_genes)

  # sum across the mappings
  mats_sum = sum(mats)

  return mats_sum, uniq_genes


if __name__ == "__main__":
  # get arguments
  parser  = argparse.ArgumentParser()
  parser.add_argument("sel_batch", type=str)
  parser.add_argument("sel_run", type=str)
  parser.add_argument("integration_f", type=str)
  parser.add_argument("h5_paths_f", type = str)
  parser.add_argument("coldata_f", type=str)
  parser.add_argument("rowdata_f", type=str)
  parser.add_argument("run_var", type = str)
  parser.add_argument("batch_var", type=str)
  parser.add_argument("clean_h5ad_f", type = str)

  # set up some locations
  args    = parser.parse_args()

  # run
  make_clean_h5ad(args.sel_batch, args.sel_run, args.integration_f, args.h5_paths_f, 
    args.coldata_f, args.rowdata_f, args.run_var, args.batch_var, args.clean_h5ad_f)