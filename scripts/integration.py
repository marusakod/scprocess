import h5py
import numpy as np
import argparse
import re
import gzip
import scipy.sparse as sp
from scipy.sparse import csc_matrix, csr_matrix, hstack
import anndata as ad
import pandas as pd
import polars as pl

# hvg_mat_f        = '/pstore/data/mus-brain-analysis/studies/test_project/output/test_hvg/top_hvgs_counts_test_project_2025-01-01.h5'
# dbl_hvg_mat_f    = '/pstore/data/mus-brain-analysis/studies/test_project/output/test_hvg/top_hvgs_doublet_counts_test_project_2025-01-01.h5'
# sample_qc_f      = '/pstore/data/mus-brain-analysis/studies/test_project/output/test_qc/qc_sample_statistics_test_project_2025-01-01.csv'
# coldata_f        = '/pstore/data/mus-brain-analysis/studies/test_project/output/test_qc/coldata_dt_all_samples_test_project_2025-01-01.txt.gz'
# demux_type       = 'none'
# exclude_mito     = 'True'
# embedding        = 'harmony'
# n_dims           = 50
# cl_method        = 'leiden'
# dbl_res          = 4
# dbl_cl_prop      = 0.5
# theta            = 0.1
# res_ls_concat    = '0.1 0.2 0.5 1 2'
# integration_f    = '/pstore/data/mus-brain-analysis/studies/test_project/output/test_integration/python_integrated_dt_test_project_2025-01-01.txt.gz'
# batch_var        = 'sample_id' 
# n_cores          = 8

def run_integration(hvg_mat_f, dbl_hvg_mat_f, sample_qc_f, coldata_f, demux_type, 
  exclude_mito, embedding, n_dims, cl_method, dbl_res, dbl_cl_prop, theta, res_ls_concat,
  integration_f, batch_var, use_gpu = False): 

  print('setting up parameters')
  exclude_mito  = bool(exclude_mito)
  res_ls        = res_ls_concat.split()
  if demux_type == "none":
    dbl_batch_var = 'sample_id'
  else:
    dbl_batch_var = 'pool_id'
  dbl_theta     = 0

  print('loading hvg matrix')
  all_hvg_mat, all_bcs   = _get_hvg_mat(hvg_mat_f, dbl_hvg_mat_f)

  print('loading relevant cell ids')
  ok_cells_df   = _get_ok_cells_df(sample_qc_f, coldata_f, all_bcs)

  print('normalizing hvg matrix')
  all_hvg_mat   = _normalize_hvg_mat(all_hvg_mat, ok_cells_df, exclude_mito)

  print('making anndata object')
  adata         = ad.AnnData(X = all_hvg_mat.T , obs = ok_cells_df.to_pandas())
  print(f"  anndata object has {adata.shape[0]} cells and {adata.shape[1]} dims")

  print('running integration to find more doublets')
  int_dbl       = _do_one_integration(adata, dbl_batch_var, cl_method, n_dims,
    dbl_res, embedding, use_gpu, theta = dbl_theta)

  print('filter to non-doublet cells')
  dbl_data      = _calc_dbl_data(int_dbl, ok_cells_df, dbl_res, dbl_cl_prop)
  adata         = _adata_filter_out_doublets(adata, dbl_data)

  print('running integration on clean data')
  int_ok        = _do_one_integration(adata, batch_var, cl_method, n_dims, res_ls, embedding, use_gpu, theta)

  print('join results')
  int_dt        = int_ok.join(dbl_data, on=["sample_id", "cell_id"], coalesce=True, how = 'full')

  print('save results')
  with gzip.open(integration_f, 'wb') as f:
    int_dt.write_csv(f)

  print('done!')

  return int_dt


def _get_ok_cells_df(sample_qc_f, coldata_f, all_bcs):
  # load files
  sample_qc     = pl.read_csv(sample_qc_f)
  all_coldata   = pl.read_csv(coldata_f)

  # checks
  if not 'sample_id' in all_coldata.columns:
    raise KeyError("column sample_id is missing from coldata file")
  if not 'sample_id' in sample_qc.columns:
    raise KeyError("column sample_id is missing from sample QC file")

  # get ok samples
  ok_samples    = sample_qc.filter(pl.col("bad_sample") == False)["sample_id"].to_list()

  # get ok cells
  ok_cells_df   = all_coldata.filter(
    ((pl.col("keep") == True) | (pl.col("dbl_class") == "doublet")) & 
    ((pl.col("sample_id").is_in(ok_samples)) | (pl.col("sample_id") == ""))
  )
  # check ok
  if not set(all_bcs) == set(ok_cells_df['cell_id'].to_list()):
    raise ValueError("barcodes from hvg mats and cell_ids don't match")

  # put cell in coldata in the order matching mat cols
  order_dt      = pl.DataFrame({
    "cell_id":    all_bcs,
    "order":      range(1, len(all_bcs) + 1)
  })
  ok_cells_df   = ok_cells_df.join(order_dt, on = 'cell_id').sort('order').drop('order')

  return ok_cells_df


def _get_hvg_mat(hvg_mat_f, dbl_hvg_mat_f):
  # get a matrix with hvgs (cells and doublets)
  all_hvg_mat = None
  all_bcs     = []
  features    = []

  # open both
  for mat_f in [hvg_mat_f, dbl_hvg_mat_f]:
    # load
    with h5py.File(mat_f, 'r') as f:
      # get components
      indptr      = f['matrix/indptr'][:]
      indices     = f['matrix/indices'][:]
      data        = f['matrix/data'][:]
      tmp_feats   = f['matrix/features/name'][:]
      barcodes    = f['matrix/barcodes'][:]
      num_rows    = f['matrix/shape'][0]
      num_cols    = f['matrix/shape'][1]
      
      # make sparse matrix
      csc_mat     = csc_matrix((data, indices, indptr), shape=(num_rows, num_cols))
      # all_hvg_mat = csc_mat if all_hvg_mat is None else hstack([all_hvg_mat, csc_mat])
      all_hvg_mat = hstack( [all_hvg_mat, csc_mat] )
      
      # sort out features and barcodes
      if len(features) == 0:
        features = tmp_feats

      # store barcodes
      barcodes    = [b.decode('utf-8') for b in barcodes]
      all_bcs.extend(barcodes)

  return all_hvg_mat, all_bcs


def _normalize_hvg_mat(hvg_mat, coldata, exclude_mito, scale_f = 10000):
  if exclude_mito:
    coldata = coldata.with_columns((pl.col("total") - pl.col("subsets_mito_sum")).alias("total_no_mito"))
    lib_sizes = coldata["total_no_mito"].to_numpy()
  else:
    lib_sizes = coldata["total"].to_numpy()

  hvg_mat = hvg_mat.tocsr() 
  hvg_mat.data /= lib_sizes[hvg_mat.indices]
  hvg_mat.data *= scale_f

  np.log1p(hvg_mat.data, out=hvg_mat.data)

  return hvg_mat


def _do_one_integration(adata, batch_var, cl_method, n_dims, res_ls, embedding, use_gpu, theta):
  # check whether we have one or more values of batch
  n_batches = len(adata.obs[batch_var].unique())
  if n_batches == 1:
    this_embedding = 'pca'
  else:
    this_embedding = embedding
  
  # move anndata to gpu if necessary
  if use_gpu:
    sc.get.anndata_to_GPU(adata)
  
  # start integration
  print(' scaling')
  sc.pp.scale(adata, max_value = 10)
  print(' PCA')
  sc.tl.pca(adata, n_comps = n_dims)
  
  sel_embed = 'X_pca'
  if this_embedding == 'harmony':
    print(' integrating with Harmony')
    if use_gpu:
      sc.pp.harmony_integrate(adata, key = batch_var, dtype=cp.float32, theta = theta)
    else:
      sce.pp.harmony_integrate(adata, key = batch_var, theta = theta)
    
    # harmony embedding instead of pca has to be used for umap and clustering
    sel_embed = 'X_pca_harmony'
  
  print(' running UMAP')
  sc.pp.neighbors(adata, n_pcs= n_dims,use_rep=sel_embed)
  sc.tl.umap(adata) # need to tell umap to use harmony

  print(' finding clusters')
  if not isinstance(res_ls, list):
    res_ls = [res_ls]
  for res in res_ls:
    if cl_method == 'leiden':
      sc.tl.leiden(
        adata, key_added=f"RNA_snn_res.{res}", resolution=float(res)
      )
    elif cl_method == 'louvain': # louvain not working in non-gpu mode
      sc.tl.louvain(
        adata, key_added=f"RNA_snn_res.{res}", resolution=float(res)
      ) 
  
  print(' recording clusters')
  clusts_dt = _get_clusts_from_adata(adata, this_embedding)
  
  print(' extracting other outputs')
  embeds_dt = _get_embeddings_from_adata(adata, this_embedding, sel_embed)

  int_dt = clusts_dt.join(embeds_dt, on = 'cell_id', how = 'inner')

  return int_dt


def _get_clusts_from_adata(adata, embedding):
  # get results
  clusts_dt = adata.obs.copy()
  clusts_dt = pl.from_pandas(clusts_dt)
  clusts_dt = clusts_dt.with_columns(pl.lit(embedding).alias('embedding'))
  
  # get clustering and id vars and subset
  cl_vs     = [col for col in clusts_dt.columns if re.match(r'RNA_snn_res.*', col)]
  all_cols  = cl_vs + ['embedding', 'cell_id', 'sample_id']
  clusts_dt = clusts_dt.select(all_cols)
  
  # get nice labels for clusters
  transform_exprs = []
  for cl_v in cl_vs:
    # count each cluster, put in order, make nice new cluster cluster names
    cl_lu = clusts_dt.select(
      cl_v
    ).filter(
      pl.col(cl_v).is_not_null()
    ).group_by(cl_v).agg(
      pl.len().alias("N")
    ).sort(
      "N", descending = True
    ).with_row_index(
      "rank", offset = 1
    ).with_columns(
      pl.format("cl{}", pl.col("rank").cast(pl.String).str.zfill(2)).alias( cl_v + ".tmp" )
    ).select([cl_v, cl_v + ".tmp"])

    # replace old values
    clusts_dt   = clusts_dt.join(
      cl_lu, on = cl_v, how = "left"
    ).drop(cl_v).rename({ cl_v + ".tmp": cl_v })

  return clusts_dt


def _get_embeddings_from_adata(adata, embedding,  sel_embed): 
  pca_array = adata.obsm[sel_embed]
  n_dims = pca_array.shape[1]

  # get pca dim names
  prefix = 'pca'
  if embedding == 'harmony': 
    prefix = 'hmny_pca'
  
  pca_col_names = [
    re.sub(r'_(\d)$', r'_0\1', f'{prefix}_{i}')
    for i in range(1, n_dims + 1)
  ]
   
  # make dt with reduced dims
  pca_dt = pl.DataFrame({
    'cell_id': adata.obs['cell_id'].to_numpy(), 
     **dict(zip(pca_col_names, pca_array.T))
  })
  
  umap_array = adata.obsm['X_umap']
  umap_col_names = ['UMAP1', 'UMAP2']
  umap_dt = pl.DataFrame({
    'cell_id': adata.obs['cell_id'].to_numpy(),
    **dict(zip(umap_col_names, umap_array.T))
  })

  # merge
  embeds_dt = pca_dt.join(umap_dt, on = 'cell_id', how = 'inner')
  return embeds_dt


def _calc_dbl_data(int_dbl, ok_cells_df, dbl_res, dbl_cl_prop):
  # get doublet cells
  dbl_ids       = ok_cells_df.filter(pl.col("dbl_class") == "doublet")["cell_id"].to_list()

  # get doublet cluster column
  dbl_res_str   = str(dbl_res)
  dbl_clust_col = f"RNA_snn_res.{dbl_res_str}"
  
  # make nice dataframe
  dbl_data      = int_dbl.select(
    'cell_id', 'sample_id',
    pl.col('UMAP1').alias('dbl_UMAP1'),  pl.col('UMAP2').alias('dbl_UMAP2'),
    pl.col(dbl_clust_col).alias('dbl_cluster')
  ).with_columns( 
    pl.col('cell_id').is_in(dbl_ids).alias('is_dbl'),
  ).with_columns(
    (pl.col('is_dbl').sum().over('dbl_cluster') /  
     pl.len().over('dbl_cluster')              
    ).alias('dbl_prop')
  ).with_columns(
    (pl.col('dbl_prop') > dbl_cl_prop).alias('in_dbl_cl')
  )

  return dbl_data


def _adata_filter_out_doublets(adata, dbl_data):
  # get non-doublet cell_ids
  ok_ids        = dbl_data.filter(
    (pl.col("is_dbl") == False) & (pl.col("in_dbl_cl") == False)
    )["cell_id"].to_list()

  # remove doublets from anndata object
  ok_ids_idx    = adata.obs['cell_id'].isin(ok_ids)
  adata         = adata[ok_ids_idx, :]

  return adata


if __name__ == "__main__":
  # get inputs
  parser = argparse.ArgumentParser()
  parser.add_argument('--hvg_mat_f',      type = str)
  parser.add_argument('--dbl_hvg_mat_f',  type = str)
  parser.add_argument('--sample_qc_f',    type = str)
  parser.add_argument('--coldata_f',      type = str)
  parser.add_argument('--demux_type',     type = str)
  parser.add_argument('--exclude_mito',   type = str)
  parser.add_argument('--embedding',      type = str)
  parser.add_argument('--n_dims',         type = int)
  parser.add_argument('--cl_method',      type = str)
  parser.add_argument('--dbl_res',        type = float)
  parser.add_argument('--dbl_cl_prop',    type = float)
  parser.add_argument('--theta',          type = float)
  parser.add_argument('--res_ls_concat',  type = str)
  parser.add_argument('--integration_f',  type = str)
  parser.add_argument('--batch_var',      type = str)
  parser.add_argument("-g", "--use_gpu",  action='store_true', 
    help='Use GPU-accelerated libraries if available.'  
  )
  args = parser.parse_args()

  # gpu vs cpu setup
  if args.use_gpu:
    print('using GPU')
    # import some GPU-specific modules
    import rapids_singlecell as sc
    import cupy as cp
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator

    # do some recommended setup
    rmm.reinitialize(
      managed_memory=False,  # Allows oversubscription (what is this?)
      pool_allocator=False,  # default is False; they in the vignette (https://rapids-singlecell.readthedocs.io/en/latest/notebooks/01_demo_gpu.html), they set this to True for harmony 
      devices=0,  # GPU device IDs to register. By default registers only GPU 0. (???)
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)
  else:
    print('not using GPU')
    # import some standard modules
    import scanpy as sc
    import scanpy.external as sce # for harmony

  # run
  run_integration(args.hvg_mat_f, args.dbl_hvg_mat_f, args.sample_qc_f, args.coldata_f,
    args.demux_type, args.exclude_mito, args.embedding, args.n_dims, args.cl_method,
    args.dbl_res, args.dbl_cl_prop, args.theta, args.res_ls_concat, args.integration_f,
    args.batch_var, args.use_gpu)

