import h5py
import numpy as np
import argparse
import re
import gzip
import gc
import scipy.sparse as sp
from scipy.sparse import csc_matrix, csr_matrix, hstack
import anndata as ad
import polars as pl


def run_integration(hvg_mat_f, dbl_hvg_mat_f, sample_qc_f, coldata_f, demux_type, 
  exclude_mito, embedding, n_dims, cl_method, dbl_res, dbl_cl_prop, theta, res_ls_concat,
  integration_f, batch_var, use_gpu = False, use_paga = False, paga_cl_res = None):
  print('setting up parameters')
  exclude_mito  = bool(exclude_mito)
  res_ls        = res_ls_concat.split()
  if demux_type == "none":
    dbl_batch_var = 'sample_id'
  else:
    dbl_batch_var = 'pool_id'
  dbl_theta     = 0

  print('loading hvg matrix')
  all_hvg_mat, bcs_passed, bcs_dbl = _get_hvg_mat(hvg_mat_f, dbl_hvg_mat_f)

  print('loading relevant cell ids')
  cells_df      = _get_cells_df(sample_qc_f, coldata_f, bcs_passed, demux_type, batch_var, bcs_dbl = bcs_dbl)

  print('normalizing hvg matrix')
  all_hvg_mat   = _normalize_hvg_mat(all_hvg_mat, cells_df, exclude_mito)

  print('making anndata object')
  adata_dbl     = ad.AnnData(X = all_hvg_mat.T , obs = cells_df.to_pandas())
  print(adata_dbl)
  print(f"  anndata object including doublets has {adata_dbl.shape[0]} cells and {adata_dbl.shape[1]} dims")

  print('running integration to find more doublets')
  int_dbl       = _do_one_integration(adata_dbl, dbl_batch_var, cl_method, n_dims,
    dbl_res, embedding, use_gpu, theta = dbl_theta, use_paga = use_paga, paga_cl = f"RNA_snn_res.{dbl_res}")
  
  del adata_dbl
  gc.collect()

  print('filter to non-doublet cells')
  dbl_data      = _calc_dbl_data(int_dbl, cells_df, dbl_res, dbl_cl_prop, demux_type, batch_var)
  adata         = _adata_filter_out_doublets(all_hvg_mat, cells_df, dbl_data)
  print(adata)
  print(f"  anndata object has {adata.shape[0]} cells and {adata.shape[1]} dims")

  del all_hvg_mat
  gc.collect()

  print('running integration on clean data')
  int_ok        = _do_one_integration(adata, batch_var, cl_method, n_dims,
    res_ls, embedding, use_gpu, theta = theta, use_paga = use_paga, paga_cl = f"RNA_snn_res.{paga_cl_res}")

  print('join results')
  int_df        = int_ok.join(dbl_data, on="cell_id", coalesce=True, how = 'full')

  print('save results')
  with gzip.open(integration_f, 'wb') as f:
    int_df.write_csv(f)

  print('done!')

  return int_df


def run_zoom_integration(hvg_mat_f, sample_qc_f, coldata_f, demux_type,
  exclude_mito, embedding, n_dims, cl_method, theta, res_ls_concat,
  integration_f, batch_var, use_gpu = False, use_paga = False, paga_cl_res = None):

  print('setting up parameters')
  exclude_mito  = bool(exclude_mito)
  res_ls        = res_ls_concat.split()

  print('loading hvg matrix')
  hvg_mat, bcs_passed, _ = _get_hvg_mat(hvg_mat_f)
  
  print('loading relevant cell ids')
  cells_df      = _get_cells_df(sample_qc_f, coldata_f, bcs_passed, demux_type, batch_var, zoom = True)

  print('normalizing hvg matrix')
  hvg_mat       = _normalize_hvg_mat(hvg_mat, cells_df, exclude_mito)

  print('making anndata object')
  adata         = ad.AnnData(X = hvg_mat.T , obs = cells_df.to_pandas())
  print(adata)
  print(f"  anndata object has {adata.shape[0]} cells and {adata.shape[1]} dims")

  print('running integration')
  int_df        = _do_one_integration(adata, batch_var, cl_method, n_dims, res_ls,
    embedding, use_gpu, theta, use_paga = use_paga, paga_cl = f"RNA_snn_res.{paga_cl_res}")

  print('save results')
  with gzip.open(integration_f, 'wb') as f:
    int_df.write_csv(f)

  print('done!')

  return int_df


def _get_hvg_mat(hvg_mat_f, dbl_hvg_mat_f = None):
  # get a matrix with hvgs (cells and doublets)
  all_hvg_mat = None
  bcs_passed  = []
  bcs_dbl     = []
  features    = []

  # open both
  for mat_f in [f for f in [hvg_mat_f, dbl_hvg_mat_f] if f is not None]:
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
      if mat_f == hvg_mat_f:
        bcs_passed  = barcodes
      elif mat_f == dbl_hvg_mat_f:
        bcs_dbl     = barcodes

  return all_hvg_mat, bcs_passed, bcs_dbl


def _get_cells_df(sample_qc_f, coldata_f, bcs_passed, demux_type, batch_var, zoom = False, bcs_dbl = []):
  # load files
  sample_qc   = pl.read_csv(sample_qc_f)
  all_coldata = pl.read_csv(coldata_f)

  # checks
  if not batch_var in all_coldata.columns:
    raise KeyError(f"column {batch_var} is missing from coldata file")
  if not batch_var in sample_qc.columns:
    raise KeyError(f"column {batch_var} is missing from sample QC file")

  # get ok samples
  bad_var     = "bad_" + batch_var
  ok_batches  = sample_qc.filter( pl.col(bad_var) == False )[batch_var].to_list()

  # subset to doublets
  if zoom:
    if not set(bcs_passed).issubset(set(all_coldata['cell_id'])):
      raise ValueError("Not all column names in hvg_mat are present in cell metadata.")
    cells_df    = all_coldata.filter( pl.col("cell_id").is_in(bcs_passed) )
  else:
    # get ok cells
    passed_idx  = all_coldata['keep']
    if not set(bcs_passed).issubset(set(all_coldata.filter(passed_idx)['cell_id'])):
      raise ValueError("qc-passed barcodes from hvg mats and cell_ids don't match")

    # get dbl cells
    if demux_type == "none":
      dbl_idx     = all_coldata["scdbl_class"] == "doublet"
    else:
      if batch_var == "sample_id":
        dbl_idx     = (all_coldata["scdbl_class"] == "doublet") | (all_coldata["demux_class"] == "doublet")
      elif batch_var == "pool_id":
        dbl_idx     = all_coldata["scdbl_class"] == "doublet"
    if not set(bcs_dbl).issubset(set(all_coldata.filter(dbl_idx)['cell_id'])):
      raise ValueError("doublet barcodes from hvg mats and cell_ids don't match")
    
    # add doublet label, filter to doublet or ok
    cells_df    = all_coldata.with_columns(
      pl.when(dbl_idx).then(True).otherwise(False).alias("is_dbl_int")
    ).filter(
      passed_idx | dbl_idx
    ).filter(
      ((pl.col(batch_var).is_in(ok_batches)) | (pl.col(batch_var).is_null()))
    )

  # put cell in coldata in the order matching mat cols
  all_bcs     = [*bcs_passed, *bcs_dbl]
  order_df    = pl.DataFrame({
    "cell_id":  all_bcs,
    "order":    range(1, len(all_bcs) + 1)
  })
  cells_df    = cells_df.join(order_df, on = 'cell_id').sort('order').drop('order')

  # check ok
  if not set(all_bcs) == set(cells_df['cell_id'].to_list()):
    raise ValueError("barcodes from hvg mats and cell_ids don't match")

  return cells_df


def _normalize_hvg_mat(hvg_mat, cells_df, exclude_mito, scale_f = 10000):
  if exclude_mito:
    cells_df    = cells_df.with_columns((pl.col("total") - pl.col("subsets_mito_sum")).alias("total_no_mito"))
    lib_sizes   = cells_df["total_no_mito"].to_numpy()
  else:
    lib_sizes   = cells_df["total"].to_numpy()

  hvg_mat     = hvg_mat.tocsr() 
  hvg_mat.data /= lib_sizes[ hvg_mat.indices ]
  hvg_mat.data *= scale_f

  np.log1p(hvg_mat.data, out=hvg_mat.data)

  return hvg_mat


def _do_one_integration(adata, batch_var, cl_method, n_dims, res_ls, embedding,
  use_gpu, theta, use_paga=False, paga_cl=None):
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
      sc.pp.harmony_integrate(adata, key = batch_var, max_iter_harmony = 5, 
        dtype=cp.float32, theta = theta)
    else:
      sce.pp.harmony_integrate(adata, key = batch_var, theta = theta)
    
    # harmony embedding instead of pca has to be used for umap and clustering
    sel_embed = 'X_pca_harmony'
  
  print(' finding neighbors')
  if np.isnan(adata.obsm[sel_embed]).any():
    raise ValueError("some NaN values in harmony output")
  sc.pp.neighbors(adata, n_pcs = n_dims, use_rep = sel_embed)

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

  # optionally run PAGA and use as init_pos for UMAP
  init_pos  = 'paga' if use_paga else 'spectral'
  if use_paga:
    import scanpy
    print(' running PAGA')
    scanpy.tl.paga(adata, groups = paga_cl)
    scanpy.pl.paga(adata)

  # run UMAP
  print(' running UMAP')
  sc.tl.umap(adata, maxiter = 750, init_pos = init_pos)

  if use_gpu:
    sc.get.anndata_to_CPU(adata)

  print(' recording clusters')
  clusts_df = _get_clusts_from_adata(adata, this_embedding, batch_var)
  
  print(' extracting other outputs')
  embeds_df = _get_embeddings_from_adata(adata, this_embedding, sel_embed)

  int_df = clusts_df.join(embeds_df, on = 'cell_id', how = 'inner')

  return int_df


def _get_clusts_from_adata(adata, embedding, batch_var):
  # get results
  clusts_df = adata.obs.copy()
  clusts_df = pl.from_pandas(clusts_df)
  clusts_df = clusts_df.with_columns(pl.lit(embedding).alias('embedding'))
  
  # get clustering and id vars and subset
  cl_vs     = [col for col in clusts_df.columns if re.match(r'RNA_snn_res.*', col)]
  sample_vs = list(set([batch_var, "sample_id"]))
  all_cols  = cl_vs + ['embedding', 'cell_id', *sample_vs]
  clusts_df = clusts_df.select(all_cols)

  # get nice labels for clusters
  for cl_v in cl_vs:
    # count each cluster, put in order, make nice new cluster cluster names
    cl_lu = clusts_df[ cl_v ].value_counts().sort(
      "count", descending = True
    ).with_row_index(
      "rank", offset = 1
    ).with_columns(
      pl.format("cl{}", pl.col("rank").cast(pl.String).str.zfill(2)).alias( cl_v + ".tmp" )
    ).select([cl_v, cl_v + ".tmp"])

    # replace old values
    clusts_df   = clusts_df.join(
      cl_lu, on = cl_v, how = "left"
    ).drop(cl_v).rename({ cl_v + ".tmp": cl_v })

  return clusts_df


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
  pca_df = pl.DataFrame({
    'cell_id': adata.obs['cell_id'].to_numpy(), 
     **dict(zip(pca_col_names, pca_array.T))
  })
  
  umap_array = adata.obsm['X_umap']
  umap_col_names = ['UMAP1', 'UMAP2']
  umap_df = pl.DataFrame({
    'cell_id': adata.obs['cell_id'].to_numpy(),
    **dict(zip(umap_col_names, umap_array.T))
  })

  # merge
  embeds_df = pca_df.join(umap_df, on = 'cell_id', how = 'inner')

  return embeds_df


def _calc_dbl_data(int_dbl, cells_df, dbl_res, dbl_cl_prop, demux_type, batch_var):
  # get doublet cells
  dbl_ids       = cells_df.filter( pl.col("is_dbl_int") == True )["cell_id"]

  # get doublet cluster column
  dbl_res_str   = str(dbl_res)
  dbl_clust_col = f"RNA_snn_res.{dbl_res_str}"
  
  # make nice dataframe
  dbl_data      = int_dbl.select(
    'cell_id',
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


def _adata_filter_out_doublets(all_hvg_mat, cells_df, dbl_data):
  # get ok ids
  ok_ids    = dbl_data.filter(
    (pl.col("is_dbl") == False) & (pl.col("in_dbl_cl") == False)
  )["cell_id"].to_list()
  
  # make adata, subset to ok ids
  adata     = ad.AnnData(X = all_hvg_mat.T, obs = cells_df.to_pandas())
  keep_idx  = adata.obs.cell_id.isin(ok_ids).to_numpy()
  adata     = adata[keep_idx, :].copy()
  adata.obs_names_make_unique()
  
  return adata


if __name__ == "__main__":
  # get inputs
  parser = argparse.ArgumentParser()
  subparsers = parser.add_subparsers(dest = "function_name", help = "Name of the function to run")

  # parser for run_integration
  parser_run_integration = subparsers.add_parser('run_integration')
  parser_run_integration.add_argument('--hvg_mat_f',      type = str)
  parser_run_integration.add_argument('--dbl_hvg_mat_f',  type = str)
  parser_run_integration.add_argument('--sample_qc_f',    type = str)
  parser_run_integration.add_argument('--coldata_f',      type = str)
  parser_run_integration.add_argument('--demux_type',     type = str)
  parser_run_integration.add_argument('--exclude_mito',   type = str)
  parser_run_integration.add_argument('--embedding',      type = str)
  parser_run_integration.add_argument('--n_dims',         type = int)
  parser_run_integration.add_argument('--cl_method',      type = str)
  parser_run_integration.add_argument('--dbl_res',        type = float)
  parser_run_integration.add_argument('--dbl_cl_prop',    type = float)
  parser_run_integration.add_argument('--theta',          type = float)
  parser_run_integration.add_argument('--res_ls_concat',  type = str)
  parser_run_integration.add_argument('--integration_f',  type = str)
  parser_run_integration.add_argument('--batch_var',      type = str)
  parser_run_integration.add_argument("-g", "--use-gpu",  action='store_true',
    help='Use GPU-accelerated libraries if available.'  
  )
  parser_run_integration.add_argument("-p", "--use-paga",  action='store_true',
    help='Use PAGA as initialization for UMAP.'
  )
  parser_run_integration.add_argument("--paga-cl-res",  type = str,
    help='The resolution of the PAGA cluster to use for initialization of UMAP.'
  )

  # parser for run_zoom_integration
  parser_run_zoom_integration = subparsers.add_parser('run_zoom_integration')
  parser_run_zoom_integration.add_argument('--hvg_mat_f',      type = str)
  parser_run_zoom_integration.add_argument('--sample_qc_f',    type = str)
  parser_run_zoom_integration.add_argument('--coldata_f',      type = str)
  parser_run_zoom_integration.add_argument('--demux_type',     type = str)
  parser_run_zoom_integration.add_argument('--exclude_mito',   type = str)
  parser_run_zoom_integration.add_argument('--embedding',      type = str)
  parser_run_zoom_integration.add_argument('--n_dims',         type = int)
  parser_run_zoom_integration.add_argument('--cl_method',      type = str)
  parser_run_zoom_integration.add_argument('--theta',          type = float)
  parser_run_zoom_integration.add_argument('--res_ls_concat',  type = str)
  parser_run_zoom_integration.add_argument('--integration_f',  type = str)
  parser_run_zoom_integration.add_argument('--batch_var',      type = str)
  parser_run_zoom_integration.add_argument("-g", "--use-gpu",  action='store_true',
    help='Use GPU-accelerated libraries if available.'  
  )
  parser_run_zoom_integration.add_argument("-p", "--use-paga",  action='store_true',
    help='Use PAGA as initialization for UMAP.'
  )
  parser_run_zoom_integration.add_argument("--paga-cl-res",  type = str,
    help='The resolution of the PAGA cluster to use for initialization of UMAP.'
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
  if args.function_name == 'run_integration':
    run_integration(args.hvg_mat_f, args.dbl_hvg_mat_f, args.sample_qc_f, args.coldata_f,
      args.demux_type, args.exclude_mito, args.embedding, args.n_dims, args.cl_method,
      args.dbl_res, args.dbl_cl_prop, args.theta, args.res_ls_concat, args.integration_f,
      args.batch_var, args.use_gpu, args.use_paga, args.paga_cl_res)
  elif args.function_name == 'run_zoom_integration':
    run_zoom_integration(args.hvg_mat_f, args.sample_qc_f, args.coldata_f,
      args.demux_type, args.exclude_mito, args.embedding, args.n_dims, args.cl_method,
      args.theta, args.res_ls_concat, args.integration_f,
      args.batch_var, args.use_gpu, args.use_paga, args.paga_cl_res)
  else:
    parser.print_help()
