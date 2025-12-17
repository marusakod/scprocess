import h5py
import numpy as np
import pandas as pd
import polars as pl
import argparse
import gzip
import os
import scanpy as sc
import scipy.sparse as sp
from scipy.sparse import csc_matrix, csr_matrix, hstack
import concurrent.futures
import re
from scanpy.preprocessing._utils import _get_mean_var
from functools import partial
import csv
from skmisc.loess import loess
import numba
import warnings


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


def get_one_csr_counts(run, hvg_paths_df, keep_df, smpl_stats_df, gene_ids, 
  RUN_VAR, batch_var, demux_type, chunk_size):
  # get input (ambient) file and output files
  filt_counts_f = hvg_paths_df.filter(pl.col(RUN_VAR) == run)["amb_filt_f"].item(0)
  out_fs        = hvg_paths_df.filter(pl.col(RUN_VAR) == run)["chunked_f"].to_list()

  # get bad samples
  bad_batches   = smpl_stats_df.filter( pl.col(f'bad_{batch_var}') == True )[batch_var].to_list()

  # get which batches we want
  if demux_type == "none":
    batches   = [run]
  else:
    batches   = hvg_paths_df.filter(pl.col(RUN_VAR) == run)[batch_var].to_list()

  # get valid barcodes
  bcs_dict = {}
  for b in batches:
    keep_cells  = keep_df.filter( pl.col(batch_var) == b )["cell_id"].to_list()
    bcs         = [cell_id.replace(run + ":", "", 1) for cell_id in keep_cells]
    bcs         = np.array(bcs)
    bcs_dict[b] = bcs

  # open input file
  with h5py.File(filt_counts_f, 'r') as f:
    indptr      = f['matrix/indptr'][:]
    indices     = f['matrix/indices'][:]
    data        = f['matrix/data'][:]
    features    = f['matrix/features/name'][:]
    barcodes    = f['matrix/barcodes'][:].astype('U16')
    num_rows    = f['matrix/shape'][0]
    num_cols    = f['matrix/shape'][1]

  # make a csc sparse matrix
  sua_csc = csc_matrix((data, indices, indptr), shape=(num_rows, num_cols))

  # loop through batches
  for b, out_f in zip(batches, out_fs):
    if b in bad_batches:
      print(f"{batch_var} {b}: No cells after QC; writing an empty file: {out_f}.")
      open(out_f, "w").close()
      continue
    
    # get indices of barcodes to keep
    bcs       = bcs_dict[b]
    keep_idx  = np.where(np.isin(barcodes, bcs))[0]
    if sum(keep_idx) == 0:
      ValueError("none of the selected barcodes found in h5 file :/")
    filt_bcs  = barcodes[keep_idx]
    
    # subset matrix
    sua_csc_qc = sua_csc[:, keep_idx]

    # merge splices, unspliced, ambiguous
    csc, uniq_features = sum_SUA(sua_csc_qc, features)
    
    # get indices of genes to keep
    gene_ids      = np.array(gene_ids)
    gs_keep_idx   = np.where(np.isin(uniq_features, gene_ids))[0]
    uniq_features = uniq_features[gs_keep_idx]
    if len(uniq_features) == 0:
      raise ValueError(f"no features selected for {batch_var} {b}")
    csc       = csc[gs_keep_idx, :]

    # convert to csr
    csr = csc.tocsr()

    # save csr matrix in chunks
    with h5py.File(out_f, 'w') as f:
      indptr_chunk_size = min(chunk_size + 1, len(csr.indptr))
      indices_chunk_size = min(len(csr.indices), chunk_size * num_cols)
      data_chunk_size = min(len(csr.data), chunk_size * num_cols)

      f.create_dataset('matrix/indptr', data=csr.indptr, compression='gzip', chunks=(indptr_chunk_size,))
      f.create_dataset('matrix/indices', data=csr.indices, compression='gzip', chunks=(indices_chunk_size,))
      f.create_dataset('matrix/data', data=csr.data, compression='gzip', chunks=(data_chunk_size,))
      f.create_dataset('matrix/shape', data=csr.shape, compression='gzip')

      # Saving the feature and barcode names
      f.create_dataset('matrix/features/name', data=np.array(uniq_features, dtype='S'), compression='gzip')
      f.create_dataset('matrix/barcodes', data=np.array(filt_bcs, dtype='S'), compression='gzip')
      
    print(f"CSR matrix for {b} successfully saved to {out_f}.")


def get_csr_counts(hvg_paths_f, cell_filter_f, keep_var, keep_vals_str, smpl_stats_f, 
  rowdata_f, RUN_VAR, batch_var, demux_type, chunk_size = 2000, n_cores = 8):
  # load things
  hvg_paths_df  = pl.read_csv(hvg_paths_f)
  smpl_stats_df = pl.read_csv(smpl_stats_f)

  # get QCed cells
  filt_df       = pl.read_csv(cell_filter_f)
  keep_vals     = keep_vals_str.split(',')
  keep_df       = (
    filt_df
    .with_columns(pl.col(keep_var).cast(pl.Utf8).alias(keep_var))  
    .filter(pl.col(keep_var).is_in(keep_vals)) 
  )
  if not keep_df.shape[0] > 0:
    raise ValueError("keep_df is empty")

  # get gene details
  rows_df       = pl.read_csv(rowdata_f)
  keep_ids      = rows_df['ensembl_id'].to_list()

  # define list of samples
  runs          = hvg_paths_df[RUN_VAR].unique().to_list()
 
  with concurrent.futures.ThreadPoolExecutor(max_workers=n_cores) as executor:
    futures = [executor.submit(get_one_csr_counts, run, hvg_paths_df, keep_df, 
      smpl_stats_df, keep_ids, RUN_VAR, batch_var, demux_type, chunk_size) for run in runs]

    # some more parallel stuff i guess
    for future in concurrent.futures.as_completed(futures):
      future.result()

  return


def read_csr_chunk(file, start_row, end_row):
  # open file
  with h5py.File(file, 'r') as f:
    num_cols  = f['matrix/shape'][1]
    indptr    = f['matrix/indptr'][start_row:end_row+1]
    start_pos = indptr[0]
    end_pos   = indptr[-1]
    data      = f['matrix/data'][start_pos:end_pos]
    indices   = f['matrix/indices'][start_pos:end_pos]

    # get features and barcodes for chunk
    features  = f['matrix/features/name'][start_row:end_row].astype(str)
    barcodes  = f['matrix/barcodes'][:].astype(str)

    # adjust indptr to start at zero for the chunk
    indptr_chunk = indptr - start_pos

  # csr for the chunk
  sparse_chonk  = csr_matrix((data, indices, indptr_chunk), shape=(end_row - start_row, num_cols))
  
  return sparse_chonk, features, barcodes


def read_full_csr(file): 
  # open file
  with h5py.File(file, 'r') as f:
    indptr    = f['matrix/indptr'][:]  
    indices   = f['matrix/indices'][:]  
    data      = f['matrix/data'][:]  
    features  = f['matrix/features/name'][:]  
    barcodes  = f['matrix/barcodes'][:]  

    num_rows  = f['matrix/shape'][0]  
    num_cols  = f['matrix/shape'][1]  

    sparse_csr  = csr_matrix((data, indices, indptr), shape=(num_rows, num_cols))

  return sparse_csr, features, barcodes


def _calculate_feature_stats(sparse_csr, features, rowdata_f):
  # calculate mean and variance for sparse matrix, put into df
  mean, var = _get_mean_var(sparse_csr, axis=1 )
  feat_stats_df = pl.DataFrame({
    'ensembl_id': features.astype(str),
    'mean':       mean,
    'variance':   var
  })

  # merge with nice labels for genes
  rows_df       = pl.read_csv(rowdata_f)
  feat_stats_df = feat_stats_df.join(rows_df, on = 'ensembl_id')
  
  return feat_stats_df


# minor modifiction of scanpy _sum_and_sum_squares_clipped
def _calculate_regularized_variance(sparse_csr, features, input_df):
  # put input in features order
  order_df  = pl.DataFrame({"ensembl_id": features.astype(str)}).with_row_index("row_ord")
  input_df  = input_df.join(order_df, on = "ensembl_id").sort("row_ord")
  clip_val  = input_df['clip_val'].to_numpy()

  # convert to column format, get values we need
  csc       = sparse_csr.tocsc()
  indices   = csc.indices
  data      = csc.data
  nnz       = csc.nnz

  # set up some variables
  n_rows        = sparse_csr.shape[0]
  sq_counts_sum = np.zeros(n_rows, dtype=np.float64)
  counts_sum    = np.zeros(n_rows, dtype=np.float64)

  # loop through non-zero values to get squared sums of clipped values
  for i in numba.prange(nnz):
    idx                 = indices[i]
    element             = min(np.float64(data[i]), clip_val[idx])
    sq_counts_sum[idx]  += element**2
    counts_sum[idx]     += element

  # make into df
  clipped_df  = pl.DataFrame({
    "ensembl_id":         np.array(features, dtype=str),
    "squared_counts_sum": sq_counts_sum,
    "counts_sum":         counts_sum
  })

  # make output
  output_df   = input_df.join(clipped_df, on = "ensembl_id")
  output_df   = output_df.with_columns(
    (((pl.col('n_cells') * pl.col('mean').pow(2)) + pl.col('squared_counts_sum') - 2 * pl.col('counts_sum') * pl.col('mean')) / 
          ((pl.col('n_cells') - 1) * pl.col('reg_std').pow(2))).alias('variances_norm')
  )

  return output_df


def calculate_estimated_vars(estim_vars_f, hvg_method, batch_var, mean_var_merged_f = None,
  mean_var_df = None, min_mean = 0.01, span: float = 0.3 ):
  # read mean var df if not specified
  if mean_var_df is None:
    assert mean_var_merged_f is not None, 'input file path or dataframe missing'
    stats_df  = pl.read_csv(mean_var_merged_f)
  else:
    stats_df = mean_var_df

  # group on relevant variable
  group_var   = batch_var if hvg_method == 'sample' else 'group'
  grouped     = stats_df.partition_by( group_var, as_dict = True )

  # go through each group
  all_results = []
  for name, group in grouped.items():
    # get values we need
    N           = group['n_cells'].unique()[0]
    num_rows    = group.shape[0]
    mean        = group['mean'].to_numpy()
    var         = group['variance'].to_numpy()
     
    # this code is a minor modification of equivalent scanpy code:
    # https://github.com/scverse/scanpy/blob/1.11.1/src/scanpy/preprocessing/_highly_variable_genes.py#L503-L717
    not_const   = var > 0
    y           = np.log10(var[not_const]) 
    x           = np.log10(mean[not_const])
    data_df     = pl.DataFrame({ "x": x, "y": y }).unique()
    model       = safe_loess( data_df["x"].to_numpy(), data_df["y"].to_numpy(), span = span )

    # store trended variance values, also get std deviation
    est_var     = np.zeros(num_rows, dtype=np.float64)
    est_var[not_const] = model.predict(x).values
    reg_std     = np.sqrt(10**est_var)

    # clip large values as in Seurat
    vmax        = np.sqrt(N)
    clip_val    = reg_std * vmax + mean
  
    # store values      
    group = group.with_columns(
      pl.lit(reg_std).alias('reg_std'),
      pl.lit(clip_val).alias('clip_val'),
      pl.lit(est_var).alias('estimated_variance')
    )
    all_results.append(group)
  
  # join into one nice dataframe
  final_df = pl.concat(all_results)

  # save or just pass back?
  if estim_vars_f is not None:
    with gzip.open(estim_vars_f, 'wb') as f:
      final_df.write_csv(f)
    return
  else:
    return final_df
  

# we have this function because sometimes there are v many points at 0, and this
# means that running loess with a span of 0.3 falls over, as it seeks to estimate 
# a curve on the nearest 30% of the datapoints. to get round this, we cheat by adding
# the minimal possible random jitter to all points.
def safe_loess(x, y, span, initial_amount=1e-16, max_attempts=5, seed=1234):
  # let's get ready
  attempts = 0
  amount   = initial_amount
  
  # try with increasing size of jitter
  while attempts < max_attempts:
    try:
      # first go is just normal data
      if attempts == 0:
        x_jitter = x
      else:
        # after that we add jitter
        np.random.seed(seed)
        jitter    = np.random.uniform(-amount, amount, len(x))
        x_jitter  = x + jitter
      
      # let's fit the model
      model     = loess(x_jitter, y, span=span, degree=2)
      model.fit()

      return model
    
    except Exception as e:
      # INCREASE JITTER!!!
      warnings.warn(f"Attempt {attempts + 1} failed: {e}. Increasing jitter amount to {amount * 10}.")
      amount *= 10
    
    attempts += 1
  
  # even jitter wasn't enough :(
  raise RuntimeError(f"Failed to fit loess model after {max_attempts} attempts with increasing jitter.")


def get_chunk_params(hvg_paths_f, rowdata_f, metadata_f, qc_smpl_stats_f, chunk_num, 
  hvg_method, batch_var, chunk_size=2000, group_var=None, group=None):
  # load rows
  total         = pl.read_csv(rowdata_f).shape[0]

  # get input paths df
  hvg_paths_df  = pl.read_csv(hvg_paths_f)
  
  # get list of good samples
  qc_df         = pl.read_csv(qc_smpl_stats_f)
  good_batches  = qc_df.filter(pl.col(f'bad_{batch_var}') == False)[batch_var].to_list()

  # get input files for selected samples
  if hvg_method == 'all':
    files = hvg_paths_df.loc[hvg_paths_df['sample_id'].isin(good_batches), 'chunked_f'].tolist()
  else:
    assert group_var is not None, "group_var must be defined."
    assert group is not None, "group must be defined."
    
    # select samples based on group
    meta = pl.read_csv(metadata_f)
    grp_samples = meta.loc[meta[group_var] == group, 'sample_id'].tolist()
    
    sel_samples = list(set(grp_samples) & set(good_batches))
    files = hvg_paths_df.loc[hvg_paths_df['sample_id'].isin(sel_samples), 'chunked_f'].tolist()

   # calculate start and end rows for the chunk
  start_row = chunk_num * chunk_size
  end_row = min(start_row + chunk_size, total)

  return files, total, start_row, end_row


def calculate_std_var_stats_for_sample(batch, batch_var, qc_smpl_stats_f, csr_f, rowdata_f, std_var_stats_f):
  # get data
  qc_df       = pl.read_csv(qc_smpl_stats_f)
  bad_batches = qc_df.filter( pl.col(f"bad_{batch_var}") == True )[batch_var].to_list()
   
  if batch in bad_batches:
    # save an empty file
    open(std_var_stats_f, 'w').close()
    return
  
  # read sparse matrix
  sparse_csr, features, _ = read_full_csr(csr_f)
  
  # calculate mean and variance
  n_cells     = sparse_csr.shape[1]
  stats_df    = _calculate_feature_stats(sparse_csr, features, rowdata_f).with_columns(
    pl.lit(batch).alias(batch_var),
    pl.lit(n_cells).alias('n_cells')
  )

  # calculate estimated variance
  estim_vars_df = calculate_estimated_vars(
    estim_vars_f      = None,
    mean_var_merged_f = None,
    hvg_method        = 'sample',
    batch_var         = batch_var,
    mean_var_df       = stats_df
    )
  
  # calculate normalised variance
  stats_df  = _calculate_regularized_variance(sparse_csr, features, estim_vars_df)

  # save
  with gzip.open(std_var_stats_f, 'wb') as f:
    stats_df.write_csv(f)

  return


def calculate_mean_var_for_chunk(hvg_paths_f, rowdata_f, metadata_f, qc_smpl_stats_f, 
  mean_var_f, chunk_num, hvg_method, batch_var, chunk_size=2000, group_var=None, 
  group=None, n_cores = 8):
  
  files, total, start_row, end_row = get_chunk_params(
    hvg_paths_f, rowdata_f, metadata_f, qc_smpl_stats_f, 
    chunk_num, hvg_method, batch_var, chunk_size, group_var, group
    )
  
  if len(files) == 0:
    if hvg_method == 'group':
      print(f"No (good) {batch_var}s found in group '{group}'. Writing an empty output file.")
      open(mean_var_f, 'w').close()
    else:
      print("Error: no files selected for calculating hvg statistics")
    return
  
  if start_row >= total:
    print(f"Start row {start_row} exceeds total number of genes ({total}). Writing an empty output file.")
    open(mean_var_f, 'w').close()
    return
  
  # read chunks for multiple samples in parallel
  merged_chunk = None
  features = []
  
  read_csr_chunk_p = partial(read_csr_chunk, start_row=start_row, end_row=end_row)
  with concurrent.futures.ThreadPoolExecutor(max_workers=n_cores) as executor:
    results = list(executor.map(read_csr_chunk_p, files))
    
    for csr_chunk, chunk_feats, _ in results:
      if len(features) == 0:
        features = chunk_feats
      
      # Merge chunks column-wise
      merged_chunk = csr_chunk if merged_chunk is None else hstack([merged_chunk, csr_chunk])


  # get number of cells in chunk
  n_cells = merged_chunk.shape[1] 

  # calculate stats
  merged_chunk_stats = _calculate_feature_stats(merged_chunk, features, rowdata_f).with_columns(
    pl.lit(n_cells).alias('n_cells'),
    pl.lit(group).alias('group')
    )

  # save
  with gzip.open(mean_var_f, 'wb') as f:
    merged_chunk_stats.write_csv(f)

  return


def calculate_std_var_stats_for_chunk(hvg_paths_f, rowdata_f, metadata_f, qc_smpl_stats_f, 
  std_var_stats_f, estim_vars_f, chunk_num, hvg_method, batch_var, chunk_size = 2000, group_var = None, 
  group = None, n_cores = 8):
  # get chunk parameters
  files, total, start_row, end_row = get_chunk_params(
    hvg_paths_f, rowdata_f, metadata_f, qc_smpl_stats_f,
    chunk_num, hvg_method, batch_var, chunk_size, group_var, group
    )
  
  if len(files) == 0:
    if hvg_method == 'group':
      print(f"No (good) {batch_var} found in group '{group}'. Writing an empty output file.")
      open(std_var_stats_f, 'w').close()
    else:
      print("Error: no files selected for calculating hvg statistics")
    return
  
  if start_row >= total:
    print(f"Start row {start_row} exceeds total number of genes ({total}). Writing an empty output file.")
    open(std_var_stats_f, 'w').close()
    return
  
  # get estimated variances
  estim_vars_df = pl.read_csv(estim_vars_f)
  
  # read chunks for multiple samples in parallel
  merged_chunk  = None
  features      = []
  
  # set up parallel stuff
  read_csr_chunk_p = partial(read_csr_chunk, start_row=start_row, end_row=end_row)
  with concurrent.futures.ThreadPoolExecutor(max_workers=n_cores) as executor:
    results = list(executor.map(read_csr_chunk_p, files))
    # loop over chunks
    for csr_chunk, chunk_feats, _ in results:
      if len(features) == 0:
        features = chunk_feats
      
      # merge chunks column-wise
      merged_chunk = csr_chunk if merged_chunk is None else hstack([merged_chunk, csr_chunk])

  # calculate stats 
  chunk_estim_vars    = estim_vars_df.filter( (pl.col('group') == group) & (pl.col('ensembl_id').is_in(features)) )
  merged_chunk_stats  = _calculate_regularized_variance(merged_chunk, features, merged_chunk_stats)

  # save
  with gzip.open(std_var_stats_f, 'wb') as f:
    merged_chunk_stats.write_csv(f)

  return


def _calculate_standardized_variance(df):
  df = df.with_columns(
    (((pl.col('n_cells') * pl.col('mean').pow(2)) + pl.col('squared_counts_sum') - 2 * pl.col('counts_sum') * pl.col('mean')) / 
          ((pl.col('n_cells') - 1) * pl.col('reg_std').pow(2))).alias('variances_norm')
  )
  return df


def _process_single_group(stats_df, empty_gs, n_hvgs, exclude_ambient):
  # get genes to exclude
  exclude_gs = []
  if exclude_ambient:
    exclude_gs = empty_gs
  else:
    exclude_gs = []

  out_df      = stats_df.select(['gene_id', 'variances_norm'])
  out_df      = out_df.with_columns(
    pl.when( pl.col('gene_id').is_in(exclude_gs) ).then(None).otherwise('variances_norm').alias('variance_tmp')
  )
  out_df      = out_df.sort('variances_norm', descending = True)
  out_df      = out_df.with_row_index('highly_variable_rank', offset = 1)
  out_df      = out_df.with_columns(
    pl.when(pl.col('highly_variable_rank') <= n_hvgs).then(1).otherwise(0).alias("highly_variable_nbatches")
  )

  return out_df


def _process_multiple_groups(stats_df, group_var, empty_gs, n_hvgs,  exclude_ambient):
  # decide whether to exclude
  exclude_gs = []
  if exclude_ambient:
    exclude_gs = empty_gs
  else:
    exclude_gs = []

  # do like for one group, but for multiple
  tmp_df      = stats_df.select([group_var, 'gene_id', 'variances_norm'])
  tmp_df      = tmp_df.with_columns(
    pl.when( pl.col('gene_id').is_in(exclude_gs) ).then(None).otherwise('variances_norm').alias('variance_tmp')
  )
  tmp_df      = tmp_df.sort([group_var, 'variances_norm'], descending = [False, True])
  tmp_df      = tmp_df.with_columns(
    highly_variable_rank  = pl.col("variance_tmp").rank(method = "average", descending = True).over(group_var)
  )

  # aggregate
  out_df      = tmp_df.group_by("gene_id").agg(
    (pl.col("highly_variable_rank") <= n_hvgs).sum().alias("highly_variable_nbatches"),
    pl.col("highly_variable_rank").median().alias("highly_variable_rank"),
    pl.col("variances_norm").mean().alias("variances_norm")
  )

  # sort
  out_df      = out_df.sort(['highly_variable_nbatches', 'highly_variable_rank'], descending = [True, False])

  return out_df


# main function to calculate highly variable genes
def calculate_hvgs(std_var_stats_f, hvg_f, empty_gs_f, hvg_method, batch_var, n_hvgs, exclude_ambient=True):
  # get empty genes
  empty_dt  = pl.read_csv(empty_gs_f)
  empty_gs  = empty_dt.filter( pl.col("is_ambient") == True )["gene_id"].to_list()

  # get stats
  group_var = batch_var if hvg_method == 'sample' else 'group'
  stats_df  = pl.read_csv(std_var_stats_f)

  # find HVGs for each
  if stats_df[ group_var ].n_unique() == 1:
    hvg_df = _process_single_group(stats_df, empty_gs, n_hvgs, exclude_ambient)
  else:
    hvg_df = _process_multiple_groups(stats_df, group_var, empty_gs, n_hvgs, exclude_ambient)

  # sort hvg_df nicely
  sort_cols = ["highly_variable_nbatches", "highly_variable_rank"]
  sort_desc = [True, False]
  hvg_df    = hvg_df.sort(sort_cols, descending = sort_desc, nulls_last = True)

  # finally add label if in top n_hvgs
  hvg_df    = hvg_df.with_row_index("rank_tmp", offset = 1).with_columns(
    pl.when( pl.col("rank_tmp") <= n_hvgs ).then(True).otherwise(False).alias("highly_variable")
  ).drop("rank_tmp")
  if hvg_df["highly_variable"].sum() != n_hvgs:
    raise ValueError(f"somehow more than %d HVGs selected")

  # save
  with gzip.open(hvg_f, 'wb') as f:
    hvg_df.write_csv(f)

  return hvg_df


# read top 2000 hvgs from each sample and save file
def create_hvg_matrix(qc_smpl_stats_f, hvg_paths_f, hvg_f, out_h5_f, demux_type, batch_var):
  # get bad samples
  qc_sample_df  = pl.read_csv(qc_smpl_stats_f)
  bad_batches   = qc_sample_df.filter( pl.col(f"bad_{batch_var}") == True)[ batch_var ].to_list()

  # get all chunked files
  hvg_paths_df  = pl.read_csv(hvg_paths_f)
  chunked_fs    = hvg_paths_df['chunked_f'].to_list()
  batches       = hvg_paths_df[ batch_var ].to_list()
  if demux_type == "none": 
    pools         = batches
  else:
    pools         = hvg_paths_df['pool_id'].to_list()

  # get all hvgs
  hvg_df        = pl.read_csv(hvg_f)
  hvg_ids       = hvg_df.filter(pl.col('highly_variable') == True)['gene_id'].to_list()

  # extract ensembl ids 
  hvg_ensembl   = []
  for gene in hvg_ids:
    parts = gene.rsplit('_', 1)
    hvg_ensembl.append(parts[-1])
  
  # initialise hvg matrix
  top_genes_mat = None
  all_barcodes  = []

  # open each file separately and extract highly variable genes
  for f, b, p in zip(chunked_fs, batches, pools):
    if b in bad_batches:
      continue
    sample_csr, features, barcodes = read_full_csr(f)

    features    = np.array(features, dtype=str)
    hvg_indices = [i for i, feature in enumerate(features) if feature in hvg_ensembl]
    csr_chunk   = sample_csr[hvg_indices, :]

    barcodes    = barcodes.astype('<U21')
    barcodes    = [f"{p}:{bc}" for bc in barcodes]  
    barcodes    = np.array(barcodes)

    # merge to other chunks column-wise
    all_barcodes.extend(barcodes)
    top_genes_mat = csr_chunk if top_genes_mat is None else hstack([top_genes_mat, csr_chunk])

  top_genes_csc = top_genes_mat.tocsc()
  
  # save to a new h5 file
  with h5py.File(out_h5_f, 'w') as f:
    f.create_dataset('matrix/data', data=top_genes_csc.data)
    f.create_dataset('matrix/indices', data=top_genes_csc.indices)
    f.create_dataset('matrix/indptr', data=top_genes_csc.indptr)
    f.create_dataset('matrix/shape', data=top_genes_csc.shape)
    f.create_dataset('matrix/features/name', data=np.array(hvg_ensembl, dtype='S'))
    f.create_dataset('matrix/barcodes', data=np.array(all_barcodes, dtype='S'))


def create_doublets_matrix(hvg_paths_f, hvg_f, qc_f, qc_smpl_stats_f, out_h5_f, RUN_VAR, demux_type, batch_var):
  # get all hvgs
  hvg_df    = pl.read_csv(hvg_f)
  hvg_ids   = hvg_df.filter( pl.col('highly_variable') == True)['gene_id'].to_list()
  
  # get qc file with all cells
  qc_df     = pl.read_csv(qc_f)

  # subset to doublets
  if demux_type == "none":
    dbl_df    = qc_df.filter( pl.col("scdbl_class") == "doublet" )
  else:
    if batch_var == "sample_id":
      dbl_df    = qc_df.filter( (pl.col("scdbl_class") == "doublet") | (pl.col("demux_class") == "doublet") )
    elif batch_var == "pool_id":
      dbl_df    = qc_df.filter( pl.col("scdbl_class") == "doublet" )

  # extract ensembl ids (maybe keep ensembl is in the df earlier)
  hvg_ensembl   = []
  for gene in hvg_ids:
    parts = gene.rsplit('_', 1)
    hvg_ensembl.append(parts[-1])
  hvg_ensembl   = np.array(hvg_ensembl)
  
  # initialise doublet matrix
  doublet_mat   = None
  all_barcodes  = []

  # get samples (or pools) that passed qc
  hvg_paths_df  = pl.read_csv(hvg_paths_f)
  qc_sample_df  = pl.read_csv(qc_smpl_stats_f)
  good_batches  = qc_sample_df.filter( pl.col(f'bad_{batch_var}') == False )[batch_var]
  keep_runs     = hvg_paths_df.filter( pl.col(batch_var).is_in(good_batches.to_list()) )[RUN_VAR].unique().to_list()

  # loop through these
  for run in keep_runs:
    # get doublets for this run
    dbls     = dbl_df.filter( pl.col(RUN_VAR) == run )["cell_id"].to_list()
    dbl_ids  = np.array(dbls)

    # get ambient output file
    filt_counts_f = hvg_paths_df.filter(pl.col(RUN_VAR) == run)["amb_filt_f"].item(0)

    # open input file
    with h5py.File(filt_counts_f, 'r') as f:
      indptr      = f['matrix/indptr'][:]
      indices     = f['matrix/indices'][:]
      data        = f['matrix/data'][:]
      features    = f['matrix/features/name'][:]
      barcodes    = f['matrix/barcodes'][:]

      num_rows    = f['matrix/shape'][0]
      num_cols    = f['matrix/shape'][1]
     
    # add sample ids to barcodes
    barcodes    = barcodes.astype('<U21')
    barcodes    = [f"{run}:{bc}" for bc in barcodes]  
    barcodes    = np.array(barcodes)

    # make a csc sparse matrix
    sua_csc     = csc_matrix((data, indices, indptr), shape=(num_rows, num_cols))
    
    # get indices of barcodes to keep
    keep_idx    = np.where(np.isin(barcodes, dbl_ids))[0]
    filt_bcs    = barcodes[keep_idx]

    # subset matrix
    sua_csc_dbl = sua_csc[:, keep_idx]

    # merge splices, unspliced, ambiguous
    csc, uniq_features = sum_SUA(sua_csc_dbl, features)

    # get indices of highly variable genes
    features    = features.astype('<U21')
    hvg_indices = [i for i, feature in enumerate(uniq_features) if feature in hvg_ensembl]
    csc         = csc[hvg_indices, :]

    # combine matrices and barcodes
    all_barcodes.extend(filt_bcs)
    doublet_mat = csc if doublet_mat is None else hstack([doublet_mat, csc])

  # save to a new h5 file
  with h5py.File(out_h5_f, 'w') as f:
    f.create_dataset('matrix/data', data=doublet_mat.data)
    f.create_dataset('matrix/indices', data=doublet_mat.indices)
    f.create_dataset('matrix/indptr', data=doublet_mat.indptr)
    f.create_dataset('matrix/shape', data=doublet_mat.shape)
    f.create_dataset('matrix/features/name', data=np.array(hvg_ensembl, dtype='S'))
    f.create_dataset('matrix/barcodes', data=np.array(all_barcodes, dtype='S'))

  return


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  subparsers = parser.add_subparsers(dest = "function_name", help = "Name of the function to run")

  # parser for get_csr_counts
  parser_makeCSR = subparsers.add_parser('get_csr_counts')
  parser_makeCSR.add_argument("hvg_paths_f", type=str)
  parser_makeCSR.add_argument("cell_filter_f", type=str)
  parser_makeCSR.add_argument("keep_var", type=str)
  parser_makeCSR.add_argument("keep_vals_str", type=str)
  parser_makeCSR.add_argument("smpl_stats_f", type=str)
  parser_makeCSR.add_argument("rowdata_f", type=str)
  parser_makeCSR.add_argument("run_var", type=str)
  parser_makeCSR.add_argument("batch_var", type=str)
  parser_makeCSR.add_argument("demux_type", type=str)
  parser_makeCSR.add_argument("-s", "--chunksize", required=False, default= 2000, type = int)
  parser_makeCSR.add_argument("-n", "--ncores", required=False, default= 8, type = int)
  

  # parser for calculate_mean_var_for_chunk
  parser_chunkCalcs = subparsers.add_parser('calculate_mean_var_for_chunk')
  parser_chunkCalcs.add_argument("hvg_paths_f", type=str)
  parser_chunkCalcs.add_argument("rowdata_f", type=str)
  parser_chunkCalcs.add_argument("metadata_f", type=str)
  parser_chunkCalcs.add_argument("qc_smpl_stats_f", type=str)
  parser_chunkCalcs.add_argument("mean_var_f", type=str)
  parser_chunkCalcs.add_argument("chunk", type=int)
  parser_chunkCalcs.add_argument("hvg_method", type=str)
  parser_chunkCalcs.add_argument("batch_var", type=str)
  parser_chunkCalcs.add_argument("-s", "--chunksize", type=int, required=False, default=2000)
  parser_chunkCalcs.add_argument("-v", "--groupvar", type=str, required=False, default=None)
  parser_chunkCalcs.add_argument("-g", "--group", type=str, required=False, default= None)
  parser_chunkCalcs.add_argument("-n", "--ncores", type=int, required=False, default = 8)

  # parser for calculate_estimated_vars
  parser_varEstim = subparsers.add_parser('calculate_estimated_vars')
  parser_varEstim.add_argument("estim_vars_f", type=str)
  parser_varEstim.add_argument("hvg_method", type=str)
  parser_varEstim.add_argument("batch_var", type=str)
  parser_varEstim.add_argument("mean_var_merged_f", type=str)

  # parser for calculate_std_var_stats_for_sample
  parser_sampleSssc = subparsers.add_parser('calculate_std_var_stats_for_sample')
  parser_sampleSssc.add_argument("sample", type=str)
  parser_sampleSssc.add_argument("batch_var", type=str)
  parser_sampleSssc.add_argument("qc_smpl_stats_f", type=str)
  parser_sampleSssc.add_argument("csr_f", type=str)
  parser_sampleSssc.add_argument("rowdata_f", type=str)
  parser_sampleSssc.add_argument("std_var_stats_f", type=str)


  # parser for calculate_std_var_stats_for_chunk
  parser_chunkSssc = subparsers.add_parser('calculate_std_var_stats_for_chunk')
  parser_chunkSssc.add_argument("hvg_paths_f", type=str)
  parser_chunkSssc.add_argument("rowdata_f", type=str)
  parser_chunkSssc.add_argument("metadata_f", type=str)
  parser_chunkSssc.add_argument("qc_smpl_stats_f", type=str)
  parser_chunkSssc.add_argument("std_var_stats_f", type=str)
  parser_chunkSssc.add_argument("estim_vars_f", type=str)
  parser_chunkSssc.add_argument("chunk_num", type=int)
  parser_chunkSssc.add_argument("hvg_method", type=str)
  parser_chunkSssc.add_argument("batch_var", type=str)
  parser_chunkSssc.add_argument("-s", "--chunksize", type=int, required=False, default=2000)
  parser_chunkSssc.add_argument("-v", "--groupvar", type=str, required=False, default=None)
  parser_chunkSssc.add_argument("-g", "--group", type=str, required=False, default= None)
  parser_chunkSssc.add_argument("-n", "--ncores", type=int, required=False, default = 8)

  
  # parser for calculate_hvgs(std_var_stats_f, hvg_f, empty_gs_fs, hvg_method, n_hvgs, exclude_empty=True
  parser_getHvgs = subparsers.add_parser('calculate_hvgs')
  parser_getHvgs.add_argument("std_var_stats_f", type=str)
  parser_getHvgs.add_argument("hvg_f", type=str)
  parser_getHvgs.add_argument("empty_gs_f", type=str)
  parser_getHvgs.add_argument("hvg_method", type=str)
  parser_getHvgs.add_argument("batch_var", type=str)
  parser_getHvgs.add_argument("n_hvgs", type=int)
  parser_getHvgs.add_argument("-e", "--noambient", action='store_true')


  # parser for create_hvg_matrix()
  parser_readHvgs = subparsers.add_parser('create_hvg_matrix')
  parser_readHvgs.add_argument("qc_smpl_stats_f", type=str)
  parser_readHvgs.add_argument("hvg_paths_f", type=str)
  parser_readHvgs.add_argument("hvg_f", type=str)
  parser_readHvgs.add_argument("out_h5_f", type=str)
  parser_readHvgs.add_argument("demux_type", type=str)
  parser_readHvgs.add_argument("batch_var", type=str)

  # parser for create_doublets_matrix()
  parser_getDoublets = subparsers.add_parser('create_doublets_matrix')
  parser_getDoublets.add_argument("hvg_paths_f", type=str)
  parser_getDoublets.add_argument("hvg_f", type=str)
  parser_getDoublets.add_argument("qc_f", type=str)
  parser_getDoublets.add_argument("qc_smpl_stats_f", type=str)
  parser_getDoublets.add_argument("out_h5_f", type=str)
  parser_getDoublets.add_argument("run_var", type=str)
  parser_getDoublets.add_argument("demux_type", type=str)
  parser_getDoublets.add_argument("batch_var", type=str)
 
  args = parser.parse_args()

  if args.function_name == 'get_csr_counts':
    get_csr_counts( 
      args.hvg_paths_f, args.cell_filter_f, args.keep_var, args.keep_vals_str, args.smpl_stats_f, args.rowdata_f,
      args.run_var, args.batch_var, args.demux_type, args.chunksize, args.ncores
    )
  elif args.function_name == 'calculate_mean_var_for_chunk':
    calculate_mean_var_for_chunk(
      args.hvg_paths_f, args.rowdata_f, args.metadata_f, args.qc_smpl_stats_f,
      args.mean_var_f, args.chunk, args.hvg_method, args.batch_var, args.chunk_size,
      args.groupvar, args.group, args.ncores
    )
  elif args.function_name == 'calculate_estimated_vars':
    calculate_estimated_vars(
      args.estim_vars_f, args.hvg_method, args.batch_var, args.mean_var_merged_f
    )
  elif args.function_name == 'calculate_std_var_stats_for_sample':
    calculate_std_var_stats_for_sample(
      args.sample, args.batch_var, args.qc_smpl_stats_f, args.csr_f, args.rowdata_f,  args.std_var_stats_f
    )
  elif args.function_name == 'calculate_std_var_stats_for_chunk':
    calculate_std_var_stats_for_chunk(
      args.hvg_paths_f, args.rowdata_f, args.metadata_f, args.qc_smpl_stats_f,
      args.std_var_stats_f, args.estim_vars_f, args.chunk_num,
      args.hvg_method, args.batch_var, args.chunk_size, args.groupvar, args.group, args.ncores
    )
  elif args.function_name == 'calculate_hvgs':
    calculate_hvgs(
      args.std_var_stats_f, args.hvg_f, args.empty_gs_f, args.hvg_method, args.batch_var, 
      args.n_hvgs, args.noambient
    )
  elif args.function_name == 'create_hvg_matrix':
    create_hvg_matrix(
      args.qc_smpl_stats_f, args.hvg_paths_f, args.hvg_f, args.out_h5_f, args.demux_type, args.batch_var
    )
  elif args.function_name == 'create_doublets_matrix': 
    create_doublets_matrix(
      args.hvg_paths_f, args.hvg_f, args.qc_f, args.qc_smpl_stats_f, 
      args.out_h5_f, args.run_var, args.demux_type, args.batch_var
    )
  else:
    parser.print_help()
  
