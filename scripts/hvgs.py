import h5py
import numpy as np
import pandas as pd
import argparse
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


def get_one_csr_counts(run, hvg_df, keep_df, qc_sample_df, gene_ids, SAMPLE_VAR, DEMUX_TYPE, chunk_size):
  # get input (ambient) file and output files
  filt_counts_f = hvg_df.loc[hvg_df[SAMPLE_VAR] == run, "amb_filt_f"].values[0]
  out_fs        = hvg_df.loc[hvg_df[SAMPLE_VAR] == run, "chunked_f"].tolist()

  # get bad samples
  bad_samples   = qc_sample_df.loc[ qc_sample_df['bad_sample'] == True, 'sample_id'].tolist()

  # get cell ids for each sample
  if DEMUX_TYPE != "none":
    samples = hvg_df.loc[hvg_df[SAMPLE_VAR] == run, "sample_id"].tolist()
  else:
    samples = [run]

  cell_ids_dict = {}
  for s in samples:
    keep_cells = keep_df.loc[keep_df['sample_id'] == s, "cell_id"].tolist()
    cell_ids   = [bc.replace(run + ":", "", 1) for bc in keep_cells]
    cell_ids   = np.array(cell_ids)
    cell_ids_dict[s]  = cell_ids

  # open input file
  with h5py.File(filt_counts_f, 'r') as f:
    indptr    = f['matrix/indptr'][:]
    indices   = f['matrix/indices'][:]
    data      = f['matrix/data'][:]
    features  = f['matrix/features/name'][:]
    barcodes  = f['matrix/barcodes'][:]

    num_rows  = f['matrix/shape'][0]
    num_cols  = f['matrix/shape'][1]

  # make a csc sparse matrix
  sua_csc = csc_matrix((data, indices, indptr), shape=(num_rows, num_cols))

  for s, out_f in zip(samples, out_fs):
    if s in bad_samples:
      # return an empty file
      open(out_f, "w").close()
      return
    
    cell_ids = cell_ids_dict[s]
  
    # get indices of barcodes to keep
    keep_idx = np.where(np.isin(barcodes, cell_ids))[0]
    filt_bcs = barcodes[keep_idx]
    
    # subset matrix
    sua_csc_qc = sua_csc[:, keep_idx]

    # merge splices, unspliced, ambiguous
    csc, uniq_features = sum_SUA(sua_csc_qc, features)
    
    # get indices of genes to keep
    gene_ids    = np.array(gene_ids)
    gs_keep_idx   = np.where(np.isin(uniq_features, gene_ids))[0]
    uniq_features = uniq_features[gs_keep_idx]
    assert len(uniq_features) != 0, "No features selected"
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
      
    print(f"CSR matrix for {s} successfully saved to {out_f}.")


def get_csr_counts(hvg_paths_f, qc_f, qc_smpl_stats_f, rowdata_f, SAMPLE_VAR, DEMUX_TYPE, chunk_size=2000, n_cores = 8):
  # load up useful things
  hvg_paths_df  = pd.read_csv(hvg_paths_f)
  qc_sample_df  = pd.read_csv(qc_smpl_stats_f)
  
  # get QCed cells
  qc_df         = pd.read_csv(qc_f, sep = '\t')
  keep_df       = qc_df[qc_df["keep"] == True]

  # get gene details
  rows_df       = pd.read_csv(rowdata_f, sep = '\t')
  keep_ids      = rows_df['ensembl_id'].tolist()

  # define list of samples
  runs          = hvg_paths_df[SAMPLE_VAR].unique()

  # do some parallel calculations!
  with concurrent.futures.ThreadPoolExecutor(max_workers=n_cores) as executor:
    # i think a parallel job for each sample?
    futures = [executor.submit(get_one_csr_counts, run, hvg_paths_df, keep_df, 
      qc_sample_df, keep_ids, SAMPLE_VAR, DEMUX_TYPE, chunk_size) for run in runs]

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
  feat_stats_df = pd.DataFrame({
    'ensembl_id': features.astype(str),
    'mean': mean,
    'variance': var
  })

  # merge with nice labels for genes
  rows_df       = pd.read_csv(rowdata_f, sep = '\t')
  feat_stats_df = pd.merge(feat_stats_df, rows_df, on = 'ensembl_id')
  
  return feat_stats_df


# minor modifiction of scanpy _sum_and_sum_squares_clipped
def calculate_sum_and_sum_squares_cilipped(sparse_csr, estim_vars_df):
  # convert to column format, get values we need
  csc       = sparse_csr.tocsc()
  indices   = csc.indices
  data      = csc.data
  nnz       = csc.nnz

  # get some more values we need
  clip_val      = estim_vars_df['clip_val'].to_numpy()

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

  # add these values to df
  estim_vars_df['squared_counts_sum'] = sq_counts_sum  
  estim_vars_df['counts_sum']         = counts_sum

  return estim_vars_df


def calculate_estimated_vars(estim_vars_f, hvg_method, mean_var_merged_f = None, 
  mean_var_df = None, span: float = 0.3 ):
  # read mean var df if not specified
  if mean_var_df is None:
    assert mean_var_merged_f is not None, 'input file path or dataframe missing'
    stats_df  = pd.read_csv(mean_var_merged_f, sep = '\t')
  else:
    stats_df = mean_var_df

  # group on relevant variable
  group_var   = 'sample_id' if hvg_method == 'sample' else 'group'
  grouped     = stats_df.groupby(group_var)

  # go through each group
  all_results = []
  for name, group in grouped:
    # get values we need
    N           = group['n_cells'].unique()[0]  
    num_rows    = group.shape[0]
    mean        = group['mean'].to_numpy()
    var         = group['variance'].to_numpy()
     
    # this code is a minor modification of equivalent scanpy code:
    # https://github.com/scverse/scanpy/blob/1.11.1/src/scanpy/preprocessing/_highly_variable_genes.py#L503-L717
    not_const   = np.array(var) > 0
    y           = np.log10(var[not_const]) 
    x           = np.log10(mean[not_const])
    model       = safe_loess(x, y, span = span)

    # store trended variance values, also get std deviation
    est_var     = np.zeros(num_rows, dtype=np.float64)
    est_var[not_const] = model.outputs.fitted_values
    reg_std     = np.sqrt(10**est_var)

    # clip large values as in Seurat
    vmax        = np.sqrt(N)
    clip_val    = reg_std * vmax + mean
  
    # store values      
    group['reg_std']            = reg_std
    group['clip_val']           = clip_val
    group['estimated_variance'] = est_var
    all_results.append(group)
  
  # join into one nice dataframe
  final_df = pd.concat(all_results)

  # save or just pass back?
  if estim_vars_f is not None:
    final_df.to_csv(estim_vars_f,  sep='\t', index=False, compression='gzip', quoting=csv.QUOTE_NONE)
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
  hvg_method, chunk_size=2000, group_var=None, group=None):
  
  rows_df = pd.read_csv(rowdata_f, sep='\t')
  total   = rows_df.shape[0]

  # get input paths df
  hvg_paths_df = pd.read_csv(hvg_paths_f)
  
  # get list of good samples
  qc_df = pd.read_csv(qc_smpl_stats_f)
  good_samples = qc_df.loc[qc_df['bad_sample'] == False, 'sample_id'].tolist()

  # get input files for selectd samples
  if hvg_method == 'all':
    files = hvg_paths_df.loc[hvg_paths_df['sample_id'].isin(good_samples), 'chunked_f'].tolist()
  else:
    assert group_var is not None, "group_var must be defined."
    assert group is not None, "group must be defined."
    
    # select samples based on group
    meta = pd.read_csv(metadata_f)
    grp_samples = meta.loc[meta[group_var] == group, 'sample_id'].tolist()
    
    sel_samples = list(set(grp_samples) & set(good_samples))
    files = hvg_paths_df.loc[hvg_paths_df['sample_id'].isin(sel_samples), 'chunked_f'].tolist()

   # calculate start and end rows for the chunk
  start_row = chunk_num * chunk_size
  end_row = min(start_row + chunk_size, total)

  return files, total, start_row, end_row


def calculate_std_var_stats_for_sample(sample, qc_smpl_stats_f, csr_f, rowdata_f, std_var_stats_f):

  qc_df = pd.read_csv(qc_smpl_stats_f)
  bad_samples = qc_df.loc[qc_df['bad_sample'] == True, 'sample_id'].tolist()
   
  if sample in bad_samples:
    # save an empty file
    open(std_var_stats_f, 'w').close()
    return
  
  # read sparse matrix
  sparse_csr, features, _ = read_full_csr(csr_f)
  

  # calculate mean and variance
  stats_df = _calculate_feature_stats(sparse_csr, features, rowdata_f)
  stats_df['sample_id'] = sample
  # count cells

  n_cells = sparse_csr.shape[1]
  stats_df['n_cells'] = n_cells

  # calculate estimated variance
  estim_vars_df = calculate_estimated_vars(
    estim_vars_f    = None,
    mean_var_merged_f   = None,
    hvg_method      = 'sample',
    mean_var_df     = stats_df
    )
  
  # get stats for standardized variance
  if sample in bad_samples:
    # save an empty file
    open(std_var_stats_f, 'w').close()
    return
  
  # read sparse matrix
  sparse_csr, features, _ = read_full_csr(csr_f)

  # calculate stats
  stats_df = calculate_sum_and_sum_squares_cilipped(sparse_csr, estim_vars_df)

  # save
  stats_df.to_csv(std_var_stats_f, sep='\t', index=False, compression='gzip', quoting=csv.QUOTE_NONE)

  return


def calculate_mean_var_for_chunk(hvg_paths_f, rowdata_f, metadata_f, qc_smpl_stats_f, mean_var_f, chunk_num, hvg_method,
                chunk_size=2000, group_var=None, group=None, n_cores = 8):
  
  files, total, start_row, end_row = get_chunk_params(
    hvg_paths_f, rowdata_f, metadata_f, qc_smpl_stats_f, 
    chunk_num, hvg_method, chunk_size, group_var, group
    )
  
  if len(files) == 0:
    if hvg_method == 'group':
      print(f"No (good) samples found in group '{group}'. Writing an empty output file.")
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
  merged_chunk_stats = _calculate_feature_stats(merged_chunk, features, rowdata_f)
  merged_chunk_stats['n_cells'] = n_cells

  # add group
  merged_chunk_stats['group'] = group

  # save
  merged_chunk_stats.to_csv(mean_var_f, sep='\t', index=False, compression='gzip', quoting=csv.QUOTE_NONE)

  return


def calculate_std_var_stats_for_chunk(hvg_paths_f, rowdata_f, metadata_f, qc_smpl_stats_f, std_var_stats_f, estim_vars_f, chunk_num, hvg_method,
              chunk_size=2000, group_var=None, group=None, n_cores = 8):
  
  files, total, start_row, end_row = get_chunk_params(
    hvg_paths_f, rowdata_f, metadata_f, qc_smpl_stats_f,
    chunk_num, hvg_method, chunk_size, group_var, group
    )
  
  if len(files) == 0:
    if hvg_method == 'group':
      print(f"No (good) samples found in group '{group}'. Writing an empty output file.")
      open(std_var_stats_f, 'w').close()
    else:
      print("Error: no files selected for calculating hvg statistics")
    return
  
  if start_row >= total:
    print(f"Start row {start_row} exceeds total number of genes ({total}). Writing an empty output file.")
    open(std_var_stats_f, 'w').close()
    return
  
  # get estimated variances
  estim_vars_df = pd.read_csv(estim_vars_f, sep = '\t')
  
  # read chunks for multiple samples in parallel
  merged_chunk = None
  features = []
  
  read_csr_chunk_p = partial(read_csr_chunk, start_row=start_row, end_row=end_row)
  with concurrent.futures.ThreadPoolExecutor(max_workers=n_cores) as executor:
    results = list(executor.map(read_csr_chunk_p, files))
    
    for csr_chunk, chunk_feats, _ in results:
      if len(features) == 0:
        features = chunk_feats
      
      # merge chunks column-wise
      merged_chunk = csr_chunk if merged_chunk is None else hstack([merged_chunk, csr_chunk])

  # calculate stats 
  chunk_estim_vars = estim_vars_df[
  (estim_vars_df['group'] == group) & 
  (estim_vars_df['ensembl_id'].isin(features))
  ]

  merged_chunk_stats = calculate_sum_and_sum_squares_cilipped(merged_chunk, chunk_estim_vars)

  # save
  merged_chunk_stats.to_csv(std_var_stats_f, sep='\t', index=False, compression='gzip', quoting=csv.QUOTE_NONE)

  return


def _calculate_standardized_variance(df):
  return df.assign(
    variances_norm=lambda d: (
      1 / ((d['n_cells'] - 1) * np.square(d['reg_std']))
    ) * (
      (d['n_cells'] * np.square(d['mean']))
      + d['squared_counts_sum']
      - 2 * d['counts_sum'] * d['mean']
    )
  )


def _rank_genes(ranking_df, n_hvgs):
  # sort and rank genes
  norm_gene_vars = ranking_df['variances_norm'].to_numpy()
  
  ranked_norm_gene_vars = np.argsort(np.argsort(-norm_gene_vars)).astype(float)
  ranked_norm_gene_vars[ranked_norm_gene_vars >= n_hvgs] = np.nan
  ma_ranked = np.ma.masked_invalid(ranked_norm_gene_vars)

  return ranked_norm_gene_vars, norm_gene_vars, ma_ranked


def _process_single_group(stats_df, empty_gs, n_hvgs, exclude_ambient):
  # get genes to exclude
  exclude_gs = []
  if exclude_ambient:
    exclude_gs = empty_gs

  ranks_df   = stats_df[~stats_df['gene_id'].isin(exclude_gs)]  
  # rank genes
  ranked_norm_gene_vars, _, _ = _rank_genes(ranks_df, n_hvgs)

  # create output dataframe
  out_df = stats_df[['gene_id', 'variances_norm']].copy()
  out_df['highly_variable_nbatches'] = 1  
  out_df['highly_variable_rank']   = np.nan

  # assign rankings to nonexcluded genes
  out_df.loc[~out_df['gene_id'].isin(exclude_gs), 'highly_variable_rank'] = ranked_norm_gene_vars

  return out_df


def _process_multiple_groups(stats_df, group_var, empty_gs, n_hvgs,  exclude_ambient):
  
  all_ranked_norm_gene_vars = []
  norm_gene_vars_list = []

  for group_name, group_df in stats_df.groupby(group_var):
    ranks_df = _process_single_group(group_df, empty_gs, n_hvgs, exclude_ambient)
    ranked_norm_gene_vars = ranks_df['highly_variable_rank'].to_numpy()
    variances       = ranks_df['variances_norm'].to_numpy()
    all_ranked_norm_gene_vars.append(ranked_norm_gene_vars)
    norm_gene_vars_list.append(variances)
    
  # merge metrics across multiple groups
  all_ranked_norm_gene_vars = np.array(all_ranked_norm_gene_vars)
  num_batches_high_var = np.sum((all_ranked_norm_gene_vars < n_hvgs).astype(int), axis=0)
  median_ranked = np.ma.median(np.ma.masked_invalid(all_ranked_norm_gene_vars), axis=0).filled(np.nan)
  variances_norm = np.mean(norm_gene_vars_list, axis=0)

  out_df = stats_df[['gene_id']].drop_duplicates().copy()
  out_df['highly_variable_nbatches'] = num_batches_high_var
  out_df['highly_variable_rank'] = median_ranked
  out_df['variances_norm'] = variances_norm

  return out_df


# main function to calculate highly variable genes
def calculate_hvgs(std_var_stats_f, hvg_f, empty_gs_f, hvg_method, n_hvgs, exclude_ambient=True):
   
  stats_df  = pd.read_csv(std_var_stats_f, sep='\t')
  group_var = 'sample_id' if hvg_method == 'sample' else 'group'

  # get empty genes
  empty_dt = pd.read_csv(empty_gs_f, sep=',', compression='gzip')
  empty_gs = empty_dt.loc[empty_dt['is_ambient'], 'gene_id'].tolist()
  
  stats_df  = _calculate_standardized_variance(stats_df)

  if stats_df[group_var].nunique() == 1:
    hvg_df = _process_single_group(stats_df, empty_gs, n_hvgs, exclude_ambient)
  else:
    hvg_df = _process_multiple_groups(stats_df, group_var, empty_gs, n_hvgs, exclude_ambient)

  sort_cols = ["highly_variable_rank", "highly_variable_nbatches"]
  sort_ascending = [True, False]

  #  label highly variable genes
  sorted_index = (
    hvg_df[sort_cols]
    .sort_values(sort_cols, ascending=sort_ascending, na_position="last")
    .index
  )

  hvg_df["highly_variable"] = False
  hvg_df.loc[sorted_index[: int(n_hvgs)], "highly_variable"] = True

  hvg_df.to_csv(hvg_f, sep='\t', index=False)
  return hvg_df


# read top 2000 hvgs from each sample and save file
def read_top_genes(qc_smpl_stats_f, hvg_paths_f, hvg_f, out_h5_f, SAMPLE_VAR):
  # get bad samples
  qc_sample_df  = pd.read_csv(qc_smpl_stats_f)
  bad_samples   = qc_sample_df.loc[ qc_sample_df['bad_sample'] == True, 'sample_id'].tolist()

  # get all chunked files
  hvg_paths_df  = pd.read_csv(hvg_paths_f)
  chunked_fs    = hvg_paths_df['chunked_f'].tolist()
  runs          = hvg_paths_df[SAMPLE_VAR].tolist()

  # get all hvgs
  hvg_df        = pd.read_csv(hvg_f, sep='\t')
  hvg_ids       = hvg_df.loc[hvg_df['highly_variable'] == True, 'gene_id'].tolist()

  # extract ensembl ids (maybe keep ensembl is in the df earlier)
  hvg_ensembl   = []
  for gene in hvg_ids:
    parts = gene.rsplit('_', 1)
    hvg_ensembl.append(parts[-1])
  
  # initialise hvg matrix
  top_genes_mat = None
  all_barcodes  = []

  # open each file separately and extract highly variable genes
  for f, run in zip(chunked_fs, runs):
    if run in bad_samples:
      continue
    sample_csr, features, barcodes = read_full_csr(f)

    features = np.array(features, dtype=str)

    hvg_indices = [i for i, feature in enumerate(features) if feature in hvg_ensembl]
    csr_chunk = sample_csr[hvg_indices, :]

    barcodes = barcodes.astype('<U21')
    barcodes = [f"{run}:{bc}" for bc in barcodes]  
    barcodes = np.array(barcodes)

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


def create_doublets_matrix(hvg_paths_f, hvg_f, qc_f, qc_smpl_stats_f, out_h5_f, SAMPLE_VAR):
  # get all hvgs
  hvg_df  = pd.read_csv(hvg_f, sep='\t')
  hvg_ids = hvg_df.loc[hvg_df['highly_variable'] == True, 'gene_id'].tolist()
  
  # get qc file with all cells
  qc_df   = pd.read_csv(qc_f, sep = '\t')

  # subset to doublets
  dbl_df  = qc_df[qc_df["dbl_class"] == "doublet"]

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
  hvg_paths_df  = pd.read_csv(hvg_paths_f)
  qc_sample_df  = pd.read_csv(qc_smpl_stats_f)
  good_samples  = qc_sample_df.loc[ qc_sample_df['bad_sample'] == False, "sample_id"].tolist()
  keep_runs     = hvg_paths_df.loc[hvg_paths_df['sample_id'].isin(good_samples), SAMPLE_VAR].unique().tolist()

  for run in keep_runs:
    # get doublets for this run
    dbls     = dbl_df.loc[dbl_df[SAMPLE_VAR] == run, "cell_id"].tolist()
    dbl_ids  = np.array(dbls)

    # get ambient output file
    filt_counts_f = hvg_paths_df.loc[hvg_paths_df[SAMPLE_VAR] == run, "amb_filt_f"].values[0]
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
  parser_makeCSR.add_argument("qc_f", type=str)
  parser_makeCSR.add_argument("qc_smpl_stats_f", type=str)
  parser_makeCSR.add_argument("rowdata_f", type=str)
  parser_makeCSR.add_argument("sample_var", type=str)
  parser_makeCSR.add_argument("demux_type", type=str)
  parser_makeCSR.add_argument("-s", "--size", required=False, default= 2000, type = int)
  parser_makeCSR.add_argument("-n", "--ncores", required=False, default= 8, type = int)
  

  # parser for calculate_mean_var_for_chunk
  parser_chunkCalcs = subparsers.add_parser('calculate_mean_var_for_chunk')
  parser_chunkCalcs.add_argument("hvg_paths_f", type=str)
  parser_chunkCalcs.add_argument("rowdata_f", type=str)
  parser_chunkCalcs.add_argument("metadata_f", type=str)
  parser_chunkCalcs.add_argument("qc_smpl_stats_f", type=str)
  parser_chunkCalcs.add_argument("mean_var_f", type=str)
  parser_chunkCalcs.add_argument("chunk", type=int)
  parser_chunkCalcs.add_argument("method", type=str)
  parser_chunkCalcs.add_argument("chunk_size", type=int)
  parser_chunkCalcs.add_argument("-v", "--groupvar", type=str, required=False, default=None)
  parser_chunkCalcs.add_argument("-g", "--group", type=str, required=False, default= None)
  parser_chunkCalcs.add_argument("-n", "--ncores", type=int, required=False, default = 8)

  # parser for calculate_estimated_vars
  parser_varEstim = subparsers.add_parser('calculate_estimated_vars')
  parser_varEstim.add_argument("estim_vars_f", type=str)
  parser_varEstim.add_argument("hvg_method", type=str)
  parser_varEstim.add_argument("mean_var_merged_f", type=str)

  # parser for calculate_std_var_stats_for_sample
  parser_sampleSssc = subparsers.add_parser('calculate_std_var_stats_for_sample')
  parser_sampleSssc.add_argument("sample", type=str)
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
  parser_chunkSssc.add_argument("-s", "--size", type=int, required=False, default=2000)
  parser_chunkSssc.add_argument("-v", "--groupvar", type=str, required=False, default=None)
  parser_chunkSssc.add_argument("-g", "--group", type=str, required=False, default= None)
  parser_chunkSssc.add_argument("-n", "--ncores", type=int, required=False, default = 8)

  
  # parser for calculate_hvgs(std_var_stats_f, hvg_f, empty_gs_fs, hvg_method, n_hvgs, exclude_empty=True
  parser_getHvgs = subparsers.add_parser('calculate_hvgs')
  parser_getHvgs.add_argument("std_var_stats_f", type=str)
  parser_getHvgs.add_argument("hvg_f", type=str)
  parser_getHvgs.add_argument("empty_gs_f", type=str)
  parser_getHvgs.add_argument("hvg_method", type=str)
  parser_getHvgs.add_argument("n_hvgs", type=int)
  parser_getHvgs.add_argument("-e", "--noambient", action='store_true')


  # parser for read_top_genes()
  parser_readHvgs = subparsers.add_parser('read_top_genes')
  parser_readHvgs.add_argument("qc_smpl_stats_f", type=str)
  parser_readHvgs.add_argument("hvg_paths_f", type=str)
  parser_readHvgs.add_argument("hvg_f", type=str)
  parser_readHvgs.add_argument("out_h5_f", type=str)
  parser_readHvgs.add_argument("sample_var", type=str)

  # parser for create_doublets_matrix()
  parser_getDoublets = subparsers.add_parser('create_doublets_matrix')
  parser_getDoublets.add_argument("hvg_paths_f", type=str)
  parser_getDoublets.add_argument("hvg_f", type=str)
  parser_getDoublets.add_argument("qc_f", type=str)
  parser_getDoublets.add_argument("qc_smpl_stats_f", type=str)
  parser_getDoublets.add_argument("out_h5_f", type=str)
  parser_getDoublets.add_argument("sample_var", type=str)

  args = parser.parse_args()

  if args.function_name == 'get_csr_counts':
    get_csr_counts(
      args.hvg_paths_f, args.qc_f, args.qc_smpl_stats_f, args.rowdata_f,
      args.sample_var, args.demux_type, args.size, args.ncores
    )
  elif args.function_name == 'calculate_mean_var_for_chunk':
    calculate_mean_var_for_chunk(
      args.hvg_paths_f, args.rowdata_f, args.metadata_f, args.qc_smpl_stats_f,
      args.mean_var_f, args.chunk, args.method, args.chunk_size,
      args.groupvar, args.group, args.ncores
    )
  elif args.function_name == 'calculate_estimated_vars':
    calculate_estimated_vars(
      args.estim_vars_f, args.hvg_method, args.mean_var_merged_f
    )
  elif args.function_name == 'calculate_std_var_stats_for_sample':
    calculate_std_var_stats_for_sample(
      args.sample, args.qc_smpl_stats_f, args.csr_f, args.rowdata_f,  args.std_var_stats_f
    )
  elif args.function_name == 'calculate_std_var_stats_for_chunk':
    calculate_std_var_stats_for_chunk(
      args.hvg_paths_f, args.rowdata_f, args.metadata_f, args.qc_smpl_stats_f,
      args.std_var_stats_f, args.estim_vars_f, args.chunk_num,
      args.hvg_method, args.size, args.groupvar, args.group, args.ncores
    )
  elif args.function_name == 'calculate_hvgs':
    calculate_hvgs(
      args.std_var_stats_f, args.hvg_f, args.empty_gs_f, args.hvg_method, 
      args.n_hvgs, args.noambient
    )
  elif args.function_name == 'read_top_genes':
    read_top_genes(
      args.qc_smpl_stats_f, args.hvg_paths_f, args.hvg_f, args.out_h5_f, args.sample_var
    )
  elif args.function_name == 'create_doublets_matrix': 
    create_doublets_matrix(
      args.hvg_paths_f, args.hvg_f, args.qc_f, args.qc_smpl_stats_f, 
      args.out_h5_f, args.sample_var
    )
  else:
    parser.print_help()
  
