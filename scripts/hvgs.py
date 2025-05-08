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


def sum_SUA(sua_mat, row_names):
    
    row_names = row_names.astype(str)
    types = ['_S$', '_U$', '_A$']
    mats = []

    for t in types:
        indices = [i for i, name in enumerate(row_names) if re.search(t,name)]
        mat = sua_mat[indices, :]
        mats.append(mat)

    genes = [re.sub(r'_[SUA]$', '', name) for name in row_names]
    uniq_genes = list(dict.fromkeys(genes))
    uniq_genes = np.array(uniq_genes)

    mats_sum = sum(mats)

    return mats_sum, uniq_genes


def get_one_csr_counts(run, hvg_df, qc_df, qc_sample_df, gene_ids, SAMPLE_VAR, DEMUX_TYPE, chunk_size):
    
    # get input (ambient) file and output files
    filt_counts_f = hvg_df.loc[hvg_df[SAMPLE_VAR] == run, "amb_filt_f"].values[0]
    out_fs = hvg_df.loc[hvg_df[SAMPLE_VAR] == run, "chunked_f"].tolist()

    # get cell ids for each sample
    if DEMUX_TYPE != "":
        samples = hvg_df.loc[hvg_df[SAMPLE_VAR] == run, "sample_id"].tolist()
    else:
        samples = [run]

    cell_ids_dict = {}

    for s in samples:
        smpl_cells = qc_df.loc[qc_df['sample_id'] == s, "cell_id"].tolist()
        cell_ids = [bc.replace(run + ":", "", 1) for bc in smpl_cells]
        cell_ids = np.array(cell_ids)
        cell_ids_dict[s] = cell_ids

    # open input file
    with h5py.File(filt_counts_f, 'r') as f:
        indptr   = f['matrix/indptr'][:]
        indices  = f['matrix/indices'][:]
        data     = f['matrix/data'][:]
        features = f['matrix/features/name'][:]
        barcodes = f['matrix/barcodes'][:]

        num_rows = f['matrix/shape'][0]
        num_cols = f['matrix/shape'][1]

    # make a csc sparse matrix
    sua_csc = csc_matrix((data, indices, indptr), shape=(num_rows, num_cols))

    # get bad samples
    bad_samples  = qc_sample_df.loc[ qc_sample_df['bad_sample'] == True, 'sample_id'].tolist()

    for s, out_f in zip(samples, out_fs):

        if s in bad_samples:
            # return an empty file
            open(out_f, "w").close()
            return
        
        cell_ids = cell_ids_dict[s]

        # get indices of barcodes to keep
        bcs_keep_idx = np.where(np.isin(barcodes, cell_ids))[0]
    
        # subset matrix
        sua_csc_qc = sua_csc[:, bcs_keep_idx]

        # merge splices, unspliced, ambiguous
        csc, uniq_features = sum_SUA(sua_csc_qc, features)
        
        # get indices of genes to keep
        gene_ids      = np.array(gene_ids)
        gs_keep_idx   = np.where(np.isin(uniq_features, gene_ids))[0]
        uniq_features = uniq_features[gs_keep_idx]
        assert len(uniq_features) != 0, "No features selected"
        csc           = csc[gs_keep_idx, :]

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
            f.create_dataset('matrix/barcodes', data=np.array(cell_ids, dtype='S'), compression='gzip')

        print(f"CSR matrix for {s} successfully saved to {out_f}.")



def get_csr_counts(hvg_paths_f, qc_f, qc_smpl_stats_f, rowdata_f, SAMPLE_VAR, DEMUX_TYPE, chunk_size=2000, n_cores = 8):
    hvg_df  = pd.read_csv(hvg_paths_f)
    qc_df   = pd.read_csv(qc_f, sep = '\t')
    qc_df   = qc_df[qc_df["keep"] == True]
    rows_df = pd.read_csv(rowdata_f, sep = '\t')
    keep_ids= rows_df['ensembl_id'].tolist()

    qc_sample_df = pd.read_csv(qc_smpl_stats_f, sep = '\t')

    runs = hvg_df[SAMPLE_VAR].unique()

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_cores) as executor:
        futures = [executor.submit(get_one_csr_counts, run, hvg_df, qc_df, qc_sample_df, keep_ids, SAMPLE_VAR, DEMUX_TYPE, chunk_size) for run in runs]

        for future in concurrent.futures.as_completed(futures):
            future.result()


def read_csr_chunk(file, start_row, end_row):
    with h5py.File(file, 'r') as f:
        num_cols = f['matrix/shape'][1]
        indptr = f['matrix/indptr'][start_row:end_row+1]
        start_pos = indptr[0]
        end_pos = indptr[-1]
        data    = f['matrix/data'][start_pos:end_pos]
        indices = f['matrix/indices'][start_pos:end_pos]

        # get features and barcodes for chunk
        features= f['matrix/features/name'][start_row:end_row].astype(str)
        barcodes= f['matrix/barcodes'][:].astype(str)

        # adjust indptr to start at zero for the chunk
        indptr_chunk = indptr - start_pos

    # csr for the chunk
    sparse_chunk = csr_matrix((data, indices, indptr_chunk), shape=(end_row - start_row, num_cols))
    
    return sparse_chunk, features, barcodes


def read_full_csr(file): 
    with h5py.File(file, 'r') as f:
        indptr   = f['matrix/indptr'][:]  
        indices  = f['matrix/indices'][:]  
        data     = f['matrix/data'][:]  
        features = f['matrix/features/name'][:]  
        barcodes = f['matrix/barcodes'][:]  

        num_rows = f['matrix/shape'][0]  
        num_cols = f['matrix/shape'][1]  

        sparse_csr = csr_matrix((data, indices, indptr), shape=(num_rows, num_cols))

    return sparse_csr, features, barcodes



def calculate_feature_stats(sparse_csr, features):

    mean, var = _get_mean_var(sparse_csr, axis=1 )
    
    feat_stats_df = pd.DataFrame({
        'feature_name': features.astype(str),
        'mean': mean,
        'variance': var
    })
    
    return feat_stats_df



# minor modifiction of scanpy _sum_and_sum_squares_clipped
def calculate_sum_and_sum_squares_cilipped(sparse_csr, estim_vars_df):
    
    n_rows  = sparse_csr.shape[0]
    csc     = sparse_csr.tocsc()

    indices = csc.indices
    data    = csc.data
    clip_val= estim_vars_df['clip_val'].to_numpy()
    nnz     = csc.nnz

    squared_counts_sum = np.zeros(n_rows, dtype=np.float64)
    counts_sum = np.zeros(n_rows, dtype=np.float64)

    for i in numba.prange(nnz):
        idx = indices[i]
        element = min(np.float64(data[i]), clip_val[idx])
        squared_counts_sum[idx] += element**2
        counts_sum[idx] += element


    estim_vars_df['squared_counts_sum'] = squared_counts_sum  
    estim_vars_df['counts_sum']         = counts_sum

    return estim_vars_df


def calculate_estimated_vars(estim_vars_f, mean_var_merged_f, hvg_method, mean_var_df = None, span: float = 0.3 ):

    if mean_var_df is None:
        # read mean var df if not specified
        assert mean_var_merged_f is not None, 'input file path or dataframe missing'
        stats_df  = pd.read_csv(mean_var_merged_f, sep = '\t')
    else:
        stats_df = mean_var_df

    group_var = 'sample_id' if hvg_method == 'sample' else 'group'
    
    grouped = stats_df.groupby(group_var)
    all_results = []

    for name, group in grouped:
        N = group['n_cells'].unique()[0]
            
        num_rows = group.shape[0]
         
        # this code is a minor modification from scanpy (https://github.com/scverse/scanpy/blob/1.11.1/src/scanpy/preprocessing/_highly_variable_genes.py#L503-L717)
        mean = group['mean'].to_numpy()
        var  = group['variance'].to_numpy()
        not_const = np.array(var) > 0
        estimat_var = np.zeros(num_rows, dtype=np.float64)

        y = np.log10(var[not_const]) 
        x = np.log10(mean[not_const])
        model = loess(x, y, span=span, degree=2)
        model.fit()
        estimat_var[not_const] = model.outputs.fitted_values
        reg_std = np.sqrt(10**estimat_var)

        # clip large values as in Seurat
        vmax = np.sqrt(N)
        clip_val = reg_std * vmax + mean
            
        group['reg_std']  = reg_std
        group['clip_val'] = clip_val
        group['estimated_variance'] = estimat_var
          
        all_results.append(group)
    
    final_df = pd.concat(all_results)

    if estim_vars_f is not None:
        final_df.to_csv(estim_vars_f,  sep='\t', index=False, compression='gzip', quoting=csv.QUOTE_NONE)
        return
    else:
        return final_df


def get_chunk_params(hvg_paths_f, rowdata_f, metadata_f, qc_smpl_stats_f, chunk_num, hvg_method,
                              chunk_size=2000, group_var=None, group=None):
    
    rows_df = pd.read_csv(rowdata_f, sep='\t')
    total   = rows_df.shape[0]

    # get input paths df
    hvg_paths_df = pd.read_csv(hvg_paths_f)
    
    # get list of good samples
    qc_df = pd.read_csv(qc_smpl_stats_f, sep = '\t')
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



def calculate_std_var_stats_for_sample(sample, qc_smpl_stats_f, csr_f, std_var_stats_f):

    qc_df = pd.read_csv(qc_smpl_stats_f, sep = '\t')
    bad_samples = qc_df.loc[qc_df['bad_sample'] == True, 'sample_id'].tolist()
   
    if sample in bad_samples:
        # save an empty file
        open(std_var_stats_f, 'w').close()
        return
    
    # read sparse matrix
    sparse_csr, features, _ = read_full_csr(csr_f)
    
    # calculate mean and variance
    stats_df = calculate_feature_stats(sparse_csr, features)
    stats_df['sample_id'] = sample

    # calculate estimated variance
    estim_vars_df = calculate_estimated_vars(
        estim_vars_f        = None,
        mean_var_merged_f   = None,
        hvg_method          = 'sample',
        mean_var_df         = stats_df
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
    merged_chunk_stats = calculate_feature_stats(merged_chunk, features)
    merged_chunk_stats['n_cells'] = n_cells

    # add group
    merged_chunk_stats['group'] = group

    # save
    merged_chunk_stats.to_csv(mean_var_f, sep='\t', index=False, compression='gzip', quoting=csv.QUOTE_NONE)

    return



# def calculate_std_var_stats_for_sample(sample, qc_smpl_stats_f, estim_vars_f, csr_f, std_var_stats_f):

#     qc_df = pd.read_csv(qc_smpl_stats_f, sep = '\t')
#     bad_samples = qc_df.loc[qc_df['bad_sample'] == True, 'sample_id'].tolist()
   
#     if sample in bad_samples:
#         # save an empty file
#         open(std_var_stats_f, 'w').close()
#         return
    
#     # read sparse matrix
#     sparse_csr, features, _ = read_full_csr(csr_f)

#     # get df with estimated variances for sample
#     estim_vars_df = pd.read_csv(estim_vars_f, sep = '\t')
#     estim_vars_smpl = estim_vars_df[estim_vars_df['sample_id'] == sample,]
    
#     # calculate stats
#     stats_df = calculate_sum_and_sum_squares_cilipped(sparse_csr, estim_vars_smpl)

#     # save
#     stats_df.to_csv(std_var_stats_f, sep='\t', index=False, compression='gzip', quoting=csv.QUOTE_NONE)

#     return



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
    (estim_vars_df['feature_name'].isin(features))
    ]

    merged_chunk_stats = calculate_sum_and_sum_squares_cilipped(merged_chunk, chunk_estim_vars)

    # save
    merged_chunk_stats.to_csv(std_var_stats_f, sep='\t', index=False, compression='gzip', quoting=csv.QUOTE_NONE)

    return


def _calculate_standardized_variance(df):
    
    return df.assign(
        norm_gene_var=lambda d: (
            1 / ((d['n_cells'] - 1) * np.square(d['reg_std']))
        ) * (
            (d['n_cells'] * np.square(d['mean']))
            + d['squared_counts_sum']
            - 2 * d['counts_sum'] * d['mean']
        )
    )


def _get_excluded_genes(empty_gs_fs, hvg_method, stats_df, group_var):

    # get list or dictionary of empty genes to exclude
    if hvg_method in ['sample', 'all']:
        empty_dt = pd.read_csv(empty_gs_fs, sep='\t')
        return empty_dt.loc[empty_dt['is_ambient'], 'gene_id'].tolist()
    else:
        empty_fs_ls = empty_gs_fs.split(',')
        uniq_groups = stats_df[group_var].nunique()
        assert uniq_groups == len(empty_fs_ls), "Number of empty gene files must match the number of groups"
        empty_dts = [pd.read_csv(f, sep='\t') for f in empty_fs_ls]

        exc_gs = {}
        for df in empty_dts:
            key = df['group'].unique()
            assert len(key) == 1, "More than one group in dataframe with empty genes"
            exc_gs[key[0]] = df.loc[df['is_ambient'], 'gene_id'].tolist()
        
        return exc_gs


def _rank_genes(ranking_df, n_hvgs):
    # sort and rank genes
    norm_gene_vars = ranking_df['norm_gene_var'].to_numpy()
    ranked_norm_gene_vars = np.argsort(np.argsort(-norm_gene_vars))
    ranked_norm_gene_vars[ranked_norm_gene_vars >= n_hvgs] = np.nan
    ma_ranked = np.ma.masked_invalid(ranked_norm_gene_vars)
    
    return ranked_norm_gene_vars, np.mean(norm_gene_vars), ma_ranked


def _process_single_group(stats_df, genes_to_exclude, n_hvgs):
    
    ranks_df = stats_df[~stats_df['gene'].isin(genes_to_exclude)]
    ranked_norm_gene_vars, variances_norm, _ = _rank_genes(ranks_df, n_hvgs)

    out_df                              = stats_df[['gene']].copy()
    out_df['gene_name']                 = out_df['gene']
    out_df['highly_variable_nbatches']  = 1  
    out_df['highly_variable_rank']      = ranked_norm_gene_vars
    out_df['variances_norm']            = variances_norm

    return out_df


def _process_multiple_groups(stats_df, group_var, exc_gs, n_hvgs, hvg_method):

    all_ranked_norm_gene_vars = []
    norm_gene_vars_list = []

    for group_name, group_df in stats_df.groupby(group_var):
        genes_to_exclude = exc_gs if hvg_method in ['sample', 'all'] else exc_gs.get(group_name, [])
        ranks_df = group_df[~group_df['gene'].isin(genes_to_exclude)]
        ranked_norm_gene_vars, variances, ma_ranked = _rank_genes(ranks_df, n_hvgs)

        all_ranked_norm_gene_vars.append(ranked_norm_gene_vars)
        norm_gene_vars_list.append(variances)

    # merge metrics across multiple groups 
    all_ranked_norm_gene_vars = np.array(all_ranked_norm_gene_vars)
    num_batches_high_var      = np.sum((all_ranked_norm_gene_vars < n_hvgs).astype(int), axis=0)
    median_ranked             = np.ma.median(np.ma.masked_invalid(all_ranked_norm_gene_vars), axis=0).filled(np.nan)
    variances_norm            = np.mean(norm_gene_vars_list, axis=0)

    out_df                             = stats_df[['gene']].drop_duplicates().copy()
    out_df['gene_name']                = out_df['gene']
    out_df['highly_variable_nbatches'] = num_batches_high_var
    out_df['highly_variable_rank']     = median_ranked
    out_df['variances_norm']        = variances_norm

    return out_df


# main function to calculate highly variable genes
def calculate_hvgs(std_var_stats_f, hvg_f, empty_gs_fs, hvg_method, n_hvgs, exclude_ambient=True):
   
    stats_df  = pd.read_csv(std_var_stats_f, sep='\t')
    group_var = 'sample_id' if hvg_method == 'sample' else 'group'

    stats_df  = _calculate_standardized_variance(stats_df)
    exc_gs    = _get_excluded_genes(empty_gs_fs, hvg_method, stats_df, group_var) if exclude_ambient else None

    if stats_df[group_var].nunique() == 1:
        hvg_df = _process_single_group(stats_df, exc_gs, n_hvgs)
    else:
        hvg_df = _process_multiple_groups(stats_df, group_var, exc_gs, n_hvgs, hvg_method)

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
    parser_varEstim.add_argument("mean_var_merged_f", type=str)
    parser_varEstim.add_argument("hvg_method", type=str)

    # parser for calculate_std_var_stats_for_sample
    parser_sampleSssc = subparsers.add_parser('calculate_std_var_stats_for_sample')
    parser_sampleSssc.add_argument("sample", type=str)
    parser_sampleSssc.add_argument("qc_smpl_stats_f", type=str)
    parser_sampleSssc.add_argument("csr_f", type=str)
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
    parser_getHvgs.add_argument("-e", "--noambient",required=False, default=True)
   
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
            args.estim_vars_f, args.mean_var_merged_f, args.hvg_method
        )
    elif args.function_name == 'calculate_std_var_stats_for_sample':
        calculate_std_var_stats_for_sample(
            args.sample, args.qc_smpl_stats_f, args.csr_f, args.std_var_stats_f
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
            args.n_hvgs, args.noempty
        )
    else:
        parser.print_help()
    



# # read top 2000 hvgs from each sample and save file
# def read_top_genes(files_dt, top_indices):

#     files = pd.read_csv(files_dt)
#     fs_ls = files['amb_out_f'].tolist()
    
#     top_genes_mx = None
    
#     for f in fs_ls:
#         with h5py.File(f, 'r') as f:
#             num_cols =  f['matrix/shape'][1]
#             indptr = f['matrix/indptr'][:]
#             indices = f['matrix/indices'][:]
#             data = f['matrix/data'][:]
#             indptr_chunk = indptr[top_indices[0]:top_indices[-1]+1] - indptr[top_indices[0]]
#             data_chunk = data[indptr_chunk[0]: indptr_chunk[-1]]
#             indices_chunk = indices[indptr_chunk[0]:indptr_chunk[-1]]
#             sparse_chunk = csr_matrix((data_chunk, indices_chunk, indptr_chunk - indptr_chunk[0]), shape = (len(top_indices, num_cols)))
#             top_genes_mx = sparse_chunk if top_genes_mx is None else csr_matrix.hstack([top_genes_mx, sparse_chunk])

#     return top_genes_mx

