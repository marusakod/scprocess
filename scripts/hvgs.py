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


def calculate_stats_for_sample(sample, qc_smpl_stats_f, csr_f, stats_f):

    qc_df = pd.read_csv(qc_smpl_stats_f, sep = '\t')
    bad_samples = qc_df.loc[qc_df['bad_sample'] == True, 'sample_id'].tolist()
   
    if sample in bad_samples:
        # save an empty file
        open(stats_f, 'w').close()
        return
    
    # read sparse matrix
    sparse_csr, features, _ = read_full_csr(csr_f)
    
    # calculate stats
    stats_df = calculate_feature_stats(sparse_csr, features)
    stats_df['sample_id'] = sample
    
    # save
    stats_df.to_csv(stats_f, sep='\t', index=False, compression='gzip')

    return



def calculate_stats_for_chunk(hvg_paths_f, rowdata_f, metadata_f, qc_smpl_stats_f, stats_f, chunk_num, hvg_method,
                              chunk_size=2000, group_var=None, group=None, n_cores = 8):
    
    # get total number of genes
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
        
        # check how many files were selected, make empty output if no files
        if len(files) == 0:
            print(f"No (good) samples found in group '{group}'. Writing an empty output file.")
            open(stats_f, 'w').close()
            return

    # calculate start and end rows for the chunk
    start_row = chunk_num * chunk_size
    end_row = min(start_row + chunk_size, total)

    # if start row exceeds total number of genes, write an empty file
    if start_row >= total:
        print(f"Start row {start_row} exceeds total number of genes ({total}). Writing an empty output file.")
        open(stats_f, 'w').close()
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

    # Calculate stats
    merged_chunk_stats = calculate_feature_stats(merged_chunk, features)

    # add group
    if hvg_method == 'group':
        merged_chunk_stats[group_var] = group

    # Save the stats to file
    merged_chunk_stats.to_csv(stats_f, sep='\t', index=False, compression='gzip')

    return



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest = "function_name", help = "Name of the function to run")
    
    # parser for get_csr_counts
    parser_makeCSR = subparsers.add_parser('get_csr_counts')
    parser_makeCSR.add_argument("paths", type=str)
    parser_makeCSR.add_argument("qc", type=str)
    parser_makeCSR.add_argument("qc_samples", type=str)
    parser_makeCSR.add_argument("rd", type=str)
    parser_makeCSR.add_argument("var", type=str)
    parser_makeCSR.add_argument("demux", type=str)
    
    # parser for calculate_stats_per_sample
    parser_sampleCalcs = subparsers.add_parser('calculate_stats_per_sample')
    parser_sampleCalcs.add_argument("sample", type=str)
    parser_sampleCalcs.add_argument("input", type=str)
    parser_sampleCalcs.add_argument("output", type=str)

    # parser for calculate_stats_per_chunk
    parser_chunkCalcs = subparsers.add_parser('calculate_stats_per_chunk')
    parser_chunkCalcs.add_argument("hvg_paths_f", type=str)
    parser_chunkCalcs.add_argument("rowdata_f", type=str)
    parser_chunkCalcs.add_argument("metadata_f", type=str)
    parser_chunkCalcs.add_argument("qc_smpl_stats_f", type=str)
    parser_chunkCalcs.add_argument("stats_f", type=str)
    parser_chunkCalcs.add_argument("chunk", type=int)
    parser_chunkCalcs.add_argument("method", type=str)
    parser_chunkCalcs.add_argument("chunk_size", type=int)
    parser_chunkCalcs.add_argument("-v", "--groupvar", type=str, required=False, default=None)
    parser_chunkCalcs.add_argument("-g", "--group", type=str, required=False, default= None)
    parser_chunkCalcs.add_argument("-n", "--ncores", type=int, required=False, default = 8)
   
    args = parser.parse_args()

    if args.function_name == 'get_csr_counts':
      get_csr_counts(args.paths, args.qc, args.qc_samples, args.rd, args.var, args.demux)
    elif args.function_name == 'calculate_stats_per_sample':
      calculate_stats_for_sample(args.sample, args.input, args.output)
    elif args.function_name == 'calculate_stats_per_chunk':
      calculate_stats_for_chunk(args.hvg_paths_f, args.rowdata_f, args.metadata_f, args.qc_smpl_stats_f,
                                args.stats_f, args.chunk, args.method, args.chunk_size, args.groupvar, args.group, args.ncores)
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

