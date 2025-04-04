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



def process_run(run, hvg_df, qc_df, SAMPLE_VAR, DEMUX_TYPE, chunk_size):
    
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
        indptr = f['matrix/indptr'][:]
        indices = f['matrix/indices'][:]
        data = f['matrix/data'][:]
        features = f['matrix/features/name'][:]
        barcodes = f['matrix/barcodes'][:]

        num_rows = f['matrix/shape'][0]
        num_cols = f['matrix/shape'][1]

    # make a csc sparse matrix
    sua_csc = csc_matrix((data, indices, indptr), shape=(num_rows, num_cols))

    for s, out_f in zip(samples, out_fs):
        
        cell_ids = cell_ids_dict[s]

        # geet indices of barcodes to keep
        keep_idx = np.where(np.isin(barcodes, cell_ids))[0]

        # subset matrix
        sua_csc_qc = sua_csc[:, keep_idx]

        # merge splices, unspliced, ambiguous
        csc, uniq_features = sum_SUA(sua_csc_qc, features)

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


def get_csr_counts(hvg_paths_f, qc_f, SAMPLE_VAR, DEMUX_TYPE, chunk_size=2000, n_cores = 8):
    hvg_df = pd.read_csv(hvg_paths_f)
    qc_df  = pd.read_csv(qc_f, sep = '\t')
    qc_df = qc_df[qc_df["keep"] == True]

    runs = hvg_df[SAMPLE_VAR].unique()

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_cores) as executor:
        futures = [executor.submit(process_run, run, hvg_df, qc_df, SAMPLE_VAR, DEMUX_TYPE, chunk_size) for run in runs]

        for future in concurrent.futures.as_completed(futures):
            future.result()



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


def calculate_stats_for_sample(sample, csr_f, stats_f):
    
    # read sparse matrix
    sparse_csr, features, _ = read_full_csr(csr_f)
    
    # calculate stats
    stats_df = calculate_feature_stats(sparse_csr, features)
    stats_df['sample_id'] = sample
    
    # save
    stats_df.to_csv(stats_f, sep='\t', index=False, compression='gzip')

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest = "function_name", help = "Name of the function to run")
    
    # parser for get_csr_counts
    parser_makeCSR = subparsers.add_parser('get_csr_counts')
    parser_makeCSR.add_argument("paths", type=str)
    parser_makeCSR.add_argument("qc", type=str)
    parser_makeCSR.add_argument("var", type=str)
    parser_makeCSR.add_argument("demux", type=str)
    
    # parser for calculate_stats_per_sample
    parser_sampleCalcs = subparsers.add_parser('calculate_stats_per_sample')
    parser_sampleCalcs.add_argument("sample", type=str)
    parser_sampleCalcs.add_argument("input", type=str)
    parser_sampleCalcs.add_argument("output", type=str)
   
    args = parser.parse_args()

    if args.function_name == 'get_csr_counts':
      get_csr_counts(args.paths, args.qc, args.var, args.demux)
    elif args.function_name == 'calculate_stats_per_sample':
      calculate_stats_for_sample(args.sample, args.input, args.output)
    else:
     parser.print_help()
    
    



# def read_csr_chunk(file, start_row, end_row, as_adata = True):
#     with h5py.File(file, 'r') as f:
#         num_cols = f['matrix/shape'][1]
        
#         # get relevant part of indptr
#         indptr = f['matrix/indptr'][start_row:end_row+1]
        
#         # get start and end positions in the data and indices arrays
#         start_pos = indptr[0]
#         end_pos = indptr[-1]
        
#         # read data and indices for the chunk
#         data    = f['matrix/data'][start_pos:end_pos]
#         indices = f['matrix/indices'][start_pos:end_pos]

#         # get features and barcodes for chunk
#         features= f['matrix/features/name'][start_row:end_row].astype(str)
#         barcodes= f['matrix/barcodes'][:].astype(str)

#         # adjust indptr to start at zero for the chunk
#         indptr_chunk = indptr - start_pos

#     # csr for the chunk
#     sparse_chunk = csr_matrix((data, indices, indptr_chunk), shape=(end_row - start_row, num_cols))
    
#     if as_adata:
#         adata_chunk = sc.AnnData(X = sparse_chunk.T, var = pd.DataFrame(index = features), obs = pd.DataFrame(index=barcodes))
#         return(adata_chunk)
#     else:
#         return sparse_chunk, features, barcodes




# def process_chunk(files, samples, start_row, end_row):
#     # initialize chunk
#     merged_chunk = None
#     features = []
#     barcodes = []

#     for f, s in zip(files, samples):
#         # read sparse chunk
#         csr_chunk, chunk_feats, chunk_bcs = read_csr_chunk(f, start_row, end_row, as_adata = False)
#         chunk_bcs = chunk_bcs.astype(str)
#         chunk_bcs = [f"{s}:{barcode}" for barcode in chunk_bcs]

#         features = chunk_feats
#         barcodes.extend(chunk_bcs)

#         # merge chunks columnwise
#         merged_chunk = csr_chunk if merged_chunk is None else hstack([merged_chunk, csr_chunk])
        
#     # convert to andata
#     merged_chunk_adata = sc.AnnData(X=merged_chunk.T, var=pd.DataFrame(index=features), obs=pd.DataFrame(index=barcodes))

#     # calculate stats
#     merged_chunk_stats = calc_stats(merged_chunk_adata)

#     return(merged_chunk_stats)




# def do_cross_smpl_calcs_parallel(files_dt, sample_var = 'sample_id', chunk_size = 2000, n_cores = 4):
#     files = pd.read_csv(files_dt)
#     fs_ls = files['norm_f'].tolist()
#     samples = files[sample_var].tolist()

#     # get number of genes
#     with h5py.File(fs_ls[0], 'r') as f:
#         num_rows = f['matrix/shape'][0]
    
#     # define args for each task
#     tasks = [(fs_ls, samples, start_row, min(start_row + chunk_size, num_rows)) for start_row in range(0, num_rows, chunk_size)]

#     with Pool(processes=n_cores) as pool:
#         all_stats = pool.starmap(process_chunk, tasks)

#     # concatenate results from all chunks
#     final_stats = pd.concat(all_stats, ignore_index=True)

#     return final_stats





# def do_cross_smpl_calcs(files_dt, sample_var='sample_id', chunk_size=2000):
#     files = pd.read_csv(files_dt)
#     fs_ls = files['norm_f'].tolist()
#     samples = files[sample_var].tolist()
    
#     # Get number of rows (genes)
#     with h5py.File(fs_ls[0], 'r') as f:
#         num_rows = f['matrix/shape'][0]
#     all_stats = []  # Initialize empty list to collect stats for all chunks
#     for start_row in range(0, num_rows, chunk_size):
#         end_row = min(start_row + chunk_size, num_rows)
#         merged_chunk = None
#         features = []  
#         barcodes = []  
                
#         for f, s in zip(fs_ls, samples):
#             # Read the chunk as a sparse matrix
#             sparse_chunk, chunk_features, chunk_barcodes = read_csr_chunk(f, start_row, end_row, as_adata=False)
            
#             chunk_barcodes = chunk_barcodes.astype(str)
#             chunk_barcodes = [f"{s}:{barcode}" for barcode in chunk_barcodes]

#             # Store feature names and barcodes
#             features = chunk_features  # Should be the same for each file
#             barcodes.extend(chunk_barcodes)
            
#             # Merge chunks as sparse matrices
#             merged_chunk = sparse_chunk if merged_chunk is None else hstack([merged_chunk, sparse_chunk])
        
#         # Convert the merged chunk to an AnnData object
#         merged_chunk_adata = sc.AnnData(X=merged_chunk.T, var=pd.DataFrame(index=features), obs=pd.DataFrame(index=barcodes))
        
#         # Get stats for merged chunk
#         merged_chunk_stats = calc_stats(merged_chunk_adata)
#         all_stats.append(merged_chunk_stats)
#     final_stats = pd.concat(all_stats, ignore_index=True)
    
#     return final_stats
    



# def do_group_wise_calcs(files_dt, split_var = None, meta_f = None, sample_var = 'sample_id', chunk_size = 2000):

#     files = pd.read_csv(files_dt)
#     #smpl_fs_ls = files['amb_out_f'].tolist()
#     #samples = files[sample_var].tolist()

#     if split_var is not None:
#         # check if meta_f available and ok
#         assert meta_f is not None, \
#             "meta_f missing for splitting samples into groups for hvg calculations"
#         assert os.path.isfile(meta_f), \
#             f"{meta_f} is not a file"
       
#         meta = pd.read_csv(meta_f)
#         assert split_var in meta.columns(), \
#             f"{split_var} variable not present in {meta_f}"
        
#         # get sample categories
#         smpl_cats = meta[split_var].unique().tolist()
#         smpl_grps = {}

#         for c in smpl_cats:
#             smpls_in_grp = meta.loc[meta[split_var] == c, 'sample_id'].tolist()
#             smpl_grps[c] = smpls_in_grp
        

#         # initalize empty dictionaty to store hvg stats per group
#         grp_stats_dict = {}

#         for grp in smpl_grps:
#             # get files
#             grp_fs = files.loc[files[sample_var] in grp, "amb_out_f"].tolist()
#             # read chunks from each file and merge
              
#               # get row number
#             with h5py.File(grp_fs[0], 'r') as f:
#                 num_rows = f['matrix/shape'][0]

#             all_grp_stats = [] # initialize empty list to collect stats for all chunks
#             for start_row in range(0, num_rows, chunk_size):
#                 end_row = min(start_row + chunk_size, num_rows)
#                 merged_chunk = None
                
#                 for f in grp_fs:
#                     chunk = read_csr_chunk(f, start_row, end_row, as_adata = True)
#                     merged_chunk = chunk if merged_chunk is None else csr_matrix.hstack([merged_chunk, chunk])

#                 # get stats for merged chunk
#                 merged_chunk_stats = calc_stats(merged_chunk)
#                 all_grp_stats.append(merged_chunk_stats) # ideally use something else instead of start_row for dictionary names
            

#             # merge group stats
#             final_grp_stats = pd.concat(all_grp_stats)
#             final_grp_stats[split_var] = grp
#             # add to dictionarty
#             grp_stats_dict[grp] = final_grp_stats

        




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

