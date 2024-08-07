# detecting 10x chemistry - requires seqkit!
# parallelisation is working but printing of messages is weird 

import os
import glob
import re
import argparse
import numpy as np
import pandas as pd
import random
import subprocess
import multiprocessing


def find_R1_files(fastq_dir, sample):
        # get all files
        all_fs = glob.glob(f"{fastq_dir}/*{sample}*")

        # get all reads
        re_R1 = re.compile('.*_R1_.*\\.fastq(\\.gz)?$')
        R1_fs = [f for f in all_fs if re_R1.match(f)]
        R1_fs = sorted(R1_fs)
        
        return R1_fs
      
def extract_raw_seqs_from_fq(fastq_all):
        entries = fastq_all.strip().split('\n@')[0:]
        seqs = []
        
        for e in entries:
            ls = e.split('\n')
            # the sequence is the second line
            seq = ls[1]
            seqs.append(seq)
        
        return seqs     


# Function to detect chemistry for one sample
def detect_sample_chemistry(fastq_dir, sample, wl_dict, r1_max=3):
   
    # Get list with R1 fastqs
    R1_f = find_R1_files(fastq_dir, sample)
    
    if len(R1_f) > r1_max:
        # Sample 3 R1 fastqs
        random.seed(12346)
        R1_f = random.sample(R1_f, r1_max)

    # Store overlap proportions 
    overlap_dict = {chem: [] for chem in wl_dict}

    for file in R1_f:
        print(f'Extracting barcodes from ' + file)
        # command = f"seqkit head -n 1000000 {file} | seqkit subseq -r 1:16 | awk '/^[NATGC]/'"
        command = f"seqkit head -n 1000000 {file} | seqkit subseq -r 1:16"
        spell_res = subprocess.run(command, shell=True, capture_output=True, text=True)
        spell_out = spell_res.stdout
        barcodes = extract_raw_seqs_from_fq(spell_out)

        # Filter unique barcodes
        barcodes = np.unique(barcodes)

        n = len(barcodes)
        print(f'Number of unique barcodes: {n}')

        # calculate fraction of barcodes in each whitelist
        for chem, wl in wl_dict.items():
            overlap = sum(np.in1d(barcodes, wl)) / n
            print(f'Overlap with whitelist {chem}: {overlap:.1%}')
            overlap_dict[chem].append(overlap)

    # check if version is consistent across fastqs
    # get the maximum overlap (and corresponding chemistry) for each fastq
    max_values = [None] * len(R1_f)
    max_keys = [None] * len(R1_f)
    
    # Loop through the indices
    for i in range(len(R1_f)):
        # initialize the maximum value and corresponding key for this index
        max_value = 0
        max_key = None
        
        # compare elements at the same index across all items
        for key, values in overlap_dict.items():
            if values[i] > max_value:
                max_value = values[i]
                max_key = key
        
        # store the maximum value and corresponding key
        max_values[i] = max_value
        max_keys[i] = max_key
    
    # check that versions are consistent
    assert len(set(max_keys)) == 1, "10x chemistry version inconsistent across fastq files"
      
    # declare chemistry
    version = max_keys[0]
    print(f'10x chemistry version for {sample} is {version}!')

    # select alevin chemistry setting and expected ori based on version detected
    if version == '3v2_5v1_5v2':
         af_chem = '10xv2'
         exp_ori = 'both'
    elif version == 'flex':
         af_chem = 'custom'
         exp_ori = 'custom'
    else:
         af_chem = '10xv3'
         exp_ori = 'fw'

    # Summarize output  
    chem_stats_df = pd.DataFrame([{
        'sample_id': sample,
        'n_R1_fastqs': len(R1_f), 
        'max_pct_overlap': round(max(max_values) * 100, 1),
        'version': version, 
        'af_chem': af_chem, 
        'exoected_ori': exp_ori
    }])

    return chem_stats_df

# function to process all samples in parallel
def process_samples(fastq_dir, meta_f, wl_dict, chem_stats_f, num_cores):
    meta = pd.read_csv(meta_f, header=0, index_col=False)
    samples = list(meta['sample_id'])
    
    with multiprocessing.Pool(processes=num_cores) as pool:
        results = pool.starmap(detect_sample_chemistry, [(fastq_dir, s, wl_dict) for s in samples])

    # concatenate all data frames in the list
    all_chem_df = pd.concat(results, ignore_index=True)
    # Save to CSV
    all_chem_df.to_csv(chem_stats_f, index=False)

    print('Done!')


def main():
    # get argument values
    argParser = argparse.ArgumentParser()
    argParser.add_argument("-f", "--fastqdir", help="project directory", required=True)
    argParser.add_argument("-m", "--metaf", help="path to metadata file", required=True)
    argParser.add_argument("-o", "--outf", help="path to output file", required=True)
    argParser.add_argument("-b", "--barcodes", help=".csv file with paths to barcode whitelists", required=True)
    argParser.add_argument("-c", "--cores", help="number of cores to use", type=int, default=multiprocessing.cpu_count())

    args = argParser.parse_args()

    fastqdir = args.fastqdir
    meta_f = args.metaf
    out_f = args.outf
    wl_csv_f = args.barcodes
    num_cores = args.cores

    # get whitelist files
    wl_df = pd.read_csv(wl_csv_f)
    
    # dictionary to hold whitelist data
    wl_dict = {}
    for index, row in wl_df.iterrows():
        chem = row['chemistry']
        barcodes_f = row['barcodes_f']
        wl_dict[chem] = np.loadtxt(barcodes_f, dtype='str')

    # process samples in parallel
    process_samples(fastq_dir=fastqdir, meta_f=meta_f, wl_dict=wl_dict, chem_stats_f=out_f, num_cores=num_cores)

if __name__ == "__main__":
    main()
