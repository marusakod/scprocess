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
         af_chem = None
         exp_ori = None
    elif version == '3v1':
         af_chem = None
         exp_ori = None
    else:
         af_chem = '10xv3'
         exp_ori = 'fw'

    # Summarize output  
    chem_stats_df = pd.DataFrame({
        'sample_id': sample,
        'n_R1_fastqs': len(R1_f), 
        'max_pct_overlap': round(max(max_values) * 100, 1),
        'version': version, 
        'af_chem': af_chem, 
        'expected_ori': exp_ori, 
        'wl_f': wl_dict[version]
    })

    return chem_stats_df

# function to process all samples in parallel
def process_samples(fastq_dir, SAMPLES, wl_dict, chem_stats_f, num_cores, custom_chem_f = None):

    # check if custom chemistry file exists 
    if custom_chem_f:
     # read file
     chem_df = pd.read_csv(custom_chem_f)
     chem_df = pd.DataFrame({
        'sample_id': ['sample1', 'sample2'],
        'version': ['3v3', '5v2']
     })
     
     # check that all samples are in the custom file
     for s in SAMPLES:
        assert s in chem_df['sample_id'].tolist(), \
         f"Sample {s} missing from {custom_chem_f}"

     # check if all chemistries are valid
     valid_chems = ['3LT', '3v2', '5v1', '5v2', '3v3', 'multiome']
     assert all([c in valid_chems for c in chem_df['version'].tolist()]), \
      f"10x chemistries specified in {custom_chem_f} are not valid. Valid values are: {', '.join(valid_chems)}"
     
     # define af_chem and expected_ori based on version
     chem_df['af_chem'] = chem_df['version'].apply(lambda x: '10xv2' if x in ['3v2', '5v1', '5v2'] else '10xv3')
     chem_df['wl_chem'] = chem_df['version'].apply(lambda x: '3v2_5v1_5v2' if x in ['3v2', '5v1', '5v2'] else x)
     chem_df['expected_ori'] = chem_df['version'].apply(lambda x: 'fw' if x in ['3LT', '3v2', '3v3', 'multiome'] else 'both')
     # get the right whitelists
     chem_df['wl_f'] = [wl_dict[c] for c in chem_df['wl_chem'].tolist()]

    else:
     with multiprocessing.Pool(processes= num_cores) as pool:
      results = pool.starmap(detect_sample_chemistry, [(fastq_dir, s, wl_dict) for s in SAMPLES])

    # concatenate all data frames in the list
      chem_df = pd.concat(results, ignore_index=True)

    # Save to CSV
    chem_df.to_csv(chem_stats_f, index=False)

    print('Done!')
    return


def list_of_strings(arg):
    return arg.split(',')


if __name__ == "__main__":
    # get argument values
    argParser = argparse.ArgumentParser()
    argParser.add_argument("-f", "--fastqdir", type=str, help="project directory", required=True)
    argParser.add_argument("-s", "--samples", type =list_of_strings, help="path to metadata file", required=True)
    argParser.add_argument("-o", "--outf", type=str,  help="path to output file", required=True)
    argParser.add_argument("-b", "--barcodes", type = str, help=".csv file with paths to barcode whitelists", required=True)
    argParser.add_argument("-c", "--cores", type = int, help="number of cores to use", default= 1)

    args = argParser.parse_args()

    wl_csv_f = args.barcodes

    # get whitelist files
    wl_df = pd.read_csv(wl_csv_f)
    
    # dictionary to hold whitelist data
    wl_dict = {}
    for index, row in wl_df.iterrows():
        chem = row['chemistry']
        barcodes_f = row['barcodes_f']
        wl_dict[chem] = np.loadtxt(barcodes_f, dtype='str')

    # process samples in parallel
    process_samples(fastq_dir=args.fastqdir, SAMPLES=args.samples, wl_dict=wl_dict, chem_stats_f=args.outf, num_cores=args.cores)

