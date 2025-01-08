import pandas as pd
import gzip
import argparse
import os
import re
import glob
import subprocess
import numpy as np




# function that gets params for simpleaf index
def parse_setup_params_for_af(genome, params_csv):
 # read params csv created in rule get_reference_genomes
 params_df = pd.read_csv(params_csv, dtype={'decoy': bool})

 filt_params_df = params_df[params_df['genome_name'] == genome] 
 filt_params_df = filt_params_df.reset_index(drop=True)
 fasta_f        = filt_params_df.loc[0, 'fasta_f']
 gtf_f          = filt_params_df.loc[0, 'gtf_f']
 dcoy           = filt_params_df.loc[0, 'decoy']
 rrna           = filt_params_df.loc[0, 'rrna']
 
 if dcoy:
  w_dcoy = 'yes'
 else:
  w_dcoy = 'no'

 # get zenodo urls for prebuild indices
 zenodo_idx_urls = { 
   'human_2020': "https://zenodo.org/records/14247195/files/alevin-idx_human_2020_rrna.tar.gz", 
   'human_2024': "https://zenodo.org/records/14247195/files/alevin-idx_human_2024_rrna.tar.gz", 
   'mouse_2020': "https://zenodo.org/records/14247195/files/alevin-idx_mouse_2020_rrna.tar.gz", 
   'mouse_2024': "https://zenodo.org/records/14247195/files/alevin-idx_mouse_2024_rrna.tar.gz"
 }

 # check if prebuild genome is available
 idx_url = None

 if (
    genome in zenodo_idx_urls and 
    dcoy and  # This works because `dcoy` is a boolean
    rrna and  # This works because `rrna` is a boolean
    (pd.isna(fasta_f) or fasta_f == "") and  # Check for NaN or empty string
    (pd.isna(gtf_f) or gtf_f == "")         # Check for NaN or empty string
 ):
    idx_url = zenodo_idx_urls[genome]
 
 return fasta_f, gtf_f, w_dcoy, idx_url



# function that makes simpleaf index
def make_af_idx(genome, params_csv, scprocess_data_dir, cores):
  # get af params for combn
  fasta_f, gtf_f, w_dcoy, idx_url = parse_setup_params_for_af(genome, params_csv)
  
  # create af home directory
  af_home = os.path.join(scprocess_data_dir, 'alevin_fry_home')
  os.makedirs(af_home, exist_ok=True)

  # specify output directory for index
  idx_out_dir = os.path.join(af_home, genome)

  if idx_url is not None:

    print('Downloading alevin index for ' + genome)
    # make output directory
    os.makedirs(idx_out_dir, exist_ok=True)
    os.chdir(idx_out_dir)
    # download index
    subprocess.run(f"wget {idx_url}", shell=True)
    idx_name = f'alevin-idx_{genome}_rrna.tar.gz'

    # untar
    subprocess.run(f'tar --strip-components=1 -xvf {idx_name}', shell=True, capture_output=False)
    # remove tar archive
    os.remove(idx_name)

  else:
    if w_dcoy == 'yes':
      print('Creating alevin index for ' + genome + ' with decoys in ' + idx_out_dir )
    else:
      print('Creating alevin index for ' + genome + ' in ' + idx_out_dir )
     
    # define whether or not to include --decoy-paths flag
    decoy_flag = f"--decoy-paths {fasta_f}" if w_dcoy == 'yes' else ""

    # code for making index
    bash_script = f"""
    #!/bin/bash
    ulimit -n 2048

    # simpleaf configuration
    export ALEVIN_FRY_HOME={af_home}
    simpleaf set-paths
  
    cd $ALEVIN_FRY_HOME

    simpleaf index \
    --output {idx_out_dir} \
    --fasta {fasta_f} \
    --gtf {gtf_f} \
    {decoy_flag} \
    --overwrite \
    --threads {cores}
    """

    # run bash script
    subprocess.run(bash_script, shell=True, executable='/bin/bash', check=True)

  print('Done')
  
  return


# make script executable from the command line
if __name__ == '__main__':
    # Define arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('genome', type=str)
    parser.add_argument('params_f', type=str)
    parser.add_argument('data_dir', type=str)
    parser.add_argument('cores', type = int)

    # Get arguments
    args = parser.parse_args()
    # run funciton that makes index
    make_af_idx(args.genome, args.params_f, args.data_dir, args.cores)

