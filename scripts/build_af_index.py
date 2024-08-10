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
 fasta_f = filt_params_df.loc[0, 'fasta_f']
 gtf_f  = filt_params_df.loc[0, 'gtf_f']
 dcoy = filt_params_df.loc[0, 'decoy']
 
 # make name for alevin index directory
 if dcoy:
  w_dcoy = 'yes'
 else:
  w_dcoy = 'no'

 return fasta_f, gtf_f, w_dcoy



# function that makes simpleaf index
def make_af_idx(genome, params_csv, scprocess_data_dir, cores):
  # get af params for combn
  fasta_f, gtf_f, w_dcoy = parse_setup_params_for_af(genome, params_csv)
  
  # create af home directory
  af_home = os.path.join(scprocess_data_dir, 'alevin_fry_home')
  os.makedirs(af_home, exist_ok=True)

  # specify output directory for index
  idx_out_dir = os.path.join(af_home, genome)

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

