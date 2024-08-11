# script for running simpleaf quant
import pandas as pd
import argparse
import os
import subprocess


def parse_params_for_af_quant(sample, chem_csv):
    # read dataframe with chemistries
 chem_df = pd.read_csv(chem_csv)
 filt_chem_df = chem_df[chem_df['sample_id'] == sample]
 filt_chem_df = filt_chem_df.reset_index(drop = True)
 chem = filt_chem_df.loc[0, 'version']

 assert chem in ['3v3', '3v2_5v1_5v2'], \
  f"{chem} chemistry not supported by scprocess"

 af_chem =  filt_chem_df.loc[0, 'af_chem']
 wl_f = filt_chem_df.loc[0, 'wl_f']
 exp_ori = filt_chem_df.loc[0,'expected_ori']
 
 return af_chem, wl_f, exp_ori



# simpleaf configuration
def get_simpleaf_counts(sample, chem_csv, AF_HOME_DIR, af_out_dir, AF_INDEX_DIR, R1_fs, R2_fs, cores):

    # get params
    AF_CHEMISTRY, whitelist_f, exp_ori = parse_params_for_af_quant(sample, chem_csv)
    #R1_fs = ','.join(R1_fs)
    #R2_fs = ','.join(R2_fs)

    bash_script = f"""
    #!/bin/bash

    # set up alevin
    export ALEVIN_FRY_HOME="{AF_HOME_DIR}"
    simpleaf set-paths
  
    # simpleaf quantfication
    simpleaf quant \
      --reads1 {R1_fs}  \
      --reads2 {R2_fs}  \
      --threads {cores} \
      --index {AF_INDEX_DIR}/index \
      --chemistry {AF_CHEMISTRY} --resolution cr-like \
      --expected-ori {exp_ori} \
      --t2g-map {AF_INDEX_DIR}/index/t2g_3col.tsv \
      --unfiltered-pl {whitelist_f} --min-reads 1 \
      --output {af_out_dir}/af_{sample}
    """
	

    print('Running simpleaf quant for ' + sample)
    # run bash script
    subprocess.run(bash_script, shell=True, executable='/bin/bash', check=True)

    print('Done')

    return

# make script work from command line

if __name__ == "__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument('sample', type = str)
   parser.add_argument('chem_csv', type=str)
   parser.add_argument('afhome', type=str)
   parser.add_argument('afout', type=str)
   parser.add_argument('afidx', type = str) # comma separated string with fastq files
   parser.add_argument('R1_fs', type=str) # comma separated string with fastq files
   parser.add_argument('R2_fs', type = str)
   parser.add_argument('cores', type=int)

   args = parser.parse_args()

   get_simpleaf_counts(args.sample,args.chem_csv, args.afhome, args.afout, args.afidx, args.R1_fs, args.R2_fs, args.cores)
   
    
