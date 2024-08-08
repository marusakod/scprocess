# functions for setting up scprocess_data

import pandas as pd
import pyranges as pr
import gzip
import argparse
import os
import yaml
import re
import glob
import subprocess
import numpy as np
import re

#configfile = '/Users/marusa/Projects/scprocess_test/configs/config-setup-template.yaml'

def get_genome_params(GENOME_NAMES, FASTA_FS, GTF_FS, MITO_STRS, DECOYS, SCPROCESS_DATA_DIR):

  # make dictionaries with all parameters that will be updated later
  fasta_dict = dict(zip(GENOME_NAMES, FASTA_FS))
  gtf_dict = dict(zip(GENOME_NAMES, GTF_FS))
  mito_str_dict = dict(zip(GENOME_NAMES, MITO_STRS))
  decoys_dict = dict(zip(GENOME_NAMES, DECOYS))


  # crete directory for reference genomes inside scprocess data
  ref_dir = os.path.join(SCPROCESS_DATA_DIR, 'reference_genomes')
  os.makedirs(ref_dir, exist_ok = True)

  # sort out predefined genomes
  pre_gnomes = {'human_2020', 'human_2024', 'mouse_2020', 'mouse_2024'}
  all_gnomes = set(GENOME_NAMES)

# which predefined genomes are in params_df
  names = list(pre_gnomes & all_gnomes)
  
  if(len(names) != 0):
    
    refs_10x_links_dict ={
      'human_2020': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz", 
      'human_2024': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz", 
      'mouse_2020': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz", 
      'mouse_2024': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
    }

    # download all required genomes from 10x
    for n in names:
      l = refs_10x_links_dict[n]
      gnome_fs = get_predefined_gnomes(ref_dir= ref_dir, name = n, link = l )
      # add fasta_f and gtf_f to dicts
      fasta_dict[n] =  gnome_fs[1]
      gtf_dict[n] = gnome_fs[0]
    # define mito strings for predefined genomes
  
    for  n in names:
      mito_str = "^MT-" if 'human_' in n else "^mt-"
      mito_str_dict[n] = mito_str

# sort out custom genomes
  cust_names = list(all_gnomes - pre_gnomes)
  if(len(cust_names) != 0):
    print('Sorting out locations for custom genomes')
    # create directories for custom genomes in ref_dir and symlinks for fasta and gtf files 
    for cn in cust_names:
      # get original fasta
      cn_fasta = fasta_dict[cn]
      # get original gtf
      cn_gtf = gtf_dict[cn]

      
      cn_ref_dir = os.path.join(ref_dir, cn)
      os.makedirs(cn_ref_dir, exist_ok = True)
      fasta_sym = os.path.join(cn_ref_dir, 'genes.fa')
      gtf_sym = os.path.join(cn_ref_dir, 'genes.gtf')

      if not os.path.exists(fasta_sym):
        os.symlink(cn_fasta, fasta_sym)
      if not os.path.exists(gtf_sym):
        os.symlink(cn_gtf, gtf_sym)

      # update paths to fasta and gtf files in dicts
      fasta_dict[cn] = fasta_sym
      gtf_dict[cn] = gtf_sym
  
  # make a txt file from all gtf files
  print('Creating txt files from gtf')
  gtf_txt_dict = {}
 
  for n, gtf in gtf_dict.items():
    # define output file
    out_txt_f = os.path.join(os.path.dirname(gtf), 'genes_gtf.txt.gz')
    # convert .gtf to .txt
    save_gtf_as_txt(gtf, out_txt_f)
    # update dictionary 
    gtf_txt_dict.update({n:out_txt_f})

  # create one dataframe from all dictionaries 
  PARAMS_DF = pd.DataFrame([fasta_dict, gtf_dict, mito_str_dict, decoys_dict])
  PARAMS_DF = PARAMS_DF.T
  PARAMS_DF.columns = ['fasta_f', 'gtf_f', 'mito_str', 'decoy']
  PARAMS_DF.index.name = 'genome_name'
  PARAMS_DF.reset_index(inplace=True)
  
  out_csv_f = os.path.join(SCPROCESS_DATA_DIR, 'setup_parameters.csv')
  PARAMS_DF.to_csv(out_csv_f, index = False)
  
  print(f'Completed writing genomes in {SCPROCESS_DATA_DIR}')
  return 
  


# get refrence genomes from 10x
def get_predefined_gnomes(ref_dir, name, link):

  print(f'Downloading {name} genome from 10x')
  
  # get tarball 
  tball = link.rsplit('/', 1)[-1]
  
  # create a new directory for the specified genome and switch to it
  gnome_dir = os.path.join(ref_dir, name)
  os.makedirs(gnome_dir, exist_ok = True)
  os.chdir(gnome_dir)
        
  # download reference into that dir
  subprocess.run(f'wget {link}', shell=True)
        
  # list tar contents
  tar_out= subprocess.run(f'tar -tzf {tball}', shell=True, capture_output=True, text=True)
  tball_all_fs = tar_out.stdout.splitlines()
        
  # look for genome.fa and genes.gtf (or genes.gtf.gz) files and extract the
  patt = re.compile(r'(genome\.fa|genes\.gtf(\.gz)?)$') 
  tball_filt_fs = [file for file in tball_all_fs if patt.search(file)]

  for t in tball_filt_fs:
    ext_spell = ['tar', '--extract', '--file=' + tball, '--strip-components=2', t]
    subprocess.run(ext_spell, check=True)

  # unzip
  gnome_fs = os.listdir()
  for f in gnome_fs:
    if f.endswith('.gtf.gz'):
      subprocess.run(f'gunzip {f}', shell=True, capture_output=False)

  # delete the tarball
  os.remove(tball)
  print('Done!')

  return [os.path.join(gnome_dir, 'genes.gtf'), os.path.join(gnome_dir, 'genome.fa')]




# download  marker genes, xgboost objects, barcode whitelists and gmt pathways from github
def get_scprocess_data(scprocess_data_dir):

  print('Downloading data from scprocessData github repo')
  # switch to scprocess data dir
  os.chdir(scprocess_data_dir)

  # download tar archive from scprocessData and extract
  subprocess.run('wget -O - https://github.com/marusakod/scprocessData/releases/download/v0.1.0/scprocess_data_archive.tar.gz | tar xvf - --strip-components=1',
  shell=True, capture_output=False)
  
  # check if all necessary directories are there
  dirs_ls = ["cellranger_ref", "gmt_pathways", "marker_genes", "xgboost"]
  for dir in dirs_ls:
    assert os.path.isdir(dir), \
      f"{dir} directory doesn't exist"
  
  print('Done!')

  return



# convert genes.gtf file to txt file
def save_gtf_as_txt(gtf_f, gtf_txt_f):
  
    gtf = pr.read_gtf(gtf_f, as_df = True)

    # restrict to just genes
    gene_annots = gtf[gtf['Feature'] == 'gene']
    # remove some columnns
    gene_annots = gene_annots[['gene_id', 'gene_name', 'gene_type', 'Chromosome', 'Start', 'End', 'Strand']]
    
    # Rename columns
    gene_annots.rename(columns={
        'gene_id': 'ensembl_id',
        'gene_name': 'symbol',
        'Chromosome': 'chromosome',
        'Start': 'start',
        'End': 'end',
        'Strand': 'strand'
    }, inplace=True)

    # calculate width and add gene_id column
    gene_annots['width'] = gene_annots['end'] - gene_annots['start'] + 1
    gene_annots['gene_id'] = gene_annots['symbol'] + '_' + gene_annots['ensembl_id']
    
    # do some sorting of chromosomes 
    chr_freq = gene_annots['chromosome'].value_counts()
    gene_annots['chr_freq'] = gene_annots['chromosome'].map(chr_freq)
    gene_annots_sorted = gene_annots.sort_values(by=['chr_freq', 'start', 'end'], ascending=[False, True, True])
    gene_annots_sorted = gene_annots_sorted.drop(columns=['chr_freq'])
    
    # make gene_id the first column
    gene_id_col = gene_annots_sorted.pop('gene_id')
    gene_annots_sorted.insert(0, 'gene_id', gene_id_col)
    
    # make sure there are no duplicated gene ids
    assert gene_annots_sorted['gene_id'].duplicated().sum() == 0

    # save
    with gzip.open(gtf_txt_f, 'wt') as f:
        gene_annots_sorted.to_csv(f, sep='\t', index=False)
        
    return


# function that gets params for simpleaf index
def parse_setup_params_for_af(genome, params_csv):
 # read params csv created in rule get_reference_genomes
 params_df = pd.read_csv(params_csv, dtype={'decoy': bool})

 filt_params_df = params_df[(params_df['genome_name'] == genome)] 
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

  if w_dcoy == 'yes':
    print('Creating alevin index for ' + genome + ' with decoys in ' + idx_out_dir )
  else:
    print('Creating alevin index for ' + genome + ' in ' + idx_out_dir )
   
  # create af home directory
  af_home = os.path.join(scprocess_data_dir, 'alevin_fry_home')
  os.makedirs(af_home, exist_ok=True)

  # specify output directory for index
  idx_out_dir = os.path.join(af_home, genome)
  
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
  subprocess.Popen(bash_script, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  print('Done')
  
  return


def list_of_strings(arg):
    return arg.split(',')

def list_of_bools(arg):
    str_ls = arg.split(',')
    bool_list = [bool(s) for s in str_ls]
    return bool_list



  # make some functions executable from the command line (make_af_idx, get_genome_params, get_scprocess_data)


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  subparsers = parser.add_subparsers(dest="function_name", help="Name of the function to run")

    # parser for get_scprocess_data
  parser_getdata = subparsers.add_parser('get_scprocess_data')
  parser_getdata.add_argument('scp_data_dir', type=str)

    # parsers for get_genome_params
  parser_getgnomes = subparsers.add_parser('get_genome_params')
  parser_getgnomes.add_argument('genomes', type=list_of_strings, help='list with all genome names')
  parser_getgnomes.add_argument('fasta_fs', type=list_of_strings, help='list with paths to all fasta files')
  parser_getgnomes.add_argument('gtf_fs', type=list_of_strings, help='list with paths to all gtf files')
  parser_getgnomes.add_argument('mito_str', type=list_of_strings, help='list with all mitochondrial gene identifiers')
  parser_getgnomes.add_argument('decoys', type=list_of_bools, help='list with bool values definiing whether or not to use decoys when building indices with simpleaf')
  parser_getgnomes.add_argument('scp_data_dir', type=str)

    # parser for make_af_idx
  parser_afidx = subparsers.add_parser('make_af_idx')
  parser_afidx.add_argument('genome', type=str, help="genome name")
  parser_afidx.add_argument('params', type=str, help='path to csv file with all parameters (output of get_genome_params)')
  parser_afidx.add_argument('scp_data_dir', type=str)
  parser_afidx.add_argument('cores', type=int)
  
  args = parser.parse_args()

  if args.function_name == 'get_scprocess_data':
      get_scprocess_data(args.scp_data_dir)
  elif args.function_name == 'get_genome_params':
      get_genome_params(args.genomes, args.fasta_fs, args.gtf_fs, args.mito_str, args.decoys, args.scp_data_dir)
  elif args.function_name == 'make_af_idx':
      make_af_idx(args.genome, args.params, args.scp_data_dir, args.cores)
  else:
    parser.print_help()


