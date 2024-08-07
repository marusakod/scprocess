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

  # make a dataframe with all params and later update
  PARAMS_DF = pd.DataFrame({
    'genome_name':GENOME_NAMES, 
    'fasta_f': FASTA_FS, 
    'gft_f':GTF_FS, 
    'mito_str': MITO_STRS, 
    'decoy':DECOYS

  })

  # crete directory for reference genomes inside scprocess data
  ref_dir = os.path.join(SCPROCESS_DATA_DIR, 'reference_genomes')
  os.makedirs(ref_dir, exist_ok = True)

  # sort out predefined genomes
  pre_gnomes = {'human_2020', 'human_2024', 'mouse_2020', 'mouse_2024'}
  all_gnomes = set(GENOME_NAMES) # only unique values bc decoys can be both True and False

# which predefined genomes are in params_df
  names = list(pre_gnomes & all_gnomes)
  
  if(len(names) != 0):
    print('Downloading 10x genomes for ' + ', '.join(names))

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
      PARAMS_DF.loc[PARAMS_DF['genome_name'] == n, 'fasta_f'] = gnome_fs[1]
      PARAMS_DF.loc[PARAMS_DF['genome_name'] == n, 'gtf_f'] = gnome_fs[0]
    # define mito strings for predefined genomes
  
    for  n in names:
      mito_str = "^MT-" if 'human_' in n else "^mt-"
      PARAMS_DF.loc[PARAMS_DF['genome_name'] == n, 'mito_str'] = mito_str

# sort out custom genomes
  cust_names = list(all_gnomes - pre_gnomes)
  if(len(cust_names) != 0):
    print('Sorting out locations for custom genomes')
    # create directories for custom genomes in ref_dir and symlinks for fasta and gtf files 
    for cn in cust_names:
      # get original fasta
      cn_fasta = PARAMS_DF[PARAMS_DF['genome_name'] == cn]['fasta_f'].unique().tolist()
      # get original gtf
      cn_gtf = PARAMS_DF[PARAMS_DF['genome_name'] == cn]['gtf_f'].unique().tolist()
      
      cn_ref_dir = os.path.join(ref_dir, cn)
      os.makedirs(cn_ref_dir, exist_ok = True)
      fasta_sym = os.path.join(cn_ref_dir, 'genes.fa')
      gtf_sym = os.path.join(cn_ref_dir, 'genes.gtf')
      os.symlink(cn_fasta, fasta_sym)
      os.symlink(cn_gtf, gtf_sym)

      # update paths to fasta and gtf files in params_df
      PARAMS_DF.loc[PARAMS_DF['genome_name'] == cn, 'fasta_f'] = fasta_sym
      PARAMS_DF.loc[PARAMS_DF['genome_name'] == cn, 'gtf_f'] = gtf_sym
  
  # make a txt file from all gtf files
  print('Creating txt files from gtf')
  gtf_df = PARAMS_DF[['genome_name', 'gtf_f']].drop_duplicates()
  gtf_dict = gtf_df.set_index('genome_name')['gtf_f'].to_dict()
  gtf_txt_dict = {}
 
  for n, gtf in gtf_dict.items():
    # define output file
    out_txt_f = os.path.join(os.path.dirname(gtf), 'genes_gtf.txt.gz')
    # convert .gtf to .txt
    save_gtf_as_txt(gtf, out_txt_f)
    # update dictionary 
    gtf_txt_dict.update({n:out_txt_f})

  
  # add gtf_txt_fs to params_df and return a csv file with all specified parameters
  PARAMS_DF['gtf_txt_f'] = PARAMS_DF['genome_name'].map(gtf_txt_dict)
  
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
  # subprocess.run(f'wget {link}', shell=True)
        
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
  #os.remove(tball)
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
  
  # remove tar archive
  os.remove('scprocess_data_archive.tar.gz')
  
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
def parse_setup_params_for_af(combn, params_csv):
 # read params csv created in rule get_reference_genomes
 params_df = pd.read_csv(params_csv, dtype={'decoy': bool})

 # split genome name and decoy
 sep = '_'
 gnome, _sep, _after = combn.rpartition(sep)
 _before, _sep, dcoy = combn.rpartition(sep)
 dcoy = bool(dcoy)
 # get all parameters for af from params_df
 filt_params_df = params_df[(params_df['genome_name'] == gnome) & (params_df['decoy'] == dcoy)]
 filt_params_df = filt_params_df.reset_index(drop=True)
 fasta_f = filt_params_df.loc[0, 'fasta_f']
 gtf_f  = filt_params_df.loc[0, 'gtf_f']

 # make name for alevin index directory
 if dcoy:
  w_dcoy = 'yes'
 else:
  w_dcoy = 'no'

 return gnome, fasta_f, gtf_f, w_dcoy



# function that makes simpleaf index
def make_af_idx(combn, params_csv, scprocess_data_dir, idx_out_dir, cores):
  # get af params for combn
  gnome, fasta_f, gtf_f, w_dcoy = parse_setup_params_for_af(combn, params_csv)

  if w_dcoy == 'yes':
    print('Creating alevin index for ' + gnome + ' with decoys in ' + idx_out_dir )
  else:
    print('Creating alevin index for ' + gnome + ' in ' + idx_out_dir )
   
  # create af home directory
  af_home = os.path.join(scprocess_data_dir, 'alevin_fry_home')
  os.makedirs(af_home, exist_ok=True)
  
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


  # make some functions executable from the command line (make_af_idx, get_genome_params, get_scprocess_data)


if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  subparsers = parser.add_subparsers(dest="function_name", help="Name of the function to run")

    # parser for get_scprocess_data
  parser_getdata = subparsers.add_parser('get_scprocess_data')
  parser_getdata.add_argument('scp_data_dir', type=str)

    # parsers for get_genome_params
  parser_getgnomes = subparsers.add_parser('get_genome_params')
  parser_getgnomes.add_argument('genome_names', type=list, help='list with all genome names')
  parser_getgnomes.add_argument('fasta_fs', type=list, help='list with paths to all fasta files')
  parser_getgnomes.add_argument('gtf_fs', type=list, help='list with paths to all gtf files')
  parser_getgnomes.add_argument('mito_str', type=list, help='list with all mitochondrial gene identifiers')
  parser_getgnomes.add_argument('decoys', type=list, help='list with bool values definiing whether or not to use decoys when building indices with simpleaf')
  parser_getgnomes.add_argument('scp_data_dir', type=str)

    # parser for make_af_idx
  parser_afidx = subparsers.add_parser('make_af_idx')
  parser_afidx.add_argument('combn', type=str, help="combination of genome name and whether or not decoys should be used i.e. 'human_2024'")
  parser_afidx.add_argument('params', type=str, help='path to csv file with all parameters (output of get_genome_params)')
  parser_afidx.add_argument('scp_data_dir', type=str)
  parser_afidx.add_argument('out_dir', type = str, help= "path to alevin index output directory")
  parser_afidx.add_argument('cores', type=int)
  
  args = parser.parse_args()

  if args.function_name == 'get_scprocess_data':
      get_scprocess_data(args.scp_data_dir)
  elif args.function_name == 'get_genome_params':
      get_genome_params(args.genome_names, args.fasta_fs, args.gtf_fs, args.mito_str, args.decoys, args.scp_data_dir)
  elif args.function_name == 'make_af_idx':
      make_af_idx(args.combn, args.params, args.scp_data_dir, args.out_dir, args.cores)
  else:
    parser.print_help()


