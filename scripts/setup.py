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


def get_setup_params(configfile):
  print('doing some checks on inputs')
  # check configfile exists
  assert os.path.exists(configfile), \
    f"Config file {configfile} does not exist"

  # open configfile
  with open(configfile, "r") as stream:
    config      = yaml.safe_load(stream)

  # check that scprocess_data_dir exists and is a directory
  assert "scprocess_data_dir" in config, "scprocess_data_dir not in config file"
  assert os.path.isdir(config["scprocess_data_dir"]), \
    f"proj_dir {config['scprocess_data_dir']} is not a directory"
  scprocess_data_dir    = config["scprocess_data_dir"]
  
  # check that genome is in config
  assert "genome" in config, "genome not defined in config file"
  
  # check that either species or custom are in the config 
  assert 'species' in config['genome'] or 'custom' in config['genome'], \
    "either 'species' or 'custom' need to be specified under genome"
    
  # crete directory for reference genomes inside scprocess data
  ref_dir = os.path.join(scprocess_data_dir, 'reference_genomes')
  os.makedirs(ref_dir, exist_ok = True)
  
  # start a dictionary with paths to fa and gtf files for all specified genomes
  fasta_dict = {}
  gtf_dict={}
  mito_str_dict={}
  decoy_dict={}
    
  if('species' in config['genome']):
  
    # crete directory for reference genomes inside scprocess data
    ref_dir = os.path.join(scprocess_data_dir, 'reference_genomes')
    os.makedirs(ref_dir, exist_ok = True)
  
    
    allowed_names = ['human_2020', 'human_2024', 'mouse_2020', 'mouse_2024']
    refs_10x_links_dict ={
      'human_2020': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz", 
      'human_2024': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz", 
      'mouse_2020': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz", 
      'mouse_2024': "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
    }
    
    names = config['genome']['species']['name']
    if isinstance(names, str):
      names = [names]
    
    assert np.isin(names, allowed_names), "unrecognized species name"
    
    # define mito strings for all names
    mito_strs = [
    "^MT-" if 'human_' in n else "^mt-"
    for n in names
    ]
    
    # update dictionary with mito strings
    for n, mito in zip(names, mito_strs):
      mito_str_dict.update({n:mito})
    
    # set defaults for decoys
    decoys = [False] * len(names)
    
    # change defaults for decoys if available
    if 'decoys' in config['genome']['species']:
      decoys = config['genome']['species']['decoys']
      if isinstance(decoys, bool):
        decoys = [decoys]
        
    # update dictionary with decoys
     for n, dcoy in zip(names, decoys):
       decoy_dict.update({n:dcoy})
      
    # download all required genomes from 10x
    for gnome in names:
      # get link
      link = refs_10x_links_dict[gnome]
      # get tarball 
      tball = link.rsplit('/', 1)[-1]
      
      gnome_dir = os.path.join(ref_dir, gnome)
      os.makedirs(gnome_dir, exist_ok = True)
      os.chdir(gnome_dir)
        
      # download reference into that dir
      subprocess.run(['wget', link], shell=True)
        
      # list tar contents
      tar_out= subprocess.run(f'tar -tzf {tball}', shell=True, capture_output=True, text=True)
      tball_fs = tar_out.stdout.splitlines()
        
      # look for genome.fa and genes.gtf files and extract them 
      for t in tball_fs:
       if t.endswith('/fasta/genome.fa') or t.endswith('/genes/genes.gtf'): # make this better (gtf file can be zipped)
        ext_spell = ['tar', '-xzf', tball, '--strip-components=1', '-C', cwd, t] # this will extract without including the name of the tarball
        subprocess.run(command, check=True)
      
      # append paths to fasta and gtf files to dictionary
      fasta_dict.update({gnome: os.path.join(gnome_dir, 'fasta/genome.fa')})
      gtf_dict.update({gnome: os.path.join(gnome_dir, 'genes/genes.gtf')})
      
       # delete the tarball
      #os.remove(tball)
    
    # check custom entries  
    if 'custom' in config['genome']:
      # check if all required entries are specified
      req = ['name', 'fasta', 'gtf', 'mito_str']
      for r in req:
        assert r in config['genome']['custom'], f"missing {r} in custom genome specification"
        
      cust_names = config['genome']['custom']['name']
      if isinstance(cust_names, str):
        cust_names =[cust_names]
        
      cust_fastas = config['genome']['custom']['fasta']
        if isinstance(cust_fastas, str):
        cust_fastas =[cust_fastas]
        
      cust_gtfs = config['genome']['custom']['gtf']
        if isinstance(cust_gtfs, str):
        cust_gtfs =[cust_gtfs]
        
      cust_mito_strs = config['genome']['custom']['mito_str']
        if isinstance(cust_mito_strs, str):
        cust_mito_strs =[cust_mito_strs]
      
      # get defaults for decoys
      cust_decoys = [False] * len(cust_names)
      
      if 'decoys' in config['genome']['custom']:
        cust_decoys = config['genome']['species']['decoys']
        if isinstance(cust_decoys), bool):
          cust_decoys = [cust_decoys]
        
      # check if all entries have the same length
      cust_ls_len = {len(lst) for lst in [cust_names, cust_fastas, cust_gtfs, cust_mito_strs, cust_decoys]}
      assert len(cust_ls_len) == 1, "custom parameters don't have the same lengths"
    
      # check if fasta and gtf files exist
      for fa in cust_fastas:
        assert os.path.isfile(fa), \
          f"file {fa} specified in configfile doesn't exist"
          
      for gtf in cust_gtfs:
        assert os.path.isfile(gtf), \
          f"file {gtf} specified in configfile doesn't exist"
      
      # update dictionaries with fasta and gtf files
      for n, f, g in zip(cust_names, cust_fastas, cust_gtfs):
        fasta_dict.update({n: f})
        gtf_dict.update({n: g})
    
  
    # make an output csv file with all specified parameters 
    # in scprocess_data dir make a file with setup params (.csv file) --> this should be the output of the first rule 
    # columns in this .csv file: genome_name, fasta_f, gft_f, decoy, mito_str
    params_df = pd.DataFrame([fasta_dict, gtf_dict, mito_str_dict, decoy_dict])
    params_df = params_df.T
    params_df.reset_index(inplace=True)
    params_df.columns = ['genome_name', 'fasta_f', 'gtf_f', 'mito_str', 'decoy']
    
    out_csv_f = os.path.join(scprocess_data_dir, 'setup_parameters.csv')
    params_df.to_csv(out_csv_f, index = False)
    
    return 

# function that downloads scprocess_data (everything else apart from alevin indices) 
def get_scprocess_data():
  # barcode whitelists
  # gmt pathway files
  # canonical marker lists
  # xgboost objects
  
    
        

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
  
  

    
