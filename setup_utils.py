import warnings
import yaml
import pandas as pd
import os
import re
import glob
import subprocess
import numpy as np


def get_setup_parameters(config):

  # open configfile
  #with open(configfile, "r") as stream:
  #  config      = yaml.safe_load(stream)

  # check that genome is in config
  assert "genome" in config, "genome not defined in config file"
  
  # check that either tenx or custom are in the config 
  assert 'tenx' in config['genome'] or 'custom' in config['genome'], \
    "either 'tenx' or 'custom' need to be specified under genome"
      
  # start lists with params for all genomes
  all_genome_names = []
  all_fasta_fs = []
  all_gtf_fs = []
  all_mito_str = []
  all_dcoys= []

  allowed_names = ['human_2020', 'human_2024', 'mouse_2020', 'mouse_2024']

  if('tenx' in config['genome']):
    
    names = config['genome']['tenx']['name']
    if isinstance(names, str):
      names = [names]

    # check that all genome names are unique
    assert len(names) == len(set(names)), "Duplicated genome names are not allowed!"
    
    assert np.isin(names, allowed_names), "unrecognized 10x genome name"
     
    # set defaults
    fasta_fs = [None] * len(names)
    gtf_fs = [None] * len(names)
    mito_strs = [None] * len(names)
    decoys = [False] * len(names)
    
    # change defaults for decoys if available
    if 'decoys' in config['genome']['tenx']:
      decoys = config['genome']['tenx']['decoys']
      if isinstance(decoys, bool):
        decoys = [decoys]
      
      #check that lenght of decoys is the same as names
      params_ls_len = {len(lst) for lst in [names, decoys]}
      assert len(params_ls_len) == 1, "10x genome parameters don't have the same lengths"
        
    # update params lists
    all_genome_names.extend(names)
    all_fasta_fs.extend(fasta_fs)
    all_gtf_fs.extend(gtf_fs)
    all_mito_str.extend(mito_strs)
    all_dcoys.extend(decoys)

    # check custom entries  
  if 'custom' in config['genome']:
      # check if all required entries are specified
    req = ['name', 'fasta', 'gtf', 'mito_str']
    for r in req:
      assert r in config['genome']['custom'], f"missing {r} in custom genome specification"
        
    cust_names = config['genome']['custom']['name']
    if isinstance(cust_names, str):
      cust_names =[cust_names]

    # check that all genome names are unique
    assert len(cust_names) == len(set(cust_names)), "Duplicated genome names are not allowed!"
    
    # check that custom names don't overlap with predefined names
    assert len(set(cust_names) & set(allowed_names)) == 0, \
      "Some names for custom genomes overlap names for 10x genomes; please define other names for custom genomes"

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
      cust_decoys = config['genome']['custom']['decoys']
      if isinstance(cust_decoys, bool):
        cust_decoys = [cust_decoys]
        
      # check if all entries have the same length
    cust_params_ls_len = {len(lst) for lst in [cust_names, cust_fastas, cust_gtfs, cust_mito_strs, cust_decoys]}
    assert len(cust_params_ls_len) == 1, "custom parameters don't have the same lengths"
    
      # check if fasta and gtf files exist
    for fa in cust_fastas:
      assert os.path.isfile(fa), \
        f"file {fa} specified in configfile doesn't exist"
          
    for gtf in cust_gtfs:
      assert os.path.isfile(gtf), \
        f"file {gtf} specified in configfile doesn't exist"
        

    # update params lists
    all_genome_names.extend(cust_names)
    all_fasta_fs.extend(cust_fastas)
    all_gtf_fs.extend(cust_gtfs)
    all_mito_str.extend(cust_mito_strs)
    all_dcoys.extend(cust_decoys)

  # return lists with all params
  return all_genome_names, all_fasta_fs, all_gtf_fs, all_mito_str, all_dcoys
  