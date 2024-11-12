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
    
    names = [entry['name'] for entry in config['genome']['tenx']]
    decoys = [entry.get('decoys', None) for entry in config['genome']['tenx']]

    
    # check that all genome names are unique
    assert len(names) == len(set(names)), "Duplicated genome names are not allowed!"
    
    assert all(name in allowed_names for name in names), "unrecognized 10x genome name"
     
    # set defaults
    fasta_fs = ['None'] * len(names)
    gtf_fs = ['None'] * len(names)
    mito_strs = ['None'] * len(names)
    
    decoys = [True if d is None else d for d in decoys]
    
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

    has_name_key      = ['name' in entry for entry in config['genome']['custom']]
    assert all(has_name_key), \
     "missing 'name' in custom genome specification"
    has_fasta_key     = ['fasta' in entry for entry in config['genome']['custom']]
    assert all(has_fasta_key), \
     "missing 'fasta' in custom genome specification"
    has_gtf_key       = ['gtf' in entry for entry in config['genome']['custom']]
    assert all(has_gtf_key), \
     "missing 'gtf' in custom genome specification"
    has_mito_str_key  = ['mito_str' in entry for entry in config['genome']['custom']]
    assert all(has_mito_str_key), \
      "missing 'mito_str' in custom genome specification"
        
    cust_names = [entry['name'] for entry in config['genome']['custom']]

    # check that all genome names are unique
    assert len(cust_names) == len(set(cust_names)), "Duplicated genome names are not allowed!"
    
    # check that custom names don't overlap with predefined names
    assert len(set(cust_names) & set(allowed_names)) == 0, \
      "Some names for custom genomes overlap names for 10x genomes; please define other names for custom genomes"

    cust_fastas =  [entry['fasta'] for entry in config['genome']['custom']]
    cust_gtfs =  [entry['gtf'] for entry in config['genome']['custom']]
    cust_mito_strs = [entry['mito_str'] for entry in config['genome']['custom']]
      
    # get decoys
    cust_decoys = [entry.get('decoys', None) for entry in config['genome']['custom']]
    cust_decoys = [True if d is None else d for d in cust_decoys]
    
        
    # check if all entries have the same length
    cust_params_ls_len = {len(lst) for lst in [cust_names, cust_fastas, cust_gtfs, cust_mito_strs, cust_decoys]}
    assert len(cust_params_ls_len) == 1, "custom parameters don't have the same lengths"
    
      # check if fasta and gtf files exist
    # for fa in cust_fastas:
    #   assert os.path.isfile(fa), \
    #     f"file {fa} specified in configfile doesn't exist"
          
    # for gtf in cust_gtfs:
    #   assert os.path.isfile(gtf), \
    #     f"file {gtf} specified in configfile doesn't exist"
        

    # update params lists
    all_genome_names.extend(cust_names)
    all_fasta_fs.extend(cust_fastas)
    all_gtf_fs.extend(cust_gtfs)
    all_mito_str.extend(cust_mito_strs)
    all_dcoys.extend(cust_decoys)

  # convert decoy list to strings
  all_dcoys = [str(b) for b in all_dcoys]

  # check if reference for tutorial needs to be added
  if 'human_2024' not in all_genome_names:
    all_genome_names.append('human_2024')
    all_dcoys.append('True')


   # convert all lists to a single string with commas between elements
  genome_names_str = ','.join(all_genome_names)
  fasta_fs_str = ','.join(all_fasta_fs)
  gtf_fs_str = ','.join(all_gtf_fs)
  mito_one_str = ','.join(all_mito_str)
  decoys_str = ','.join(all_dcoys)


  # return lists with all params
  return genome_names_str, fasta_fs_str, gtf_fs_str, mito_one_str, decoys_str



  