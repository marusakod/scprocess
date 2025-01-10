import warnings
import yaml
import pandas as pd
import os
import re
import glob
import subprocess
import numpy as np


def get_setup_parameters(config):

  # # open configfile
  # with open(configfile, "r") as stream:
  #  config      = yaml.safe_load(stream)

  # check that genome is in config
  assert "genome" in config, "genome not defined in config file"
  
  # check that either tenx or custom are in the config 
  assert 'tenx' in config['genome'] or 'custom' in config['genome'], \
    "either 'tenx' or 'custom' need to be specified under genome"
      
  # start lists with params for all genomes
  all_genome_names = []
  all_fasta_fs     = []
  all_gtf_fs       = []
  all_idx_dirs     = []
  all_mito_str     = []
  all_dcoys        = []
  all_rrnas        = []

  allowed_names = ['human_2020', 'human_2024', 'mouse_2020', 'mouse_2024']

  if('tenx' in config['genome']):
    
    names      = [entry['name'] for entry in config['genome']['tenx']]
    decoys     = [entry.get('decoys', None) for entry in config['genome']['tenx']]
    rrnas      = [entry.get('rrnas', None) for entry in config['genome']['tenx']]
    
    # check that all genome names are unique
    assert len(names) == len(set(names)), "Duplicated genome names are not allowed!"
    
    assert all(name in allowed_names for name in names), "unrecognized 10x genome name"
     
    # set defaults
    fasta_fs = ['None'] * len(names)
    gtf_fs = ['None'] * len(names)
    mito_strs = ['None'] * len(names)
    index_dirs = ['None'] * len(names)
    
    decoys = [True if d is None else d for d in decoys]
    rrnas  = [True if r is None else r for r in rrnas]
    
    #check that lenght of decoys is the same as names
    params_ls_len = {len(lst) for lst in [names, decoys, rrnas]}
    assert len(params_ls_len) == 1, "10x genome parameters don't have the same lengths"
        
    # update params lists
    all_genome_names.extend(names)
    all_fasta_fs.extend(fasta_fs)
    all_gtf_fs.extend(gtf_fs)
    all_idx_dirs.extend(index_dirs)
    all_mito_str.extend(mito_strs)
    all_dcoys.extend(decoys)
    all_rrnas.extend(rrnas)

    # check custom entries  
  if 'custom' in config['genome']:
    # check if all required entries are specified

    has_name_key      = ['name' in entry for entry in config['genome']['custom']]
    assert all(has_name_key), \
     "missing 'name' in custom genome specification"
    has_gtf_key       = ['gtf' in entry for entry in config['genome']['custom']]
    assert all(has_gtf_key), \
     "missing 'gtf' in custom genome specification"
    has_mito_str_key  = ['mito_str' in entry for entry in config['genome']['custom']]
    assert all(has_mito_str_key), \
      "missing 'mito_str' in custom genome specification"
    
    # check that either 'fasta' or 'index_dir' is present in each entry
    has_fasta_or_index_dir =  [('fasta' in entry or 'index dir' in entry) for entry in config['genome']['custom']]
    assert all(has_fasta_or_index_dir), \
        "Each custom genome requires either 'fasta' or 'index_dir'"
    
    # check for mutual exclusivity (don't allow both 'fasta' and 'index_dir')
    no_both_fasta_and_index_dir = [not('fasta' in entry and 'index_dir' in entry) for entry in config['genome']['custom']]
    assert all(no_both_fasta_and_index_dir), \
        "Custom genome cannot have both 'fasta' and 'index_dir' defined at the same time"

        
    cust_names     = [entry['name'] for entry in config['genome']['custom']]

    # check that all genome names are unique
    assert len(cust_names) == len(set(cust_names)), "Duplicated genome names are not allowed!"
    
    # check that custom names don't overlap with predefined names
    assert len(set(cust_names) & set(allowed_names)) == 0, \
      "Some names for custom genomes overlap names for 10x genomes; please define other names for custom genomes"

    cust_fastas =  [entry['fasta'] for entry in config['genome']['custom']]
    cust_gtfs =  [entry['gtf'] for entry in config['genome']['custom']]
    cust_index_dirs = [entry.get('index_dir', None) for entry in config['genome']['custom']]
    cust_mito_strs = [entry['mito_str'] for entry in config['genome']['custom']]
      
    # get decoys
    cust_decoys = [entry.get('decoys', None) for entry in config['genome']['custom']]
    cust_decoys = [True if d is None else d for d in cust_decoys]
    cust_rrnas  = ['None'] * len(cust_names)
    
        
    # check if all entries have the same length
    cust_params_ls_len = {len(lst) for lst in [cust_names, cust_fastas, cust_gtfs, cust_mito_strs, cust_decoys]}
    assert len(cust_params_ls_len) == 1, "custom parameters don't have the same lengths"
    
    # check if all specified fastas exist
    for fa in cust_fastas:
      if fa is not None:
        assert os.path.isfile(fa), \
          f"file {fa} specified in configfile doesn't exist"
      
    # check if all specified index directories are ok
    for idx_dir in cust_index_dirs:
      if idx_dir is not None:
         assert check_valid_index(idx_dir), \
          "Alevin index incomplete"

    # check if gtf files exist
    for gtf in cust_gtfs:
      assert os.path.isfile(gtf), \
        f"file {gtf} specified in configfile doesn't exist"
        
    # update params lists
    all_genome_names.extend(cust_names)
    all_fasta_fs.extend(cust_fastas)
    all_gtf_fs.extend(cust_gtfs)
    all_idx_dirs.extend(cust_index_dirs)
    all_mito_str.extend(cust_mito_strs)
    all_dcoys.extend(cust_decoys)
    all_rrnas.extend(cust_rrnas)

  # convert decoy list to strings
  all_dcoys = [str(d) for d in all_dcoys]
  # convert rrnas list to strings
  all_rrnas = [str(r) for r in all_rrnas]

  # check if reference for tutorial needs to be added
  if 'human_2024' not in all_genome_names:
    all_genome_names.append('human_2024')
    all_dcoys.append('True')
    all_rrnas.append('True')
    all_fasta_fs.append('None')
    all_gtf_fs.append('None')
    all_idx_dirs.append('None')
    all_mito_str.append('None')


   # convert all lists to a single string with commas between elements
  genome_names_str = ','.join(all_genome_names)
  fasta_fs_str     = ','.join(all_fasta_fs)
  idx_dirs_str     = ','.join(all_idx_dirs)
  gtf_fs_str       = ','.join(all_gtf_fs)
  mito_one_str     = ','.join(all_mito_str)
  decoys_str       = ','.join(all_dcoys)
  rrnas_str        = ','.join(all_rrnas)


  # return lists with all params
  return genome_names_str, fasta_fs_str, gtf_fs_str, idx_dirs_str, mito_one_str, decoys_str, rrnas_str

def check_valid_index(idx_path):
   
    # check if main directory exists
    if not os.path.isdir(idx_path):
        raise ValueError(f"The provided path '{idx_path}' is not a valid directory.")

    # define expected directories and files
    req_dirs = {
        "index_dir": os.path.join(idx_path, 'index'),
        "ref_dir": os.path.join(idx_path, 'ref')
    }
    req_fs = {
        "index_info.json": os.path.join(idx_path, 'index_info.json'),
        "simpleaf_index_log.json": os.path.join(idx_path, 'simpleaf_index_log.json'),
        "gene_id_to_name.tsv": os.path.join(req_dirs["ref_dir"], 'gene_id_to_name.tsv'),
        "roers_make-ref.json": os.path.join(req_dirs["ref_dir"], 'roers_make-ref.json'),
        "roers_ref.fa": os.path.join(req_dirs["ref_dir"], 'roers_ref.fa'),
        "t2g_3col.tsv": os.path.join(req_dirs["ref_dir"], 't2g_3col.tsv')
    }

    # Check required directories
    for dir_name, dir_path in req_dirs.items():
        if not os.path.isdir(dir_path):
            raise ValueError(f"The required directory '{dir_name}' is missing at '{idx_path}'.")

    # Check required files
    missing_fs = [file_name for file_name, file_path in req_fs.items() if not os.path.isfile(file_path)]
    if missing_fs:
        raise ValueError(f"The following required files are missing: {', '.join(missing_fs)}")

    return True

  
