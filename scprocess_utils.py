# load modules
import warnings
import yaml
import pandas as pd
import os
import re
import glob
import datetime
import subprocess



# find fastq files for a sample
def find_fastq_files(fastqs_dir, sample, read):
  # get all files
  all_fs    = glob.glob(f"{fastqs_dir}/*{sample}*")

  # get all reads
  re_R1     = re.compile('.*_R1.*\\.fastq(\\.gz)?$')
  R1_fs     = [ f for f in all_fs if re_R1.match(f) ]
  R1_fs     = sorted(R1_fs)
  re_R2     = re.compile('.*_R2.*\\.fastq(\\.gz)?$')
  R2_fs     = [ f for f in all_fs if re_R2.match(f) ]
  R2_fs     = sorted(R2_fs)

  # check they match
  R1_chk    = [ re.sub("_R1", "", f) for f in R1_fs ]
  R2_chk    = [ re.sub("_R2", "", f) for f in R2_fs ]
  assert R1_chk == R2_chk, "R1 and R2 fastq files do not match for " + sample

  if read == "R1":
    sel_fs = R1_fs
  elif read == "R2":
    sel_fs = R2_fs

  return sel_fs


# function to exclude samples without valid fastq files
def exclude_samples_without_fastq_files(FASTQ_DIR, PROJ_DIR, SAMPLES):
 

  # get fastq files for each sample
  chk_samples   = []
  for sample in SAMPLES:
    R1_fs = find_fastq_files(FASTQ_DIR, sample, "R1")
    if len(R1_fs) > 0:
      chk_samples.append(sample)
    else:
      print(f"WARNING: no fastq files found for sample {sample}; excluded")

  return chk_samples


# get list of samples
def get_project_parameters(config):
  # check expected variables are in the config file
  for v in ["proj_dir", "fastq_dir" "short_tag", "full_tag", "date_stamp", "your_name", "affiliation", "sample_metadata"]:
    assert v in config, f"{v} not in config file"

  ## what is specified in config directory?
  PROJ_DIR      = config["proj_dir"]
  FASTQ_DIR     = config["fastq_dir"]
  SHORT_TAG     = config["short_tag"]
  FULL_TAG      = config["full_tag"]
  YOUR_NAME     = config["your_name"]
  AFFILIATION   = config["affiliation"]
  DATE_STAMP    = config["date_stamp"]

  # check whether date is given as datetime object
  if isinstance(DATE_STAMP, datetime.date):
    # if so, convert to string
    DATE_STAMP  = DATE_STAMP.strftime("%Y-%m-%d")

  date_regex    = re.compile("^20[0-9]{2}-[0-9]{2}-[0-9]{2}$")
  assert date_regex.match(DATE_STAMP), f"{DATE_STAMP} does not match date format YYYY-MM-DD"
  
  ## get samples
  METADATA_F    = config["sample_metadata"]
  if not os.path.isfile(METADATA_F):
    METADATA_F    = os.path.join(PROJ_DIR, METADATA_F)
  assert os.path.isfile(METADATA_F), f"sample metadata file {METADATA_F} does not exist"

  samples_df    = pd.read_csv( METADATA_F )
  assert "sample_id" in samples_df.columns, "sample_id not in sample metadata"
  SAMPLES       = samples_df["sample_id"].dropna().tolist()

  # remove some samples
  EXC_SAMPLES   = None
  if ('exclude_samples' in config) and (config['exclude_samples'] is not None):
    EXC_SAMPLES = config["exclude_samples"]
    for s in EXC_SAMPLES:
      if s not in SAMPLES:
        warnings.warn(f"sample {s} specified in exclude_samples but not in metadata file", UserWarning)
    to_keep       = set(SAMPLES) - set(EXC_SAMPLES)
    SAMPLES       = [s for s in SAMPLES if s in to_keep]

  # sort out metadata variables
  METADATA_VARS = []
  if ('metadata_vars' in config) and (config['metadata_vars'] is not None):
    METADATA_VARS  = config["metadata_vars"]
    for var in METADATA_VARS:
      assert var in samples_df.columns, f"{var} not in sample metadata"

  return PROJ_DIR, FASTQ_DIR, SHORT_TAG, FULL_TAG, YOUR_NAME, AFFILIATION, METADATA_F, METADATA_VARS, EXC_SAMPLES, SAMPLES, DATE_STAMP


# define alevin parameters
def get_alevin_parameters(config, scprocess_data_dir):

  # check that alevin is in config
  assert "alevin" in config, \
    "alevin not defined in config file"

  setup_params_f  = os.path.join(scprocess_data_dir, 'setup_parameters.csv')

  # from setup_parameters.csv get valid values for species
  setup_params= pd.read_csv(setup_params_f)
  valid_species = setup_params['genome_name'].tolist()
  
  # check that species is defined in alevin and is valid
  assert "species" in config['alevin'], \
    "species not defined in configfile"
  
  SPECIES = config['alevin']['species']
  
  assert SPECIES in valid_species, \
   f"species {SPECIES} not defined"
  
  # get chemistry file if defined
  CHEMISTRY_F = None
  if 'chemistry' in config['alevin']:
    CHEMISTRY_F = config['alevin']['chemistry']
    assert os.path.isfile(CHEMISTRY_F), \
     "Custom chemistry file doesn't exist"

  # get mito strings from setup params
  AF_MITO_STR = setup_params[setup_params['genome_name'] == SPECIES, 'mito_str']

  # get af index directory and check if exists
  AF_HOME_DIR = os.path.join(scprocess_data_dir, 'alevin_fry_home') # check if this exists in scprocess script
  AF_INDEX_DIR = os.path.join(AF_HOME_DIR, SPECIES)
  assert os.path.isdir(AF_INDEX_DIR), \
    f"alevin index for {SPECIES} doesn't exist"
  
  # get gtf txt file, check that exists
  AF_GTF_DT_F = setup_params[setup_params['genome_name'] == SPECIES, 'gtf_txt_f']

  return SPECIES, AF_MITO_STR, AF_HOME_DIR, AF_INDEX_DIR, AF_GTF_DT_F, CHEMISTRY_F

