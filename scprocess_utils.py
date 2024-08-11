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
def exclude_samples_without_fastq_files(FASTQ_DIR, SAMPLES):
 

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
  for v in ["proj_dir", "fastq_dir", "short_tag", "full_tag", "date_stamp", "your_name", "affiliation", "sample_metadata"]:
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
  AF_MITO_STR = setup_params.loc[setup_params['genome_name'] == SPECIES, 'mito_str'].values[0]

  # get af index directory and check if exists
  AF_HOME_DIR = os.path.join(scprocess_data_dir, 'alevin_fry_home') # check if this exists in scprocess script
  AF_INDEX_DIR = os.path.join(AF_HOME_DIR, SPECIES)
  assert os.path.isdir(AF_INDEX_DIR), \
    f"alevin index for {SPECIES} doesn't exist"
  
  # get gtf txt file, check that exists
  AF_GTF_DT_F = setup_params.loc[setup_params['genome_name'] == SPECIES, 'gtf_txt_f'].values[0]

  return SPECIES, AF_MITO_STR, AF_HOME_DIR, AF_INDEX_DIR, AF_GTF_DT_F, CHEMISTRY_F

# define cellbender parameters
def get_cellbender_parameters(config):
  # set default values
  DO_CELLBENDER                 = True
  CELLBENDER_VERSION            = 'v0.3.0'
  CELLBENDER_PROP_MAX_KEPT      = 0.9
  CUSTOM_CELLBENDER_PARAMS_F    = None
  FORCE_EXPECTED_CELLS          = None
  FORCE_TOTAL_DROPLETS_INCLUDED = None
  FORCE_LOW_COUNT_THRESHOLD     = None
  CELLBENDER_LEARNING_RATE      = 1e-4

  # change defaults if specified
  if ('cellbender' in config) and (config['cellbender'] is not None):
    if 'do_cellbender' in config['cellbender']:
      DO_CELLBENDER                 = config['cellbender']['do_cellbender']
    if 'version' in config['cellbender']:
      CELLBENDER_VERSION            = config['cellbender']['version']
    if 'cb_max_prop_kept' in config['cellbender']:
      CELLBENDER_PROP_MAX_KEPT      = config['cellbender']['cb_max_prop_kept']
    if 'custom_cellbender_params' in config['cellbender']:
      CUSTOM_CELLBENDER_PARAMS_F    = config['cellbender']['custom_cellbender_params']
    if 'force_expected_cells' in config['cellbender']:
      FORCE_EXPECTED_CELLS          = config['cellbender']['force_expected_cells']
    if 'force_total_droplets_included' in config['cellbender']:
      FORCE_TOTAL_DROPLETS_INCLUDED = config['cellbender']['force_total_droplets_included']
    if 'force_low_count_threshold' in config['cellbender']:
      FORCE_LOW_COUNT_THRESHOLD     = config['cellbender']['force_low_count_threshold']
    if 'learning_rate' in config['cellbender']:
      CELLBENDER_LEARNING_RATE      = config['cellbender']['learning_rate']

  # get cellbender image
  if CELLBENDER_VERSION == 'v0.3.0':
    CELLBENDER_IMAGE              = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.0'
  elif CELLBENDER_VERSION == 'v0.2.0':
    CELLBENDER_IMAGE              = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.2.0'
  else:
    raise ValueError(f"selected cellbender version {CELLBENDER_VERSION} not supported")

  # check for a dumb combination of parameter settings
  if not DO_CELLBENDER and CUSTOM_CELLBENDER_PARAMS_F is not None:
    print("WARNING: do_cellbender is set to False, and a custom cellbender parameters file is specified. Custom file will be ignored.")
    CUSTOM_CELLBENDER_PARAMS_F = None

  # some checks on custom parameters for cellbender
  if CUSTOM_CELLBENDER_PARAMS_F is not None: 
    # if a custom file is specified, and do_cellbender is False, give a warning then skip this bit
    # does the file exist?
    assert os.path.exists(CUSTOM_CELLBENDER_PARAMS_F), \
      f"specified path {CUSTOM_CELLBENDER_PARAMS_F} does not exist"

    # do the samples match the metadata?
    params_df   = pd.read_csv(CUSTOM_CELLBENDER_PARAMS_F)
    meta_df     = pd.read_csv(config["sample_metadata"])
    assert all(params_df.columns.values == \
      ['sample_id','total_droplets_included','expected_cells','low_count_threshold','learning_rate']), \
      f'column names in {CUSTOM_CELLBENDER_PARAMS_F} are not correct'
    assert set(params_df.sample_id.values) <= set(meta_df.sample_id.values)
  
  return CELLBENDER_IMAGE, CELLBENDER_PROP_MAX_KEPT, DO_CELLBENDER, CUSTOM_CELLBENDER_PARAMS_F, \
    FORCE_EXPECTED_CELLS, FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE


