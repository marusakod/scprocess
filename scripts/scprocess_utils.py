# load modules
import os
import sys
import re
import pathlib
import warnings
import yaml
import pandas as pd
import csv
import math
import glob
import gzip
import datetime
import subprocess
import json
import jsonschema

### much setup

MB_PER_GB   = 1024

# do some checks of setup before running scprocess
def check_setup_before_running_scprocess(scprocess_dir, extraargs):
  # check that SCPROCESS_DATA_DIR exists
  scprocess_data_dir = os.getenv('SCPROCESS_DATA_DIR')
  if not scprocess_data_dir:
    raise ValueError('SCPROCESS_DATA_DIR is not defined in .bashrc')
  if not os.path.isdir(scprocess_data_dir):
    raise FileNotFoundError("SCPROCESS_DATA_DIR is not a directory")
  
  # check that spcrocess_data_dir has some files
  scsetup_dirs = ['cellranger_ref', 'gmt_pathways', 'marker_genes', 'xgboost', 'alevin_fry_home']
  scsetup_full_dirs = [os.path.join(scprocess_data_dir, d) for d in scsetup_dirs]
  for d in scsetup_full_dirs:
    if not os.path.isdir(d):
      raise FileNotFoundError(f"A directory is missing in {scprocess_data_dir}; consider (re)running setup.\nMissing directory:\n{d}")
  
  # check that setup csv exists
  scsetup_csv = os.path.join(scprocess_data_dir,'index_parameters.csv')
  if not os.path.isfile(scsetup_csv):
    raise FileNotFoundError(f"{scsetup_csv} is missing; consider (re)running setup.")
  
  # check if cluster profile is defined
  setup_configfile = os.path.join(scprocess_data_dir, 'scprocess_setup.yaml')
  if not os.path.exists(setup_configfile):
    raise FileNotFoundError(f"scprocess_setup.yaml does not exist in {scprocess_data_dir}")

  # load setup config file
  with open(setup_configfile, "r") as stream:
    setup_config      = yaml.safe_load(stream)

  # if there is a cluster profile check it
  if ('profile' in setup_config) and setup_config['profile'] is not None:
    # check if profile exists
    profile     = setup_config['profile']
    profile_dir = os.path.join(scprocess_dir, 'profiles', profile)
    profile_f   = os.path.join(profile_dir, 'config.yaml')
    if not os.path.isfile(profile_f):
      raise FileNotFoundError(f"cluster configuration file {profile_f} does not exist")
    
    # if ok, add profile to snakemake call
    extraargs = extraargs + ' --workflow-profile ' + profile_dir

  return extraargs


# get validated config yaml list
def get_config_validated_against_schema(config_path, schema_f):
  try:
    # open schema file
    with open(schema_f, "r") as f:
      schema  = json.load(f)

    # open and parse the yaml file
    with open(config_path, "r") as f:
      CONFIG  = yaml.safe_load(f)
    
    # validate the parsed yaml config against the json schema
    jsonschema.validate( instance = CONFIG, schema = schema )
    print("  config file has correct format")

  except json.decoder.JSONDecodeError as e:
    print(f"problem with scprocess schema file:\n  {str(e)}")
    sys.exit(1)
  except jsonschema.ValidationError as e:
    print(f"problem with your config file:\n  {str(e)}")
    sys.exit(1)
  except jsonschema.SchemaError as e:
    print(f"schema error: {str(e)}")
    sys.exit(1)
  except yaml.YAMLError as e:
    print(f"YAML Parsing error: {e}")
    sys.exit(1)
  except Exception as e:
    print(f"Unexpected error: {e}")
    sys.exit(1)

  return CONFIG


def get_default_config_from_schema(schema_f):
  # mess about w schema path
  schema_p  = pathlib.Path(schema_f)
  path_bits = list(schema_p.parts)
  if path_bits[0] == '..':
    path_bits = path_bits[1:]
  schema_p  = str(pathlib.Path(*path_bits))
  with open(schema_p, 'r') as f:
    schema = json.load(f)

  # extract defaults from the schema
  default_config = {}
  for key, props in schema.get('properties', {}).items():
    # if the key has a top-level default (uncommon for object types)
    if 'default' in props:
      default_config[key] = props['default']
    # if the key is an object, recursively extract its property defaults
    elif props.get('type') == 'object' and 'properties' in props:
      # create a nested dictionary of defaults for this section
      section_defaults = {}
      for sub_key, sub_props in props['properties'].items():
        if 'default' in sub_props:
          section_defaults[sub_key] = sub_props['default']
      
      if section_defaults:
        default_config[key] = section_defaults

  return default_config


# check proj dir is wflowr
def check_proj_dir_is_wflowr(proj_dir):
  # check that proj_dir is a workflowr directory 
  wflowr_fs_ls = ['_workflowr.yml', '.gitignore', '.Rprofile', '.gitattributes',
    'analysis/_site.yml', 'analysis/about.Rmd', 'analysis/index.Rmd', 'analysis/license.Rmd', 
    'public/.nojekyll']
  wflowr_fs_full_ls = [os.path.join(proj_dir, f) for f in wflowr_fs_ls]
  for f in wflowr_fs_full_ls:
    if not os.path.isfile(f):
      raise FileNotFoundError(f"proj_dir {CONFIG['proj_dir']} is not a workflowr project; you can create a workflowr project using `scprocess newproj`")


# check sample metadata
def load_and_check_sample_metadata(CONFIG):
  # unpack
  proj_dir        = CONFIG['project']['proj_dir']
  sample_metadata = CONFIG['project']['sample_metadata']
  if 'multiplexing' in CONFIG:
    demux_type      = CONFIG['multiplexing']['demux_type']
  else:
    demux_type      = None

  # check sample metadata file exists
  if not os.path.isabs(sample_metadata):
    sample_metadata = os.path.join(proj_dir, sample_metadata)
  if not os.path.isfile(sample_metadata):
    raise FileNotFoundError(f"sample metadata {sample_metadata} is not a valid file")

  # load it
  sample_df   = pd.read_csv(sample_metadata)
  
  # define correct sample variable, check it's present
  if demux_type is None:
    sample_var  = "sample_id"
  else:
    sample_var  = "pool_id"
  if sample_var not in sample_df.columns:
   raise KeyError(f"'{sample_var}' not present in sample metadata")

  # some checks for multiplexing
  if demux_type == "hto":
    # check that hto fastq path is present
    if "hto_id" not in sample_df.columns:
      raise KeyError("'hto_id' not present in sample metadata")
    
    # check that hto fastq path is present
    hto_fastq_dir = CONFIG['multiplexing']['fastq_dir'] 
    if not os.path.isabs(hto_fastq_dir):
      hto_fastq_dir = os.path.join(proj_dir, hto_fastq_dir)
    if not os.path.isdir(hto_fastq_dir):
      raise FileNotFoundError(f"{hto_fastq_dir} is not a valid directory")

  # check that sample_id values are unique
  sample_ns   = sample_df[ "sample_id" ].value_counts()
  if not all(sample_ns == 1):
    raise ValueError("'sample_id' values in metadata csv not unique")

  # check columns of sample_df
  if any(' ' in col for col in sample_df.columns):
    raise ValueError("some column names in metadata csv contain spaces.")

  # check whether metadata_vars are present in sample metadata
  if "metadata_vars" in CONFIG["project"]:
    for var in CONFIG["project"]["metadata_vars"]:
      if var not in sample_df.columns:
        raise KeyError(f"metadata_var '{var}' not present in sample metadata")

  Warning("add check for 'exclude' here")
  return sample_df, sample_var


# get ok list of samples
def get_and_check_samples(sample_df, sample_var, CONFIG):
  # get samples
  SAMPLES  = sample_df[sample_var].tolist()

  # check that samples don't have '_R1' or '_R2' in their names 
  if any('_R1' in s or '_R2' in s for s in SAMPLES):
    raise ValueError(f"One or more {sample_var} values contain '_R1' or '_R2'. Please ensure all values exclude these substrings.")

  # compare every sample id to every other sample id
  smpl_overlaps = []
  for i, sample in enumerate(SAMPLES):
    for j, other_sample in enumerate(SAMPLES):
      if i != j and sample in other_sample:
        smpl_overlaps.append((sample, other_sample))

  # if samples overlap print error
  if smpl_overlaps:
    msg = "The following sample_id values are problematic (one is a subset of the other):\n"
    for sample, other_sample in smpl_overlaps:
      msg += f"  - '{sample}' is a subset of '{other_sample}'\n"
    raise ValueError(msg)

  # check that fastq files exist for samples
  sample_fastqs = get_sample_fastqs(CONFIG, SAMPLES, is_hto = False)
  SAMPLES       = [ s for s in SAMPLES if s in sample_fastqs ]

  # check that hto fastq files exist for samples
  demux_type    = "none"
  if ('multiplexing' in CONFIG) and CONFIG['multiplexing'] is not None:
    demux_type      = CONFIG['multiplexing']['demux_type']
  if demux_type == "hto":
    sample_htos   = get_sample_fastqs(CONFIG, is_hto = True )
    SAMPLES       = [ s for s in SAMPLES if s in sample_htos ]

  # check that there's still something left
  if not len(SAMPLES) > 0:
    raise ValueError("No matching fastq files found.")

  return SAMPLES


# get all fastq files
def get_sample_fastqs(CONFIG, SAMPLES, is_hto = False):
  # get place to look for fastq files
  if is_hto:
    tmp_ls      = CONFIG['multiplexing']
  else:
    tmp_ls      = CONFIG['project']
  if "fastq_dir" in tmp_ls:
    fastq_dir   = tmp_ls['fastq_dir']
    arv_uuid    = None
  else:
    fastq_dir   = None
    arv_uuid    = tmp_ls['arv_uuid']

  # get 
  if fastq_dir is not None:
    fastq_dict  = _list_fastq_files(fastq_dir, CONFIG)
  elif arv_uuid is not None:
    fastq_dict  = _list_fastq_files_arvados(arv_uuid)

  # get fastq files for each sample
  fastq_fs      = fastq_dict["fastqs"]
  sample_fastqs = {}
  for sample in SAMPLES:
    # get R1 and R2 files matching each sample
    R1_regex      = rf".*{sample}.*(_|\.)R1.*\.fastq\.gz"
    R1_fs         = [f for f in fastq_fs if re.match(R1_regex, f) ]
    R2_regex      = rf".*{sample}.*(_|\.)R2.*\.fastq\.gz"
    R2_fs         = [f for f in fastq_fs if re.match(R2_regex, f) ]

    # check have full set of files
    check_R1      = [re.sub(r'(?<=(_|\.))R1', 'R0', f) for f in R1_fs]
    check_R2      = [re.sub(r'(?<=(_|\.))R2', 'R0', f) for f in R2_fs]
    if len(R1_fs) == 0:
      import pdb; pdb.set_trace()
      print(f"  WARNING: no {[ "hto " if is_hto else ""]}fastq files found for sample {sample}; excluded.")
    elif set(check_R1) != set(check_R2):
      print(f"  WARNING: {[ "hto " if is_hto else ""]}fastq files found for sample {sample} but R1 and R2 don't match; excluded.")
    else:
      sample_fastqs[sample] = {"where": fastq_dict["where"], "R1_fs": R1_fs, "R2_fs": R2_fs}

  return sample_fastqs


# get all fastq files in directory
def _list_fastq_files(fastq_dir, CONFIG):
  # check directory exists
  fastq_dir   = _check_directory_within_project(fastq_dir, CONFIG['project']['proj_dir'])

  # get all files
  all_fs      = os.listdir(fastq_dir)

  # filter to just fastqs
  fastq_fs    = [ f for f in all_fs if re.match(r".+\.fastq\.gz", f) ]

  return { "where": fastq_dir, "fastqs": fastq_fs }


# helper function to check for existence of directory
def _check_directory_within_project(input_dir, proj_dir):
  # if boring do nothing
  if input_dir is None:
    return None

  # otherwise check if relative or absolute path, and check exists
  if not os.path.isabs(input_dir):
    input_dir = os.path.join(proj_dir, input_dir)
  if not os.path.isdir(input_dir):
    raise FileNotFoundError(f"{input_dir} if not a valid directory")

  return input_dir


# get all fastq files in arvados uuid
def _list_fastq_files_arvados(arv_uuid):
  # import relevant packages
  import arvados
  import collections
  import pathlib

  # set up
  arv_token   = os.environ["ARVADOS_API_TOKEN"]

  # set arvados API info
  arv_client  = arvados.api('v1', host = 'api.arkau.roche.com',
    token = arv_token, insecure = True, num_retries = 2 )
  arv_colln   = arvados.collection.Collection(arv_uuid, arv_client)

  # get all files within this uuid
  stream_q    = collections.deque([pathlib.PurePosixPath('.')])
  arv_files   = []
  while stream_q:
    stream_path = stream_q.popleft()
    tmp_colln   = arv_colln.find(str(stream_path))
    for item_name in tmp_colln:
      try:
        my_file = tmp_colln.open(item_name)
        arv_files.append( os.path.join(str(stream_path), item_name) )
      except IsADirectoryError:
        # item_name refers to a stream. Queue it to walk later.
        stream_q.append(stream_path / item_name)
        continue

  # filter to just fastqs
  fastq_fs    = [ f for f in arv_files if re.match(r".+\.fastq\.gz", f) ]

  return { "where": arv_uuid, "fastqs": fastq_fs }


### parameter definitions

def get_config(configfile, schema_f):
  # check configfile exists
  if not os.path.exists(configfile):
    raise FileNotFoundError(f"Config file {configfile} does not exist")
  
  # open configfile
  with open(configfile, "r") as stream:
    config      = yaml.safe_load(stream)

  # check against schema
  with open(schema_f, "r") as f:
    schema      = json.load(f)

  # validate config against schema
  try:
    jsonschema.validate(instance=config, schema=schema)
    print("YAML file is valid!")
  except jsonschema.exceptions.ValidationError as e:
    print(f"Validation error: {e.message}")

  # check each set of parameters
  _check_project_parameters(config)
  os.path.isdir(config["proj_dir"]), \
    f"proj_dir {config['proj_dir']} is not a directory"


def _check_project_parameters(config):
  # check fastq vs arvados
  has_fastq     = "fastq_dir" in config
  has_arv_uuid  = "arv_uuid" in config
  if has_fastq + has_arv_uuid != 1:
    KeyError('config file must contain exactly one of "fastq_dir" and "arv_uuid"')

  # define these
  if has_fastq and not has_arv_uuid:
    FASTQ_DIR     = config["fastq_dir"]
    ARV_UUID      = None
  if not has_fastq and has_arv_uuid:
    FASTQ_DIR     = None
    ARV_UUID      = config["arv_uuid"]

  # check if selected species is valid
  index_params_f  = os.path.join(scprocess_data_dir, 'index_parameters.csv')

  # from index_parameters.csv get valid values for species
  index_params      = pd.read_csv(index_params_f)
  valid_species     = index_params['genome_name'].tolist()
  valid_species_str = ', '.join(valid_species)
  if not SPECIES in valid_species:
    raise ValueError(f"species {SPECIES} not defined. Valid values are {valid_species_str}")

  # check whether date is given as datetime object
  if isinstance(DATE_STAMP, datetime.date):
    # if so, convert to string
    DATE_STAMP  = DATE_STAMP.strftime("%Y-%m-%d")

  date_regex    = re.compile("^20[0-9]{2}-[0-9]{2}-[0-9]{2}$")
  if not date_regex.match(DATE_STAMP):
    raise ValueError(f"{DATE_STAMP} does not match date format YYYY-MM-DD")

  # get fastqs
  if FASTQ_DIR is not None:
    if not os.path.isabs(FASTQ_DIR):
      FASTQ_DIR = os.path.join(PROJ_DIR, FASTQ_DIR)
  
  # get samples
  METADATA_F    = config["sample_metadata"]
  if not os.path.isabs(METADATA_F):
    METADATA_F = os.path.join(PROJ_DIR, METADATA_F)

  samples_df    = pd.read_csv( METADATA_F )
  SAMPLES       = samples_df["sample_id"].dropna().tolist()

  # get multiplexing params
  DEMUX_TYPE, HTO_FASTQ_DIR, FEATURE_REF, DEMUX_F, BATCH_VAR, EXC_POOLS, POOL_IDS, \
    POOL_SAMPLES, SAMPLE_MAPPING = get_multiplexing_parameters(config, PROJ_DIR, samples_df)

  # define sample variable for alevin, ambient, doublets and sample mapping dictionary
  if DEMUX_TYPE != "none":
    SAMPLE_VAR = "pool_id"
    # exclude samples in excluded pools
    SAMPLES = POOL_SAMPLES
  else:
    SAMPLE_VAR = "sample_id"
  
  # remove some samples
  EXC_SAMPLES   = None
  if ('exclude' in config) and (config['exclude'] is not None):
    if ('sample_id' in config['exclude']) and (config['exclude']['sample_id'] is not None):
      EXC_SAMPLES = config['exclude']["sample_id"]
      for s in EXC_SAMPLES:
        if s not in samples_df["sample_id"].dropna().tolist():
          warnings.warn(f"sample {s} specified in exclude_samples but not in metadata file", UserWarning)
      to_keep       = set(SAMPLES) - set(EXC_SAMPLES)
      SAMPLES       = [s for s in SAMPLES if s in to_keep]


  # check custom sample config file if exists
  CUSTOM_SAMPLE_PARAMS_F = None
  if ('custom_sample_params' in config) and (config['custom_sample_params'] is not None):
    CUSTOM_SAMPLE_PARAMS_F = config["custom_sample_params"]
    # check if exists
    if not os.path.isabs(CUSTOM_SAMPLE_PARAMS_F):
      CUSTOM_SAMPLE_PARAMS_F = os.path.join(PROJ_DIR, CUSTOM_SAMPLE_PARAMS_F)
    assert os.path.isfile(CUSTOM_SAMPLE_PARAMS_F)
    # open the file and check if all samples can be found in SAMPLES or POOL_IDS
    with open(CUSTOM_SAMPLE_PARAMS_F) as f:
      custom_smpl_params = yaml.load(f, Loader=yaml.FullLoader)
      custom_smpls = list(custom_smpl_params.keys())
      if DEMUX_TYPE != "none":
        for s in custom_smpls:
          assert s in POOL_IDS, f"{s} in custom_sample_params file doesn't match any pool_id values in sample_metadata"
      else:
        for s in custom_smpls:
          assert s in SAMPLES, f"{s} in custom_sample_params file doesn't match any sample_id values in sample_metadata"


  # sort out metadata variables
  METADATA_VARS = []
  if ('metadata_vars' in config) and (config['metadata_vars'] is not None):
    METADATA_VARS  = config["metadata_vars"]
    for var in METADATA_VARS:
      assert var in samples_df.columns, f"{var} not in sample metadata"
      # check that there are less than 10 unique values (otherwise probably not a categorical variable)
      var_vals  = samples_df[var].tolist()
      if all(isinstance(x, float) and math.isnan(x) for x in var_vals):
        continue
      uniq_vals = len(set(var_vals))
      assert uniq_vals <= 10, \
        f"{var} variable has more than 10 unique values"
  

  return PROJ_DIR, FASTQ_DIR, SHORT_TAG, FULL_TAG, YOUR_NAME, AFFILIATION, METADATA_F, \
    METADATA_VARS, EXC_SAMPLES, SAMPLES, DATE_STAMP, CUSTOM_SAMPLE_PARAMS_F, SPECIES, \
    DEMUX_TYPE, HTO_FASTQ_DIR, FEATURE_REF, DEMUX_F, BATCH_VAR, EXC_POOLS, POOL_IDS, \
    SAMPLE_VAR, SAMPLE_MAPPING


# check parameters for qc
def check_qc_parameters(config):
  # define some hard values; maybe move these to schema?
  QC_HARD_MIN_COUNTS  = 200
  QC_HARD_MIN_FEATS   = 100
  QC_HARD_MAX_MITO    = 0.5

  # some checks
  config['qc']['exclude_mito']        = _safe_boolean(config['qc']['exclude_mito'])

  # make sure they're consistent
  config['qc']['qc_hard_min_counts']  = min(QC_HARD_MIN_COUNTS, config['qc']['qc_min_counts'])
  config['qc']['qc_hard_min_feats']   = min(QC_HARD_MIN_FEATS, config['qc']['qc_min_feats'])
  config['qc']['qc_hard_max_mito']    = max(QC_HARD_MAX_MITO, config['qc']['qc_max_mito'])

  return config


# get parameters for hvgs
def check_hvg_parameters(config, METADATA_F, AF_GTF_DT_F):
  # define dummy group names for all
  if config['hvg']['hvg_method'] == 'all':
    config['hvg']['hvg_group_names'] = ['all_samples']
  # if groups, check that the values are ok
  elif config['hvg']['hvg_method'] == 'groups':
    # check that value of metadata_split_var matches a column in sample metadata
    hvg_split_var = config['hvg']['hvg_metadata_split_var']
    meta          = pd.read_csv(METADATA_F)
    if not hvg_split_var in meta.columns():
      raise KeyError(f"{hvg_split_var} is not a column in the sample metadata file.")
    
    # check number of unique group values
    uniq_groups = meta[ hvg_split_var ].unique().tolist()
    if len(uniq_groups) == meta.shape[0]:
      raise ValueError(f"Number of unique values in '{hvg_split_var}' is the same as the number of samples.")

    # store nice names
    config['hvg']['hvg_group_names'] = [n.replace(" ", "_") for n in uniq_groups]

  # get number of gene chunks if method is 'groups' or 'all'
  if config['hvg']['hvg_method'] in ['groups', 'all']:
    # get total number of genes
    gtf_df      = pd.read_csv(AF_GTF_DT_F,  sep = '\t')
    num_genes   = gtf_df.shape[0]

    # chunk them up and name them
    num_chunks  = (num_genes + config['hvg']['hvg_chunk_size'] - 1) // config['hvg']['hvg_chunk_size']
    chunk_names = [f"chunk_{i+1}" for i in range(num_chunks)]
    
    # add to config
    config['hvg']['hvg_num_chunks'] = num_chunks
    config['hvg']['hvg_chunk_names'] = chunk_names

  return config


# check parameters for integration
def check_integration_parameters(config):
  # nothing to do here at the moment; leaving in case it's useful later

  return config


# get parameters for marker genes
def check_marker_genes_parameters(config, PROJ_DIR, SCPROCESS_DATA_DIR):
  # set some more default values
  config['marker_genes']['mkr_gsea_dir'] = os.path.join(SCPROCESS_DATA_DIR, 'gmt_pathways')

  # get custom marker files
  custom_mkr_names, custom_mkr_paths = _get_custom_marker_genes_specs(config, PROJ_DIR, SCPROCESS_DATA_DIR)
  config['marker_genes']['custom_mkr_names'] = custom_mkr_names
  config['marker_genes']['custom_mkr_paths'] = custom_mkr_paths

  return config


# get specified custom marker genes
def _get_custom_marker_genes_specs(config, PROJ_DIR, SCPROCESS_DATA_DIR):
  # set defaults
  custom_mkr_names = ""
  custom_mkr_paths = ""

  # populate with custom sets
  if 'custom_sets' in config["marker_genes"]:
    custom_sets = config["marker_genes"]["custom_sets"]
    mkr_names = []
    mkr_paths = []
    for i, gene_set in enumerate(custom_sets):
      # get name and file
      name      = gene_set["name"]
      file_path = gene_set.get("file", os.path.join(SCPROCESS_DATA_DIR, 'marker_genes', f"{name}.csv"))

      # check whether it exists
      if not os.path.isabs(file_path):
        file_path = os.path.join(PROJ_DIR, file_path)
      if not os.path.isfile(file_path):
        raise FileNotFoundError(f"File not found for marker set '{name}'")
      if not file_path.endswith(".csv"):
        raise ValueError(f"File for custom marker set '{name}' is not a csv file")

      # check csv file contents
      mkrs_df   = pd.read_csv(file_path)
      req_col   = "label"
      opt_cols  = ["symbol", "ensembl_id"]
      if not req_col in mkrs_df.columns:
        raise KeyError(f"File '{file_path}' is missing the mandatory column 'label'.")
      if not any(col in mkrs_df.columns for col in opt_cols):
        raise KeyError(f"File '{file_path}' must contain at least one of 'symbol' or 'ensembl_id' column.")

      # Store validated values
      mkr_names.append(name)
      mkr_paths.append(file_path)

    custom_mkr_names = ",".join(mkr_names)
    custom_mkr_paths = ",".join(mkr_paths)
  
  return custom_mkr_names, custom_mkr_paths


# get parameters for multiplexing
def get_multiplexing_parameters(config, PROJ_DIR, sample_metadata):
  #set defaults
  DEMUX_TYPE     = "none"
  HTO_FASTQ_DIR  = None
  FEATURE_REF    = None
  DEMUX_F        = ""
  BATCH_VAR      = "sample_id"
  POOL_IDS       = ""
  SAMPLE_IDS     = ""
  EXC_POOLS      = ""
  SAMPLE_MAPPING = None
  
  if ('multiplexing' in config) and (config['multiplexing'] is not None):
    POOL_IDS = list(set(sample_metadata["pool_id"].tolist()))
    SAMPLE_IDS = list(sample_metadata["sample_id"].tolist())
    DEMUX_TYPE = config['multiplexing']['demux_type']
    SAMPLE_MAPPING = sample_metadata.groupby("pool_id")["sample_id"].apply(list).to_dict()

    if DEMUX_TYPE == 'hto':
      # check feature ref specified and valid
      assert 'feature_ref' in config['multiplexing'], \
       "feature_ref not specified in the config file"
      FEATURE_REF     = config['multiplexing']['feature_ref']
      if not os.path.isabs(FEATURE_REF):
        FEATURE_REF = os.path.join(PROJ_DIR, FEATURE_REF)
      assert os.path.isfile(FEATURE_REF), \
        f"feature reference {FEATURE_REF} is not a valid file"
      
      # check for columns in feature ref
      feature_ref = pd.read_csv(FEATURE_REF)
      assert all(col in feature_ref.columns for col in ["hto_id", "sequence"])

      # check that all hto_ids in feature ref file are unique
      assert all(feature_ref["hto_id"].value_counts() == 1), \
       "hto_id values in feature reference file not unique"
      hto_ids = feature_ref["hto_id"].tolist()

      # check if all hto_id values in sample metadata match the ones in feature reference
      assert all(hto in hto_ids for hto in list(set(sample_metadata["hto_id"]))), \
        "One or more hto_id values in sample metadata don't match hto_id values in the feature reference file"

      # get fastqs
      HTO_FASTQ_DIR = config['multiplexing']['fastq_dir']
      if not os.path.isabs(HTO_FASTQ_DIR):
        HTO_FASTQ_DIR = os.path.join(PROJ_DIR, HTO_FASTQ_DIR)
        
    else:
      assert 'demux_output' in config['multiplexing'], \
       "demux_output not specified in the config file"
      DEMUX_F = config['multiplexing']['demux_output']
      # check if valid
      if not os.path.isabs(DEMUX_F):
       DEMUX_F = os.path.join(PROJ_DIR, DEMUX_F)
      assert os.path.isfile(DEMUX_F), \
       f"file {DEMUX_F} doesn't exist"

      # check if looks ok 
      demux_dt = pd.read_csv(DEMUX_F)
      for col in ["pool_id", "sample_id", "cell_id"]:
        assert col in demux_dt.columns, \
        f"{col} not present in demux_output"

      # check if samples in metadata and demux_dt match
      assert {x for x in set(demux_dt['sample_id']) if pd.notna(x)} == set(sample_metadata['sample_id']), \
        "Unique values for sample_id don't match in demux_output and sample_metadata"
      
      assert set(demux_dt['pool_id']) == set(sample_metadata['pool_id']), \
        "Unique values for pool_id don't match in demux_output and sample_metadata"

    if 'batch_var' in config['multiplexing']:
      BATCH_VAR = config['multiplexing']['batch_var']
      assert BATCH_VAR in ["sample_id", "pool_id"], \
       "Invalid value of batch_var parameter. Must be either sample_id or pool_id"
  
    if ('exclude' in config) and (config['exclude'] is not None):
      if ('pool_id' in config['exclude']) and (config['exclude']['pool_id'] is not None):
        EXC_POOLS = config['exclude']["pool_id"]
        EXC_SAMPLES = [] # exclude all samples in excluded pools
        for p in EXC_POOLS:
          if p not in POOL_IDS:
            warnings.warn(f"sample {p} specified in exclude_pools but not in sample_metadata file", UserWarning)
          EXC_SAMPLES.extend(SAMPLE_MAPPING[p])

        pools_to_keep  = set(POOL_IDS) - set(EXC_POOLS)
        smpls_to_keep  = set(SAMPLE_IDS) - set(EXC_SAMPLES)

        POOL_IDS = [p for p in POOL_IDS if p in pools_to_keep]
        SAMPLE_IDS = [s for s in SAMPLE_IDS if s in smpls_to_keep]
      
  return DEMUX_TYPE, HTO_FASTQ_DIR, FEATURE_REF, DEMUX_F, BATCH_VAR, EXC_POOLS, POOL_IDS, SAMPLE_IDS, SAMPLE_MAPPING


# get list of samples
def get_project_parameters(config, scprocess_data_dir):
  ## what is specified in config directory?
  proj_dict     = config['project']
  PROJ_DIR      = proj_dict["proj_dir"]
  SHORT_TAG     = proj_dict["short_tag"]
  FULL_TAG      = proj_dict["full_tag"]
  YOUR_NAME     = proj_dict["your_name"]
  AFFILIATION   = proj_dict["affiliation"]
  DATE_STAMP    = proj_dict["date_stamp"]
  SPECIES       = proj_dict["species"]

  # check fastq vs arvados
  has_fastq     = "fastq_dir" in proj_dict
  has_arv_uuid  = "arv_uuid" in proj_dict
  if has_fastq + has_arv_uuid != 1:
    KeyError('config file must contain exactly one of "fastq_dir" and "arv_uuid"')

  # define these
  if has_fastq and not has_arv_uuid:
    FASTQ_DIR     = proj_dict["fastq_dir"]
    ARV_UUID      = None
  if not has_fastq and has_arv_uuid:
    FASTQ_DIR     = None
    ARV_UUID      = proj_dict["arv_uuid"]

  # check if selected species is valid
  index_params_f  = os.path.join(scprocess_data_dir, 'index_parameters.csv')

  # from index_parameters.csv get valid values for species
  index_params      = pd.read_csv(index_params_f)
  valid_species     = index_params['genome_name'].tolist()
  valid_species_str = ', '.join(valid_species)

  assert SPECIES in valid_species, f"species {SPECIES} not defined. Valid values are {valid_species_str}"

  # check whether date is given as datetime object
  if isinstance(DATE_STAMP, datetime.date):
    # if so, convert to string
    DATE_STAMP  = DATE_STAMP.strftime("%Y-%m-%d")

  date_regex    = re.compile("^20[0-9]{2}-[0-9]{2}-[0-9]{2}$")
  assert date_regex.match(DATE_STAMP), f"{DATE_STAMP} does not match date format YYYY-MM-DD"

  # get fastqs
  if FASTQ_DIR is not None:
    if not os.path.isabs(FASTQ_DIR):
      FASTQ_DIR = os.path.join(PROJ_DIR, FASTQ_DIR)
  
  # get samples
  METADATA_F    = proj_dict["sample_metadata"]
  if not os.path.isabs(METADATA_F):
    METADATA_F = os.path.join(PROJ_DIR, METADATA_F)

  samples_df    = pd.read_csv( METADATA_F )
  SAMPLES       = samples_df["sample_id"].dropna().tolist()

  # get multiplexing params
  DEMUX_TYPE, HTO_FASTQ_DIR, FEATURE_REF, DEMUX_F, BATCH_VAR, EXC_POOLS, POOL_IDS, \
    POOL_SAMPLES, SAMPLE_MAPPING = get_multiplexing_parameters(proj_dict, PROJ_DIR, samples_df)

  # define sample variable for alevin, ambient, doublets and sample mapping dictionary
  if DEMUX_TYPE != "none":
    SAMPLE_VAR  = "pool_id"
    # exclude samples in excluded pools
    SAMPLES     = POOL_SAMPLES
  else:
    SAMPLE_VAR  = "sample_id"
  
  # remove some samples
  EXC_SAMPLES   = None
  if ('exclude' in proj_dict) and (proj_dict['exclude'] is not None):
    if ('sample_id' in proj_dict['exclude']) and (proj_dict['exclude']['sample_id'] is not None):
      EXC_SAMPLES = proj_dict['exclude']["sample_id"]
      for s in EXC_SAMPLES:
        if s not in samples_df["sample_id"].dropna().tolist():
          warnings.warn(f"sample {s} specified in exclude_samples but not in metadata file", UserWarning)
      to_keep       = set(SAMPLES) - set(EXC_SAMPLES)
      SAMPLES       = [s for s in SAMPLES if s in to_keep]

  # check custom sample proj_dict file if exists
  CUSTOM_SAMPLE_PARAMS_F = None
  if ('custom_sample_params' in proj_dict) and (proj_dict['custom_sample_params'] is not None):
    CUSTOM_SAMPLE_PARAMS_F = proj_dict["custom_sample_params"]
    # check if exists
    if not os.path.isabs(CUSTOM_SAMPLE_PARAMS_F):
      CUSTOM_SAMPLE_PARAMS_F = os.path.join(PROJ_DIR, CUSTOM_SAMPLE_PARAMS_F)
    assert os.path.isfile(CUSTOM_SAMPLE_PARAMS_F)
    # open the file and check if all samples can be found in SAMPLES or POOL_IDS
    with open(CUSTOM_SAMPLE_PARAMS_F) as f:
      custom_smpl_params = yaml.load(f, Loader=yaml.FullLoader)
      custom_smpls = list(custom_smpl_params.keys())
      if DEMUX_TYPE != "none":
        for s in custom_smpls:
          assert s in POOL_IDS, f"{s} in custom_sample_params file doesn't match any pool_id values in sample_metadata"
      else:
        for s in custom_smpls:
          assert s in SAMPLES, f"{s} in custom_sample_params file doesn't match any sample_id values in sample_metadata"

  # sort out metadata variables
  METADATA_VARS = []
  if ('metadata_vars' in proj_dict) and (proj_dict['metadata_vars'] is not None):
    METADATA_VARS  = proj_dict["metadata_vars"]
    for var in METADATA_VARS:
      assert var in samples_df.columns, f"{var} not in sample metadata"
      # check that there are less than 10 unique values (otherwise probably not a categorical variable)
      var_vals  = samples_df[var].tolist()
      if all(isinstance(x, float) and math.isnan(x) for x in var_vals):
        continue
      uniq_vals = len(set(var_vals))
      assert uniq_vals <= 10, \
        f"{var} variable has more than 10 unique values"

  return PROJ_DIR, FASTQ_DIR, SHORT_TAG, FULL_TAG, YOUR_NAME, AFFILIATION, METADATA_F, \
    METADATA_VARS, EXC_SAMPLES, SAMPLES, DATE_STAMP, CUSTOM_SAMPLE_PARAMS_F, SPECIES, \
    DEMUX_TYPE, HTO_FASTQ_DIR, FEATURE_REF, DEMUX_F, BATCH_VAR, EXC_POOLS, POOL_IDS, \
    SAMPLE_VAR, SAMPLE_MAPPING


# get parameters for mapping
def get_mapping_parameters(config, scprocess_data_dir, SPECIES):
  # get chemisty
  CHEMISTRY     = config['project']['tenx_chemistry']
  
  # from index_parameters.csv get valid values for species
  idx_params_f  = os.path.join(scprocess_data_dir, 'index_parameters.csv')
  index_params  = pd.read_csv(idx_params_f)
     
  # get mito strings from setup params
  AF_MITO_STR   = index_params.loc[index_params['genome_name'] == SPECIES, 'mito_str'].values[0]

  # get af index directory and check if exists
  AF_HOME_DIR   = os.path.join(scprocess_data_dir, 'alevin_fry_home') # check if this exists in scprocess script
  AF_INDEX_DIR  = os.path.join(AF_HOME_DIR, SPECIES)
  if not os.path.isdir(AF_INDEX_DIR):
    raise FileNotFoundError(f"alevin index for {SPECIES} doesn't exist")
  
  # get gtf txt file, check that exists
  AF_GTF_DT_F = index_params.loc[index_params['genome_name'] == SPECIES, 'gtf_txt_f'].values[0]

  return AF_MITO_STR, AF_HOME_DIR, AF_INDEX_DIR, AF_GTF_DT_F, CHEMISTRY


# get parameters for ambient
def get_ambient_parameters(config):
  # set default values
  AMBIENT_METHOD                  = 'decontx'
  CELLBENDER_VERSION              = 'v0.3.2'
  CELLBENDER_PROP_MAX_KEPT        = 0.9
  FORCE_EXPECTED_CELLS            = None
  FORCE_TOTAL_DROPLETS_INCLUDED   = None
  FORCE_LOW_COUNT_THRESHOLD       = None
  CELLBENDER_LEARNING_RATE        = 1e-4
  CELLBENDER_POSTERIOR_BATCH_SIZE = 128 # parameter only available in v0.3.2. Smaller values for lower GPU memory
  CELL_CALLS_METHOD               = 'barcodeRanks'

  # change defaults if specified
  if ('ambient' in config) and (config['ambient'] is not None):
    if 'ambient_method' in config['ambient']:
      AMBIENT_METHOD                 = config['ambient']['ambient_method']
    if 'cell_calling' in config['ambient']:
      CELL_CALLS_METHOD             = config['ambient']['cell_calling']
    if 'cellbender_version' in config['ambient']:
      CELLBENDER_VERSION            = config['ambient']['cellbender_version']
    if 'cb_max_prop_kept' in config['ambient']:
      CELLBENDER_PROP_MAX_KEPT      = config['ambient']['cb_max_prop_kept']
    if 'cb_force_expected_cells' in config['ambient']:
      FORCE_EXPECTED_CELLS          = config['ambient']['cb_force_expected_cells']
    if 'cb_force_total_droplets_included' in config['ambient']:
      FORCE_TOTAL_DROPLETS_INCLUDED = config['ambient']['cb_force_total_droplets_included']
    if 'cb_force_low_count_threshold' in config['ambient']:
      FORCE_LOW_COUNT_THRESHOLD     = config['ambient']['cb_force_low_count_threshold']
    if 'cb_force_learning_rate' in config['ambient']:
      CELLBENDER_LEARNING_RATE      = config['ambient']['cb_force_learning_rate']
    if 'cb_posterior_batch_size' in config['ambient']:
      CELLBENDER_POSTERIOR_BATCH_SIZE =config['ambient']['cb_posterior_batch_size']
      if CELLBENDER_VERSION != 'v0.3.2':
        warnings.warn(f"'cb_posterior_batch_size' is only supported in CellBender v0.3.2. Ignoring for CellBender {CELLBENDER_VERSION}.")

  # get cellbender image (maybe skip this if cellbender is not selected?)
  if CELLBENDER_VERSION   == 'v0.3.2':
    CELLBENDER_IMAGE  = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.2'
  elif CELLBENDER_VERSION == 'v0.3.0':
    CELLBENDER_IMAGE  = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.0'
  elif CELLBENDER_VERSION == 'v0.2.0':
    CELLBENDER_IMAGE  = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.2.0'
  else:
    raise ValueError(f"selected cellbender version {CELLBENDER_VERSION} not supported")

  # some checks on custom parameters for cellbender
      
  return CELLBENDER_IMAGE, CELLBENDER_VERSION, CELLBENDER_PROP_MAX_KEPT, AMBIENT_METHOD, CELL_CALLS_METHOD, \
    FORCE_EXPECTED_CELLS, FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE, CELLBENDER_POSTERIOR_BATCH_SIZE 


# get parameters for labelling celltypes
def get_label_celltypes_parameters(config, SPECIES, SCPROCESS_DATA_DIR): 
  # set some more default values
  LBL_XGB_F       = None
  LBL_XGB_CLS_F   = None
  LBL_TISSUE      = ""
  LBL_GENE_VAR    = "gene_id"
  LBL_SEL_RES_CL  = "RNA_snn_res.2"
  LBL_MIN_PRED    = 0.8
  LBL_MIN_CL_PROP = 0.5
  LBL_MIN_CL_SIZE = 100

  # change defaults if specified
  if ('label_celltypes' in config) and (config['label_celltypes'] is not None):

    assert 'lbl_tissue' in config['label_celltypes'], \
      "lbl_tissue parameter missing in config file"
    LBL_TISSUE      = config['label_celltypes']['lbl_tissue']
    
    if 'lbl_sel_res_cl' in config['label_celltypes']:
      LBL_SEL_RES_CL  = config['label_celltypes']['lbl_sel_res_cl']
    if 'lbl_min_pred' in config['label_celltypes']:
      LBL_MIN_PRED    = config['label_celltypes']['lbl_min_pred']
    if 'lbl_min_cl_prop' in config['label_celltypes']:
      LBL_MIN_CL_PROP = config['label_celltypes']['lbl_min_cl_prop']
    if 'lbl_min_cl_size' in config['label_celltypes']:
      LBL_MIN_CL_SIZE = config['label_celltypes']['lbl_min_cl_size']

    # check that classifier name is valid
    valid_boosts = ['', 'human_cns', 'mouse_cns']
    assert LBL_TISSUE in valid_boosts, \
      f"value {LBL_TISSUE} for 'lbl_tissue' parameter is not valid"
 
    # pick labeller
    xgb_dir  = os.path.join(SCPROCESS_DATA_DIR, 'xgboost')
    assert os.path.isdir(xgb_dir)
    
    if LBL_TISSUE == '':
      LBL_XGB_F       = ""
      LBL_XGB_CLS_F   = ""
    elif LBL_TISSUE == 'human_cns':
      LBL_XGB_F       = os.path.join(xgb_dir, "Siletti_Macnair-2025-07-23/xgboost_obj_hvgs_Siletti_Macnair_2025-07-23.rds")
      LBL_XGB_CLS_F   = os.path.join(xgb_dir, "Siletti_Macnair-2025-07-23/allowed_cls_Siletti_Macnair_2025-07-23.csv")
    else: 
      raise ValueError(f"{LBL_TISSUE} classifier is unfortunately not available yet")

    if LBL_TISSUE != '':
      # check these are ok
      assert os.path.isfile(LBL_XGB_F), \
        f"file {LBL_XGB_F} doesn't exist; consider (re)runnning scprocess setup"
      assert os.path.isfile(LBL_XGB_CLS_F), \
        f"file {LBL_XGB_CLS_F} doesn't exist; consider (re)runnning scprocess setup"
 
  return LBL_XGB_F, LBL_XGB_CLS_F, LBL_GENE_VAR, LBL_SEL_RES_CL, LBL_MIN_PRED, LBL_MIN_CL_PROP, LBL_MIN_CL_SIZE, LBL_TISSUE


# get parameters for pb and empties
def get_pb_empties_parameters(config):
  # get parameters for filtering edger results
  AMBIENT_GENES_LOGFC_THR = 0
  AMBIENT_GENES_FDR_THR   = 0.01
  
  if ('pb_empties' in config) and (config['pb_empties'] is not None):
    if 'ambient_genes_logfc_thr' in config['pb_empties']:
      AMBIENT_GENES_LOGFC_THR = config['pb_empties']['ambient_genes_logfc_thr']
    if 'ambient_genes_fdr_thr'   in config['pb_empties']:
      AMBIENT_GENES_FDR_THR   = config['pb_empties']['ambient_genes_fdr_thr']

  return AMBIENT_GENES_LOGFC_THR, AMBIENT_GENES_FDR_THR


# get parameters for zoom
def get_zoom_parameters(config, LBL_TISSUE, LBL_XGB_CLS_F, METADATA_F, 
  AF_GTF_DT_F, PROJ_DIR, SHORT_TAG, FULL_TAG, DATE_STAMP, SCPROCESS_DATA_DIR):
  # if (rule_name != 'zoom') or ('zoom' not in config) or (config['zoom'] is None):
  if ('zoom' not in config) or (config['zoom'] is None):
    ZOOM_NAMES        = []
    ZOOM_PARAMS_DICT  = []
    ZOOM_NAMES_SUBSET = []
  else:
    ZOOM_NAMES    = list(config['zoom'].keys())
    assert len(ZOOM_NAMES) == len(set(ZOOM_NAMES)), \
     "all subset labels for zoom must be unique"
    
    ZOOM_YAMLS = list(config['zoom'].values())
    for i, zoom_f in enumerate(ZOOM_YAMLS):
      if not os.path.isabs(zoom_f):
        ZOOM_YAMLS[i] = os.path.join(PROJ_DIR, zoom_f)
      else:
        ZOOM_YAMLS[i] = zoom_f
      assert os.path.isfile(ZOOM_YAMLS[i]), \
        f"file {ZOOM_YAMLS[i]} doesn't exist"
    
    # make dictionary of zoom params
    ZOOM_PARAMS_DICT = dict(zip(
      ZOOM_NAMES,
      [ _get_one_zoom_parameters(zoom_f, LBL_TISSUE, LBL_XGB_CLS_F, METADATA_F, AF_GTF_DT_F,
        PROJ_DIR, SHORT_TAG, FULL_TAG, DATE_STAMP, SCPROCESS_DATA_DIR
      ) for zoom_f in ZOOM_YAMLS ]
      ))

    # get all zoom names for which subset sces should be created
    ZOOM_NAMES_SUBSET = [zoom_name for zoom_name in ZOOM_NAMES if ZOOM_PARAMS_DICT[zoom_name]["MAKE_SUBSET_SCES"]]

  return ZOOM_NAMES, ZOOM_PARAMS_DICT, ZOOM_NAMES_SUBSET


### other functions

# get list of clusters to check zoom specification
def _get_cl_ls(PROJ_DIR, SHORT_TAG, FULL_TAG, DATE_STAMP, SEL_RES):
  # specify harmony outputs
  int_dir     = f"{PROJ_DIR}/output/{SHORT_TAG}_integration"
  int_f       = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'

  # get list of clusters
  int_dt      = pd.read_csv(int_f)
  cl_col      = f"RNA_snn_res.{SEL_RES}"
  cl_ls       = list(int_dt[cl_col].unique())
  cl_ls       = [cl for cl in cl_ls if str(cl) != "nan"]
  cl_ls       = sorted(cl_ls)

  return int_f, cl_ls


# get parameters for one zoom specification
def _get_one_zoom_parameters(zoom_yaml_f, LBL_TISSUE, LBL_XGB_CLS_F, METADATA_F, 
  AF_GTF_DT_F, PROJ_DIR, SHORT_TAG, FULL_TAG, DATE_STAMP, SCPROCESS_DATA_DIR):
  # set defaults
  LABELS           = ""
  LABELS_F         = ""
  LABELS_SOURCE    = ""
  LBL_SEL_RES_CL   = "RNA_snn_res.2"
  CLUSTER_RES      = None
  CUSTOM_LABELS_F  = ""
  MIN_N_SAMPLE     = 10
  MAKE_SUBSET_SCES = True

  # unpack
  with open(zoom_yaml_f, "r") as stream:
    this_zoom = yaml.safe_load(stream)

  # check required params
  for key in ['labels', 'labels_source']:
    assert key in this_zoom and this_zoom[key] is not None, \
      f"{key} parameter missing from file {zoom_yaml_f}"

  LABELS        = this_zoom['labels']
  LABELS_SOURCE = this_zoom['labels_source']

  valid_sources = ['xgboost', 'clusters', 'custom']
  assert LABELS_SOURCE in valid_sources, \
    f'labels_source must be one of {valid_sources}'

  if LABELS_SOURCE == 'custom':
    assert 'custom_labels_f' in this_zoom and this_zoom['custom_labels_f'] is not None, \
      f"custom_labels_f parameter missing from file {zoom_yaml_f}"
    CUSTOM_LABELS_F = this_zoom['custom_labels_f']

    # check that exists
    if not os.path.isabs(CUSTOM_LABELS_F):
      CUSTOM_LABELS_F = os.path.join(PROJ_DIR, CUSTOM_LABELS_F)
    assert os.path.isfile(CUSTOM_LABELS_F), \
      f"file {CUSTOM_LABELS_F} doesn't exist"

    # check that columns are ok
    custom_lbls_dt = pd.read_csv(CUSTOM_LABELS_F)
    for col in ['sample_id', 'cell_id', 'label']:
      assert col in custom_lbls_dt.columns, \
        f"column {col} not present in {CUSTOM_LABELS_F}"
        
    # check that all labels are in the labels column of the custom file
    assert all([lbl in set(custom_lbls_dt['label'].tolist()) for lbl in LABELS])
  
    LABELS_F   = CUSTOM_LABELS_F
    LABELS_VAR = 'label'

  if LABELS_SOURCE == 'xgboost':
    assert LBL_TISSUE != "", \
      "lbl_tissue parameter is not defined"
    # check that labels are ok
    xgb_allow_lbls = pd.read_csv(LBL_XGB_CLS_F)['cluster'].tolist()
    for lbl in LABELS:
      assert lbl in xgb_allow_lbls, \
        f"{lbl} is not a valid label name for the {LBL_TISSUE} classifier"
    
    lbl_dir     = f"{PROJ_DIR}/output/{SHORT_TAG}_label_celltypes"
    LABELS_F    = lbl_dir + '/cell_annotations_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
    assert os.path.exists(LABELS_F), \
      f"{LABELS_F} doesn't exist; consider (re)running rule label_celltypes"
    
    if 'lbl_sel_res_cl' in this_zoom:
      LBL_SEL_RES_CL  = this_zoom['lbl_sel_res_cl']
    
    LABELS_VAR = "cl_pred_" + LBL_SEL_RES_CL

  if LABELS_SOURCE == 'clusters':
    assert 'cluster_res' in this_zoom and this_zoom['cluster_res'] is not None, \
      f"cluster_res parameter missing from file {zoom_yaml_f}"
    CLUSTER_RES = this_zoom['cluster_res']

    # get list of all clusters to check if cluster names are valid
    LABELS_F, cl_ls = _get_cl_ls(PROJ_DIR, SHORT_TAG, FULL_TAG, DATE_STAMP, CLUSTER_RES)
    assert set(LABELS).issubset(cl_ls)
    LABELS_VAR = f"RNA_snn_res.{CLUSTER_RES}"

  # check optional parameters
  if 'min_n_sample' in this_zoom:
    MIN_N_SAMPLE = this_zoom['min_n_sample']
  if 'make_subset_sces' in this_zoom:
    MAKE_SUBSET_SCES = this_zoom['make_subset_sces']
    MAKE_SUBSET_SCES = int(_safe_boolean(MAKE_SUBSET_SCES))

  # hvg params
  HVG_PARAMS  = get_hvg_parameters(this_zoom, METADATA_F, AF_GTF_DT_F)
  HVG_KEYS    = ["HVG_METHOD", "HVG_SPLIT_VAR", "HVG_CHUNK_SIZE", "HVG_NUM_CHUNKS", 
    "HVG_GROUP_NAMES", "HVG_CHUNK_NAMES", "N_HVGS", "EXCLUDE_AMBIENT_GENES"]
  HVG_DICT    = dict(zip(HVG_KEYS, HVG_PARAMS))

  # pb_empties params
  PB_EMPTIES_PARAMS = get_pb_empties_parameters(this_zoom) 
  PB_EMPTIES_KEYS   = ["AMBIENT_GENES_LOGFC_THR", "AMBIENT_GENES_FDR_THR"]
  PB_EMPTIES_DICT   = dict(zip(PB_EMPTIES_KEYS, PB_EMPTIES_PARAMS))

  # integration params
  INT_PARAMS  = get_integration_parameters(this_zoom)
  INT_KEYS    = ["INT_CL_METHOD", "INT_REDUCTION", "INT_N_DIMS", "INT_THETA", "INT_RES_LS"]
  INT_DICT    = dict(zip(INT_KEYS, INT_PARAMS[:5]))   
  
  # marker gene params
  MKR_PARAMS  = get_marker_genes_parameters(this_zoom, PROJ_DIR, SCPROCESS_DATA_DIR)
  MKR_KEYS    = ["MKR_SEL_RES", "MKR_GSEA_DIR", "MKR_MIN_CL_SIZE", "MKR_MIN_CELLS", 
    "MKR_NOT_OK_RE", "MKR_MIN_CPM_MKR", "MKR_MIN_CPM_GO", "MKR_MAX_ZERO_P", "MKR_GSEA_CUT", 
    "CUSTOM_MKR_NAMES", "CUSTOM_MKR_PATHS"]
  MKR_DICT    = dict(zip(MKR_KEYS, MKR_PARAMS))

  # combine all parameters into a single dictionary
  params = {
    "LABELS": LABELS,
    "LABELS_F": LABELS_F, 
    "LABELS_VAR": LABELS_VAR, 
    "LABELS_SOURCE": LABELS_SOURCE,
    "CLUSTER_RES": CLUSTER_RES,
    "CUSTOM_LABELS_F": CUSTOM_LABELS_F,
    "LBL_SEL_RES_CL": LBL_SEL_RES_CL, 
    "MIN_N_SAMPLE": MIN_N_SAMPLE,
    "MAKE_SUBSET_SCES": MAKE_SUBSET_SCES,
    **HVG_DICT,
    **PB_EMPTIES_DICT, 
    **INT_DICT,
    **MKR_DICT
  }

  return params


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
def exclude_samples_without_fastq_files(FASTQS_DIR, SAMPLES, HTO = False):
  # get parameters

  # get fastq files for each sample
  chk_samples   = []
  for sample in SAMPLES:
    R1_fs = find_fastq_files(FASTQS_DIR, sample, "R1")
    if len(R1_fs) > 0:
      chk_samples.append(sample)
    else:
      if HTO:
        print(f"WARNING: no hto fastq files found for sample {sample}; excluded")
      else:
        print(f"WARNING: no fastq files found for sample {sample}; excluded")

  return chk_samples


# remove samples and pools from sample mapping
def update_sample_mapping_after_exclusions(SAMPLE_MAPPING, POOL_IDS, SAMPLES):

  if SAMPLE_MAPPING is None:
    return None
  # filter to keep specific pool ids
  SAMPLE_MAPPING = {pool_id: sample_ids for pool_id, sample_ids in SAMPLE_MAPPING.items() if pool_id in POOL_IDS}
  # filter to keep specific sample_ids
  for pool_id in SAMPLE_MAPPING:
    SAMPLE_MAPPING[pool_id] = [sample_id for sample_id in SAMPLE_MAPPING[pool_id] if sample_id in SAMPLES]

  return SAMPLE_MAPPING


# nice boolean values from yaml inputs
def _safe_boolean(val):
  if type(val) is bool:
    res = val
  elif val in ["True", "true"]:
    res = True
  elif val in ["False", "false"]:
    res = False
  else:
    raise ValueError('{val} is not a boolean')

  return res


# helper function to merge multiple zipped csv/txt files
def merge_tmp_files(in_files, out_file):
  df_ls     = [pd.read_csv(f, compression='gzip', sep='\t') for f in in_files if gzip.open(f, 'rb').read(1)]
  df_merged = pd.concat(df_ls, ignore_index=True)
  df_merged.to_csv(out_file, sep='\t', index=False, compression='gzip', quoting=csv.QUOTE_NONE)


# HVGs function: make df with list of chunked counts files
def make_hvgs_input_df(DEMUX_TYPE, SAMPLE_VAR, runs, ambient_outs_yamls, SAMPLE_MAPPING, FULL_TAG, DATE_STAMP, hvg_dir):
  # define output variable
  df_list = []

  # populate list
  for r, yaml_file in zip(runs, ambient_outs_yamls):
    # get filtered ambient outputs
    with open(yaml_file) as f:
      amb_outs = yaml.load(f, Loader=yaml.FullLoader)

    amb_filt_f = amb_outs['filt_counts_f']

    if DEMUX_TYPE != "none":
      # get sample ids for pool
      sample_ids = SAMPLE_MAPPING.get(r, [])

      for sample_id in sample_ids:
        hvg_df = pd.DataFrame({
          SAMPLE_VAR: [r],
          'amb_filt_f': [amb_filt_f],
          'sample_id': [sample_id]
        })

        df_list.append(hvg_df)
    else:
      hvg_df = pd.DataFrame({
        SAMPLE_VAR: [r],
        'amb_filt_f': [amb_filt_f]
      })
      df_list.append(hvg_df)

  # merge dfs for all runs
  hvg_df_full = pd.concat(df_list, ignore_index=True)

  # add path to chunked file
  hvg_df_full['chunked_f'] = hvg_df_full['sample_id'].apply(lambda s: f"{hvg_dir}/chunked_counts_{s}_{FULL_TAG}_{DATE_STAMP}.h5")

  return hvg_df_full


# zoom function: get list of all mean var files for zooms
def get_zoom_std_var_stats_files(zoom_name, zoom_dir, ZOOM_PARAMS_DICT, FULL_TAG, DATE_STAMP, SAMPLES):
  hvg_method = ZOOM_PARAMS_DICT[zoom_name]['HVG_METHOD']

  if hvg_method == "sample":
    return [
      zoom_dir + f'/{zoom_name}/tmp_std_var_stats_{sample}_sample_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
      for sample in SAMPLES
    ]
  else:
    group_names = ZOOM_PARAMS_DICT[zoom_name]['HVG_GROUP_NAMES']
    num_chunks = ZOOM_PARAMS_DICT[zoom_name]['HVG_NUM_CHUNKS']

    return [
      zoom_dir + f'/{zoom_name}/tmp_std_var_stats_{group}_group_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
      for group in group_names
      for chunk in range(num_chunks)
    ]


# zoom function: make df with good / bad sample labels for a specific zoom
def extract_zoom_sample_statistics(qc_stats_f, SAMPLES, LABELS_F, LABELS_VAR, LABELS, MIN_N_SAMPLE, AMBIENT_METHOD):
  # load inputs
  qc_df     = pd.read_csv(qc_stats_f)
  qc_df     = qc_df.drop('n_cells', axis=1)
  lbls_dt   = pd.read_csv(LABELS_F, compression='gzip')

  # keep selected labels
  lbls_dt   = lbls_dt[ lbls_dt[LABELS_VAR].isin(LABELS) ]
  
  # count the number of cells per sample
  zoom_sample_stats = (
    lbls_dt.groupby('sample_id')
    .size()
    .reset_index(name='n_cells')
  )
  
  # add empty samples
  empty_ss  = list(set(SAMPLES) - set(zoom_sample_stats["sample_id"].tolist()))
  empty_df  = pd.DataFrame({ "sample_id": empty_ss, "n_cells": 0 })
  zoom_sample_stats = pd.concat([zoom_sample_stats, empty_df])

  # identify samples that do not meet the minimum cell threshold
  zoom_sample_stats['bad_zoom_qc'] = zoom_sample_stats['n_cells'] < MIN_N_SAMPLE
  
  # merge new and existing sample stats
  sample_df = qc_df.merge(zoom_sample_stats, on='sample_id',how='left')
  
  # update 'bad_sample' column
  if AMBIENT_METHOD == 'cellbender':
    sample_df['bad_sample'] = (
      sample_df['bad_bender'] | sample_df['bad_qc'] | sample_df['bad_zoom_qc']
    )
  else:
    sample_df['bad_sample'] = (
      sample_df['bad_qc'] | sample_df['bad_zoom_qc']
    )

  # check that at least 2 good samples remain
  good_smpls_count = (sample_df['bad_sample'] == False).sum()
  assert good_smpls_count >= 2, \
    "Fewer than 2 samples available for this zoom."
  
  return sample_df


# zoom function: get list of all mean var files for zooms
def get_zoom_raw_mean_var_files(zoom_name, ZOOM_PARAMS_DICT, FULL_TAG, DATE_STAMP):
  group_names = ZOOM_PARAMS_DICT[zoom_name]['HVG_GROUP_NAMES']
  num_chunks = ZOOM_PARAMS_DICT[zoom_name]['HVG_NUM_CHUNKS']

  return [
    zoom_dir + f'/{zoom_name}/tmp_mean_var_{group}_group_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
    for group in group_names
    for chunk in range(num_chunks)
  ]


# zoom function: specify some optional outputs for zoom (at the moment only FGSEA outputs)
def get_zoom_conditional_outputs(species, zoom_dir, FULL_TAG, DATE_STAMP):
  if species in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020']:
    return {
      'fgsea_go_bp_f': zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_go_bp_' + DATE_STAMP + '.txt.gz', 
      'fgsea_go_cc_f': zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_go_cc_' + DATE_STAMP + '.txt.gz',
      'fgsea_go_mf_f': zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_go_mf_' + DATE_STAMP + '.txt.gz',
      'fgsea_paths_f': zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_paths_' + DATE_STAMP + '.txt.gz',
      'fgsea_hlmk_f':  zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_hlmk_' + DATE_STAMP + '.txt.gz'
    }
  else:
    return {}


# end
