# load modules
import os
import sys
import re
import copy
import pathlib
import warnings
import yaml
import polars as pl
import csv
import math
import glob
import gzip
import datetime
import subprocess
import shutil
try:
  import snakemake
except ImportError:
  snakemake = None
import json
import jsonschema
from jsonschema.exceptions import best_match

### not much setup

# do some checks of setup before running scprocess
def check_setup_before_running_scprocess(scprocess_dir, extraargs):
  # check that SCPROCESS_DATA_DIR exists
  scdata_dir  = pathlib.Path(os.getenv('SCPROCESS_DATA_DIR'))
  if not scdata_dir:
    raise ValueError('SCPROCESS_DATA_DIR is not defined an environment variable')
  if not scdata_dir.is_dir():
    raise FileNotFoundError("SCPROCESS_DATA_DIR is not a directory")
  
  # check that spcrocess_data_dir has some files
  scsetup_dirs = ['cellranger_ref', 'gmt_pathways', 'marker_genes', 'xgboost', 'alevin_fry_home']
  scsetup_full_dirs = [ scdata_dir / d for d in scsetup_dirs]
  for d in scsetup_full_dirs:
    if not os.path.isdir(d):
      raise FileNotFoundError(f"A directory is missing in {scdata_dir}; consider (re)running setup.\nMissing directory:\n{d}")
  
  # check that setup csv exists
  scsetup_csv = scdata_dir / 'index_parameters.csv'
  if not os.path.isfile(scsetup_csv):
    raise FileNotFoundError(f"{scsetup_csv} is missing; consider (re)running setup.")
  
  # check if cluster profile is defined
  setup_configfile  = scdata_dir / 'scprocess_setup.yaml'
  if not os.path.exists(setup_configfile):
    raise FileNotFoundError(f"scprocess_setup.yaml does not exist in {scdata_dir}")

  # load and validate setup config file
  with open(setup_configfile, "r") as stream:
    setup_cfg     = yaml.safe_load(stream)
  setup_schema_f  = scprocess_dir / "resources/schemas/setup.schema.json"
  setup_cfg       = check_setup_config(setup_cfg, setup_schema_f, scprocess_dir)

  # add profile or local_cores to snakemake call
  if _uses_cluster_profile(setup_cfg):
    profile_dir = get_cluster_profile_dir(scprocess_dir, scdata_dir, setup_cfg)
    extraargs.append('--workflow-profile'),
    extraargs.append(str(profile_dir))

  else:
    extraargs.append('--cores')
    extraargs.append(str(setup_cfg['user']['local_cores']))

  return scdata_dir, extraargs, setup_cfg


def _is_plain_profile_name(profile_name):
  path = pathlib.PurePath(profile_name)
  return profile_name not in ["", ".", ".."] and path.name == profile_name


def _normalise_profile_config(setup_cfg):
  user = setup_cfg.get('user', {})
  if 'profile_template' in user:
    profile_template = user['profile_template']
    profile_name     = user.get('profile_name', profile_template)
  elif 'profile' in user:
    profile_template = user['profile']
    profile_name     = user.get('profile_name', profile_template)
  else:
    return setup_cfg

  for key, value in {
      'profile_template': profile_template,
      'profile_name': profile_name,
  }.items():
    if not _is_plain_profile_name(value):
      raise ValueError(f"user.{key} must be a profile directory name, not a path: {value}")

  user['profile_template'] = profile_template
  user['profile_name']     = profile_name
  return setup_cfg


def _uses_cluster_profile(setup_cfg):
  user = setup_cfg.get('user', {})
  return any(key in user for key in ['profile', 'profile_template'])


def get_cluster_profile_dir(scprocess_dir, scdata_dir, setup_cfg):
  setup_cfg   = _normalise_profile_config(setup_cfg)
  profile     = setup_cfg['user']['profile_name']
  profile_dir = scdata_dir / 'profiles' / profile
  profile_f   = profile_dir / 'config.yaml'
  if not profile_f.is_file():
    raise FileNotFoundError(
      f"cluster configuration file {profile_f} does not exist. "
      "Run `scprocess setup` to copy the profile template into SCPROCESS_DATA_DIR, "
      "then edit the local profile if needed."
    )

  setup_cfg['user']['profile_dir'] = profile_dir

  return profile_dir


def initialise_cluster_profile(scprocess_dir, scdata_dir, setup_cfg, dryrun=False):
  setup_cfg = _normalise_profile_config(setup_cfg)
  if not _uses_cluster_profile(setup_cfg):
    return False

  template_name = setup_cfg['user']['profile_template']
  profile_name  = setup_cfg['user']['profile_name']
  template_dir  = scprocess_dir / 'profiles' / template_name
  template_f    = template_dir / 'config.yaml'
  profile_dir   = scdata_dir / 'profiles' / profile_name
  profile_f     = profile_dir / 'config.yaml'

  if profile_f.is_file():
    setup_cfg['user']['profile_dir'] = profile_dir
    return False
  if profile_dir.exists():
    raise FileNotFoundError(
      f"local profile directory {profile_dir} exists but does not contain config.yaml"
    )
  if not template_f.is_file():
    raise FileNotFoundError(f"profile template {template_f} does not exist")

  if dryrun:
    print(f"Would create profile template copy at {profile_dir}")
    return True

  profile_dir.parent.mkdir(parents=True, exist_ok=True)
  shutil.copytree(template_dir, profile_dir)
  setup_cfg['user']['profile_dir'] = profile_dir
  print(f"Created profile template copy at {profile_dir}")
  print(f"Edit {profile_f} if needed, then rerun `scprocess setup`.")
  return True


def get_conda_prefix(setup_cfg):
  """Get conda-prefix from the Snakemake profile, if configured."""
  if 'profile_dir' not in setup_cfg.get('user', {}):
    return None
  profile_f = setup_cfg['user']['profile_dir'] / 'config.yaml'
  if not profile_f.is_file():
    return None
  with open(profile_f) as f:
    profile = yaml.safe_load(f)
  return profile.get('conda-prefix')


def get_filtered_counts_file(run, ambient_method, amb_dir, date_stamp):
  if ambient_method == "cellbender":
    return f'{amb_dir}/ambient_{run}/bender_{run}_{date_stamp}_filtered.h5'
  if ambient_method in ["decontx_background", "decontx_cluster"]:
    return f'{amb_dir}/ambient_{run}/decontx_{run}_{date_stamp}_filtered.h5'
  if ambient_method == "none":
    return f'{amb_dir}/ambient_{run}/uncorrected_{run}_{date_stamp}_filtered.h5'
  raise ValueError(f"Unknown ambient_method '{ambient_method}' for filtered counts file.")


### much checking

# wrapper for checking setup
def check_setup_config(setup_cfg, schema_f, scprocess_dir):
  # start with defaults, overwrite with setup_cfg values
  schema      = _load_schema_file(schema_f)
  defaults    = _get_default_config_from_schema(schema)
  snakemake.utils.update_config(defaults, setup_cfg)
  setup_cfg   = defaults

  # check file is ok
  _validate_object_against_schema(setup_cfg, schema_f, "setup config")

  setup_cfg = _normalise_profile_config(setup_cfg)

  # check that all genome names are unique
  if ('ref_txomes' in setup_cfg) and ('custom' in setup_cfg['ref_txomes']):
    gen_names     = [ spec['name'] for spec in setup_cfg['ref_txomes']['custom'] ]
    not_unique    = len(set(gen_names)) != len(gen_names)
    if not_unique:
      raise KeyError("custom reference transcriptomes do not have unique names")

  return setup_cfg


# wrapper for checking
def check_config(config, schema_f, scdata_dir, scprocess_dir):
  # start with defaults, overwrite with config values
  schema      = _load_schema_file(schema_f)
  defaults    = _get_default_config_from_schema(schema, config)
  snakemake.utils.update_config(defaults, config)
  config      = defaults

  # check file is ok
  _validate_object_against_schema(config, schema_f, "config")

  # get parameters
  config      = _check_project_parameters(config, scdata_dir, scprocess_dir)
  config      = _check_arvados_parameters(config, scdata_dir)
  config      = _check_multiplexing_parameters(config)
  config      = _check_mapping_parameters(config, scdata_dir)
  config      = _check_ambient_parameters(config)
  config      = _check_qc_parameters(config)
  config      = _check_hvg_parameters(config)
  config      = _check_integration_parameters(config)
  config      = _check_marker_genes_parameters(config, scdata_dir)
  config      = _check_pb_empties_parameters(config)
  config      = _check_shiny_parameters(config)
  config      = _check_train_xgboost_parameters(config)

  return config


# get all default values from scheme file
def _load_schema_file(schema_f):
  # mess about w schema path
  schema_p  = pathlib.Path(schema_f)
  path_bits = list(schema_p.parts)
  if path_bits[0] == '..':
    path_bits = path_bits[1:]
  schema_p  = str(pathlib.Path(*path_bits))
  with open(schema_p, 'r') as f:
    schema = json.load(f)

  return schema


def _get_default_config_from_schema(schema, config=None):
  # extract defaults from the schema
  default_config = {}
  for key, props in schema.get('properties', {}).items():
    # if the key has a top-level default (uncommon for object types)
    if 'default' in props:
      default_config[key] = props['default']
    # if the key is an object, recursively extract its property defaults
    elif props.get('type') == 'object' and 'properties' in props:
      # skip sections absent from the user config whose required fields lack
      # defaults — injecting partial defaults for those causes validation errors
      if config is not None and key not in config and props.get('required'):
        sub_props = props.get('properties', {})
        if any('default' not in sub_props.get(r, {}) for r in props['required']):
          continue
      # create a nested dictionary of defaults for this section
      section_defaults = {}
      for sub_key, sub_props in props['properties'].items():
        if 'default' in sub_props:
          section_defaults[sub_key] = sub_props['default']

      if section_defaults:
        default_config[key] = section_defaults

  return default_config


# get validated config yaml list
def _validate_object_against_schema(config, schema_f, file_desc):
  try:
    # open schema file
    with open(schema_f, "r") as f:
      schema  = json.load(f)

    # validate the parsed yaml config against the json schema
    jsonschema.validate( instance = config, schema = schema )
    print(f"  {file_desc} file has correct format")

  except json.decoder.JSONDecodeError as e:
    print(f"problem with schema file:\n  {str(e)}")
    sys.exit(1)
  except jsonschema.ValidationError as e:
    best = best_match([e])
    path = ".".join(str(p) for p in best.absolute_path) if best.absolute_path else "(root)"
    print(f"problem with your {file_desc} file:\n  At '{path}': {best.message}")
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

  return config


# check parameters for project
def _check_project_parameters(config, scdata_dir, scprocess_dir):
  # do some path stuff
  config["project"]['proj_dir'] = pathlib.Path(config["project"]['proj_dir'])
  project_dc  = config["project"]

  # check fastq vs arvados
  if not project_dc['proj_dir'].is_dir():
    raise FileNotFoundError(f"proj_dir {project_dc['proj_dir']} is not a directory")

  # check project directory is ok
  _check_proj_dir_is_wflowr(config)

  # check fastq vs arvados
  has_fastq     = "fastq_dir" in project_dc
  has_arv_uuids = "arv_uuids" in project_dc
  if has_fastq + has_arv_uuids != 1:
    raise KeyError('"project" part of config file must contain exactly one of "fastq_dir" and "arv_uuids"')

  # do some checks if fastq_dir is specified
  if has_fastq and not has_arv_uuids:
    config["project"]["fastq_dir"] = _check_fastq_dir_or_cache(config["project"]["fastq_dir"], config, cache_glob="fastq_metadata_*.yaml")

  # check if selected ref_txome or probe_set is valid
  index_params_f    = scdata_dir / 'index_parameters.csv'
  index_params      = pl.read_csv(index_params_f)
  is_flex                      = config['project'].get('tenx_assay_type', 'poly_a') == 'flex'
  config['project']['is_flex'] = is_flex

  if is_flex:
    valid_probe_sets     = index_params.filter(pl.col('reference_type') == 'probe_set')['reference'].to_list()
    valid_probe_sets_str = ', '.join(valid_probe_sets)
    if not config['project']['probe_set'] in valid_probe_sets:
      raise ValueError(f"probe_set {config['project']['probe_set']} not defined. Valid values are {valid_probe_sets_str}")
  else:
    valid_ref_txomes    = index_params.filter(pl.col('reference_type') == 'ref_txome')['reference'].to_list()
    valid_ref_txome_str = ', '.join(valid_ref_txomes)
    if not config['project']['ref_txome'] in valid_ref_txomes:
      raise ValueError(f"ref_txome {config['project']['ref_txome']} not defined. Valid values are {valid_ref_txome_str}")

  # check whether date is given as datetime object
  date_regex    = re.compile("^20[0-9]{2}-[0-9]{2}-[0-9]{2}$")
  if not date_regex.match(config['project']['date_stamp']):
    raise ValueError(f"{config['project']['date_stamp']} does not match date format YYYY-MM-DD")

  # check samples
  config["project"]["sample_metadata"] = _check_path_exists_in_project(config["project"]["sample_metadata"], config, what = "file")

  # load it, do some checks
  samples_df  = pl.read_csv(config["project"]["sample_metadata"])
  _check_samples_df(samples_df, config)

  # check custom parameters file
  if 'custom_sample_params' in config['project']:
    # check path
    custom_f    = config["project"]["custom_sample_params"]
    custom_f    = _check_path_exists_in_project(custom_f, config, what = "file")

    # open and parse the yaml file
    with open(custom_f, "r") as f:
      custom_sample_params = yaml.safe_load(f)

    # load, check against schema
    schema_f    = scprocess_dir / "resources/schemas/custom_sample_params.schema.json"
    _validate_object_against_schema(custom_sample_params, schema_f, "custom sample parameter")

    # file was fine so store the output
    config["project"]["custom_sample_params"] = custom_f

  return config


def _check_arvados_parameters(config, scdata_dir):
  # read setup cfg to get arvados_setup (do not write it into project config)
  scdata_setup_f  = scdata_dir / 'scprocess_setup.yaml'
  arvados_setup   = None
  if scdata_setup_f.is_file():
    try:
      with open(scdata_setup_f, 'r') as sf:
        setup_cfg   = yaml.safe_load(sf)
        arv_dict    = setup_cfg.get('arvados', {})
    except Exception:
      arv_dict  = {}

  # check whether consistent
  if 'arv_uuids' in config['project'] and len(arv_dict) == 0:
    raise ValueError("arv_uuids specified in project config but no arvados section found in scprocess_setup.yaml")

  # add arvados parameters to config if present in setup
  if len(arv_dict) > 0:
    # check that the arv_uuids match the start of the arv_instance
    arv_instance = arv_dict.get('arv_instance', None)
    if 'arv_uuids' in config['project'] and arv_instance is not None:
      for uuid in config['project']['arv_uuids']:
        if not uuid.startswith(arv_instance):
          raise ValueError(f"arv_uuids specified in project config do not match the prefix of arv_instance specified in scprocess_setup.yaml") 

    # if ok then store
    config['arvados'] = {
      'arv_instance': arv_instance
    }

  return config


# check proj dir is wflowr
def _check_proj_dir_is_wflowr(config):
  # check that proj_dir is a workflowr directory 
  wflowr_fs_ls = ['_workflowr.yml', '.gitignore', '.Rprofile', '.gitattributes',
    'analysis/_site.yml', 'analysis/about.Rmd', 'analysis/index.Rmd', 'analysis/license.Rmd', 
    'public/.nojekyll']
  proj_dir  = config['project']['proj_dir']
  wflowr_fs_full_ls = [os.path.join(proj_dir, f) for f in wflowr_fs_ls]
  for f in wflowr_fs_full_ls:
    if not os.path.isfile(f):
      raise FileNotFoundError(f"proj_dir {config['project']['proj_dir']} has a missing file and isn't a workflowr project:\n  {f}\nYou can create a workflowr project using `scprocess newproj`")


# helper function
def _check_path_exists_in_project(path_to_check, config, what):
  # boring case
  if path_to_check is None:
    return path_to_check

  # store for error reporting
  tmp           = path_to_check
  path_to_check = pathlib.Path(path_to_check)

  # if not an absolute path, add project directory to it
  if not path_to_check.is_absolute():
    if "join" in config:
      path_to_check = pathlib.Path(config["join"]["proj_dir"]) / path_to_check
    else:
      path_to_check = config["project"]["proj_dir"] / path_to_check

  # check if directory or file
  if what == "dir":
    if not path_to_check.is_dir():
      raise FileNotFoundError(f"the directory {tmp} does not exist")
  elif what == "file":
    if not path_to_check.is_file():
      raise FileNotFoundError(f"the file {tmp} does not exist")
  else:
    raise ValueError()

  return path_to_check


def _check_fastq_dir_or_cache(fastq_dir_raw, config, cache_glob):
  """Check that a FASTQ directory exists, falling back to cached metadata from a prior mapping run."""
  try:
    return _check_path_exists_in_project(fastq_dir_raw, config, what = "dir")
  except FileNotFoundError:
    af_dir = pathlib.Path(config["project"]["proj_dir"]) / "output" / f"{config['project']['short_tag']}_mapping"
    if af_dir.is_dir() and any(af_dir.glob(cache_glob)):
      fastq_path = pathlib.Path(fastq_dir_raw)
      if not fastq_path.is_absolute():
        fastq_path = config["project"]["proj_dir"] / fastq_path
      warnings.warn(f"FASTQ directory {fastq_path} does not exist, but cached FASTQ metadata found.")
      return fastq_path
    raise


def _check_samples_df(samples_df, config):
  # check for sample_id
  if "sample_id" not in samples_df.columns:
    raise KeyError(f"'sample_id' not present in sample metadata file")

  # for flex data, check probe_id column is present
  if config['project'].get('tenx_assay_type', 'poly_a') == 'flex':
    if 'probe_id' not in samples_df.columns:
      raise KeyError("'probe_id' not present in sample metadata file; required for flex data")

  # check for pool_id
  if not config['multiplexing']['demux_type'] == "none":
    if "pool_id" not in samples_df.columns:
      raise KeyError(f"'pool_id' not present in sample metadata file")
    library_var = 'pool_id'
  else:
    library_var = 'sample_id'
  
  lib_ids = samples_df[library_var].unique().to_list()

  # check that libraries don't have '_R1/.R1' or '_R2/.R2' in their names 
  forbidden    = ["_R1", "_R2", ".R1", ".R2"]
  invalid_libs = [r for r in set(lib_ids) if any(sub in r for sub in forbidden)]

  if invalid_libs:
    raise ValueError(f"One or more {library_var} values contain '_R1/.R1' or '_R2/.R2'. Please ensure all elements exclude these substrings")

  # some checks for multiplexing
  if config['multiplexing']['demux_type'] == "hto":
    if "hto_id" not in samples_df.columns:
      raise KeyError("'hto_id' not present in sample metadata")

  if config['multiplexing']['demux_type'] == "ocm":
    if "ocm_id" not in samples_df.columns:
      raise KeyError("'ocm_id' not present in sample metadata when demux_type is 'ocm'")

  # check that sample_id values are unique
  if not samples_df[ "sample_id" ].n_unique() == samples_df.shape[0]:
    raise ValueError("'sample_id' values in metadata csv not unique")

  # check minimum sample counts per pool (or globally for demux_type=none)
  demux_type = config['multiplexing']['demux_type']
  if demux_type == 'none':
    n_samples = samples_df['sample_id'].n_unique()
    if n_samples < 2:
      raise ValueError(f"At least 2 sample_ids are required; found {n_samples}.")
  else:
    pool_counts = samples_df.group_by('pool_id').agg(pl.col('sample_id').n_unique().alias('n_samples'))
    small_pools = pool_counts.filter(pl.col('n_samples') < 2)
    if small_pools.shape[0] > 0:
      pool_list = ', '.join(small_pools['pool_id'].to_list())
      raise ValueError(f"Each pool must contain at least 2 sample_ids. The following pools have fewer: {pool_list}")

  # check columns of samples_df
  if any(' ' in col for col in samples_df.columns):
    raise ValueError("some column names in metadata csv contain spaces.")
  
  # check that sample_ids or pool_ids are not overlapping
  _check_lib_ids(lib_ids, library_var)

  # sort out metadata variables
  if 'metadata_vars' in config["project"]:
    # load up sample file
    for var in config["project"]["metadata_vars"]:
      # check variable exists
      if not var in samples_df.columns:
        raise KeyError(f"{var} not a column in sample metadata file")

      # some checks on metadata variables
      var_col     = samples_df[var]
      if var_col.dtype == pl.String:
        if all(var_col == "NA"):
          print(f"  {var} variable has all values 'NA'; probably not a useful variable")
        continue
      # check that there are less than 10 unique values (otherwise probably not a categorical variable)
      if var_col.unique().shape[0] > 10:
        raise ValueError(f"{var} variable has more than 10 unique values; prob not a categorical variable")
  else:
    config['project']['metadata_vars'] = []

  return


def _check_lib_ids(lib_ids, lib_var):
  lib_overlaps = []

  # compare every library id to every other library id
  for i, lib in enumerate(lib_ids):
    for j, other_lib in enumerate(lib_ids):
      if i != j and lib in other_lib:
        lib_overlaps.append((lib, other_lib))

  # if libraries overlap print error
  if lib_overlaps:
    msg = f"The following {lib_var} values are problematic (one is a subset of the other):\n"
    for lib, other_lib in lib_overlaps:
      msg += f"  - '{lib}' is a subset of '{other_lib}'\n"
    raise ValueError(msg)

  return


# check parameters for multiplexing
def _check_multiplexing_parameters(config):
  # load up samples
  samples_df  = pl.read_csv(config['project']['sample_metadata'])

  # do some things if demux_type is hto
  if config['multiplexing']['demux_type'] == 'none':
    return config

  elif config['multiplexing']['demux_type'] == 'hto':
    # check feature ref specified and valid
    config["multiplexing"]["feature_ref"] = _check_path_exists_in_project(config["multiplexing"]["feature_ref"], config, what = "file")

    # check fastq vs arvados
    has_fastq     = "fastq_dir" in config["multiplexing"]
    has_arv_uuids = "arv_uuids" in config["multiplexing"]
    if has_fastq + has_arv_uuids != 1:
      raise KeyError('"multiplexing" part of config file must contain exactly one of "fastq_dir" and "arv_uuids"')

    # do some checks if fastq_dir is specified
    if has_fastq and not has_arv_uuids:
      config["multiplexing"]["fastq_dir"] = _check_fastq_dir_or_cache(config["multiplexing"]["fastq_dir"], config, cache_glob="fastq_metadata_hto_*.yaml")

    # check for columns in feature ref
    feat_ref_df   = pl.read_csv(config["multiplexing"]["feature_ref"])
    if not all(col in feat_ref_df.columns for col in ["hto_id", "sequence"]):
      raise KeyError("'hto_id' and 'sequence' must both be columns in the feature_ref file")

    # check that all hto_ids in feature ref file are unique
    if not feat_ref_df["hto_id"].n_unique() == feat_ref_df.shape[0]:
      raise ValueError("hto_id values in feature reference file not unique")
    hto_ids = feat_ref_df["hto_id"].to_list()

    # check if all hto_id values in sample metadata match the ones in feature reference
    if not all(hto in hto_ids for hto in list(set(samples_df["hto_id"]))):
      raise ValueError("One or more hto_id values in sample_metadata don't match hto_id values in the feature reference file")

  elif config['multiplexing']['demux_type'] == 'custom':
    # check specified file is ok
    config['multiplexing']['demux_output'] = _check_path_exists_in_project(config['multiplexing']['demux_output'], config, what = "file")

    # check columns look ok (value-matching checks are done in the check_demux_ids rule)
    demux_df    = pl.read_csv(config['multiplexing']['demux_output'], n_rows = 10)
    for col in ["pool_id", "sample_id", "cell_id"]:
      if not col in demux_df.columns:
        raise KeyError(f"{col} not present in demux_output")

    # check if samples in metadata and demux_df match
    if set(demux_df['sample_id']) > set(samples_df['sample_id']):
      raise ValueError("Some values for 'sample_id' in demux_output don't have a match in sample_metadata")
    # if set(samples_df['sample_id']) > set(demux_df['sample_id']):
    #   raise ValueError("Some values for 'sample_id' in sample_metadata don't have a match in demux_output")
    if set(demux_df['pool_id']) != set(samples_df['pool_id']):
      raise ValueError("Values for pool_id don't match across demux_output and sample_metadata")

  elif config['multiplexing']['demux_type'] == 'flex':
    pass

  elif config['multiplexing']['demux_type'] == 'ocm':
    # validate ocm_id values
    valid_ocm_ids = {"OB1", "OB2", "OB3", "OB4"}
    ocm_ids = set(samples_df["ocm_id"].to_list())
    invalid = ocm_ids - valid_ocm_ids
    if invalid:
      raise ValueError(f"Invalid ocm_id values: {invalid}. Valid values are {valid_ocm_ids}")

    # check uniqueness of ocm_id within each pool
    pool_ocm = samples_df.select("pool_id", "ocm_id")
    dupes = pool_ocm.group_by("pool_id", "ocm_id").len().filter(pl.col("len") > 1)
    if dupes.shape[0] > 0:
      raise ValueError(f"Duplicate ocm_id values within pools:\n{dupes}")

    # check max 4 samples per pool
    pool_counts = samples_df.group_by("pool_id").len()
    over_4 = pool_counts.filter(pl.col("len") > 4)
    if over_4.shape[0] > 0:
      raise ValueError(f"OCM supports max 4 samples per pool, but these pools have more: {over_4['pool_id'].to_list()}")

    # check OCM overhang map file exists
    scdata_dir = pathlib.Path(os.getenv('SCPROCESS_DATA_DIR'))
    ocm_overhang_f = scdata_dir / 'cellranger_ref' / 'ocm_overhang_map.txt'
    if not ocm_overhang_f.is_file():
      raise FileNotFoundError(
        f"OCM overhang map file not found: {ocm_overhang_f}\n"
        f"This file should be created during scprocess setup.")

    config['multiplexing']['ocm_overhang_f'] = str(ocm_overhang_f)

  return config


def check_demux_ids(demux_f, metadata_f, check_f):
  demux_ids    = pl.scan_csv(demux_f).select('pool_id', 'sample_id').unique().collect()
  metadata_ids = pl.read_csv(metadata_f).select('pool_id', 'sample_id')

  demux_pools = set(demux_ids['pool_id'])
  meta_pools  = set(metadata_ids['pool_id'])
  if demux_pools != meta_pools:
    missing_in_demux = meta_pools - demux_pools
    missing_in_meta  = demux_pools - meta_pools
    msg = 'Values for pool_id do not match across demux_output and sample_metadata.'
    if missing_in_demux:
      msg += f' In metadata but not demux: {missing_in_demux}'
    if missing_in_meta:
      msg += f' In demux but not metadata: {missing_in_meta}'
    raise ValueError(msg)

  demux_samples = set(demux_ids['sample_id'])
  meta_samples  = set(metadata_ids['sample_id'])
  if demux_samples > meta_samples:
    missing = demux_samples - meta_samples
    raise ValueError(f"Some values for 'sample_id' in demux_output do not have a match in sample_metadata: {missing}")

  open(check_f, 'w').close()


# check parameters for mapping
def _check_mapping_parameters(config, scdata_dir):
  # load index parameters
  idx_params_f  = scdata_dir / 'index_parameters.csv'
  index_params  = pl.read_csv(idx_params_f)

  config['mapping_af'] = {}
  config['mapping_af']['alevin_fry_home'] = scdata_dir / 'alevin_fry_home'
  config['mapping_af']['wl_lu_f']         = scdata_dir / 'cellranger_ref/cellranger_whitelists.csv'
  
  tenx_chemistry = config['project']['tenx_chemistry']
  
  if config['project']['tenx_assay_type'] == 'flex':
    probe_set      = config['project']['probe_set']
    if tenx_chemistry == "none":
      tenx_chemistry = 'flexv1' if probe_set in ['human_v1', 'mouse_v1'] else 'flexv2'
    else:
      #check that specified chemistry matches the probe set, if not raise warning and change to correct one
      if probe_set in ['human_v1', 'mouse_v1'] and tenx_chemistry != 'flexv1':
        warnings.warn(f"Specified tenx_chemistry '{tenx_chemistry}' does not match the probe set '{probe_set}'. Changing tenx_chemistry to 'flexv1'.")
        tenx_chemistry = 'flexv1'
      elif probe_set not in ['human_v1', 'mouse_v1'] and tenx_chemistry != 'flexv2':
        warnings.warn(f"Specified tenx_chemistry '{tenx_chemistry}' does not match the probe set '{probe_set}'. Changing tenx_chemistry to 'flexv2'.")
        tenx_chemistry = 'flexv2'
    
    # look up mito string and gene info from index_parameters.csv
    idx_row     = index_params.filter((pl.col('reference') == probe_set))
    mito_str    = idx_row['mito_str'][0] # this might not be necessary for flex
    gene_info_f = idx_row['gene_info_f'][0]

    # probe set specific paths
    probe_idx_dir  = scdata_dir / 'alevin_fry_home' / 'probe_sets' / probe_set
    probe_set_f    = scdata_dir / 'probe_sets' / probe_set / f'{probe_set}_probe_set.csv'
    probe_bcs_f    = scdata_dir / 'cellranger_ref' / f'cellranger_probe_barcodes_{tenx_chemistry}.tsv' # this has to be downloaded in the setup step!!!

    if not probe_idx_dir.is_dir():
      raise FileNotFoundError(f"alevin index for probe set '{probe_set}' doesn't exist at {probe_idx_dir}")
    if not probe_set_f.is_file():
      raise FileNotFoundError(f"probe set CSV '{probe_set_f}' doesn't exist")
    if not probe_bcs_f.is_file():
      raise FileNotFoundError(f"probe barcodes file '{probe_bcs_f}' doesn't exist.")

    config['mapping_af']['is_flex']        = True
    config['mapping_af']['tenx_chemistry'] = tenx_chemistry
    config['mapping_af']['af_index_dir']   = probe_idx_dir
    config['mapping_af']['probeset_f']    = probe_set_f
    config['mapping_af']['probe_bcs_f']    = probe_bcs_f
    config['mapping_af']['af_mito_str']    = mito_str
    config['mapping_af']['gene_info_f']    = gene_info_f
    config['mapping_af']['geometry']       = '1{b[16]u[12]x[0-3]hamming(f[TTGCTAGGACCG],1)s[10]x:}2{r:}' if tenx_chemistry == 'flexv2' else ''

  else:
    ref_txome = config['project']['ref_txome']
  
    config['mapping_af']['is_flex']        = False
    config['mapping_af']['tenx_chemistry'] = tenx_chemistry
    config['mapping_af']['af_index_dir'] = scdata_dir / 'alevin_fry_home' / 'ref_txomes' / ref_txome
    if not pathlib.Path(config['mapping_af']['af_index_dir']).is_dir():
      raise FileNotFoundError(f"alevin index for '{ref_txome}' doesn't exist at: {config['mapping_af']['af_index_dir']}")

    idx_row = index_params.filter( (pl.col('reference') == ref_txome))
    config['mapping_af']['af_mito_str'] = idx_row['mito_str'][0]
    config['mapping_af']['gene_info_f'] = idx_row['gene_info_f'][0]

  return config


# check parameters for ambient
def _check_ambient_parameters(config):
  # get cellbender image (maybe skip this if cellbender is not selected?)
  if config['ambient']['cb_version']   == 'v0.3.2':
    cellbender_image  = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.2'
  elif config['ambient']['cb_version'] == 'v0.3.0':
    cellbender_image  = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.0'
  elif config['ambient']['cb_version'] == 'v0.2.0':
    cellbender_image  = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.2.0'
  else:
    raise ValueError(f"selected cellbender version {config['ambient']['cb_version']} not supported")
  config['ambient']['cellbender_image'] = cellbender_image

  # check posterior batch size
  if config['ambient']['cb_version'] == 'v0.3.2':
    if not 'cb_posterior_batch_size' in config['ambient']:
      config['ambient']['cb_posterior_batch_size'] = 128

  return config


# check parameters for qc
def _validate_qc_bounds(qc, context="qc"):
  for min_key, max_key in (
      ('qc_min_counts', 'qc_max_counts'),
      ('qc_min_feats', 'qc_max_feats')):
    min_val = qc.get(min_key)
    max_val = qc.get(max_key)
    if min_val is not None and max_val is not None and max_val < min_val:
      raise ValueError(
        f"{context}: {max_key} ({max_val}) must be greater than or equal to "
        f"{min_key} ({min_val})"
      )


def _check_qc_parameters(config):
  # define some hard values; maybe move these to schema?
  QC_HARD_MIN_COUNTS  = 200
  QC_HARD_MIN_FEATS   = 100
  QC_HARD_MAX_MITO    = 0.5

  # some checks
  config['qc']['exclude_mito']        = _safe_boolean(config['qc']['exclude_mito'])
  _validate_qc_bounds(config['qc'])

  # make sure they're consistent
  config['qc']['qc_hard_min_counts']  = min(QC_HARD_MIN_COUNTS, config['qc']['qc_min_counts'])
  config['qc']['qc_hard_min_feats']   = min(QC_HARD_MIN_FEATS, config['qc']['qc_min_feats'])
  config['qc']['qc_hard_max_mito']    = max(QC_HARD_MAX_MITO, config['qc']['qc_max_mito'])

  return config


# check parameters for pb and empties
def _check_pb_empties_parameters(config):
  # nothing to do here at the moment; leaving in case it's useful later

  return config


# load valid palette names from the shared resource file (mirrors shiny.R)
_VALID_PALETTES_FILE = pathlib.Path(__file__).parent.parent / "resources" / "valid_palettes.json"

def _load_valid_palette_names():
  with open(_VALID_PALETTES_FILE) as f:
    groups = json.load(f)
  return {name for key, names in groups.items() if not key.startswith('_') for name in names}


def _check_palette_name(name, context, valid_names):
  if name not in valid_names:
    raise ValueError(
      f"Unknown palette '{name}' in {context}. "
      f"See the scprocess documentation for valid palette names."
    )


def _load_sample_metadata(config):
  """Load sample metadata from project config or all source projects in a join config."""
  if 'project' in config:
    meta_f = config['project'].get('sample_metadata')
    if meta_f and pathlib.Path(meta_f).is_file():
      return pl.read_csv(meta_f)
  elif 'join' in config and 'projects' in config:
    meta_dfs = []
    for pid, pcfg in config['projects'].items():
      pcfg_f = pcfg.get('config')
      if pcfg_f and pathlib.Path(pcfg_f).is_file():
        with open(pcfg_f) as f:
          proj_cfg = yaml.safe_load(f)
        proj_meta_f = proj_cfg.get('project', {}).get('sample_metadata')
        if proj_meta_f:
          proj_dir = pathlib.Path(proj_cfg['project']['proj_dir'])
          proj_meta_path = pathlib.Path(proj_meta_f)
          if not proj_meta_path.is_absolute():
            proj_meta_path = proj_dir / proj_meta_path
          if proj_meta_path.is_file():
            meta_dfs.append(pl.read_csv(proj_meta_path))
    if meta_dfs:
      return pl.concat(meta_dfs, how="diagonal_relaxed")
  return None


def _check_metadata_palette_values(shiny_cfg, config, metadata_vars, valid_palette_names):
  meta_df = _load_sample_metadata(config)

  for var, spec in shiny_cfg['metadata_palettes'].items():
    if metadata_vars and var not in metadata_vars:
      warnings.warn(
        f"shiny.metadata_palettes key '{var}' is not in metadata_vars; it will be ignored")
    if isinstance(spec, str):
      _check_palette_name(spec, f'shiny.metadata_palettes.{var}', valid_palette_names)
    elif isinstance(spec, dict):
      if spec.get('palette'):
        _check_palette_name(spec['palette'], f'shiny.metadata_palettes.{var}.palette', valid_palette_names)
      if spec.get('values') and meta_df is not None and var in meta_df.columns:
        actual_values = set(meta_df[var].drop_nulls().unique().to_list())
        specified_values = set(spec['values'])
        missing = actual_values - specified_values
        if missing:
          raise ValueError(
            f"shiny.metadata_palettes.{var}.values is missing values found in sample metadata: "
            f"{sorted(missing)}")


# check parameters for shiny app build
def _check_shiny_parameters(config):

  if 'shiny' not in config:
    return config

  shiny_cfg = config['shiny']

  # resolve metadata_vars from whichever config level is present
  if 'project' in config:
    base_metadata_vars = config['project'].get('metadata_vars', [])
  elif 'join' in config:
    base_metadata_vars = config['join'].get('metadata_vars', [])
  else:
    base_metadata_vars = []
  metadata_vars = shiny_cfg.get('metadata_vars', base_metadata_vars)

  valid_palette_names = _load_valid_palette_names()

  # check var_names length matches metadata_vars
  if shiny_cfg.get('var_names') is not None:
    if len(shiny_cfg['var_names']) != len(metadata_vars):
      raise ValueError(
        f"shiny.var_names has {len(shiny_cfg['var_names'])} entries but "
        f"metadata_vars has {len(metadata_vars)}; they must be the same length"
      )

  # check metadata_combns values are in metadata_vars
  if shiny_cfg.get('metadata_combns'):
    for pair in shiny_cfg['metadata_combns']:
      for v in pair:
        if v not in metadata_vars:
          raise ValueError(
            f"shiny.metadata_combns references '{v}' which is not in metadata_vars: "
            f"{metadata_vars}"
          )

  # check cluster_palette name
  if shiny_cfg.get('cluster_palette'):
    _check_palette_name(shiny_cfg['cluster_palette'], 'shiny.cluster_palette', valid_palette_names)

  # check metadata_palettes
  if shiny_cfg.get('metadata_palettes'):
    _check_metadata_palette_values(shiny_cfg, config, metadata_vars, valid_palette_names)

  # check home_md if specified
  if shiny_cfg.get('home_md'):
    home_md_f = _check_path_exists_in_project(shiny_cfg['home_md'], config, what='file')
    # verify file is readable as UTF-8 and non-empty
    try:
      content = pathlib.Path(home_md_f).read_text(encoding='utf-8')
    except UnicodeDecodeError:
      raise ValueError(f"home_md file '{shiny_cfg['home_md']}' is not valid UTF-8 text")
    if not content.strip():
      raise ValueError(f"home_md file '{shiny_cfg['home_md']}' is empty")
    config['shiny']['home_md'] = home_md_f

  # check annotation_csv if specified
  if shiny_cfg.get('annotation_csv'):
    annot_f  = _check_path_exists_in_project(shiny_cfg['annotation_csv'], config, what='file')
    annot_df = pl.read_csv(annot_f)

    # check required columns
    required_cols = ['cluster', 'cluster_name']
    missing_cols  = [c for c in required_cols if c not in annot_df.columns]
    if missing_cols:
      raise KeyError(
        f"annotation_csv must contain columns: {', '.join(missing_cols)}. "
        f"Found: {', '.join(annot_df.columns)}"
      )

    # check no duplicate cluster values
    if annot_df['cluster'].n_unique() < annot_df.shape[0]:
      raise ValueError("annotation_csv has duplicate values in the 'cluster' column")

    # warn on non-hex values in optional colour column
    if 'colour' in annot_df.columns:
      hex_re        = re.compile(r'^#[0-9A-Fa-f]{6}$')
      invalid_cols  = [
        str(c) for c in annot_df['colour'].drop_nulls().to_list()
        if not hex_re.match(str(c))
      ]
      if invalid_cols:
        warnings.warn(
          f"annotation_csv 'colour' column contains values that are not hex colours (#RRGGBB): "
          f"{', '.join(invalid_cols[:5])}"
          + (f" (and {len(invalid_cols) - 5} more)" if len(invalid_cols) > 5 else "")
        )

    config['shiny']['annotation_csv'] = annot_f

  return config


def get_shiny_app_tag(config):
  """Return the configured main-app tag or its context-specific default."""
  shiny_cfg = config.get('shiny', {})
  if 'join' in config:
    default_tag = config['join']['name']
  elif 'project' in config:
    default_tag = config['project']['short_tag']
  else:
    raise KeyError("Config must contain either a 'project' or 'join' section")
  return shiny_cfg.get('app_tag', default_tag)


# check parameters for hvgs
def _check_hvg_parameters(config):
  # for hto/custom with only 1 pool, ambient gene detection is not possible
  # (ambient genes are calculated across pseudobulks, which are per pool for these demux types)
  demux_type = config['multiplexing']['demux_type']
  if demux_type in ['hto', 'custom'] and config['hvg']['hvg_exclude_ambient_genes']:
    meta_df = pl.read_csv(config['project']['sample_metadata'])
    n_pools = meta_df['pool_id'].n_unique()
    if n_pools < 2:
      config['hvg']['hvg_exclude_ambient_genes'] = False
      print(
        f"  WARNING: hvg_exclude_ambient_genes is True but ambient gene detection requires "
        f">=2 pools; only {n_pools} pool found with demux_type='{demux_type}'. "
        f"Setting hvg_exclude_ambient_genes=False. To suppress this warning, set "
        f"hvg_exclude_ambient_genes: false in your config."
      )

  # check if any genes were specified to be excluded
  if 'hvg_exclude_from_file' in config['hvg']:
    # check file exists
    exc_f       = pathlib.Path(config['hvg']['hvg_exclude_from_file'])
    if not exc_f.is_file():
      raise FileNotFoundError("file specified in 'hvg_exclude_from_file' does not exist")

    # check file has correct columns
    exc_df      = pl.read_csv(exc_f)
    if not ((exc_df.columns == ['gene_id']) | (exc_df.columns == ['symbol'])):
      raise KeyError("file specified in 'hvg_exclude_from_file' must have exactly one column, called either 'gene_id' or 'symbol'")

    # check for duplicates
    gene_col    = exc_df.columns[0]
    exc_vals    = exc_df[ gene_col ]
    if exc_vals.n_unique() < len(exc_vals):
      raise ValueError("duplicated values found in file specified in 'hvg_exclude_from_file'")

    # check values are in relevant ref genome
    gtf_df      = pl.read_csv(config['mapping_af']['gene_info_f'], separator = "\t")
    all_vals    = gtf_df[ gene_col ]
    absent_vals = set(exc_vals) - set(all_vals)
    if len(absent_vals) > 0:
      raise ValueError(f"the following genes were specified in 'hvg_exclude_from_file' but were not found in the reference transcriptome: {', '.join(absent_vals)}")
  else:
    config['hvg']['hvg_exclude_from_file'] = None

  # define dummy group names for all
  if config['hvg']['hvg_method'] == 'all':
    config['hvg']['hvg_group_names']        = ['all_samples']
    config['hvg']['hvg_metadata_split_var'] = None

  # if groups, check that the values are ok
  elif config['hvg']['hvg_method'] == 'groups':
    # check that value of metadata_split_var matches a column in sample metadata
    hvg_split_var = config['hvg']['hvg_metadata_split_var']
    meta_df       = pl.read_csv(config['project']['sample_metadata'])
    if not hvg_split_var in meta_df.columns:
      raise KeyError(f"{hvg_split_var} is not a column in the sample metadata file.")
    
    # check number of unique group values
    uniq_groups = meta_df[ hvg_split_var ].unique().to_list()
    if len(uniq_groups) == meta_df.shape[0]:
      raise ValueError(f"Number of unique values in '{hvg_split_var}' is the same as the number of samples.")

    # store nice names
    config['hvg']['hvg_group_names'] = [g.replace(" ", "_") for g in uniq_groups]

  # get number of gene chunks if method is 'groups' or 'all'
  if config['hvg']['hvg_method'] in ['groups', 'all']:
    # get total number of genes
    gtf_df      = pl.read_csv(config['mapping_af']['gene_info_f'], separator = "\t")
    num_genes   = gtf_df.shape[0]

    # chunk them up and name them
    num_chunks  = (num_genes + config['hvg']['hvg_chunk_size'] - 1) // config['hvg']['hvg_chunk_size']
    chunk_names = [f"chunk_{i+1}" for i in range(num_chunks)]
    
    # add to config
    config['hvg']['hvg_num_chunks']   = num_chunks
    config['hvg']['hvg_chunk_names']  = chunk_names

  return config


# check parameters for integration
def _check_integration_parameters(config):
  # if paga is specified, check that the specified leiden resolution is in the list of valid values
  if config['integration']['int_use_paga']:
    valid_res_ls  = config['integration']['int_res_ls']
    if not config['integration']['int_paga_cl_res'] in valid_res_ls:
      raise ValueError(f"leiden/louvain resolution specified for paga integration must be one of {valid_res_ls}")

  return config


# check parameters for marker genes
def _check_marker_genes_parameters(config, scdata_dir):
  # set some more default values
  config['marker_genes']['mkr_gsea_dir'] = scdata_dir / 'gmt_pathways'

  # get custom marker files
  proj_dir = config['project']['proj_dir']
  custom_mkr_names, custom_mkr_paths = _get_custom_marker_genes_specs(
    config, scdata_dir, proj_dir)
  config['marker_genes']['custom_mkr_names'] = custom_mkr_names
  config['marker_genes']['custom_mkr_paths'] = custom_mkr_paths

  return config


# check specified custom marker genes
def _get_custom_marker_genes_specs(config, scdata_dir, proj_dir):
  custom_mkr_names = ""
  custom_mkr_paths = ""

  if 'mkr_custom_genesets' in config["marker_genes"]:
    mkr_names = []
    mkr_paths = []
    for i, gene_set in enumerate(config["marker_genes"]["mkr_custom_genesets"]):
      name      = gene_set["name"]
      file_path = pathlib.Path(gene_set.get("file", scdata_dir / 'marker_genes' / f"{name}.csv"))

      if not file_path.is_absolute():
        file_path = proj_dir / file_path
      if not file_path.is_file():
        raise FileNotFoundError(f"File not found for marker set '{name}'")
      if not file_path.suffix == ".csv":
        raise ValueError(f"File for custom marker set '{name}' is not a csv file")

      # check csv file contents
      mkrs_df   = pl.read_csv(file_path)
      req_col   = "label"
      opt_cols  = ["symbol", "ensembl_id"]
      if not req_col in mkrs_df.columns:
        raise KeyError(f"File '{file_path}' is missing the mandatory column 'label'.")
      if not any(col in mkrs_df.columns for col in opt_cols):
        raise KeyError(f"File '{file_path}' must contain at least one of 'symbol' or 'ensembl_id' column.")
      if "symbol" in mkrs_df.columns and any(mkrs_df["symbol"].is_duplicated()):
        raise KeyError(f"File '{file_path}' cannot have any duplicated values in the 'symbol' column")
      if "ensembl_id" in mkrs_df.columns and any(mkrs_df["ensembl_id"].is_duplicated()):
        raise KeyError(f"File '{file_path}' cannot have any duplicated values in the 'ensembl_id' column")

      # Store validated values
      mkr_names.append(name)
      mkr_paths.append(str(file_path))

    custom_mkr_names = ",".join(mkr_names)
    custom_mkr_paths = ",".join(mkr_paths)
  
  return custom_mkr_names, custom_mkr_paths


# get parameters for zoom
def get_zoom_parameters(config, zoom_schema_f, scdata_dir):
  # if (rule_name != 'zoom') or ('zoom' not in config) or (config['zoom'] is None):
  if 'zoom' not in config:
    ZOOM_PARAMS   = {}
  else:
    # get names and files
    zoom_yamls    = [ pathlib.Path(f) for f in config['zoom']]

    # make dictionary of zoom params from yamls
    zoom_ls       = [_get_one_zoom_parameters(zoom_yaml_f, zoom_schema_f, copy.deepcopy(config)) 
                      for zoom_yaml_f in zoom_yamls]
    zoom_ns       = [z['zoom']['name'] for z in zoom_ls]
    if len(zoom_ns) != len(set(zoom_ns)):
      raise ValueError("names in specified zoom parameter yaml files are not unique")
    ZOOM_PARAMS   = {z['zoom']['name']: z for z in zoom_ls}

  return ZOOM_PARAMS


def _warn_zoom_qc_thresholds(zoom_qc, main_qc, zoom_name):
  checks = [
    ('qc_min_counts', 'higher'), ('qc_min_feats', 'higher'),
    ('qc_max_counts', 'lower'), ('qc_max_feats', 'lower'),
    ('qc_max_mito', 'lower'), ('qc_min_mito', 'higher'),
    ('qc_max_splice', 'lower'), ('qc_min_splice', 'higher'),
  ]
  for key, direction in checks:
    zoom_val = zoom_qc.get(key)
    main_val = main_qc.get(key)
    if zoom_val is None or main_val is None:
      continue
    is_less_strict = (
      (direction == 'higher' and zoom_val < main_val) or
      (direction == 'lower' and zoom_val > main_val)
    )
    if is_less_strict:
      warnings.warn(
        f"zoom '{zoom_name}': {key}={zoom_val} is less strict than main config ({main_val}). "
        f"This will have no effect since cells already passed main QC."
      )


# get parameters for one zoom specification
def _get_one_zoom_parameters(zoom_yaml_f, zoom_schema_f, config):
  # check file exists
  zoom_yaml_f   = _check_path_exists_in_project(zoom_yaml_f, config, what = "file")

  # load things specified by zoom
  with open(zoom_yaml_f, "r") as stream:
    zoom_config   = yaml.safe_load(stream)

  # update with zoom defaults if not specified
  zoom_schema   = _load_schema_file(zoom_schema_f)
  zoom_defaults = _get_default_config_from_schema(zoom_schema, zoom_config)
  snakemake.utils.update_config(zoom_defaults, zoom_config)
  zoom_config   = zoom_defaults

  # check file is ok
  _validate_object_against_schema(zoom_config, zoom_schema_f, "zoom config")

  # start with defaults, overwrite with config values
  defaults      = config.copy()
  defaults.pop('hvg', None)
  snakemake.utils.update_config(defaults, zoom_config)
  zoom_config   = defaults

  # check hvgs option
  zoom_config   = _check_hvg_parameters(zoom_config)
  zoom_config   = _check_integration_parameters(zoom_config)
  zoom_config   = _check_shiny_parameters(zoom_config)
  zoom_config   = _check_train_xgboost_parameters(zoom_config)

  if 'train_xgboost' in zoom_config:
    _zxgb = zoom_config['train_xgboost']
    if 'ref_tag' not in _zxgb:
      _zxgb['ref_tag'] = f"xgboost_{config['project']['full_tag']}_zoom_{zoom_config['zoom']['name']}"
    else:
      _zxgb['ref_tag'] = f"xgboost_{_zxgb['ref_tag']}"

  # get useful things
  SHORT_TAG     = config['project']['short_tag']
  FULL_TAG      = config['project']['full_tag']
  DATE_STAMP    = config['project']['date_stamp']

  # find file for each option
  if zoom_config['zoom']['labels_source'] == 'clusters':
    labels_f      = f"output/{SHORT_TAG}_integration/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz"

  # if using xgboost or celltypist, check those things
  elif zoom_config['zoom']['labels_source'] in ['celltypist', 'scprocess']:
    labeller      = zoom_config['zoom']['labels_source']
    model         = zoom_config['zoom']['model']
    labels_f      = f"output/{SHORT_TAG}_label_celltypes/labels_{labeller}_model_{model}_{FULL_TAG}_{DATE_STAMP}.csv.gz"

  # unpack
  elif zoom_config['zoom']['labels_source'] == 'custom':
    labels_f      = pathlib.Path(zoom_config['zoom']['custom_labels_f'])
  
  # check file exists
  labels_f      = _check_path_exists_in_project(labels_f, config, what = "file")

  # get selected labels from config (cluster validation deferred to check rule)
  sel_labels    = zoom_config['zoom']['sel_labels']

  # store original labels path for the filtering rule's input
  zoom_config['zoom']['_original_labels_f'] = labels_f
  zoom_config['zoom']['sel_labels'] = sel_labels

  # point labels_f to the filtered labels output (produced by zoom_filter_cells_qc rule)
  zoom_name = zoom_config['zoom']['name']
  zoom_config['zoom']['labels_f'] = config['project']['proj_dir'] / \
    f"output/{SHORT_TAG}_zoom/{zoom_name}/filtered_labels_{FULL_TAG}_{zoom_name}_{DATE_STAMP}.csv.gz"

  # warn if any zoom QC thresholds are less strict than main config (they'd be no-ops)
  _validate_qc_bounds(zoom_config.get('qc', {}), f"zoom '{zoom_name}' qc")
  _warn_zoom_qc_thresholds(zoom_config.get('qc', {}), config.get('qc', {}), zoom_name)

  return zoom_config


def check_config_ok_for_rule(config, rule):
  if rule == 'label_celltypes':
    if 'label_celltypes' not in config:
      raise KeyError(
        "no 'label_celltypes' section found in config file. "
        "Add a 'label_celltypes' block to your config to use -r label_celltypes.")
  if rule == 'train_xgboost':
    if 'train_xgboost' not in config:
      raise KeyError(
        "no 'train_xgboost' section found in config file. "
        "Add a 'train_xgboost' block (requires at least 'annots_f') to your config to use -r train_xgboost.")
  if rule == 'zoom':
    if 'zoom' not in config:
      raise KeyError(
        "no 'zoom' section found in config file. "
        "Add a 'zoom' block to your config to use -r zoom.")


# get variables for each library
def get_lib_parameters(config, scprocess_data_dir):
  # define run variable
  if config['multiplexing']['demux_type'] == "none":
    LIB_VAR       = "sample_id"
  else:
    LIB_VAR       = "pool_id"
  
  # get all libraries
  metadata_f  = config["project"]["sample_metadata"]
  samples_df  = pl.read_csv( metadata_f )
  LIBS        = samples_df[ LIB_VAR ].drop_nulls().unique().to_list()

  # get any custom parameters
  custom_lib_params = _get_custom_parameters(config, LIB_VAR)

  # should we exclude any libraries?
  LIBS        = _do_exclusions(LIBS, config, LIB_VAR)

  # get fastq files
  RNA_FQS     = _get_fastqs(config, LIBS, is_hto = False)
  missing_libs = [l for l in LIBS if l not in RNA_FQS]
  if missing_libs:
    raise ValueError(f"no FASTQ files found for the following libraries: {missing_libs}")
  LIBS        = list(RNA_FQS.keys())
  if len(LIBS) == 0:
    raise ValueError("no libraries with FASTQs")

  # get HTO files
  HTO_FQS     = {}
  if config['multiplexing']['demux_type'] == "hto":
    HTO_FQS     = _get_fastqs(config, LIBS, is_hto = True)
    LIBS        = [r for r in LIBS if r in HTO_FQS]

  # load sample file, populate everything from config
  LIB_PARAMS  = {
    lib_name: _get_lib_parameters_one_lib(lib_name, config, RNA_FQS, HTO_FQS, scprocess_data_dir, custom_lib_params)
    for lib_name in sorted(LIBS)
  }

  return LIB_PARAMS, LIB_VAR


# get variables for each run
def get_run_parameters(config, scprocess_data_dir, LIB_VAR, LIBS):
  # define run variable
  if config['multiplexing']['demux_type'] in ["none", "flex", "ocm"]:
    RUN_VAR       = "sample_id"
  else:
    RUN_VAR       = "pool_id"

  if RUN_VAR == LIB_VAR:
    RUNS = LIBS
  else:
    # get all runs from metadata
    metadata_f  = config["project"]["sample_metadata"]
    samples_df  = pl.read_csv( metadata_f )
    RUNS        = samples_df[ RUN_VAR ].drop_nulls().unique().to_list()
    # should we exclude runs?
    RUNS        = _do_exclusions(RUNS, config, RUN_VAR)

  # get any custom parameters
  custom_run_params = _get_custom_parameters(config, RUN_VAR)

  # for flex, build a probe_id lookup from sample metadata
  if config['project']['is_flex']:
    metadata_f  = config["project"]["sample_metadata"]
    samples_df  = pl.read_csv(metadata_f)
    probe_id_map = dict(zip(samples_df["sample_id"].to_list(), samples_df["probe_id"].to_list()))
  else:
    probe_id_map = {}

  # for OCM, build an ocm_id lookup from sample metadata
  ocm_id_map = {}
  if config['multiplexing']['demux_type'] == 'ocm':
    metadata_f  = config["project"]["sample_metadata"]
    samples_df  = pl.read_csv(metadata_f)
    ocm_id_map  = dict(zip(samples_df["sample_id"].to_list(), samples_df["ocm_id"].to_list()))

  # load sample file, populate everything from config
  RUN_PARAMS  = {
    run_name: _get_run_parameters_one_run(
      run_name, config, scprocess_data_dir, custom_run_params,
      probe_id=probe_id_map.get(run_name),
      ocm_id=ocm_id_map.get(run_name)
    )
    for run_name in sorted(RUNS)
  }

  return RUN_PARAMS, RUN_VAR


# get all fastq files
def _get_fastqs(config, RUNS, is_hto = False):
  # get place to look for fastq files
  if is_hto:
    tmp_ls        = config['multiplexing']
    tmp_ls['arv_instance'] = config['project'].get('arv_instance', [])
  else:
    tmp_ls        = config['project']

  # compute cache directory
  af_dir          = pathlib.Path(config['project']['proj_dir']) / "output" / f"{config['project']['short_tag']}_mapping"
  cache_prefix    = "fastq_metadata_hto" if is_hto else "fastq_metadata"

  # try loading all runs from cache
  cached_fastqs   = {}
  uncached_runs   = []
  for run in RUNS:
    cache_f       = af_dir / f"{cache_prefix}_{run}.yaml"
    if cache_f.exists():
      cached_fastqs[run] = _read_fastq_metadata_cache(cache_f)
    else:
      uncached_runs.append(run)

  if not uncached_runs:
    return cached_fastqs

  # scan for uncached runs
  if "fastq_dir" in tmp_ls:
    fastq_dir     = tmp_ls['fastq_dir']
    try:
      fastq_dict  = _list_fastq_files_dir(fastq_dir)
    except (FileNotFoundError, OSError) as e:
      if cached_fastqs:
        warnings.warn(
          f"FASTQ directory {fastq_dir} is not accessible ({e}), but cached metadata "
          f"exists for {len(cached_fastqs)}/{len(RUNS)} libraries. Missing caches for: "
          f"{uncached_runs}"
        )
      raise FileNotFoundError(
        f"FASTQ files are not accessible at {fastq_dir} and no cached metadata exists "
        f"for: {uncached_runs}. If the FASTQs were deleted after mapping, re-download "
        f"them and re-run: scprocess run <config> -r mapping"
      ) from e
  else:
    arv_uuids     = tmp_ls['arv_uuids']
    arv_instance  = tmp_ls['arv_instance']
    fastq_dict    = _list_fastq_files_arvados(arv_uuids, arv_instance)

  scanned_fastqs  = _match_fastqs_to_runs(fastq_dict, RUNS, is_hto)

  # write cache for all successfully resolved runs
  af_dir.mkdir(parents=True, exist_ok=True)
  for run, data in scanned_fastqs.items():
    cache_f       = af_dir / f"{cache_prefix}_{run}.yaml"
    _write_fastq_metadata_cache(cache_f, run, data)

  return scanned_fastqs


def _match_fastqs_to_runs(fastq_dict, RUNS, is_hto):
  wheres        = fastq_dict["wheres"]
  fastq_fs      = fastq_dict["fastqs"]
  fq_sizes_gb   = fastq_dict["fastq_sizes"]

  fastqs        = {}
  for run in RUNS:
    # get R1 and R2 files matching each run
    R1_regex      = rf".*{run}.*(_|\.)R1.*\.fastq\.gz"
    R1_fs         = [f for f in fastq_fs if re.match(R1_regex, f) ]
    R1_fs         = sorted(R1_fs)
    R2_regex      = rf".*{run}.*(_|\.)R2.*\.fastq\.gz"
    R2_fs         = [f for f in fastq_fs if re.match(R2_regex, f) ]
    R2_fs         = sorted(R2_fs)

    # find where
    wheres_R1     = [ wheres[i] for i,f in enumerate(fastq_fs) if re.match(R1_regex, f) ]
    wheres_R2     = [ wheres[i] for i,f in enumerate(fastq_fs) if re.match(R2_regex, f) ]
    this_where    = list(set(wheres_R1 + wheres_R2))
    if len(this_where) > 1:
      raise ValueError(f"FASTQ files for run {run} seem to be found in more than location")

    # get file sizes
    R1_fs_size_gb = [ fq_sizes_gb[i] for i,f in enumerate(fastq_fs) if re.match(R1_regex, f) ]
    R2_fs_size_gb = [ fq_sizes_gb[i] for i,f in enumerate(fastq_fs) if re.match(R2_regex, f) ]

    # check have full set of files
    check_R1      = [re.sub(r'(?<=(_|\.))R1', 'R0', f) for f in R1_fs]
    check_R2      = [re.sub(r'(?<=(_|\.))R2', 'R0', f) for f in R2_fs]
    if len(R1_fs) == 0:
      print(f"  WARNING: no {'hto ' if is_hto else ''}fastq files found for run {run}; excluded.")
    elif set(check_R1) != set(check_R2):
      print(f"  WARNING: {'hto ' if is_hto else ''}fastq files found for run {run} but R1 and R2 don't match; excluded.")
    else:
      fastqs[run] = {
        "where":          this_where[0],
        "R1_fs":          R1_fs,
        "R2_fs":          R2_fs,
        "R1_fs_size_gb":  round(sum(R1_fs_size_gb), 1),
        "R2_fs_size_gb":  round(sum(R2_fs_size_gb), 1)
      }

  return fastqs


def _write_fastq_metadata_cache(cache_f, run, data):
  cache_data = {
    "lib":            run,
    "where":          str(data["where"]),
    "R1_fs":          data["R1_fs"],
    "R2_fs":          data["R2_fs"],
    "R1_fs_size_gb":  data["R1_fs_size_gb"],
    "R2_fs_size_gb":  data["R2_fs_size_gb"]
  }
  with open(cache_f, 'w') as f:
    yaml.dump(cache_data, f, default_flow_style=False)


def _read_fastq_metadata_cache(cache_f):
  with open(cache_f, 'r') as f:
    cache_data = yaml.safe_load(f)
  return {
    "where":          pathlib.Path(cache_data["where"]),
    "R1_fs":          cache_data["R1_fs"],
    "R2_fs":          cache_data["R2_fs"],
    "R1_fs_size_gb":  cache_data["R1_fs_size_gb"],
    "R2_fs_size_gb":  cache_data.get("R2_fs_size_gb", 0.0)
  }


# get all fastq files in directory
def _list_fastq_files_dir(fastq_dir):
  # get all files
  all_fs          = os.listdir(fastq_dir)

  # filter to just fastqs
  fastq_fs        = [ f for f in all_fs if re.match(r".+\.fastq\.gz", f) ]
  wheres          = [ fastq_dir for f in fastq_fs]
  fastq_sizes_gb  = [ (fastq_dir / f).stat().st_size  / BYTES_PER_GB for f in fastq_fs ]
  fastq_fs_strs   = [ f"\"{f}\"" for f in fastq_fs ]

  return { "wheres": wheres, "fastqs": fastq_fs_strs, "fastq_sizes": fastq_sizes_gb }


# get all fastq files in all arvados uuids
def _list_fastq_files_arvados(arv_uuids, arv_instance):
  # get for each UUID
  wheres      = []
  fastqs      = []
  fastq_sizes = []

  # import relevant packages
  import arvados
  import collections
  import pathlib

  # set up arvados access
  arv_token   = os.environ["ARVADOS_API_TOKEN"]
  arv_client  = arvados.api('v1', host = f'api.{arv_instance}.roche.com',
    token = arv_token, insecure = True, num_retries = 2 )

  # check it worked
  try:
    user_info = arv_client.users().current().execute()
    print(f"  Arvados token is valid: logged in as: {user_info['full_name']} ({user_info['uuid']})")
  except Exception as e:
    print(f"  Arvados token is invalid or expired. Error: {e}")

  # get all fastq files in given arvados uuid
  def _list_fastq_files_arvados_one_uuid(arv_uuid, arv_instance):
    # define variables
    arv_files   = []
    wheres      = {}
    file_sizes  = {}

    # access this collection
    arv_colln   = arvados.collection.Collection(arv_uuid, arv_client)

    # get all files within this uuid
    stream_q    = collections.deque([pathlib.PurePosixPath('.')])
    while stream_q:
      stream_path = stream_q.popleft()
      tmp_colln   = arv_colln.find(str(stream_path))
      for item_name in tmp_colln:
        try:
          # open file
          my_file     = tmp_colln.open(item_name)

          # store results
          f             = os.path.join(str(stream_path), item_name)
          arv_files.append( f )
          if f in wheres:
            raise ValueError(f"file {item_name}")
          wheres[f]     = arv_uuid
          file_sizes[f] = tmp_colln[item_name].size()

        except IsADirectoryError:
          # item_name refers to a stream. Queue it to walk later.
          stream_q.append(stream_path / item_name)
          continue

    # filter to just fastqs
    fastq_re        = r".+\.fastq\.gz"
    fastq_fs        = [ f for f in arv_files if re.match(fastq_re, f) ]
    fastq_sizes_gb  = [ round(file_sizes[f] / BYTES_PER_GB, 1) for f in fastq_fs ]
    wheres          = [ wheres[f] for f in fastq_fs ]
    fastq_fs_strs   = [ f"\"{f}\"" for f in fastq_fs ]

    return { "wheres": wheres, "fastqs": fastq_fs_strs, "fastq_sizes": fastq_sizes_gb }

  # Iterate through each UUID in the list
  for arv_uuid in arv_uuids:
    # Get the dictionary result for one UUID
    result = _list_fastq_files_arvados_one_uuid(arv_uuid, arv_instance)

    # Extend the combined lists with the data from the current result
    # Note: We use .extend() for efficient list concatenation
    wheres.extend(result["wheres"])
    fastqs.extend(result["fastqs"])
    fastq_sizes.extend(result["fastq_sizes"])

  return {"wheres": wheres, "fastqs": fastqs, "fastq_sizes": fastq_sizes }


# load custom parameters if defined
def _get_custom_parameters(config, sel_var):
  # sort out custom sample parameters
  if not 'custom_sample_params' in config['project']: 
    return {}

  # open file
  with open(config['project']['custom_sample_params']) as f:
    # get all samples with custom params
    custom_params = yaml.load(f, Loader=yaml.FullLoader)

  # load files
  if not sel_var in custom_params:
    return {}

  return custom_params[ sel_var ]


# exclude some 
def _do_exclusions(LIST, config, var):
  EXC_LIST    = []
  if 'exclude' in config['project']:
    if var in config['project']['exclude']:
      # check if we need to exclude anything
      EXC_LIST    = config['project']['exclude'][var]
      for r in EXC_LIST:
        if r not in LIST:
          warnings.warn(f"{var} {r} specified in 'exclude' but not in sample_metadata file", UserWarning)
  # use this to get list to keep
  to_keep     = set(LIST) - set(EXC_LIST)
  LIST        = [r for r in LIST if r in to_keep]

  return LIST


# get parameters for one library
def _get_lib_parameters_one_lib(lib_name, config, RNA_FQS, HTO_FQS, scdata_dir, custom_lib_params):
  # get file with valid whitelists
  wl_df_f         = scdata_dir / 'cellranger_ref/cellranger_whitelists.csv'
  wl_df           = pl.read_csv(wl_df_f)
  
  # set up chemistry etc for flex
  if config['mapping_af']['is_flex']:
    # chemistry and whitelist are determined by the probe set
    flex_version    = config['mapping_af']['tenx_chemistry']
    tenx_chemistry  = flex_version
    af_chemistry    = '10x-flexv1-gex-3p' if flex_version == 'flexv1' else '10x-flexv2-gex-3p'
    expected_ori    = 'fw'  # flex is always forward orientation
    wl_row          = wl_df.filter(pl.col('chemistry') == flex_version)
    gex_whitelist_f = scdata_dir / 'cellranger_ref' / wl_row['gex_barcodes_f'].item()

  else:
    # allow per-library chemistry override via custom_sample_params
    tenx_chemistry = config['mapping_af']['tenx_chemistry']
    if lib_name in custom_lib_params:
      if 'project' in custom_lib_params[lib_name] and 'tenx_chemistry' in custom_lib_params[lib_name]['project']:
        tenx_chemistry = custom_lib_params[lib_name]['project']['tenx_chemistry']

    # map 10x kit name to simpleaf --chemistry flag (10xv2 uses 16bp barcodes, 10xv3 uses 18bp); default is "none" then auto-detect at mapping time
    af_chemistry = 'none'
    if tenx_chemistry in ['3v2', '5v1', '5v2']:
      af_chemistry = '10xv2'
    elif tenx_chemistry in ['3v3', '3v4', '5v3', 'multiome']:
      af_chemistry = '10xv3'

    # 3' chemistries read forward, 5' chemistries read reverse complement; default is "none" then auto-detect at mapping time
    expected_ori = 'none'
    if tenx_chemistry in ['5v1', '5v2', '5v3']:
      expected_ori = 'rc'
    elif tenx_chemistry in ['3v2', '3v3', '3v4', 'multiome']:
      expected_ori = 'fw'

    # look up barcode whitelist files; 'none' defers to auto-detection at mapping time
    if tenx_chemistry == 'none':
      gex_whitelist_f = 'none'
      hto_whitelist_f = 'none'
    else:
      wl_gex_f        = wl_df.filter(pl.col('chemistry') == tenx_chemistry)['gex_barcodes_f'].item()
      wl_hto_f        = wl_df.filter(pl.col('chemistry') == tenx_chemistry)['hto_barcodes_f'].item()
      gex_whitelist_f = scdata_dir / 'cellranger_ref' / wl_gex_f
      hto_whitelist_f = scdata_dir / 'cellranger_ref' / wl_hto_f

  # make dictionary for RNA mapping
  mapping_dc  = {
    "where":              RNA_FQS[lib_name]["where"],
    "R1_fs":              RNA_FQS[lib_name]["R1_fs"],
    "R2_fs":              RNA_FQS[lib_name]["R2_fs"],
    "R1_fs_size_gb":      RNA_FQS[lib_name]["R1_fs_size_gb"],
    "tenx_chemistry":     tenx_chemistry,
    "af_chemistry":       af_chemistry,
    "expected_ori":       expected_ori,
    "gex_whitelist_f":    gex_whitelist_f,
    "geometry":           config['mapping_af'].get('geometry', ''),
  }

  # make dictionary for HTO mapping
  if config['multiplexing']['demux_type'] == "hto":
    multiplexing_dc  = {
      "where":              HTO_FQS[lib_name]["where"],
      "R1_fs":              HTO_FQS[lib_name]["R1_fs"],
      "R2_fs":              HTO_FQS[lib_name]["R2_fs"],
      "R1_fs_size_gb":      HTO_FQS[lib_name]["R1_fs_size_gb"],
      "af_chemistry":       af_chemistry, 
      "hto_whitelist_f":    hto_whitelist_f
    }
  else:
    multiplexing_dc   = {}

  # make dict of dicts
  out_dc      = {
    "mapping_af":   mapping_dc,
    "multiplexing": multiplexing_dc,
  }

  return out_dc


# get parameters for one run
def _get_run_parameters_one_run(run_name, config, scdata_dir, custom_run_params, probe_id=None, ocm_id=None, lib_name=None, LIB_PARAMS=None):

  knee1 = ""
  shin1 = ""
  knee2 = ""
  shin2 = ""

  if run_name in custom_run_params:
    if 'mapping' in custom_run_params[run_name]:
      if 'knee1' in custom_run_params[run_name]['mapping']:
        knee1 = custom_run_params[run_name]['mapping']['knee1']
      if 'shin1' in custom_run_params[run_name]['mapping']:
        shin1 = custom_run_params[run_name]['mapping']['shin1']
      if 'knee2' in custom_run_params[run_name]['mapping']:
        knee2 = custom_run_params[run_name]['mapping']['knee2']
      if 'shin2' in custom_run_params[run_name]['mapping']:
        shin2 = custom_run_params[run_name]['mapping']['shin2']

  # get R1 size from the corresponding library
  if LIB_PARAMS is not None and lib_name is not None and lib_name in LIB_PARAMS:
    r1_size_gb = LIB_PARAMS[lib_name]["mapping_af"]["R1_fs_size_gb"]
  else:
    r1_size_gb = 0

  # make dictionary for mapping
  mapping_dc  = {
    "knee1":              knee1,
    "shin1":              shin1,
    "knee2":              knee2,
    "shin2":              shin2,
    "R1_fs_size_gb":      r1_size_gb
  }
  if probe_id is not None:
    mapping_dc["probe_id"] = probe_id
  if ocm_id is not None:
    mapping_dc["ocm_id"] = ocm_id
  # make dictionary for ambient
  ambient_dc = {
    "cb_expected_cells":          "",
    "cb_total_droplets_included": "",
    "cb_low_count_threshold":     "",
    "cb_learning_rate":           "",
    "cb_empty_training_fraction": "",
    "cb_posterior_batch_size":    ""
  }

  # check for global cellbender params defined in config
  if config['ambient']['ambient_method'] == 'cellbender':
    for v in ambient_dc:
      if v in config['ambient']:
        ambient_dc[v] = config['ambient'][v]

  # check custom_sample_params_f for sample specific params
  if (run_name in custom_run_params) and ('ambient' in custom_run_params[run_name]):
    amb_run_params = custom_run_params[run_name]['ambient']
    for v in ambient_dc:
      if v in amb_run_params:
        ambient_dc[v] = amb_run_params[v]

  # make dict of dicts
  out_dc      = {
    "mapping": mapping_dc,
    "ambient": ambient_dc
  }

  return out_dc


# get variables for each run
def get_batch_parameters(config, RUNS, scprocess_data_dir):
  # define batch variable
  BATCH_VAR   = config['integration']['int_batch_var']

  # get samples
  metadata_f  = config["project"]["sample_metadata"]
  samples_df  = pl.read_csv( metadata_f )
  SAMPLES     = samples_df[ "sample_id" ].drop_nulls().to_list()

  # get parameters if batch_var is sample_id
  if BATCH_VAR == "sample_id":
    # should we exclude any runs?
    SAMPLES     = _do_exclusions(SAMPLES, config, "sample_id")

    # use to define batches
    BATCHES     = SAMPLES

  # get parameters if batch_var is pool_id
  elif BATCH_VAR == "pool_id":
    BATCHES     = RUNS

  # get any custom parameters
  custom_batch_params = _get_custom_parameters(config, BATCH_VAR)

  # load sample file, populate everything from config
  BATCH_PARAMS = {
    batch_name: _get_batch_parameters_one_batch(batch_name, config, custom_batch_params)
    for batch_name in sorted(BATCHES)
  }

  return BATCH_PARAMS, BATCH_VAR, SAMPLES


# get parameters for one sample
def _get_batch_parameters_one_batch(batch_name, config, custom_batch_params):
  # make dictionary for mapping
  qc_dc   = {
    # set defaults
    "qc_min_counts":  config['qc']['qc_min_counts'],
    "qc_max_counts":  config['qc']['qc_max_counts'],
    "qc_min_feats":   config['qc']['qc_min_feats'],
    "qc_max_feats":   config['qc']['qc_max_feats'],
    "qc_min_mito":    config['qc']['qc_min_mito'],
    "qc_max_mito":    config['qc']['qc_max_mito'],
    "qc_min_splice":  config['qc']['qc_min_splice'],
    "qc_max_splice":  config['qc']['qc_max_splice'],
    "qc_min_cells":   config['qc']['qc_min_cells']
  }
  # add sample-specific QC parameters
  if batch_name in custom_batch_params:
    if 'qc' in custom_batch_params[batch_name]:
      for v in qc_dc:
        if v in custom_batch_params[batch_name]['qc']:
          qc_dc[v]    = custom_batch_params[batch_name]['qc'][v]

  _validate_qc_bounds(qc_dc, f"qc for {batch_name}")

  # make dict of dicts
  out_dc  = {
    "qc":   qc_dc
  }

  return out_dc


def _get_pool_to_sample_map(filter_ids, sample_metadata):
  df = pl.read_csv(sample_metadata)
  grouped = (
    df.filter(pl.col("pool_id").is_in(filter_ids))
    .sort(["pool_id", "sample_id"])
    .group_by("pool_id", maintain_order=True)
    .agg(pl.col("sample_id"))
  )
  # returns a dict: {pool_id: [sample_ids]}
  return dict(zip(grouped["pool_id"], grouped["sample_id"]))


def get_runs_to_batches(config, RUNS, BATCHES, BATCH_VAR, LIBS):

  demux_type      = config['multiplexing']['demux_type']
  sample_metadata = config['project']['sample_metadata']

  if demux_type in ["none", "flex", "ocm"]:
    if RUNS != BATCHES:
      raise ValueError(f"RUNS and BATCHES must match for demux_type '{demux_type}'")

    RUNS_TO_SAMPLES = {s: [s] for s in BATCHES}
    RUNS_TO_BATCHES = {s: [s] for s in BATCHES}

    if demux_type in ["flex", "ocm"]:
      # build sample_id -> pool_id map from metadata
      df           = pl.read_csv(sample_metadata)
      all_run_lib  = dict(zip(df["sample_id"].to_list(), df["pool_id"].to_list()))
      RUNS_TO_LIBS = {r: all_run_lib[r] for r in RUNS}
    else:
      RUNS_TO_LIBS = {r: r for r in RUNS}

    return RUNS_TO_BATCHES, RUNS_TO_SAMPLES, RUNS_TO_LIBS

  else:
    RUNS_TO_SAMPLES = _get_pool_to_sample_map(RUNS, sample_metadata)

    if BATCH_VAR == "pool_id":
      if RUNS != BATCHES:
        raise ValueError("RUNS and BATCHES must match if BATCH_VAR is 'pool_id'")
      RUNS_TO_BATCHES = {s: [s] for s in sorted(BATCHES)}

    elif BATCH_VAR == "sample_id":
      RUNS_TO_BATCHES = {
        pool: sorted([s for s in samples if s in BATCHES])
        for pool, samples in RUNS_TO_SAMPLES.items()
      }
      RUNS_TO_SAMPLES = RUNS_TO_BATCHES

    if LIBS != RUNS:
      raise ValueError(f"LIBS and RUNS must match for demux_type '{demux_type}'")
    RUNS_TO_LIBS = {r: r for r in RUNS}

    return RUNS_TO_BATCHES, RUNS_TO_SAMPLES, RUNS_TO_LIBS


# get parameters for labelling celltypes
def get_labeller_parameters(config, schema_f, scdata_dir):
  # if none, done
  if not 'label_celltypes' in config:
    return []

  # get defaults for label parameters
  schema          = _load_schema_file(schema_f)
  label_schema    = schema["properties"]["label_celltypes"]["items"]
  label_defaults  = _get_default_config_from_schema(label_schema)

  # get things we need for checks
  typist_ls_f     = scdata_dir / 'celltypist/celltypist_models.csv'
  xgboost_ls_f    = scdata_dir / 'xgboost/xgboost_models.csv'
  mdls_typist     = pl.read_csv(typist_ls_f)['model'].to_list()
  xgboost_df      = pl.read_csv(xgboost_ls_f)
  mdls_xgboost    = xgboost_df['model'].to_list()

  # check that selected models are valid
  def _check_one_label_celltypes_parameters(entry):
    # check that parameters for celltypist are ok
    if entry['labeller'] == 'celltypist':
      if not entry['model'] in mdls_typist:
        raise KeyError(
          f"The value {entry['model']} specified in label_celltypes is not a valid celltypist model.\n"
          f"The following are valid models:\n{', '.join(mdls_typist)}"
          )

    # check that parameters for scprocess are ok
    elif entry['labeller'] == 'scprocess':
      if not entry['model']  in mdls_xgboost:
        raise KeyError(
          f"the value {entry['model']} specified in label_celltypes is not a valid scprocess model. "
          f"These models are currently available: {', '.join(mdls_xgboost)}"
        )

      row = xgboost_df.filter(pl.col('model') == entry['model']).row(0, named=True)
      entry['model_f']     = row['model_f']
      entry['cls_f']       = row['cls_f']
      entry['genes_f']     = row['genes_f']
      entry['label_map_f'] = row.get('label_map_f', '')

      for key in ['model_f', 'cls_f', 'genes_f']:
        if not pathlib.Path(entry[key]).is_file():
          raise FileNotFoundError(f"file {entry[key]} doesn't exist; consider (re)running scprocess setup")

    # resolve label_map_f
    if 'label_map' in entry:
      lm_path = _check_path_exists_in_project(entry['label_map'], config, what = "file")
      entry['label_map_f'] = str(lm_path)
    elif entry['labeller'] == 'celltypist':
      default_lm = scdata_dir / 'celltypist' / f"{entry['model']}_label_map.csv"
      entry['label_map_f'] = str(default_lm) if default_lm.is_file() else ''

    if 'label_map_f' not in entry:
      entry['label_map_f'] = ''

    if entry['label_map_f'] != '':
      lm_df = pl.read_csv(entry['label_map_f'], n_rows=5)
      if 'fine_label' not in lm_df.columns or 'coarse_label' not in lm_df.columns:
        raise ValueError(f"label_map file must have 'fine_label' and 'coarse_label' columns: {entry['label_map_f']}")

    # add defaults
    for v in label_defaults:
      if not v in entry:
        entry[v] = label_defaults[v]

    # validate save_cluster_names_file
    if entry.get('save_cluster_names_file', False) and 'marker_genes' not in config:
      warnings.warn(
        f"save_cluster_names_file is true for {entry['labeller']}/{entry['model']} "
        f"but no marker_genes block is configured; skipping cluster names file.")
      entry['save_cluster_names_file'] = False

    return entry

  # apply this to each specified model
  LABELLER_PARAMS = [ _check_one_label_celltypes_parameters(entry) for entry in config['label_celltypes'] ]

  return LABELLER_PARAMS


# --- join parameter helpers ---

def get_join_project_parameters(config):
  """Load source project configs and derive per-project paths and batch keys."""
  project_ids = list(config['projects'].keys())
  project_cfgs = {}

  for pid in project_ids:
    cfg_f = config['projects'][pid]['config']
    with open(cfg_f) as f:
      project_cfgs[pid] = yaml.safe_load(f)

  def _dir(pid):
    return pathlib.Path(project_cfgs[pid]['project']['proj_dir'])

  def _short_tag(pid):
    return project_cfgs[pid]['project']['short_tag']

  def _full_tag(pid):
    return project_cfgs[pid]['project']['full_tag']

  def _date(pid):
    return project_cfgs[pid]['project']['date_stamp']

  def _var_stats_f(pid):
    zoom_name = config['projects'][pid].get('zoom_name')
    if zoom_name:
      zoom_dir = _dir(pid) / f"output/{_short_tag(pid)}_zoom"
      return zoom_dir / zoom_name / f"standardized_variance_stats_{_full_tag(pid)}_{zoom_name}_{_date(pid)}.csv.gz"
    return _dir(pid) / f"output/{_short_tag(pid)}_hvg" / f"standardized_variance_stats_{_full_tag(pid)}_{_date(pid)}.csv.gz"

  def _h5ads_yaml_f(pid):
    return _dir(pid) / f"output/{_short_tag(pid)}_integration" / f"h5ads_clean_paths_{_full_tag(pid)}_{_date(pid)}.yaml"

  def _integrated_dt_f(pid):
    zoom_name = config['projects'][pid].get('zoom_name')
    if zoom_name:
      zoom_dir = _dir(pid) / f"output/{_short_tag(pid)}_zoom"
      return zoom_dir / zoom_name / f"integrated_dt_{_full_tag(pid)}_{zoom_name}_{_date(pid)}.csv.gz"
    return _dir(pid) / f"output/{_short_tag(pid)}_integration" / f"integrated_dt_{_full_tag(pid)}_{_date(pid)}.csv.gz"

  def _sample_meta_f(pid):
    meta_f = pathlib.Path(project_cfgs[pid]['project']['sample_metadata'])
    if not meta_f.is_absolute():
      meta_f = _dir(pid) / meta_f
    return meta_f

  def _qc_all_f(pid):
    return (
      _dir(pid) / f"output/{_short_tag(pid)}_qc" /
      f"qc_all_samples_{_full_tag(pid)}_{_date(pid)}.csv.gz"
    )

  # build batch keys from per-project h5ads YAMLs
  h5ads_yaml_fs = [str(_h5ads_yaml_f(pid)) for pid in project_ids]
  batch_keys = []
  for pid, h5yaml in zip(project_ids, h5ads_yaml_fs):
    with open(h5yaml) as fh:
      h5paths = yaml.safe_load(fh)
    for bk in h5paths:
      batch_keys.append(f"{pid}_{bk}")

  return {
    'project_ids':    project_ids,
    'project_cfgs':   project_cfgs,
    'var_stats_fs':   [str(_var_stats_f(pid)) for pid in project_ids],
    'h5ads_yaml_fs':  h5ads_yaml_fs,
    'integrated_fs':  [str(_integrated_dt_f(pid)) for pid in project_ids],
    'sample_meta_fs': [str(_sample_meta_f(pid)) for pid in project_ids],
    'qc_all_fs':      [str(_qc_all_f(pid)) for pid in project_ids],
    'batch_keys':     batch_keys,
  }



def get_join_source_labels_f(project_cfgs, pid, labeller, model):
  """Return path to a source project's aggregated labels file, or None if unavailable."""
  pcfg = project_cfgs[pid]
  if 'label_celltypes' not in pcfg:
    return None
  matches = [e for e in pcfg['label_celltypes']
    if e.get('labeller') == labeller and e.get('model') == model]
  if not matches:
    return None
  proj_dir  = pathlib.Path(pcfg['project']['proj_dir'])
  short_tag = pcfg['project']['short_tag']
  full_tag  = pcfg['project']['full_tag']
  date      = pcfg['project']['date_stamp']
  lbl_dir   = proj_dir / f"output/{short_tag}_label_celltypes"
  f = lbl_dir / f"labels_{labeller}_model_{model}_{full_tag}_{date}.csv.gz"
  return str(f) if f.is_file() else None


def get_join_batch_sources(project_cfgs, project_ids, batch_keys, labeller, model):
  """Split batch keys into those with reusable labels and those needing fresh runs."""
  reuse_label_fs = []
  fresh_batches = []
  for pid in project_ids:
    proj_label_f = get_join_source_labels_f(project_cfgs, pid, labeller, model)
    pid_batches = [bk for bk in batch_keys if bk.startswith(f"{pid}_")]
    if proj_label_f:
      reuse_label_fs.append(proj_label_f)
    else:
      fresh_batches.extend(pid_batches)
  return fresh_batches, reuse_label_fs


def get_shiny_targets(config, scprocess_dir, scdata_dir, dryrun, zoom_name=None):
  """Determine shiny build targets and create output directories.

  Returns (targets, proj_dir) where targets is a list of Snakemake target
  rule names or file paths.
  """
  scprocess_dir = pathlib.Path(scprocess_dir)
  _is_join = 'join' in config
  if _is_join:
    join_schema_f = scprocess_dir / "resources/schemas/join.schema.json"
    config = check_join_config(config, join_schema_f, scdata_dir)
  else:
    schema_f = scprocess_dir / "resources/schemas/config.schema.json"
    config = check_config(config, schema_f, scdata_dir, scprocess_dir)

  proj_cfg = config['join'] if _is_join else config['project']
  proj_dir = pathlib.Path(proj_cfg['proj_dir'])
  date_stamp = proj_cfg['date_stamp']
  docs_dir = proj_dir / "public"

  if _is_join and zoom_name is not None:
    raise ValueError("--zoom is not supported for join configs.")

  if zoom_name is None:
    app_tag = get_shiny_app_tag(config)
    if not dryrun:
      os.makedirs(docs_dir / f"shiny_{app_tag}", exist_ok=True)
    return ["build_shiny_app"], proj_dir

  if zoom_name == "all":
    zoom_schema_f = pathlib.Path(scprocess_dir) / "resources/schemas/zoom.schema.json"
    ZOOM_PARAMS = get_zoom_parameters(config, zoom_schema_f, scdata_dir)
    zoom_names = list(ZOOM_PARAMS.keys())
    if not zoom_names:
      raise ValueError("No zoom specs found in config. Add a 'zoom:' section to your config file.")
    if not dryrun:
      for zn in zoom_names:
        os.makedirs(docs_dir / f"shiny_zoom_{zn}", exist_ok=True)
    return ["build_all_zoom_shiny_apps"], proj_dir

  zoom_schema_f = pathlib.Path(scprocess_dir) / "resources/schemas/zoom.schema.json"
  ZOOM_PARAMS = get_zoom_parameters(config, zoom_schema_f, scdata_dir)
  zoom_names = list(ZOOM_PARAMS.keys())
  if zoom_name not in zoom_names:
    raise ValueError(
      f"Zoom '{zoom_name}' not found in config. "
      f"Available zooms: {', '.join(zoom_names) if zoom_names else '(none)'}")
  if not dryrun:
    os.makedirs(docs_dir / f"shiny_zoom_{zoom_name}", exist_ok=True)
  return [str(docs_dir / f"shiny_zoom_{zoom_name}" / f".shiny_built_{date_stamp}")], proj_dir


def get_train_xgboost_targets(config, scprocess_dir, scdata_dir, zoom_name=None):
  """Determine train_xgboost targets, optionally for a specific zoom."""
  scprocess_dir = pathlib.Path(scprocess_dir)

  schema_f = scprocess_dir / "resources/schemas/config.schema.json"
  config = check_config(config, schema_f, scdata_dir, scprocess_dir)
  proj_dir = pathlib.Path(config['project']['proj_dir'])

  if zoom_name is None:
    return ["train_xgboost"], proj_dir

  zoom_schema_f = scprocess_dir / "resources/schemas/zoom.schema.json"
  ZOOM_PARAMS = get_zoom_parameters(config, zoom_schema_f, scdata_dir)
  zoom_names = list(ZOOM_PARAMS.keys())

  if zoom_name == "all":
    if not zoom_names:
      raise ValueError("No zoom specs found in config.")
    zooms_with_xgb = [zn for zn in zoom_names if 'train_xgboost' in ZOOM_PARAMS[zn]]
    if not zooms_with_xgb:
      raise ValueError("No zooms have a train_xgboost section.")
    return ["zoom_train_xgboost_all"], proj_dir

  if zoom_name not in zoom_names:
    raise ValueError(
      f"Zoom '{zoom_name}' not found in config. "
      f"Available zooms: {', '.join(zoom_names) if zoom_names else '(none)'}")
  if 'train_xgboost' not in ZOOM_PARAMS[zoom_name]:
    raise ValueError(f"Zoom '{zoom_name}' does not have a train_xgboost section.")

  ref_tag = ZOOM_PARAMS[zoom_name]['train_xgboost']['ref_tag']
  short_tag = config['project']['short_tag']
  mkr_sel_res = ZOOM_PARAMS[zoom_name]['marker_genes']['mkr_sel_res']
  model_f = str(proj_dir / f"output/{short_tag}_zoom/{zoom_name}/train_xgboost/{ref_tag}_xgboost_model.json")
  html_f = str(proj_dir / f"public/{short_tag}_zoom_{zoom_name}_{mkr_sel_res}.html")
  return [model_f, html_f], proj_dir


def prep_resource_params(config, schema_f, scprocess_dir, LIB_PARAMS=None, BATCHES=None,
                         RUNS_TO_LIBS=None, chem_stats_paths=None):
  # add default resource values
  schema      = _load_schema_file(schema_f)
  defaults    = _get_default_config_from_schema(schema)
  defaults    = defaults.get('resources', {})

  # get user resource values
  user_vals   = config.get('resources', {}).copy()

  # if same as default, remove from user vals
  for n in list(user_vals):
    if n in defaults:
      if defaults[n] == user_vals[n]:
        del user_vals[n]

  # load resource model CSV
  lm_f        = scprocess_dir / "resources/snakemake" / LM_PARAMS_FILENAME
  lm_df       = pl.read_csv(lm_f)

  # get sizes, n-batches (optional — not available for shiny/join standalone workflows)
  R1_sizes    = { lib: vals["mapping_af"]["R1_fs_size_gb"] for lib, vals in LIB_PARAMS.items() } if LIB_PARAMS else {}
  n_batches   = len(BATCHES) if BATCHES else 0

  # build per-run chemistry lookup from LIB_PARAMS
  run_chemistries = {}
  if RUNS_TO_LIBS and LIB_PARAMS:
    for run, lib in RUNS_TO_LIBS.items():
      run_chemistries[run] = LIB_PARAMS[lib]["mapping_af"].get("tenx_chemistry", "none")

  # make full dict of useful things
  RESOURCE_PARAMS = {}
  RESOURCE_PARAMS["defaults"]               = defaults
  RESOURCE_PARAMS["user_vals"]              = user_vals
  RESOURCE_PARAMS["lm_df"]                  = lm_df
  RESOURCE_PARAMS["R1_sizes"]               = R1_sizes
  RESOURCE_PARAMS["n_batches"]              = n_batches
  RESOURCE_PARAMS["n_run_mapping"]          = config.get('resources', {}).get('n_run_mapping', 8)
  RESOURCE_PARAMS["run_chemistries"]        = run_chemistries
  RESOURCE_PARAMS["RUNS_TO_LIBS"]           = RUNS_TO_LIBS or {}
  RESOURCE_PARAMS["chem_stats_paths"]       = chem_stats_paths or {}
  RESOURCE_PARAMS["_resolved_chemistries"]  = {}

  return RESOURCE_PARAMS


def get_resources(RESOURCE_PARAMS, rules, input, rule, param, attempt, run = None, csv_rule = None):
  attempt_exp = 1.5
  if not hasattr(rules, rule):
    raise ValueError(f'rule {rule} is not defined.')
  param_name  = f'mins_{rule}' if param == 'time' else f'gb_{rule}'
  lookup_rule = csv_rule if csv_rule is not None else rule

  defaults    = RESOURCE_PARAMS['defaults']
  user_vals   = RESOURCE_PARAMS['user_vals']

  # 1. user override takes priority (keyed by actual rule name)
  if param_name in user_vals:
    if param == 'memory':
      mem_mb    = user_vals[param_name] * MB_PER_GB
    else:
      time_min  = user_vals[param_name]

  # 2. resource model CSV (keyed by lookup_rule to allow model reuse)
  elif _rule_in_csv(RESOURCE_PARAMS, lookup_rule, param):
    filt_df   = _get_csv_row(RESOURCE_PARAMS, lookup_rule, param, run)
    mem_mb, time_min = _predict_from_csv(filt_df, RESOURCE_PARAMS, input, lookup_rule, param, run)

  # 3. schema default fallback (for rules not in CSV)
  elif param_name in defaults:
    if param == 'memory':
      mem_mb    = defaults[param_name] * MB_PER_GB
    else:
      time_min  = defaults[param_name]

  else:
    raise KeyError(f'No resource model or schema default for {param_name}.')

  if param == 'memory':
    mem_mb     *= attempt_exp**(attempt - 1)
    param_val   = mem_mb
  else:
    param_val   = time_min
  param_val   = pl.Series("dummy", [param_val]).ceil().cast(pl.Int32).item()

  return param_val


def _rule_in_csv(RESOURCE_PARAMS, rule, param):
  lm_df = RESOURCE_PARAMS['lm_df']
  return lm_df.filter((pl.col("param") == param) & (pl.col("rule") == rule)).shape[0] > 0


def _get_csv_row(RESOURCE_PARAMS, rule, param, run):
  lm_df     = RESOURCE_PARAMS['lm_df']
  rule_rows = lm_df.filter((pl.col("param") == param) & (pl.col("rule") == rule))
  has_chem  = (rule_rows["chemistry"] != "").any()

  if has_chem and run is not None:
    chemistry     = _resolve_run_chemistry(run, RESOURCE_PARAMS)
    chem_group    = _get_chemistry_group(chemistry)
    chem_row      = rule_rows.filter(pl.col("chemistry") == chem_group)
    if chem_row.shape[0] == 1:
      return chem_row
    generic_row   = rule_rows.filter(pl.col("chemistry") == "")
    if generic_row.shape[0] == 1:
      return generic_row
    raise ValueError(f"No matching chemistry row for rule '{rule}', param '{param}', chemistry '{chem_group}'.")

  generic_row = rule_rows.filter(pl.col("chemistry") == "")
  if generic_row.shape[0] == 1:
    return generic_row
  if rule_rows.shape[0] == 1:
    return rule_rows
  raise ValueError(f"Ambiguous CSV rows for rule '{rule}', param '{param}'.")


def _predict_from_csv(filt_df, RESOURCE_PARAMS, input, rule, param, run):
  has_model = not filt_df["rq_slope"].is_null().all()
  floor_val = filt_df["floor"].item()
  mem_mb    = None
  time_min  = None

  if has_model:
    param_est = _estimate_resource_parameter(filt_df, RESOURCE_PARAMS, input, rule, run)
    param_est = max(floor_val, param_est)
    if param == "memory":
      mem_mb    = param_est
    else:
      time_s    = param_est
      time_min  = time_s / 60
  else:
    if param == "memory":
      mem_mb    = floor_val
    else:
      time_min  = floor_val / 60

  return mem_mb, time_min


def _resolve_input_attr(input, attr):
  aliases = {
    'af_h5_f': ['clean_h5_f', 'filt_counts_f'],
  }
  if hasattr(input, attr):
    return getattr(input, attr)
  for alias in aliases.get(attr, []):
    if hasattr(input, alias):
      return getattr(input, alias)
  return None


def _estimate_resource_parameter(filt_lm_df, RESOURCE_PARAMS, input, rule, run):
  x_rq      = filt_lm_df['model_var'].item()
  intercept = filt_lm_df['rq_intercept'].item()
  slope     = filt_lm_df['rq_slope'].item()

  if x_rq.startswith('input.'):
    input_attr  = x_rq.replace("input.", "")
    resolved    = _resolve_input_attr(input, input_attr)
    if resolved is not None:
      x_val     = os.path.getsize(resolved) / BYTES_PER_GB
    else:
      raise ValueError(f"'{input_attr}' is not a valid input attribute for rule '{rule}'.")

  elif x_rq == 'raw_data_size':
    if run is None:
      raise ValueError(f'run argument should be defined')
    x_val       = RESOURCE_PARAMS['R1_sizes'][run]

  elif x_rq == 'n_smpls_pre_qc':
    x_val       = RESOURCE_PARAMS['n_batches']

  elif x_rq == 'n_run_mapping':
    x_val       = RESOURCE_PARAMS['n_run_mapping']

  else:
    raise ValueError(f"Unknown variable '{x_rq}' for scaling resources.")

  return intercept + (slope * x_val)


def _resolve_run_chemistry(run, RESOURCE_PARAMS):
  chemistry = RESOURCE_PARAMS['run_chemistries'].get(run)
  if chemistry is not None and chemistry != "none":
    return chemistry

  lib = RESOURCE_PARAMS['RUNS_TO_LIBS'].get(run, run)
  cache = RESOURCE_PARAMS['_resolved_chemistries']
  if lib not in cache:
    chem_stats_f = RESOURCE_PARAMS['chem_stats_paths'].get(lib)
    if chem_stats_f:
      try:
        with open(chem_stats_f) as f:
          stats = yaml.safe_load(f)
        cache[lib] = stats['selected_tenx_chemistry']
      except (FileNotFoundError, KeyError, TypeError):
        cache[lib] = "3v3"
    else:
      cache[lib] = "3v3"
  return cache[lib]


def _get_chemistry_group(tenx_chemistry):
  if tenx_chemistry == "3v4":
    return "3v4"
  return "3v3"


### helpers

# some useful global variables
BYTES_PER_GB  = 1024**3
MB_PER_GB     = 1024
LM_PARAMS_FILENAME = "resources_lm_params_2026-07-07.csv"

# nice boolean values from yaml inputs
def _safe_boolean(val):
  if type(val) is bool:
    res = val
  elif val in ["True", "true"]:
    res = True
  elif val in ["False", "false"]:
    res = False
  else:
    raise ValueError(f'{val} is not a boolean')

  return res


# helper function to merge multiple zipped csv/txt files
def merge_tmp_files(in_files, out_file):
  df_ls     = [ pl.read_csv(f) for f in in_files if gzip.open(f, 'rb').read(1) ]
  df_merged = pl.concat(df_ls)
  with gzip.open(out_file, 'wb') as f:
    df_merged.write_csv(f)


def check_ranger_url(ranger_url):
  
  if not ranger_url.startswith("https://"):
    # check that link if valid
    raise ValueError(f"Invalid URL: Link must start with 'https://'")
    
  # extract version
  ranger_version = re.search(r"cellranger-(\d+\.\d+\.\d+)", ranger_url).group(1)
    
  # check that version if > 9.0.0
  version_parts = [int(p) for p in ranger_version.split(".")]
  min_version = [9, 0, 0]
  if version_parts <= min_version:
    raise ValueError(f"Provided download link for CellRanger version {ranger_version}. scprocess requires version > 9.0.0.")

  return ranger_version
  

# HVGs function: make df with list of chunked counts files
def make_hvgs_input_df(runs, ambient_outs_yamls, RUN_VAR, BATCH_VAR, BATCHES_TO_RUNS, 
  DEMUX_TYPE, FULL_TAG, DATE_STAMP, hvg_dir):
  # loop through ambient yaml files to populate listf
  df_list = []
  for r, yaml_file in zip(runs, ambient_outs_yamls):
    # get filtered ambient outputs
    with open(yaml_file) as f:
      amb_outs = yaml.load(f, Loader=yaml.FullLoader)
    amb_filt_f = amb_outs['filt_counts_f']

    # if no multiplexing or flex (samples already extracted per-run), simple
    if DEMUX_TYPE in ("none", "flex"):
      tmp_df = pl.DataFrame({
        BATCH_VAR:    r,
        'amb_filt_f': amb_filt_f
      })
      df_list.append(tmp_df)
    else:
      # get sample ids for pool
      sel_batches = BATCHES_TO_RUNS.get(r, [])
      tmp_df      = pl.DataFrame({
          RUN_VAR:      r,
          'amb_filt_f': amb_filt_f,
          BATCH_VAR:    sel_batches
        })
      df_list.append(tmp_df)

  # merge dfs for all runs
  chunk_pat   = f"{hvg_dir}/chunked_counts_{{}}_{FULL_TAG}_{DATE_STAMP}.h5"
  hvg_df_full = pl.concat(df_list).with_columns(
    pl.format(chunk_pat, BATCH_VAR).alias('chunked_f')
  )

  return hvg_df_full



def _check_train_xgboost_parameters(config):
  if 'train_xgboost' not in config:
    return config

  xgb = config['train_xgboost']

  # resolve annots_f relative to proj_dir
  annots_f = pathlib.Path(xgb['annots_f'])
  if not annots_f.is_absolute():
    annots_f = pathlib.Path(config['project']['proj_dir']) / annots_f
  xgb['annots_f'] = str(annots_f)

  # check annotations file
  if not annots_f.is_file():
    raise FileNotFoundError(f"Annotations file not found: {annots_f}")
  annots_df = pl.read_csv(str(annots_f), n_rows=5)
  if "cell_id" not in annots_df.columns:
    raise ValueError("annots_f must have 'cell_id' column")
  if "annotation" not in annots_df.columns:
    raise ValueError("annots_f must have 'annotation' column")

  # check label map if provided
  if xgb.get("label_map_f") is not None:
    label_map_f = pathlib.Path(xgb["label_map_f"])
    if not label_map_f.is_absolute():
      label_map_f = pathlib.Path(config['project']['proj_dir']) / label_map_f
    xgb['label_map_f'] = str(label_map_f)
    if not label_map_f.is_file():
      raise FileNotFoundError(f"Label map file not found: {label_map_f}")
    map_df = pl.read_csv(str(label_map_f), n_rows=5)
    if "annotation" not in map_df.columns:
      raise ValueError("label_map_f must have 'annotation' column")
    if "coarse_label" not in map_df.columns:
      raise ValueError("label_map_f must have 'coarse_label' column")

  return config


def get_train_xgboost_parameters(config, schema_f):
  """Return validated train_xgboost params dict, or None if section absent."""
  if 'train_xgboost' not in config:
    return None

  schema     = _load_schema_file(schema_f)
  xgb_schema = schema["properties"]["train_xgboost"]
  xgb_defaults = _get_default_config_from_schema(xgb_schema)

  xgb = config['train_xgboost']
  for v in xgb_defaults:
    if v not in xgb:
      xgb[v] = xgb_defaults[v]

  if 'ref_tag' not in xgb:
    xgb['ref_tag'] = f"xgboost_{config['project']['full_tag']}"
  else:
    xgb['ref_tag'] = f"xgboost_{xgb['ref_tag']}"

  return xgb


# ---------------------------------------------------------------------------
# Join config validation
# ---------------------------------------------------------------------------

def _apply_join_defaults(cfg, schema_props):
  """Recursively apply JSON Schema defaults to cfg in-place."""
  for key, prop in schema_props.items():
    if key not in cfg and 'default' in prop:
      cfg[key] = prop['default']
    elif key in cfg and prop.get('type') == 'object' and 'properties' in prop:
      _apply_join_defaults(cfg[key], prop['properties'])


def check_join_config(config, join_schema_f, scdata_dir):
  """Validate and augment a join.yaml config dict.

  Mirrors check_config() for project configs:
    1. Validates against join.schema.json
    2. Validates each referenced project config against config.schema.json
    3. Applies schema defaults for hvg, integration, marker_genes, and shiny sections
    4. Calls _check_shiny_parameters for shiny-specific cross-checks

  Returns the (possibly modified) config dict.
  """
  join_schema = _load_schema_file(join_schema_f)
  defaults    = _get_default_config_from_schema(join_schema, config)
  snakemake.utils.update_config(defaults, config)
  config      = defaults

  _validate_object_against_schema(config, join_schema_f, "join config")

  # validate each referenced project config against the project schema
  proj_schema_f = pathlib.Path(join_schema_f).parent / "config.schema.json"
  if proj_schema_f.is_file():
    proj_schema   = _load_schema_file(proj_schema_f)
    for pid, proj_entry in config.get('projects', {}).items():
      cfg_f         = proj_entry.get('config')
      if cfg_f and os.path.isfile(cfg_f):
        with open(cfg_f) as f:
          proj_cfg      = yaml.safe_load(f)
        proj_defaults = _get_default_config_from_schema(proj_schema, proj_cfg)
        defaults_copy = copy.deepcopy(proj_defaults)
        snakemake.utils.update_config(defaults_copy, proj_cfg)
        proj_cfg      = defaults_copy
        proj_errors   = sorted(
          jsonschema.Draft202012Validator(proj_schema).iter_errors(proj_cfg),
          key=lambda e: e.path)
        if proj_errors:
          raise ValueError(f"Validation errors in project '{pid}' config ({cfg_f}):\n" +
            "\n".join(f"  {list(e.path)}: {e.message}" for e in proj_errors))

  for section in ['hvg', 'integration', 'marker_genes', 'shiny', 'resources']:
    if section not in config:
      config[section] = {}
    _apply_join_defaults(config[section],
      join_schema['properties'].get(section, {}).get('properties', {}))

  config = _check_shiny_parameters(config)

  # validate train_xgboost section if present
  if 'train_xgboost' in config:
    xgb = config['train_xgboost']
    _apply_join_defaults(xgb,
      join_schema['properties'].get('train_xgboost', {}).get('properties', {}))

    join_name = config['join']['name']
    if 'ref_tag' not in xgb:
      xgb['ref_tag'] = f"xgboost_{join_name}"
    else:
      xgb['ref_tag'] = f"xgboost_{xgb['ref_tag']}"

    # resolve annots_f
    join_dir = pathlib.Path(config['join']['proj_dir'])
    annots_f = pathlib.Path(xgb['annots_f'])
    if not annots_f.is_absolute():
      annots_f = join_dir / annots_f
    xgb['annots_f'] = str(annots_f)
    if not annots_f.is_file():
      raise FileNotFoundError(f"Annotations file not found: {annots_f}")
    annots_df = pl.read_csv(str(annots_f), n_rows=5)
    if "cell_id" not in annots_df.columns:
      raise ValueError("train_xgboost.annots_f must have 'cell_id' column")
    if "annotation" not in annots_df.columns:
      raise ValueError("train_xgboost.annots_f must have 'annotation' column")

    # check label_map_f if provided
    if xgb.get("label_map_f") is not None:
      label_map_f = pathlib.Path(xgb["label_map_f"])
      if not label_map_f.is_absolute():
        label_map_f = join_dir / label_map_f
      xgb['label_map_f'] = str(label_map_f)
      if not label_map_f.is_file():
        raise FileNotFoundError(f"Label map file not found: {label_map_f}")
      map_df = pl.read_csv(str(label_map_f), n_rows=5)
      if "annotation" not in map_df.columns:
        raise ValueError("label_map_f must have 'annotation' column")
      if "coarse_label" not in map_df.columns:
        raise ValueError("label_map_f must have 'coarse_label' column")

    # inject output_dir
    join_name = config['join']['name']
    xgb['output_dir'] = str(join_dir / f"output/{join_name}_train_xgboost")

  # derive marker gene paths (mirrors _check_marker_genes_parameters for non-join)
  scdata_dir = pathlib.Path(scdata_dir)
  config['marker_genes']['mkr_gsea_dir'] = str(scdata_dir / 'gmt_pathways')
  join_dir = pathlib.Path(config['join']['proj_dir'])
  custom_mkr_names, custom_mkr_paths = _get_custom_marker_genes_specs(
    config, scdata_dir, join_dir)
  config['marker_genes']['custom_mkr_names'] = custom_mkr_names
  config['marker_genes']['custom_mkr_paths'] = custom_mkr_paths

  # look up gene_info_f from index_parameters.csv
  is_flex   = config['join'].get('tenx_assay_type', 'poly_a') == 'flex'
  ref_field = 'probe_set' if is_flex else 'ref_txome'
  join_ref  = config['join'].get(ref_field)
  idx_params_f = scdata_dir / 'index_parameters.csv'
  idx_params = pl.read_csv(idx_params_f)
  config['marker_genes']['gene_info_f'] = idx_params.filter(
    pl.col('reference') == join_ref)['gene_info_f'][0]

  # check ref_txome / probe_set consistency across all projects
  for pid, proj_entry in config.get('projects', {}).items():
    cfg_f = proj_entry.get('config')
    if cfg_f and os.path.isfile(cfg_f):
      with open(cfg_f) as f:
        proj_cfg = yaml.safe_load(f)
      proj_is_flex = proj_cfg['project'].get('tenx_assay_type', 'poly_a') == 'flex'
      proj_field   = 'probe_set' if proj_is_flex else 'ref_txome'
      if proj_field != ref_field:
        raise ValueError(
          f"Cannot join projects with different reference types: "
          f"join uses {ref_field}, but project '{pid}' uses {proj_field}")
      proj_ref = proj_cfg['project'].get(ref_field)
      if proj_ref and proj_ref != join_ref:
        raise ValueError(
          f"Project '{pid}' {ref_field}={proj_ref!r} does not match "
          f"join {ref_field}={join_ref!r}")

  # check h5ads YAML files exist (integration must be complete)
  # h5ads are always in the integration directory, even when using zoom outputs
  for pid, proj_entry in config.get('projects', {}).items():
    cfg_f = proj_entry.get('config')
    if cfg_f and os.path.isfile(cfg_f):
      with open(cfg_f) as f:
        proj_cfg = yaml.safe_load(f)
      proj_dir    = pathlib.Path(proj_cfg['project']['proj_dir'])
      short_tag   = proj_cfg['project']['short_tag']
      full_tag    = proj_cfg['project']['full_tag']
      date_stamp  = proj_cfg['project']['date_stamp']
      int_dir     = proj_dir / f"output/{short_tag}_integration"
      h5ads_f     = int_dir / f"h5ads_clean_paths_{full_tag}_{date_stamp}.yaml"
      if not h5ads_f.is_file():
        raise FileNotFoundError(
          f"h5ads YAML not found for project '{pid}': {h5ads_f}\n"
          f"  scprocess integration must be completed for this project before running join.")

  # validate label_celltypes models and resolve paths
  lbl_cfg     = config.get('label_celltypes', [])
  if lbl_cfg:
    lbl_schema_props = join_schema['properties']['label_celltypes']['items']['properties']
    for entry in lbl_cfg:
      for key, prop in lbl_schema_props.items():
        if key not in entry and 'default' in prop:
          entry[key] = prop['default']

    typist_ls_f = scdata_dir / 'celltypist/celltypist_models.csv'
    mdls_typist = pl.read_csv(typist_ls_f)['model'].to_list() if typist_ls_f.is_file() else []
    mdls_scproc = ['human_cns']
    for entry in lbl_cfg:
      if entry['labeller'] == 'celltypist':
        if entry['model'] not in mdls_typist:
          raise ValueError(f"CellTypist model '{entry['model']}' not found. Valid: {', '.join(mdls_typist)}")
      elif entry['labeller'] == 'scprocess':
        if entry['model'] not in mdls_scproc:
          raise ValueError(f"scprocess model '{entry['model']}' not found. Valid: {', '.join(mdls_scproc)}")
        xgb_dir = scdata_dir / 'xgboost'
        if entry['model'] == 'human_cns':
          entry['xgb_f']     = str(xgb_dir / 'Siletti_Macnair-2025-07-23/xgboost_obj_hvgs_Siletti_Macnair_2025-07-23.rds')
          entry['xgb_cls_f'] = str(xgb_dir / 'Siletti_Macnair-2025-07-23/allowed_cls_Siletti_Macnair_2025-07-23.csv')
        if not pathlib.Path(entry['xgb_f']).is_file():
          raise FileNotFoundError(f"XGBoost model file not found: {entry['xgb_f']}")
        if not pathlib.Path(entry['xgb_cls_f']).is_file():
          raise FileNotFoundError(f"XGBoost classes file not found: {entry['xgb_cls_f']}")

  # check that each metadata_var is present in at least one project's sample
  # metadata file. Unlike the standard pipeline (which requires all vars in a
  # single file), join allows vars to be present in only a subset of projects
  # (e.g. when projects have different experimental designs).
  metadata_vars = config['join'].get('metadata_vars', [])
  if metadata_vars:
    for var in metadata_vars:
      found = False
      for pid, proj_entry in config.get('projects', {}).items():
        cfg_f = proj_entry.get('config')
        if cfg_f and os.path.isfile(cfg_f):
          with open(cfg_f) as f:
            proj_cfg = yaml.safe_load(f)
          # resolve sample metadata path (may be relative to project dir)
          meta_f = proj_cfg['project'].get('sample_metadata', '')
          if meta_f:
            proj_dir = pathlib.Path(proj_cfg['project']['proj_dir'])
            meta_p   = pathlib.Path(meta_f)
            if not meta_p.is_absolute():
              meta_p = proj_dir / meta_p
            # read header only to check column names
            if meta_p.is_file():
              samples_df = pl.read_csv(meta_p, n_rows=0)
              if var in samples_df.columns:
                found = True
                break
      if not found:
        raise KeyError(
          f"metadata_var '{var}' not found in any project's sample metadata file")

  return config
