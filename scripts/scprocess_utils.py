# load modules
import os
import sys
import re
import pathlib
import warnings
import yaml
import pandas as pd
import polars as pl
import csv
import math
import glob
import gzip
import datetime
import subprocess
import snakemake
import json
import jsonschema

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
  setup_configfile = scdata_dir / 'scprocess_setup.yaml'
  if not os.path.exists(setup_configfile):
    raise FileNotFoundError(f"scprocess_setup.yaml does not exist in {scdata_dir}")

  # load setup config file
  with open(setup_configfile, "r") as stream:
    setup_config      = yaml.safe_load(stream)

  # if there is a cluster profile check it
  if ('profile' in setup_config) and setup_config['profile'] is not None:
    # check if profile exists
    profile     = setup_config['profile']
    profile_dir = scprocess_dir / 'profiles' / profile
    profile_f   = profile_dir / 'config.yaml'
    if not os.path.isfile(profile_f):
      raise FileNotFoundError(f"cluster configuration file {profile_f} does not exist")
    
    # if ok, add profile to snakemake call
    extraargs = extraargs + ' --workflow-profile ' + str(profile_dir)

  return scdata_dir, extraargs


### much checking

# wrapper for checking
def check_config(config, schema_f, scdata_dir, scprocess_dir):
  # start with defaults, overwrite with config values
  schema      = _load_schema_file(schema_f)
  defaults    = _get_default_config_from_schema(schema)
  snakemake.utils.update_config(defaults, config)
  config      = defaults

  # check file is ok
  _validate_object_against_schema(config, schema_f, "config")

  # get parameters
  config      = _check_project_parameters(config, scdata_dir, scprocess_dir)
  config      = _check_multiplexing_parameters(config)
  config      = _check_mapping_parameters(config, scdata_dir)
  config      = _check_ambient_parameters(config)
  config      = _check_qc_parameters(config)
  config      = _check_hvg_parameters(config)
  config      = _check_integration_parameters(config)
  config      = _check_marker_genes_parameters(config, scdata_dir)
  config      = _check_pb_empties_parameters(config)

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


def _get_default_config_from_schema(schema):
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
    print(f"problem with your {file_desc} file:\n  {str(e)}")
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
def _check_project_parameters(config, scprocess_data_dir, scprocess_dir):
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
    KeyError('"project" part of config file must contain exactly one of "fastq_dir" and "arv_uuids"')

  # do some checks if fastq_dir is specified
  if has_fastq and not has_arv_uuids:
    config["project"]["fastq_dir"] = _check_path_exists_in_project(config["project"]["fastq_dir"], config, what = "dir")

  # check if selected species is valid
  index_params_f    = scprocess_data_dir / 'index_parameters.csv'

  # from index_parameters.csv get valid values for species
  index_params      = pd.read_csv(index_params_f)
  valid_species     = index_params['genome_name'].tolist()
  valid_species_str = ', '.join(valid_species)
  if not config['project']['species'] in valid_species:
    raise ValueError(f"species {config['project']['species']} not defined. Valid values are {valid_species_str}")

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


def _check_samples_df(samples_df, config):
  # check for sample_id
  if "sample_id" not in samples_df.columns:
    raise KeyError(f"'sample_id' not present in sample metadata file")
  
  # check for pool_id
  if not config['multiplexing']['demux_type'] == "none":
    if "pool_id" not in samples_df.columns:
      raise KeyError(f"'pool_id' not present in sample metadata file")

  # some checks for multiplexing
  if config['multiplexing']['demux_type'] == "hto":
    # check that hto fastq path is present
    if "hto_id" not in samples_df.columns:
      raise KeyError("'hto_id' not present in sample metadata")

  # check that sample_id values are unique
  if not samples_df[ "sample_id" ].n_unique() == samples_df.shape[0]:
    raise ValueError("'sample_id' values in metadata csv not unique")

  # check columns of samples_df
  if any(' ' in col for col in samples_df.columns):
    raise ValueError("some column names in metadata csv contain spaces.")

  # sort out metadata variables
  if 'metadata_vars' in config["project"]:
    # load up sample file
    for var in config["project"]["metadata_vars"]:
      # check variable exists
      if not var in samples_df.columns:
        raise KeyError(f"{var} not a column in sample metadata file")

      # check that there are less than 10 unique values (otherwise probably not a categorical variable)
      var_col     = samples_df[var]
      if var_col.dtype == pl.String:
        continue
      if var_col.unique().shape[0] > 10:
        raise ValueError(f"{var} variable has more than 10 unique values; prob not a categorical variable")
  else:
    config['project']['metadata_vars'] = []

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
      KeyError('"multiplexing" part of config file must contain exactly one of "fastq_dir" and "arv_uuids"')

    # do some checks if fastq_dir is specified
    if has_fastq and not has_arv_uuids:
      config["multiplexing"]["fastq_dir"] = _check_path_exists_in_project(config["multiplexing"]["fastq_dir"], config, what = "dir")

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

    # check if looks ok 
    demux_df    = pl.read_csv(config['multiplexing']['demux_output'])
    for col in ["pool_id", "sample_id", "cell_id"]:
      if not col in demux_df.columns:
        raise KeyError(f"{col} not present in demux_output")

    # check if samples in metadata and demux_df match
    if set(demux_df['sample_id']) > set(samples_df['sample_id']):
      raise ValueError("Some values for 'sample_id' in demux_output don't have a match in sample_metadata")
    if set(samples_df['sample_id']) > set(demux_df['sample_id']):
      raise ValueError("Some values for 'sample_id' in sample_metadata don't have a match in demux_output")    
    if set(demux_df['pool_id']) != set(samples_df['pool_id']):
      raise ValueError("Values for pool_id don't match across demux_output and sample_metadata")
      
  return config


# check parameters for mapping
def _check_mapping_parameters(config, scdata_dir):
  # from index_parameters.csv get valid values for species
  idx_params_f  = scdata_dir / 'index_parameters.csv'
  index_params  = pl.read_csv(idx_params_f)
     
  # get mito strings from setup params
  species           = config['project']['species']
  config['mapping'] = {}
  config['mapping']['af_mito_str'] = index_params.filter(pl.col('genome_name') == species)['mito_str'][0]

  # get af index directory and check if exists
  config['mapping']['alevin_fry_home']  = scdata_dir / 'alevin_fry_home'
  config['mapping']['af_index_dir']     = scdata_dir / 'alevin_fry_home' / species
  if not pathlib.Path(config['mapping']['af_index_dir']).is_dir():
    raise FileNotFoundError(f"alevin index for {species} doesn't exist")
  
  # get gtf txt file, check that exists
  config['mapping']['af_gtf_dt_f'] = index_params.filter(pl.col('genome_name') == species)['gtf_txt_f'][0]

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
def _check_qc_parameters(config):
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


# check parameters for pb and empties
def _check_pb_empties_parameters(config):
  # nothing to do here at the moment; leaving in case it's useful later

  return config


# check parameters for hvgs
def _check_hvg_parameters(config):
  # define dummy group names for all
  if config['hvg']['hvg_method'] == 'all':
    config['hvg']['hvg_group_names'] = ['all_samples']
  # if groups, check that the values are ok
  elif config['hvg']['hvg_method'] == 'groups':
    # check that value of metadata_split_var matches a column in sample metadata
    hvg_split_var = config['hvg']['hvg_metadata_split_var']
    meta          = pd.read_csv(config['project']['sample_metadata'])
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
    gtf_df      = pd.read_csv(config['project']['af_gtf_dt_f'],  sep = '\t')
    num_genes   = gtf_df.shape[0]

    # chunk them up and name them
    num_chunks  = (num_genes + config['hvg']['hvg_chunk_size'] - 1) // config['hvg']['hvg_chunk_size']
    chunk_names = [f"chunk_{i+1}" for i in range(num_chunks)]
    
    # add to config
    config['hvg']['hvg_num_chunks'] = num_chunks
    config['hvg']['hvg_chunk_names'] = chunk_names

  return config


# check parameters for integration
def _check_integration_parameters(config):
  # nothing to do here at the moment; leaving in case it's useful later

  return config


# check parameters for marker genes
def _check_marker_genes_parameters(config, scdata_dir):
  # set some more default values
  config['marker_genes']['mkr_gsea_dir'] = scdata_dir / 'gmt_pathways'

  # get custom marker files
  custom_mkr_names, custom_mkr_paths = _get_custom_marker_genes_specs(config, scdata_dir)
  config['marker_genes']['custom_mkr_names'] = custom_mkr_names
  config['marker_genes']['custom_mkr_paths'] = custom_mkr_paths

  return config


# check specified custom marker genes
def _get_custom_marker_genes_specs(config, scdata_dir):
  # set defaults
  custom_mkr_names = ""
  custom_mkr_paths = ""

  # populate with custom sets
  if 'mkr_custom_genesets' in config["marker_genes"]:
    mkr_names = []
    mkr_paths = []
    for i, gene_set in enumerate(config["marker_genes"]["mkr_custom_genesets"]):
      # get name and file
      name      = gene_set["name"]
      # second argument is default in case file is missing
      file_path = pathlib.Path(gene_set.get("file", scdata_dir / 'marker_genes' / f"{name}.csv"))

      # check whether it exists
      if not file_path.is_absolute():
        file_path = config['project']['proj_dir'] / file_path
      if not file_path.is_file():
        raise FileNotFoundError(f"File not found for marker set '{name}'")
      if not file_path.suffix == ".csv":
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
    zoom_ls       = [_get_one_zoom_parameters(zoom_yaml_f, zoom_schema_f, config, scdata_dir) for zoom_yaml_f in zoom_yamls]
    ZOOM_PARAMS   = {z['zoom']['name']: z for z in zoom_ls}

  return ZOOM_PARAMS


# get parameters for one zoom specification
def _get_one_zoom_parameters(zoom_yaml_f, zoom_schema_f, config, scdata_dir):
  # check file exists
  zoom_yaml_f   = _check_path_exists_in_project(zoom_yaml_f, config, what = "file")

  # load things specified by zoom
  with open(zoom_yaml_f, "r") as stream:
    zoom_config   = yaml.safe_load(stream)

  # update with zoom defaults if not specified
  zoom_schema   = _load_schema_file(zoom_schema_f)
  zoom_defaults = _get_default_config_from_schema(zoom_schema)
  snakemake.utils.update_config(zoom_defaults, zoom_config)
  zoom_config   = zoom_defaults

  # check file is ok
  _validate_object_against_schema(zoom_config, zoom_schema_f, "zoom config")

  # start with defaults, overwrite with config values
  defaults    = config.copy()
  snakemake.utils.update_config(defaults, zoom_config)
  zoom_config = defaults

  # get useful things
  SHORT_TAG   = config['project']['short_tag']
  FULL_TAG    = config['project']['full_tag']
  DATE_STAMP  = config['project']['date_stamp']

  # find file for each option
  if zoom_config['zoom']['labels_source'] == 'clusters':
    labels_f    = f"output/{SHORT_TAG}_integration/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz"

  # if using xgboost or celltypist, check those things
  elif zoom_config['zoom']['labels_source'] in ['celltypist', 'xgboost']:
    labeller    = zoom_config['zoom']['labels_source']
    model       = zoom_config['zoom']['model']
    labels_f    = f"output/{SHORT_TAG}_label_celltypes/labels_{labeller}_model_{model}_{FULL_TAG}_{DATE_STAMP}.csv.gz"

  # unpack
  elif zoom_config['zoom']['labels_source'] == 'custom':
    labels_f    = pathlib.Path(zoom_config['zoom']['custom_labels_f'])
  
  # check file exists
  labels_f    = _check_path_exists_in_project(labels_f, config, what = "file")
  
  # get list of all clusters to check if cluster names are valid
  sel_labels  = _check_zoom_clusters_in_file(labels_f, zoom_config)

  # add this file to params list
  zoom_config['zoom']['labels_f']   = labels_f
  zoom_config['zoom']['sel_labels'] = sel_labels

  return zoom_config


# get list of clusters to check zoom specification
def _check_zoom_clusters_in_file(labels_f, zoom_config):
  # get list of clusters
  labels_var  = zoom_config['zoom']['labels_col']
  lbl_dt      = pl.read_csv(labels_f)
  labels      = lbl_dt[ labels_var ].unique().to_list()
  labels      = [cl for cl in labels if cl is not None]
  # labels      = sorted(labels)

  # check that the zoom clusters match these
  sel_labels  = zoom_config['zoom']['sel_labels']
  missing_cls = set(sel_labels) - set(labels)
  if len(missing_cls) > 0:
    raise ValueError(f"the following labels were specified in the zoom params yaml but are not present in the file:\n  {", ".join(missing_cls)}")

  return sel_labels


### much input

# get variables for each run
def get_run_parameters(config, scprocess_data_dir):
  # define run variable
  if config['multiplexing']['demux_type'] == "none":
    RUN_VAR       = "sample_id"
  else:
    RUN_VAR       = "pool_id"

  # move to another function
  metadata_f  = config["project"]["sample_metadata"]
  samples_df  = pl.read_csv( metadata_f )
  RUNS        = samples_df[ RUN_VAR ].drop_nulls().to_list()

  # get any custom parameters
  custom_run_params = _get_custom_parameters(config, RUN_VAR)

  # should we exclude any runs?
  RUNS        = _do_exclusions(RUNS, config, RUN_VAR)

  # get fastq files
  RNA_FQS     = _get_fastqs(config, RUNS, is_hto = False)
  RUNS        = list(RNA_FQS.keys())

  # get HTO files
  HTO_FQS     = {}
  if config['multiplexing']['demux_type'] == "hto":
    HTO_FQS     = _get_fastqs(config, RUNS, is_hto = True)
    RUNS        = [r for r in RUNS if r in HTO_FQS]

  # load sample file, populate everything from config
  RUN_PARAMS  = {
    run_name: _get_run_parameters_one_run(run_name, config, RNA_FQS, HTO_FQS, scprocess_data_dir, custom_run_params)
    for run_name in RUNS
  }

  return RUN_PARAMS, RUN_VAR


# get all fastq files
def _get_fastqs(config, RUNS, is_hto = False):
  # get place to look for fastq files
  if is_hto:
    tmp_ls      = config['multiplexing']
  else:
    tmp_ls      = config['project']
  if "fastq_dir" in tmp_ls:
    fastq_dir   = tmp_ls['fastq_dir']
    arv_uuids   = None
  else:
    fastq_dir   = None
    arv_uuids   = tmp_ls['arv_uuids']

  # get 
  if fastq_dir is not None:
    fastq_dict    = _list_fastq_files_dir(fastq_dir)
  elif arv_uuids is not None:
    fastq_dict    = _list_fastq_files_arvados(arv_uuids)

  # get fastq files for each sample
  wheres        = fastq_dict["wheres"]
  fastq_fs      = fastq_dict["fastqs"]
  fq_sizes_gb   = fastq_dict["fastq_sizes"]

  # loop through samples
  fastqs        = {}
  for run in RUNS:
    # get R1 and R2 files matching each run
    R1_regex      = rf".*{run}.*(_|\.)R1.*\.fastq\.gz"
    R1_fs         = [f for f in fastq_fs if re.match(R1_regex, f) ]
    R2_regex      = rf".*{run}.*(_|\.)R2.*\.fastq\.gz"
    R2_fs         = [f for f in fastq_fs if re.match(R2_regex, f) ]

    # find where
    wheres_R1     = [ wheres[i] for i,f in enumerate(fastq_fs) if re.match(R1_regex, f) ]
    wheres_R2     = [ wheres[i] for i,f in enumerate(fastq_fs) if re.match(R2_regex, f) ]
    this_where    = list(set(wheres_R1 + wheres_R2))
    if len(this_where) > 1:
      raise ValueError(f"FASTQ files for run {run} seem to be found in more than location")

    # get R1 filesize
    R1_fs_size_gb = [ fq_sizes_gb[i] for i,f in enumerate(fastq_fs) if re.match(R1_regex, f) ]

    # check have full set of files
    check_R1      = [re.sub(r'(?<=(_|\.))R1', 'R0', f) for f in R1_fs]
    check_R2      = [re.sub(r'(?<=(_|\.))R2', 'R0', f) for f in R2_fs]
    if len(R1_fs) == 0:
      print(f"  WARNING: no {"hto " if is_hto else ""}fastq files found for run {run}; excluded.")
    elif set(check_R1) != set(check_R2):
      print(f"  WARNING: {"hto " if is_hto else ""}fastq files found for run {run} but R1 and R2 don't match; excluded.")
    else:
      fastqs[run] = {
        "where":          this_where[0], 
        "R1_fs":          R1_fs, 
        "R2_fs":          R2_fs, 
        "R1_fs_size_gb":  round(sum(R1_fs_size_gb), 1)
      }

  return fastqs


# get all fastq files in directory
def _list_fastq_files_dir(fastq_dir):
  # get all files
  all_fs      = os.listdir(fastq_dir)

  # filter to just fastqs
  fastq_fs        = [ f for f in all_fs if re.match(r".+\.fastq\.gz", f) ]
  wheres          = [ fastq_dir for f in all_fs]
  fastq_sizes_gb  = [ (fastq_dir / f).stat().st_size  / BYTES_PER_GB for f in fastq_fs ]

  return { "wheres": wheres, "fastqs": fastq_fs, "fastq_sizes": fastq_sizes_gb }


# get all fastq files in all arvados uuids
def _list_fastq_files_arvados(arv_uuids):
  # get for each UUID
  wheres      = []
  fastqs      = []
  fastq_sizes = []

  # get all fastq files in given arvados uuid
  def _list_fastq_files_arvados_one_uuid(arv_uuid):
    # import relevant packages
    import arvados
    import collections
    import pathlib

    # define variables
    arv_files   = []
    wheres      = {}
    file_sizes  = {}

    # set up arvados access
    arv_token   = os.environ["ARVADOS_API_TOKEN"]
    arv_client  = arvados.api('v1', host = 'api.arkau.roche.com',
      token = arv_token, insecure = True, num_retries = 2 )

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

    return { "wheres": wheres, "fastqs": fastq_fs, "fastq_sizes": fastq_sizes_gb }

  # Iterate through each UUID in the list
  for arv_uuid in arv_uuids:
    # Get the dictionary result for one UUID
    result = _list_fastq_files_arvados_one_uuid(arv_uuid)

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


# get parameters for one 10x run
def _get_run_parameters_one_run(run_name, config, RNA_FQS, HTO_FQS, scdata_dir, custom_run_params):
  # set defaults
  sample_chem = config['project']['tenx_chemistry']
  knee1       = ""
  shin1       = ""
  knee2       = ""
  shin2       = ""
  if run_name in custom_run_params:
    if 'tenx_chemistry' in custom_run_params[run_name]:
      sample_chem = custom_run_params[run_name]['tenx_chemistry']
    if 'mapping' in custom_run_params[run_name]:
      if 'knee1' in custom_run_params[run_name]['mapping']:
        knee1       = custom_run_params[run_name]['mapping']['knee1']
      if 'shin1' in custom_run_params[run_name]['mapping']:
        shin1       = custom_run_params[run_name]['mapping']['shin1']
      if 'knee2' in custom_run_params[run_name]['mapping']:
        knee2       = custom_run_params[run_name]['mapping']['knee2']
      if 'shin2' in custom_run_params[run_name]['mapping']:
        shin2       = custom_run_params[run_name]['mapping']['shin2']

  # get af chemistry and expected orientation
  if sample_chem in ['3v2', '5v1', '5v2']:
    af_chemistry = '10xv2' 
  else: 
    af_chemistry = '10xv3'

  # get expected orientation
  if sample_chem in ['5v1', '5v2', '5v3']:
    expected_ori = 'rc'
  else:
    expected_ori = 'fw'

  # sort out whitelist file
  wl_df_f     = scdata_dir / 'cellranger_ref/cellranger_whitelists.csv'
  wl_df       = pl.read_csv(wl_df_f)
  wl_f        = wl_df.filter( pl.col('chemistry') == sample_chem )['barcodes_f'].item()
  wl_trans_f  = wl_df.filter( pl.col('chemistry') == sample_chem )['translation_f'].item()
  if type(wl_trans_f) == str:
    whitelist_trans_f = scdata_dir / 'cellranger_ref' / wl_trans_f
  else:
    whitelist_trans_f = None
  whitelist_f = scdata_dir / 'cellranger_ref' / wl_f

  # make dictionary for mapping
  mapping_dc  = {
    "where":              RNA_FQS[run_name]["where"],
    "R1_fs":              RNA_FQS[run_name]["R1_fs"],
    "R2_fs":              RNA_FQS[run_name]["R2_fs"],
    "R1_fs_size_gb":      RNA_FQS[run_name]["R1_fs_size_gb"],
    "af_chemistry":       af_chemistry, 
    "expected_ori":       expected_ori, 
    "whitelist_f":        whitelist_f, 
    "whitelist_trans_f":  whitelist_trans_f,
    "knee1":              knee1,
    "shin1":              shin1,
    "knee2":              knee2,
    "shin2":              shin2
  }

  # make dictionary for mapping
  if config['multiplexing']['demux_type'] == "hto":
    multiplexing_dc  = {
      "where":              HTO_FQS[run_name]["where"],
      "R1_fs":              HTO_FQS[run_name]["R1_fs"],
      "R2_fs":              HTO_FQS[run_name]["R2_fs"],
      "R1_fs_size_gb":      HTO_FQS[run_name]["R1_fs_size_gb"],
      "af_chemistry":       af_chemistry, 
      "whitelist_f":        whitelist_f,
      "whitelist_trans_f":  whitelist_trans_f
    }
  else:
    multiplexing_dc   = {}

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
    "mapping":      mapping_dc,
    "multiplexing": multiplexing_dc,
    "ambient":      ambient_dc
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
    for batch_name in BATCHES
  }

  return BATCH_PARAMS, BATCH_VAR, SAMPLES


# get parameters for one sample
def _get_batch_parameters_one_batch(batch_name, config, custom_batch_params):
  # make dictionary for mapping
  qc_dc   = {
    # set defaults
    "qc_min_counts":  config['qc']['qc_min_counts'],
    "qc_min_feats":   config['qc']['qc_min_feats'],
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

  # make dict of dicts
  out_dc  = {
    "qc":   qc_dc
  }

  return out_dc


# get mapping from samples to runs
def get_runs_to_batches(config, RUNS, BATCHES, BATCH_VAR):
  # if nothing to do return none
  if config['multiplexing']['demux_type'] == "none":
    # make dictionary if we can
    if not RUNS == BATCHES:
      raise ValueError("RUNS and BATCHES should be identical when demux_type is 'none'")
    RUNS_TO_BATCHES = { s: [s] for s in BATCHES }
    RUNS_TO_SAMPLES = { s: [s] for s in BATCHES }

  else:
    # get sample_metadata, convert to dictionary
    sample_metadata = pl.read_csv(config['project']['sample_metadata'])
    tmp_df          = sample_metadata.group_by("pool_id").agg(pl.col("sample_id").alias("sample_id_list"))
    pool_vals       = tmp_df["pool_id"]
    RUNS_TO_SAMPLES = { pool_id: tmp_df.filter(pl.col("pool_id") == pool_id)["sample_id_list"].to_list()[0] for pool_id in pool_vals }
    
    # filter out any RUNS that shouldn't be there
    RUNS_TO_SAMPLES = { pool_id: sample_ids for pool_id, sample_ids in RUNS_TO_SAMPLES.items() if pool_id in RUNS }

    # now choose for runs to batches
    if BATCH_VAR == "pool_id":
      # make dictionary if we can
      if not RUNS == BATCHES:
        raise ValueError("RUNS and BATCHES should be identical if 'int_batch_var' is 'pool_id'")
      RUNS_TO_BATCHES = { s: [s] for s in BATCHES }

    elif BATCH_VAR == "sample_id":
      # filter out any samples that shouldn't be there
      for pool_id in RUNS_TO_SAMPLES:
        RUNS_TO_SAMPLES[pool_id] = [ sample_id for sample_id in RUNS_TO_SAMPLES[pool_id] if sample_id in BATCHES ]
      # duplicate for batches
      RUNS_TO_BATCHES = RUNS_TO_SAMPLES

  return RUNS_TO_BATCHES, RUNS_TO_SAMPLES


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
  mdls_typist     = pl.read_csv(typist_ls_f)['model'].to_list()
  mdls_scprocess  = ['human_cns']

  # check that selected models are valid
  def _check_one_label_celltypes_parameters(entry):
    # check that parameters for celltypist are ok
    if entry['labeller'] == 'celltypist':
      if not entry['model'] in mdls_typist:
        raise KeyError(
          f"The value {entry['model']} specified in label_celltypes is not a valid celltypist model.\n"
          f"The following are valid models:\n{", ".join(mdls_typist)}"
          )

    # check that parameters for scprocess are ok
    elif entry['labeller'] == 'scprocess':
      if not entry['model'] in mdls_scprocess:
        raise KeyError(
          f"the value {entry['model']} specified in label_celltypes is not a valid scprocess model"
          f"These models are currently available: {", ".join(mdls_scprocess)}"
          )
    
      # pick labeller
      xgb_dir  = os.path.join(scdata_dir, 'xgboost')
      if not pathlib.Path(xgb_dir).is_dir():
       raise FileNotFoundError(f"xgboost directory '{xgb_dir}' not found")
    
      # get relevant values for labeller
      if entry['model'] == 'human_cns':
        entry['xgb_f']      = os.path.join(xgb_dir, "Siletti_Macnair-2025-07-23/xgboost_obj_hvgs_Siletti_Macnair_2025-07-23.rds")
        entry['xgb_cls_f']  = os.path.join(xgb_dir, "Siletti_Macnair-2025-07-23/allowed_cls_Siletti_Macnair_2025-07-23.csv")

      # check these are ok
      if not pathlib.Path(entry['xgb_f']).is_file():
        raise FileNotFoundError(f"file {entry['xgb_f']} doesn't exist; consider (re)runnning scprocess setup")
      if not pathlib.Path(entry['xgb_cls_f']).is_file():
        raise FileNotFoundError(f"file {entry['xgb_cls_f']} doesn't exist; consider (re)runnning scprocess setup")

    # check that parameters for scprocess are ok
    elif entry['labeller'] == 'custom':
      entry['custom_f']   = _check_path_exists_in_project(entry['custom_f'], config, what = "file")
  
    # add defaults
    for v in label_defaults:
      if not v in entry:
        entry[v] = label_defaults[v]

    return entry

  # apply this to each specified model
  LABELLER_PARAMS = [ _check_one_label_celltypes_parameters(entry) for entry in config['label_celltypes'] ]

  return LABELLER_PARAMS



def get_resources(rule, all_rules, param, lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, run = None):
  
  if not hasattr(all_rules, rule):
    raise ValueError(f'rule {rule} is not defined.')
  
  # get lm params
  lm_df = pl.read_csv(lm_f)

  # get resource name in the config 
  if param == 'time':
    config_param_name = f'mins_{rule}'
  else: 
    config_param_name = f'gb_{rule}'

  # get lm params for rule and param (memory or time)
  filt_lm_df    = lm_df.filter((pl.col("param") == param) & (pl.col("rule") == rule))

  # Load schema and extract default resource values
  schema        = _load_schema_file(schema_f)
  defaults      = _get_default_config_from_schema(schema)
  res_defaults  = defaults['resources']

  # if no lm params are defined
  if filt_lm_df["rq_slope"].is_null().all():
    # make sure default is specified in schema
    if config_param_name not in res_defaults.keys():
      raise ValueError(f'Default value for {config_param_name} is missing from JSON schema.')
    # use values from config file
    param_val = config['resources'].get(config_param_name, None)
    if param == 'memory':
      param_val *= MB_PER_GB

  else: # if lm params are defined
    # make sure no default value defined in schema
    if config_param_name in res_defaults.keys():
      raise ValueError(f'Default value for {config_param_name} should not be specified in JSON schema.')
    
    param_val = config['resources'].get(config_param_name, None)
    if param_val is not None:
      if param == 'memory':
        param_val *= MB_PER_GB
    else:
      # get rq params
      x_rq      = filt_lm_df['model_var'].unique().to_list()[0]
      intercept = filt_lm_df['rq_intercept'].unique().to_list()[0] 
      slope     = filt_lm_df['rq_slope'].unique().to_list()[0]
      buffer    = filt_lm_df['buffer'].unique().to_list()[0]

      # get the name of x var
      if x_rq.startswith('input.'):
        input_attr = x_rq.replace("input.", "")
        if hasattr(input, input_attr):
          x_val = os.path.getsize(getattr(input, input_attr)) // MB_PER_GB**2 
        else:
          raise ValueError(f"'{input_attr}' is not a valid input attribute for rule '{rule}'.")
      elif x_rq == 'raw_data_size':
        # use raw data size
        if run is None:
          raise ValueError(f'run argument should be defined')
        x_val = RUN_PARAMS[run]["mapping"]["R1_fs_size_gb"]
      elif x_rq == 'n_smpls_pre_qc':
        # use the number of samples
        x_val = len(BATCHES)
      else:
        raise ValueError(f"Unknown variable '{x_rq}' for scaling resources.")

      param_val = intercept + slope * x_val

      if param == 'time':
        param_val /= 60  # Convert minutes to hours
      param_val += buffer

  return param_val


### helpers

# some useful global variables
BYTES_PER_GB  = 1024**3
MB_PER_GB     = 1024


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
  df_ls     = [ pl.read_csv(f) for f in in_files if gzip.open(f, 'rb').read(1) ]
  df_merged = pl.concat(df_ls)
  with gzip.open(out_file, 'wb') as f:
    df_merged.write_csv(f)


# HVGs function: make df with list of chunked counts files
def make_hvgs_input_df(runs, ambient_outs_yamls, RUN_VAR, BATCH_VAR, BATCHES_TO_RUNS, 
  DEMUX_TYPE, FULL_TAG, DATE_STAMP, hvg_dir):
  # loop through ambient yaml files to populate list
  df_list = []
  for r, yaml_file in zip(runs, ambient_outs_yamls):
    # get filtered ambient outputs
    with open(yaml_file) as f:
      amb_outs = yaml.load(f, Loader=yaml.FullLoader)
    amb_filt_f = amb_outs['filt_counts_f']

    # if no multiplexing, simple
    if DEMUX_TYPE == "none":
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


# end
