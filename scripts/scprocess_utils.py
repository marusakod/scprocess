# load modules
import warnings
import yaml
import pandas as pd
import math
import os
import re
import glob
import datetime
import subprocess


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

    LABELS = this_zoom['labels']
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
    HVG_PARAMS = get_hvg_parameters(this_zoom, METADATA_F, AF_GTF_DT_F)
    HVG_KEYS = ["HVG_METHOD", "HVG_SPLIT_VAR", "HVG_CHUNK_SIZE", "HVG_NUM_CHUNKS", "HVG_GROUP_NAMES", "HVG_CHUNK_NAMES", "N_HVGS", "EXCLUDE_AMBIENT_GENES"]
    HVG_DICT = dict(zip(HVG_KEYS, HVG_PARAMS))

    # pb_empties params
    PB_EMPTIES_PARAMS = get_pb_empties_parameters(this_zoom) 
    PB_EMPTIES_KEYS = ["AMBIENT_GENES_LOGFC_THR", "AMBIENT_GENES_FDR_THR"]
    PB_EMPTIES_DICT = dict(zip(PB_EMPTIES_KEYS, PB_EMPTIES_PARAMS))

    # integration params
    INT_PARAMS = get_integration_parameters(this_zoom)
    INT_KEYS = ["INT_CL_METHOD", "INT_REDUCTION", "INT_N_DIMS", "INT_THETA", "INT_RES_LS"]
    INT_DICT = dict(zip(INT_KEYS, INT_PARAMS[:5]))   
    
    # marker gene params
    MKR_PARAMS = get_marker_genes_parameters(this_zoom, PROJ_DIR, SCPROCESS_DATA_DIR)
    MKR_KEYS = ["MKR_SEL_RES", "MKR_GSEA_DIR", "MKR_MIN_CL_SIZE", "MKR_MIN_CELLS", "MKR_NOT_OK_RE", "MKR_MIN_CPM_MKR", "MKR_MIN_CPM_GO", 
                "MKR_MAX_ZERO_P", "MKR_GSEA_CUT", "CUSTOM_MKR_NAMES", "CUSTOM_MKR_PATHS"]
    MKR_DICT = dict(zip(MKR_KEYS, MKR_PARAMS))

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

    if DEMUX_TYPE == 'af':
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
  # check expected variables are in the config file
  for v in ["proj_dir", "fastq_dir", "short_tag", "full_tag", "date_stamp", "your_name", "affiliation", "sample_metadata", "species"]:
    assert v in config, f"{v} needs to be an entry at the top level of the config file"

  ## what is specified in config directory?
  PROJ_DIR      = config["proj_dir"]
  FASTQ_DIR     = config["fastq_dir"]
  SHORT_TAG     = config["short_tag"]
  FULL_TAG      = config["full_tag"]
  YOUR_NAME     = config["your_name"]
  AFFILIATION   = config["affiliation"]
  DATE_STAMP    = config["date_stamp"]
  SPECIES       = config["species"]

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
  if not os.path.isabs(FASTQ_DIR):
    FASTQ_DIR = os.path.join(PROJ_DIR, FASTQ_DIR)
  
  # get samples
  METADATA_F    = config["sample_metadata"]
  if not os.path.isabs(METADATA_F):
    METADATA_F = os.path.join(PROJ_DIR, METADATA_F)

  samples_df    = pd.read_csv( METADATA_F )
  SAMPLES       = samples_df["sample_id"].dropna().tolist()

  # get multiplexing params
  DEMUX_TYPE, HTO_FASTQ_DIR, FEATURE_REF, DEMUX_F, BATCH_VAR, EXC_POOLS, POOL_IDS, POOL_SAMPLES, SAMPLE_MAPPING = \
  get_multiplexing_parameters(config, PROJ_DIR, samples_df)

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
    DEMUX_TYPE, HTO_FASTQ_DIR, FEATURE_REF, DEMUX_F, BATCH_VAR, EXC_POOLS, POOL_IDS, SAMPLE_VAR, SAMPLE_MAPPING


# remove samples and pools from sample mapping
def filter_sample_mapping(SAMPLE_MAPPING, POOL_IDS, SAMPLES):

  if SAMPLE_MAPPING is None:
    return None
  # filter to keep specific pool ids
  SAMPLE_MAPPING = {pool_id: sample_ids for pool_id, sample_ids in SAMPLE_MAPPING.items() if pool_id in POOL_IDS}
  # filter to keep specific sample_ids
  for pool_id in SAMPLE_MAPPING:
    SAMPLE_MAPPING[pool_id] = [sample_id for sample_id in SAMPLE_MAPPING[pool_id] if sample_id in SAMPLES]

  return SAMPLE_MAPPING


# define alevin parameters
def get_alevin_parameters(config, scprocess_data_dir, SPECIES):

  # check that alevin is in config
  assert "alevin" in config, \
    "alevin not defined in config file"
  
  # check that chemistry is defined in the config
  assert 'chemistry' in config['alevin'], \
    "chemistry not defined in config file"
  
  # get chemisty
  CHEMISTRY = config['alevin']['chemistry']
  # check if valid
  valid_chems = ['3LT', '3v2', '5v1', '5v2', '3v3', 'multiome', '3v4', '5v3']
  assert CHEMISTRY in valid_chems, \
    "chemistry not valid"
  
  # get setup params
  index_params_f  = os.path.join(scprocess_data_dir, 'index_parameters.csv')

  # from index_parameters.csv get valid values for species
  index_params= pd.read_csv(index_params_f)
     
  # get mito strings from setup params
  AF_MITO_STR = index_params.loc[index_params['genome_name'] == SPECIES, 'mito_str'].values[0]

  # get af index directory and check if exists
  AF_HOME_DIR = os.path.join(scprocess_data_dir, 'alevin_fry_home') # check if this exists in scprocess script
  AF_INDEX_DIR = os.path.join(AF_HOME_DIR, SPECIES)
  assert os.path.isdir(AF_INDEX_DIR), \
    f"alevin index for {SPECIES} doesn't exist"
  
  # get gtf txt file, check that exists
  AF_GTF_DT_F = index_params.loc[index_params['genome_name'] == SPECIES, 'gtf_txt_f'].values[0]

  return AF_MITO_STR, AF_HOME_DIR, AF_INDEX_DIR, AF_GTF_DT_F, CHEMISTRY


# ambient
def get_ambient_parameters(config):
  # set default values
  AMBIENT_METHOD                  = 'decontx'
  CELLBENDER_VERSION              = 'v0.3.0'
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
    CELLBENDER_IMAGE              = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.2'
  elif CELLBENDER_VERSION   == 'v0.3.0':
    CELLBENDER_IMAGE              = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.0'
  elif CELLBENDER_VERSION == 'v0.2.0':
    CELLBENDER_IMAGE              = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.2.0'
  else:
    raise ValueError(f"selected cellbender version {CELLBENDER_VERSION} not supported")

  # some checks on custom parameters for cellbender
      
  return CELLBENDER_IMAGE, CELLBENDER_VERSION, CELLBENDER_PROP_MAX_KEPT, AMBIENT_METHOD, CELL_CALLS_METHOD, \
    FORCE_EXPECTED_CELLS, FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE, CELLBENDER_POSTERIOR_BATCH_SIZE 




def get_qc_parameters(config):
  # set default values
  QC_HARD_MIN_COUNTS  = 200
  QC_HARD_MIN_FEATS   = 100
  QC_HARD_MAX_MITO    = 0.5
  QC_MIN_COUNTS       = 500
  QC_MIN_FEATS        = 300
  QC_MIN_MITO         = 0
  QC_MAX_MITO         = 0.1
  QC_MIN_SPLICE       = 0
  QC_MAX_SPLICE       = 0.75
  QC_MIN_CELLS        = 500
  DBL_MIN_FEATS       = 100
  EXCLUDE_MITO        = True

  # change defaults if specified
  if ('qc' in config) and (config['qc'] is not None):
    if 'qc_hard_min_counts' in config['qc']:
      QC_HARD_MIN_COUNTS  = config['qc']['qc_hard_min_counts']
    if 'qc_hard_min_feats'  in config['qc']:
      QC_HARD_MIN_FEATS   = config['qc']['qc_hard_min_feats']
    if 'qc_hard_max_mito'   in config['qc']:
      QC_HARD_MAX_MITO    = config['qc']['qc_hard_max_mito']
    if 'qc_min_counts'      in config['qc']:
      QC_MIN_COUNTS       = config['qc']['qc_min_counts']
    if 'qc_min_feats'       in config['qc']:
      QC_MIN_FEATS        = config['qc']['qc_min_feats']
    if 'qc_min_mito'        in config['qc']:
      QC_MIN_MITO         = config['qc']['qc_min_mito']
    if 'qc_max_mito'        in config['qc']:
      QC_MAX_MITO         = config['qc']['qc_max_mito']
    if 'qc_min_splice'      in config['qc']:
      QC_MIN_SPLICE       = config['qc']['qc_min_splice']
    if 'qc_max_splice'      in config['qc']:
      QC_MAX_SPLICE       = config['qc']['qc_max_splice']
    if 'qc_min_cells'       in config['qc']:
      QC_MIN_CELLS        = config['qc']['qc_min_cells']
    if 'dbl_min_feats'      in config['qc']:
      DBL_MIN_FEATS       = config['qc']['dbl_min_feats']
    if 'exclude_mito'       in config['qc']:
      EXCLUDE_MITO        = config['qc']['exclude_mito']
      EXCLUDE_MITO        = int(_safe_boolean(EXCLUDE_MITO))


  # make sure they're consistent
  QC_HARD_MIN_COUNTS  = min(QC_HARD_MIN_COUNTS, QC_MIN_COUNTS)
  QC_HARD_MIN_FEATS   = min(QC_HARD_MIN_FEATS, QC_MIN_FEATS)
  QC_HARD_MAX_MITO    = max(QC_HARD_MAX_MITO, QC_MAX_MITO)

  return QC_HARD_MIN_COUNTS, QC_HARD_MIN_FEATS, QC_HARD_MAX_MITO, QC_MIN_COUNTS, QC_MIN_FEATS, \
    QC_MIN_MITO, QC_MAX_MITO, QC_MIN_SPLICE, QC_MAX_SPLICE, QC_MIN_CELLS, DBL_MIN_FEATS, EXCLUDE_MITO


# define hvg parameters 
def get_hvg_parameters(config, METADATA_F, AF_GTF_DT_F): 

  # set defaults
  HVG_METHOD     = 'sample'
  HVG_SPLIT_VAR  = None
  HVG_CHUNK_SIZE = 2000
  N_HVGS         = 2000
  NUM_CHUNKS     = None
  EXCLUDE_AMBIENT_GENES = True
  GROUP_NAMES    = []
  CHUNK_NAMES    = []


  if ('hvg' in config) and (config['hvg'] is not None):

    if 'hvg_n_hvgs' in config['hvg']:
      N_HVGS          = config['hvg']['hvg_n_hvgs']
    
    if 'exclude_ambient_genes' in config['hvg']:
      EXCLUDE_AMBIENT_GENES = config['hvg']['hvg_exclude_ambient_genes']
      EXCLUDE_AMBIENT_GENES = _safe_boolean(EXCLUDE_AMBIENT_GENES)
      
    if 'hvg_method' in config['hvg']:
      HVG_METHOD      = config['hvg']['hvg_method']
      # check if valid
      valid_methods = ['sample', 'all', 'groups']
      assert HVG_METHOD in valid_methods, \
        f"Invalid hvg method '{HVG_METHOD}'. Must be one of {valid_methods}."
      
      if HVG_METHOD == 'all':
        GROUP_NAMES = ['all_samples']
      
      # if method is groups check that group variable is specified
      if HVG_METHOD == 'groups':
        assert ('hvg_metadata_split_var' in config['hvg']) and (config['hvg']['hvg_metadata_split_var'] is not None), \
          "The 'hvg_metadata_split_var' parameter must be defined when the hvg method is 'groups'."
        HVG_SPLIT_VAR = config['hvg']['hvg_metadata_split_var']

        # check that value of metadata_split_var matches a column in sample metadata
        meta = pd.read_csv(METADATA_F)
        assert HVG_SPLIT_VAR in meta.columns(), \
          f"{HVG_SPLIT_VAR} is not a column in the sample metadata file."
        
        # check number of unique group values
        uniq_groups = meta[HVG_SPLIT_VAR].unique().tolist()
        if len(uniq_groups) == meta.shape[0]:
          print(f"Number of unique values in '{HVG_SPLIT_VAR}' is the same as the number of samples; switching to 'sample' method for calculating hvgs.")
          HVG_METHOD = 'sample'
          HVG_SPLIT_VAR = None

        GROUP_NAMES = uniq_groups
        # replace spaces with underscores
        GROUP_NAMES = [n.replace(" ", "_") for n in GROUP_NAMES]

    # get number of gene chunks if method is 'groups' or 'all'
      if HVG_METHOD in ['groups', 'all']:
        gtf_df = pd.read_csv(AF_GTF_DT_F,  sep = '\t')
        num_genes = gtf_df.shape[0]

        NUM_CHUNKS  = (num_genes + HVG_CHUNK_SIZE - 1) // HVG_CHUNK_SIZE
        
        # make a list of chunk names
        CHUNK_NAMES = [f"chunk_{i+1}" for i in range(NUM_CHUNKS)]

  return HVG_METHOD, HVG_SPLIT_VAR, HVG_CHUNK_SIZE, NUM_CHUNKS, GROUP_NAMES, CHUNK_NAMES, N_HVGS, EXCLUDE_AMBIENT_GENES


# define integration parameters
def get_integration_parameters(config): 
  # set some more default values
  INT_CL_METHOD   = 'louvain'
  INT_REDUCTION   = 'harmony'
  INT_N_DIMS      = 50
  INT_THETA       = 0.1
  INT_RES_LS      = [0.1, 0.2, 0.5, 1, 2]
  INT_DBL_RES     = 4
  INT_DBL_CL_PROP = 0.5

  # change defaults if specified
  if ('integration' in config) and (config['integration'] is not None):
    if 'cl_method' in config['integration']:
      INT_CL_METHOD    = config['integration']['cl_method']
      valid_methods = ['leiden', 'louvain']
      assert INT_CL_METHOD in valid_methods, \
        f"Invalid clustering algorithm '{INT_CL_METHOD}'. Must be one of {valid_methods}."
    if 'reduction' in config['integration']:
      INT_REDUCTION    = config['integration']['reduction']
      valid_reductions = ['pca', 'harmony']
      assert INT_REDUCTION in valid_reductions, \
        f"Invalid reduction option '{INT_REDUCTION}'. Must be one of {valid_reductions}."
    if 'int_n_dims' in config['integration']:
      INT_N_DIMS      = config['integration']['int_n_dims']
    if 'int_dbl_res' in config['integration']:
      INT_DBL_RES     = config['integration']['int_dbl_res']
    if 'int_dbl_cl_prop' in config['integration']:
      INT_DBL_CL_PROP = config['integration']['int_dbl_cl_prop']
    if 'int_theta' in config['integration']:
      INT_THETA       = config['integration']['int_theta']
    if 'int_res_ls' in config['integration']:
      INT_RES_LS      = config['integration']['int_res_ls']

  return INT_CL_METHOD, INT_REDUCTION, INT_N_DIMS, INT_THETA, INT_RES_LS, INT_DBL_RES, INT_DBL_CL_PROP


def get_custom_marker_genes_parameters(config, PROJ_DIR, SCPROCESS_DATA_DIR):
  # set defaults
  CUSTOM_MKR_NAMES = ""
  CUSTOM_MKR_PATHS = ""
  
  if ('marker_genes' in config) and (config['marker_genes'] is not None):
    if 'custom_sets' in config["marker_genes"]:
      custom_sets = config["marker_genes"]["custom_sets"]
      mkr_names = []
      mkr_paths = []
      for i, gene_set in enumerate(custom_sets):
        assert "name" in gene_set, \
          f"Entry {i+1} in 'custom_sets' is missing a 'name' field."

        name = gene_set["name"]

        # check for commas in names
        assert "," not in name, \
          f"Custom marker set name '{name}' contains a comma, which is not allowed."

        file_path = gene_set.get("file", os.path.join(SCPROCESS_DATA_DIR, 'marker_genes', f"{name}.csv"))

        if not os.path.isabs(file_path):
          file_path = os.path.join(PROJ_DIR, file_path)

        assert os.path.isfile(file_path), \
          f"File not found for marker set '{name}'"

        assert file_path.endswith(".csv"), \
          f"File for custom marker set '{name}' is not a CSV file"

        # check csv file contents
        mkrs_df = pd.read_csv(file_path)

        req_col = "label"
        opt_cols = ["symbol", "ensembl_id"]
      
        assert req_col in mkrs_df.columns, \
          f"File '{file_path}' is missing the mandatory column 'label'."
      
        assert any(col in mkrs_df.columns for col in opt_cols), \
          f"File '{file_path}' must contain either 'symbol' or 'ensembl_id' column."
          

        # Store validated values
        mkr_names.append(name)
        mkr_paths.append(file_path)

      CUSTOM_MKR_NAMES = ",".join(mkr_names)
      CUSTOM_MKR_PATHS = ",".join(mkr_paths)
  
  return CUSTOM_MKR_NAMES, CUSTOM_MKR_PATHS


# define marker_genes parameters
def get_marker_genes_parameters(config, PROJ_DIR, SCPROCESS_DATA_DIR): 
  # set some more default values
  MKR_SEL_RES     = 0.2
  MKR_GSEA_DIR    = os.path.join(SCPROCESS_DATA_DIR, 'gmt_pathways')
  MKR_MIN_CL_SIZE = 1e2
  MKR_MIN_CELLS   = 10
  MKR_NOT_OK_RE   = "(lincRNA|lncRNA|pseudogene|antisense)"
  MKR_MIN_CPM_MKR = 50
  MKR_MIN_CPM_GO  = 1
  MKR_MAX_ZERO_P  = 0.5
  MKR_GSEA_CUT    = 0.1


  # change defaults if specified
  if ('marker_genes' in config) and (config['marker_genes'] is not None):
    if 'mkr_sel_res' in config['marker_genes']:
      MKR_SEL_RES     = config['marker_genes']['mkr_sel_res']
    if 'mkr_min_cl_size' in config['marker_genes']:
      MKR_MIN_CL_SIZE = config['marker_genes']['mkr_min_cl_size']
    if 'mkr_min_cells' in config['marker_genes']:
      MKR_MIN_CELLS   = config['marker_genes']['mkr_min_cells']
    if 'mkr_not_ok_re' in config['marker_genes']:
      MKR_NOT_OK_RE   = config['marker_genes']['mkr_not_ok_re']
    if 'mkr_min_cpm_mkr' in config['marker_genes']:
      MKR_MIN_CPM_MKR = config['marker_genes']['mkr_min_cpm_mkr']
    if 'mkr_min_cpm_go' in config['marker_genes']:
      MKR_MIN_CPM_GO  = config['marker_genes']['mkr_min_cpm_go']
    if 'mkr_max_zero_p' in config['marker_genes']:
      MKR_MAX_ZERO_P  = config['marker_genes']['mkr_max_zero_p']
    if 'mkr_gsea_cut' in config['marker_genes']:
      MKR_GSEA_CUT    = config['marker_genes']['mkr_gsea_cut']

  # get custom marker files
  CUSTOM_MKR_NAMES, CUSTOM_MKR_PATHS = get_custom_marker_genes_parameters(config, PROJ_DIR, SCPROCESS_DATA_DIR)

  return MKR_SEL_RES, MKR_GSEA_DIR, MKR_MIN_CL_SIZE, MKR_MIN_CELLS, MKR_NOT_OK_RE, MKR_MIN_CPM_MKR, MKR_MIN_CPM_GO, MKR_MAX_ZERO_P, MKR_GSEA_CUT, CUSTOM_MKR_NAMES, CUSTOM_MKR_PATHS


# define marker_genes parameters
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
    valid_boosts = ['human_cns', 'mouse_cns']
    assert LBL_TISSUE in valid_boosts, \
      f"value {LBL_TISSUE} for 'lbl_tissue' parameter is not valid"
 
    # pick labeller
    xgb_dir  = os.path.join(SCPROCESS_DATA_DIR, 'xgboost')
    assert os.path.isdir(xgb_dir)
    
    if LBL_TISSUE == 'human_cns':    
      LBL_XGB_F       = os.path.join(xgb_dir, "Siletti_Macnair-2025-07-23/xgboost_obj_hvgs_Siletti_Macnair_2025-07-23.rds")
      LBL_XGB_CLS_F   = os.path.join(xgb_dir, "Siletti_Macnair-2025-07-23/allowed_cls_Siletti_Macnair_2025-07-23.csv")
    else: 
      raise ValueError(f"{LBL_TISSUE} classifier is unfortunatelly not available yet")

    # check these are ok
    assert os.path.isfile(LBL_XGB_F), \
      f"file {LBL_XGB_F} doesn't exist; consider (re)runnning scprocess setup"
    assert os.path.isfile(LBL_XGB_CLS_F), \
      f"file {LBL_XGB_CLS_F} doesn't exist; consider (re)runnning scprocess setup"
 
  return LBL_XGB_F, LBL_XGB_CLS_F, LBL_GENE_VAR, LBL_SEL_RES_CL, LBL_MIN_PRED, LBL_MIN_CL_PROP, LBL_MIN_CL_SIZE, LBL_TISSUE


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


# define marker_genes parameters
def get_zoom_parameters(config, LBL_TISSUE, LBL_XGB_CLS_F, METADATA_F, AF_GTF_DT_F,
   PROJ_DIR, SHORT_TAG, FULL_TAG, DATE_STAMP, SCPROCESS_DATA_DIR):   
  if ('zoom' not in config) or (config['zoom'] is None):
    ZOOM_NAMES       = []
    ZOOM_PARAMS_DICT = []
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
    
    ZOOM_PARAMS_DICT = dict(zip(
      ZOOM_NAMES,
      [ _get_one_zoom_parameters(zoom_f, LBL_TISSUE, LBL_XGB_CLS_F, METADATA_F, AF_GTF_DT_F,
        PROJ_DIR, SHORT_TAG, FULL_TAG, DATE_STAMP, SCPROCESS_DATA_DIR
      ) for zoom_f in ZOOM_YAMLS ]
      ))

  return ZOOM_NAMES, ZOOM_PARAMS_DICT

  
# get rule resource parameters
def get_resource_parameters(config):
  # set default values
  RETRIES                         = 3
  MB_RUN_MAPPING                  = 8192
  MB_SAVE_ALEVIN_TO_H5            = 8192
  MB_RUN_AMBIENT                  = 8192
  MB_GET_BARCODE_QC_METRICS       = 8192
  MB_MAKE_HTO_SCE_OBJECTS         = 8192
  MB_COMBINE_SCDBLFINDER_OUTPUTS  = 8192
  MB_RUN_QC                       = 8192
  MB_RUN_HVGS                     = 8192
  MB_PB_MAKE_PBS                  = 8192
  MB_PB_CALC_EMPTY_GENES          = 8192
  MB_RUN_INTEGRATION              = 8192
  MB_MAKE_CLEAN_SCES              = 8192
  MB_RUN_MARKER_GENES             = 16384
  MB_RENDER_HTMLS                 = 8192
  MB_LABEL_CELLTYPES              = 8192
  MB_ZOOM_RUN_ZOOM                = 8192
  MB_MAKE_SUBSET_SCES             = 8192
  MB_ZOOM_RENDER_TEMPLATE_RMD     = 4096

  # change defaults if specified
  if ('resources' in config) and (config['resources'] is not None):
    if 'retries' in config['resources']:
      RETRIES                         = config['resources']['retries']
    if 'mb_run_mapping' in config['resources']:
      MB_RUN_MAPPING                  = config['resources']['mb_run_mapping']
    if 'mb_save_alevin_to_h5' in config['resources']:
      MB_SAVE_ALEVIN_TO_H5            = config['resources']['mb_save_alevin_to_h5']
    if 'mb_run_ambient' in config['resources']:
      MB_RUN_AMBIENT               = config['resources']['mb_run_ambient']
    if 'mb_get_barcode_qc_metrics' in config['resources']:
      MB_GET_BARCODE_QC_METRICS       = config['resources']['mb_get_barcode_qc_metrics']
    if 'mb_run_qc' in config['resources']:
      MB_RUN_QC                       = config['resources']['mb_run_qc']
    if 'mb_run_hvgs' in config['resources']:
      MB_RUN_HVGS                     = config['resources']['mb_run_hvgs']
    if 'mb_run_integration' in config['resources']:
      MB_RUN_INTEGRATION              = config['resources']['mb_run_integration']
    if 'mb_make_clean_sces' in config['resources']:
      MB_MAKE_CLEAN_SCES              = config['resources']['mb_make_clean_sces']
    if 'mb_run_marker_genes' in config['resources']:
      MB_RUN_MARKER_GENES             = config['resources']['mb_run_marker_genes']
    if 'mb_render_htmls'     in config['resources']:
      MB_RENDER_HTMLS                = config['resources']['mb_render_htmls']
    if 'mb_label_celltypes' in config['resources']:
      MB_LABEL_CELLTYPES          = config['resources']['mb_label_celltypes']
    if 'mb_pb_make_pbs' in config['resources']:
      MB_PB_MAKE_PBS                  = config['resources']['mb_pb_make_pbs']
    if 'mb_pb_calc_empty_genes' in config['resources']:
      MB_PB_CALC_EMPTY_GENES          = config['resources']['mb_pb_calc_empty_genes']
    if 'mb_zoom_run_zoom' in config['resources']:
      MB_ZOOM_RUN_ZOOM                = config['resources']['mb_zoom_run_zoom']
    if 'mb_zoom_render_template_rmd' in config['resources']:
      MB_ZOOM_RENDER_TEMPLATE_RMD     = config['resources']['mb_zoom_render_template_rmd']
    if 'mb_make_subset_sces'         in config['resources']:
      MB_MAKE_SUBSET_SCES             = config['resources']['mb_make_subset_sces']
    if 'mb_make_hto_sce_objects' in config['resources']:
      MB_MAKE_HTO_SCE_OBJECTS         = config['resources']['mb_make_hto_sce_objects']

  return RETRIES, MB_RUN_MAPPING, MB_SAVE_ALEVIN_TO_H5, \
    MB_RUN_AMBIENT, MB_GET_BARCODE_QC_METRICS, \
    MB_RUN_QC, MB_RUN_HVGS, \
    MB_RUN_INTEGRATION, MB_MAKE_CLEAN_SCES, \
    MB_RUN_MARKER_GENES, MB_RENDER_HTMLS, \
    MB_LABEL_CELLTYPES, \
    MB_PB_MAKE_PBS, MB_PB_CALC_EMPTY_GENES, MB_MAKE_HTO_SCE_OBJECTS, \
    MB_ZOOM_RUN_ZOOM, MB_ZOOM_RENDER_TEMPLATE_RMD, MB_MAKE_SUBSET_SCES

