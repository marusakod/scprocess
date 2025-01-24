# load modules
import warnings
import yaml
import pandas as pd
import os
import re
import glob
import datetime
import subprocess


def __get_cl_ls(config, mito_str, scprocess_data_dir):
  # get parameters
  PROJ_DIR, FASTQ_DIR, SHORT_TAG, FULL_TAG, _, _, _, _, _, _, DATE_STAMP, _, _ = \
    get_project_parameters(config, scprocess_data_dir)
  _, _, _, _, _, _, _, INT_SEL_RES = \
    get_integration_parameters(config, mito_str)

  # specify harmony outputs
  int_dir     = f"{PROJ_DIR}/output/{SHORT_TAG}_integration"
  hmny_f      = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'

  # get list of clusters
  hmny_dt     = pd.read_csv(hmny_f)
  cl_col      = f"RNA_snn_res.{INT_SEL_RES}"
  cl_ls       = list(hmny_dt[cl_col].unique())
  cl_ls       = [cl for cl in cl_ls if str(cl) != "nan"]
  cl_ls       = sorted(cl_ls)

  return(cl_ls)


def __get_one_zoom_parameters(config, zoom_name, cl_ls):
  # unpack
  this_zoom   = config['zoom'][ zoom_name ]

  # check parameters
  p_names     = this_zoom.keys()
  req_names   = [ "sel_cls", "n_hvgs", "n_dims", "zoom_res", "min_n_sample", "min_n_cl", "n_train" ]
  missing_ns  = [ p_name for p_name in req_names if p_name not in p_names ]
  assert len(missing_ns) == 0, f"the following parameters are missing from zoom {zoom_name}:\n{'_'.join(missing_ns)}"

  # check one by one
  assert set( this_zoom[ "sel_cls" ] ).issubset( cl_ls )
  assert this_zoom[ "zoom_res" ] > 0
  assert isinstance(this_zoom[ "n_hvgs" ], int)
  assert isinstance(this_zoom[ "n_dims" ], int)
  assert isinstance(this_zoom[ "min_n_cl" ], int)
  assert isinstance(this_zoom[ "n_train" ], int)

  # add name to dictionary
  this_zoom[ 'zoom_name' ] = zoom_name

  return this_zoom


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
  DEMUX_TYPE     = None
  HTO_FASTQ_DIR  = None
  FEATURE_REF    = None
  DEMUX_F        = ""
  BATCH_VAR      = "sample_id"
  POOL_IDS       = ""
  EXC_POOLS      = ""
  
  if ('multiplexing' in config) and (config['multiplexing'] is not None):
    POOL_IDS = list(set(sample_metadata["pool_id"].tolist()))
    DEMUX_TYPE = config['multiplexing']['demux_type']

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
      assert set(demux_dt['sample_id']) == set(sample_metadata['sample_id']), \
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
        for p in EXC_POOLS:
          if p not in POOL_IDS:
            warnings.warn(f"sample {p} specified in exclude_pools but not in sample_metadata file", UserWarning)
        to_keep  = set(POOL_IDS) - set(EXC_POOLS)
        POOL_IDS = [p for p in POOL_IDS if p in to_keep]
      
  return DEMUX_TYPE, HTO_FASTQ_DIR, FEATURE_REF, DEMUX_F, BATCH_VAR, EXC_POOLS, POOL_IDS

      


# get list of samples
def get_project_parameters(config, scprocess_data_dir):
  # check expected variables are in the config file
  for v in ["proj_dir", "fastq_dir", "short_tag", "full_tag", "date_stamp", "your_name", "affiliation", "sample_metadata", "species"]:
    assert v in config, f"{v} not in config file"

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
  setup_params_f  = os.path.join(scprocess_data_dir, 'setup_parameters.csv')

    # from setup_parameters.csv get valid values for species
  setup_params= pd.read_csv(setup_params_f)
  valid_species = setup_params['genome_name'].tolist()

  assert SPECIES in valid_species, \
   f"species {SPECIES} not defined"

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
  DEMUX_TYPE, HTO_FASTQ_DIR, FEATURE_REF, DEMUX_F, BATCH_VAR, EXC_POOLS, POOL_IDS = \
  get_multiplexing_parameters(config, PROJ_DIR, samples_df)

  # define sample variable for alevin, ambient, doublets
  if DEMUX_TYPE is not None:
    sample_var = "pool_id"
  else:
    sample_var = "sample_id"

  # remove some samples
  EXC_SAMPLES   = None
  if ('exclude' in config) and (config['exclude'] is not None):
    if ('sample_id' in config['exclude']) and (config['exclude']['sample_id'] is not None):
      EXC_SAMPLES = config["exclude_samples"]
      for s in EXC_SAMPLES:
        if s not in SAMPLES:
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
      if DEMUX_TYPE is not None:
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
      uniq_vals = len(set(samples_df[var].tolist()))
      assert uniq_vals <= 10, \
        f"{var} variable has more than 10 unique values"
  

  return PROJ_DIR, FASTQ_DIR, SHORT_TAG, FULL_TAG, YOUR_NAME, AFFILIATION, METADATA_F, \
    METADATA_VARS, EXC_SAMPLES, SAMPLES, DATE_STAMP, CUSTOM_SAMPLE_PARAMS_F, SPECIES, \
    DEMUX_TYPE, HTO_FASTQ_DIR, FEATURE_REF, DEMUX_F, BATCH_VAR, EXC_POOLS, POOL_IDS


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
  setup_params_f  = os.path.join(scprocess_data_dir, 'setup_parameters.csv')

    # from setup_parameters.csv get valid values for species
  setup_params= pd.read_csv(setup_params_f)
     
  # get mito strings from setup params
  AF_MITO_STR = setup_params.loc[setup_params['genome_name'] == SPECIES, 'mito_str'].values[0]

  # get af index directory and check if exists
  AF_HOME_DIR = os.path.join(scprocess_data_dir, 'alevin_fry_home') # check if this exists in scprocess script
  AF_INDEX_DIR = os.path.join(AF_HOME_DIR, SPECIES)
  assert os.path.isdir(AF_INDEX_DIR), \
    f"alevin index for {SPECIES} doesn't exist"
  
  # get gtf txt file, check that exists
  AF_GTF_DT_F = setup_params.loc[setup_params['genome_name'] == SPECIES, 'gtf_txt_f'].values[0]

  return AF_MITO_STR, AF_HOME_DIR, AF_INDEX_DIR, AF_GTF_DT_F, CHEMISTRY


# ambient
def get_ambient_parameters(config):
  # set default values
  AMBIENT_METHOD                = 'decontx'
  CELLBENDER_VERSION            = 'v0.3.0'
  CELLBENDER_PROP_MAX_KEPT      = 0.9
  FORCE_EXPECTED_CELLS          = None
  FORCE_TOTAL_DROPLETS_INCLUDED = None
  FORCE_LOW_COUNT_THRESHOLD     = None
  CELLBENDER_LEARNING_RATE      = 1e-4
  CELL_CALLS_METHOD             = 'barcodeRanks'

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
      FORCE_EXPECTED_CELLS          = config['ambient']['force_expected_cells']
    if 'cb_force_total_droplets_included' in config['ambient']:
      FORCE_TOTAL_DROPLETS_INCLUDED = config['ambient']['force_total_droplets_included']
    if 'cb_force_low_count_threshold' in config['ambient']:
      FORCE_LOW_COUNT_THRESHOLD     = config['ambient']['force_low_count_threshold']
    if 'cb_force_learning_rate' in config['ambient']:
      CELLBENDER_LEARNING_RATE      = config['ambient']['cb_force_learning_rate']

  # get cellbender image (maybe skip this if cellbender is not selected?)
  if CELLBENDER_VERSION == 'v0.3.0':
    CELLBENDER_IMAGE              = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.0'
  elif CELLBENDER_VERSION == 'v0.2.0':
    CELLBENDER_IMAGE              = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.2.0'
  else:
    raise ValueError(f"selected cellbender version {CELLBENDER_VERSION} not supported")

  # some checks on custom parameters for cellbender
      
  return CELLBENDER_IMAGE, CELLBENDER_PROP_MAX_KEPT, AMBIENT_METHOD, CELL_CALLS_METHOD, \
    FORCE_EXPECTED_CELLS, FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE




def get_make_sce_parameters(config):  
  # set default values
  SCE_BENDER_PROB = 0.5

  # change defaults if specified
  if ('make_sce' in config) and (config['make_sce'] is not None):
    if 'sce_bender_prob' in config['make_sce']:
      SCE_BENDER_PROB = config['make_sce']['sce_bender_prob']

  return SCE_BENDER_PROB


# define doublet_id parameters
def get_doublet_id_parameters(config):  
  # set default values
  DBL_MIN_FEATS = 100

  # change defaults if specified
  if ('doublet_id' in config) and (config['doublet_id'] is not None):
    if 'dbl_min_feats' in config['doublet_id']:
      DBL_MIN_FEATS = config['doublet_id']['dbl_min_feats']

  return DBL_MIN_FEATS

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
  QC_FILTER_BENDER    = False

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
    if 'qc_filter_bender'   in config['qc']:
      QC_FILTER_BENDER    = config['qc']['qc_filter_bender']

  # make sure they're consistent
  QC_HARD_MIN_COUNTS  = min(QC_HARD_MIN_COUNTS, QC_MIN_COUNTS)
  QC_HARD_MIN_FEATS   = min(QC_HARD_MIN_FEATS, QC_MIN_FEATS)
  QC_HARD_MAX_MITO    = max(QC_HARD_MAX_MITO, QC_MAX_MITO)

  return QC_HARD_MIN_COUNTS, QC_HARD_MIN_FEATS, QC_HARD_MAX_MITO, QC_MIN_COUNTS, QC_MIN_FEATS, \
    QC_MIN_MITO, QC_MAX_MITO, QC_MIN_SPLICE, QC_MAX_SPLICE, QC_MIN_CELLS, QC_FILTER_BENDER


# define integration parameters
def get_integration_parameters(config, mito_str): 
  # set default values
  INT_EXC_REGEX   = mito_str


  # set some more default values
  INT_N_DIMS      = 50
  INT_N_HVGS      = 2000
  INT_DBL_RES     = 4
  INT_DBL_CL_PROP = 0.5
  INT_THETA       = 0.1
  INT_RES_LS      = [0.1, 0.2, 0.5, 1, 2]
  INT_SEL_RES     = 0.2

  # change defaults if specified
  if ('integration' in config) and (config['integration'] is not None):
    if 'int_n_hvgs' in config['integration']:
      INT_N_HVGS      = config['integration']['int_n_hvgs']
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
    if 'int_sel_res' in config['integration']:
      INT_SEL_RES     = config['integration']['int_sel_res']

  return INT_EXC_REGEX, INT_N_HVGS, INT_N_DIMS, INT_DBL_RES, INT_DBL_CL_PROP, INT_THETA, INT_RES_LS, INT_SEL_RES


# define marker_genes parameters
def get_marker_genes_parameters(config, SPECIES, SCPROCESS_DATA_DIR): 
  # set some more default values
  MKR_GSEA_DIR    = os.path.join(SCPROCESS_DATA_DIR, 'gmt_pathways')
  MKR_MIN_CL_SIZE = 1e2
  MKR_MIN_CELLS   = 10
  MKR_NOT_OK_RE   = "(lincRNA|lncRNA|pseudogene|antisense)"
  MKR_MIN_CPM_MKR = 50
  MKR_MIN_CPM_GO  = 1
  MKR_MAX_ZERO_P  = 0.5
  MKR_GSEA_CUT    = 0.1

  # specify canonical marker file; this should depend on both species and organ; could be added to config as input param
  if SPECIES in ["human_2020", "human_2024"]:
    MKR_CANON_F     = os.path.join(SCPROCESS_DATA_DIR,"marker_genes/canonical_brain_celltype_markers_human_2023-10-17.txt")
  elif SPECIES in ["mouse_2020", "mouse_2024"]:
    MKR_CANON_F     = os.path.join(SCPROCESS_DATA_DIR,"marker_genes/canonical_brain_celltype_markers_mouse_2023-10-17.txt")
  else:
    MKR_CANON_F = ""

  # change defaults if specified
  if ('marker_genes' in config) and (config['marker_genes'] is not None):
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

  return MKR_GSEA_DIR, MKR_MIN_CL_SIZE, MKR_MIN_CELLS, MKR_NOT_OK_RE, MKR_MIN_CPM_MKR, MKR_MIN_CPM_GO, MKR_MAX_ZERO_P, MKR_GSEA_CUT, MKR_CANON_F


# define marker_genes parameters
def get_label_celltypes_parameters(config, SPECIES, SCPROCESS_DATA_DIR): 
  # set some more default values
  CUSTOM_LABELS_F = ""
  LBL_XGB_F       = None
  LBL_XGB_CLS_F   = None
  LBL_TISSUE      = ""
  LBL_GENE_VAR    = "gene_id"
  LBL_SEL_RES_CL  = "RNA_snn_res.2"
  LBL_MIN_PRED    = 0.8
  LBL_MIN_CL_PROP = 0.5
  LBL_MIN_CL_SIZE = 100
  LBL_SCE_SUBSETS = None

  # change defaults if specified
  if ('label_celltypes' in config) and (config['label_celltypes'] is not None):
    if 'custom_labels' in config['label_celltypes']:
      CUSTOM_LABELS_F = config['label_celltypes']['custom_labels']
    if 'lbl_tissue' in config['label_celltypes']:
      LBL_TISSUE      = config['label_celltypes']['lbl_tissue']
    if 'lbl_sel_res_cl' in config['label_celltypes']:
      LBL_SEL_RES_CL  = config['label_celltypes']['lbl_sel_res_cl']
    if 'lbl_min_pred' in config['label_celltypes']:
      LBL_MIN_PRED    = config['label_celltypes']['lbl_min_pred']
    if 'lbl_min_cl_prop' in config['label_celltypes']:
      LBL_MIN_CL_PROP = config['label_celltypes']['lbl_min_cl_prop']
    if 'lbl_min_cl_size' in config['label_celltypes']:
      LBL_MIN_CL_SIZE = config['label_celltypes']['lbl_min_cl_size']
    if 'sce_subsets' in config['label_celltypes']:
      # get sce subset definitions
      LBL_SCE_SUBSETS = config['label_celltypes']['sce_subsets']
      assert len(LBL_SCE_SUBSETS) == len(set(LBL_SCE_SUBSETS)), "sce_subsets contains duplicate entries"

    # pick either a classifier or provide a file with custom labels
    assert any(lbl != "" for lbl in [CUSTOM_LABELS_F, LBL_TISSUE]), \
     "'custom_labels' or 'lbl_tissue' parameter have to be defined in the config file"
    
    # if both are specified, use only custom label
    if all(lbl != "" for lbl in [CUSTOM_LABELS_F, LBL_TISSUE]):
      print("both 'custom_labels' and 'lbl_tissue' parameters are defined; 'custom_labels' will be used")
      labels = 'custom'
    elif CUSTOM_LABELS_F != "":
      labels = 'custom'
    else:
      labels = 'xgboost'

    
    if labels == 'xgboost':
      # check that classifier name is valid
      valid_boosts = ['human_cns', 'mouse_cns', 'human_pmbc', 'mouse_pbmc']
      assert LBL_TISSUE in valid_boosts, \
        f"value {LBL_TISSUE} for 'lbl_tissue' parameter is not valid"
 
      # pick labeller
      xgb_dir  = os.path.join(SCPROCESS_DATA_DIR, 'xgboost')
      assert os.path.isdir(xgb_dir)
    
      if LBL_TISSUE == 'human_cns':    
        LBL_XGB_F       = os.path.join(xgb_dir, "Siletti_Macnair-2024-03-11/xgboost_obj_hvgs_Siletti_Macnair_2024-03-11.rds")
        LBL_XGB_CLS_F   = os.path.join(xgb_dir, "Siletti_Macnair-2024-03-11/allowed_cls_Siletti_Macnair_2024-03-11.csv")
      else: 
        raise ValueError(f"{LBL_TISSUE} classifier is unfortunatelly not available yet")

      # check these are ok
      assert os.path.isfile(LBL_XGB_F)
      assert os.path.isfile(LBL_XGB_CLS_F)
      
    else:
      # check custom labels file
      assert os.path.isfile(CUSTOM_LABELS_F), \
       f"{CUSTOM_LABELS_F} file doesn't exist"
      
      custom_labels = pd.read_csv(CUSTOM_LABELS_F)
      cols = custom_labels.columns
      valid_cols = ['cell_id', 'label']
      assert all(c in valid_cols for c in cols), \
        f"Column names in {CUSTOM_LABELS_F} are incorrect."
 
  return LBL_XGB_F, LBL_XGB_CLS_F, LBL_GENE_VAR, LBL_SEL_RES_CL, LBL_MIN_PRED, LBL_MIN_CL_PROP, LBL_MIN_CL_SIZE, LBL_SCE_SUBSETS, LBL_TISSUE, CUSTOM_LABELS_F


# define metacells parameters
def get_metacells_parameters(config): 
  # set some more default values
  META_SUBSETS    = []
  META_MAX_CELLS  = [100]

  # change defaults if specified
  if ('metacells' in config) and (config['metacells'] is not None):
    if 'celltypes' in config['metacells']:
      META_SUBSETS  = config['metacells']['celltypes']
    if 'max_cells' in config['metacells']:
      META_MAX_CELLS  = config['metacells']['max_cells']
 
  return META_SUBSETS, META_MAX_CELLS


# define pseudobulk parameters
def get_pb_empties_parameters(config): 
  # set some more default values
  PB_SUBSETS          = []
  PB_DO_ALL           = False

  # change defaults if specified
  if ('pb_empties' in config) and (config['pb_empties'] is not None):
    if 'subsets' in config['pb_empties']:
      PB_SUBSETS          = config['pb_empties']['subsets']


  # parse out empties and all
  if "all" in PB_SUBSETS:
    PB_DO_ALL     = True
    PB_SUBSETS    = [x for x in PB_SUBSETS if x != "all"]

  return PB_SUBSETS, PB_DO_ALL


# define marker_genes parameters
def get_zoom_parameters(config, MITO_STR, scprocess_data_dir): 
  if ('zoom' not in config) or (config['zoom'] is None):
    ZOOM_NAMES    = []
    ZOOM_SPEC_LS  = []
  else:
    cl_ls         = __get_cl_ls(config, MITO_STR, scprocess_data_dir)
    ZOOM_NAMES    = list(config['zoom'].keys())
    ZOOM_SPEC_LS  = dict(zip(
      ZOOM_NAMES,
      [ __get_one_zoom_parameters(config, zoom_name, cl_ls) for zoom_name in config['zoom'] ]
      ))

  return ZOOM_NAMES, ZOOM_SPEC_LS



  
# get rule resource parameters
def get_resource_parameters(config):
  # set default values
  RETRIES                         = 1
  MB_RUN_ALEVIN_FRY               = 8192
  MB_SAVE_ALEVIN_TO_H5            = 8192
  MB_RUN_AMBIENT                  = 8192
  MB_GET_BARCODE_QC_METRICS       = 4096
  MB_RUN_SCDBLFINDER              = 4096
  MB_COMBINE_SCDBLFINDER_OUTPUTS  = 8192
  MB_RUN_QC                       = 8192
  MB_MAKE_SCE_OBJECT              = 8192
  MB_RUN_HARMONY                  = 8192
  MB_RUN_MARKER_GENES             = 8192
  MB_HTML_MARKER_GENES            = 8192
  MB_LBL_LABEL_CELLTYPES          = 8192
  MB_LBL_SAVE_SUBSET_SCES         = 8192
  MB_LBL_RENDER_TEMPLATE_RMD      = 4096
  MB_META_SAVE_METACELLS          = 8192
  MB_PB_MAKE_PBS                  = 8192
  MB_PB_CALC_EMPTY_GENES          = 8192
  MB_ZOOM_RUN_ZOOM                = 8192
  MB_ZOOM_RENDER_TEMPLATE_RMD     = 4096

  # change defaults if specified
  if ('resources' in config) and (config['resources'] is not None):
    if 'retries' in config['resources']:
      RETRIES                         = config['resources']['retries']
    if 'mb_run_alevin_fry' in config['resources']:
      MB_RUN_ALEVIN_FRY               = config['resources']['mb_run_alevin_fry']
    if 'mb_save_alevin_to_h5' in config['resources']:
      MB_SAVE_ALEVIN_TO_H5            = config['resources']['mb_save_alevin_to_h5']
    if 'mb_run_ambient' in config['resources']:
      MB_RUN_AMBIENT               = config['resources']['mb_run_ambient']
    if 'mb_get_barode_qc_metrics' in config['resources']:
      MB_GET_BARCODE_QC_METRICS    = config['resources']['mb_get_barcode_qc_metrics']
    if 'mb_run_scdblfinder' in config['resources']:
      MB_RUN_SCDBLFINDER              = config['resources']['mb_run_scdblfinder']
    if 'mb_combine_scdblfinder_outputs' in config['resources']:
      MB_COMBINE_SCDBLFINDER_OUTPUTS  = config['resources']['mb_combine_scdblfinder_outputs']
    if 'mb_run_qc' in config['resources']:
      MB_RUN_QC                       = config['resources']['mb_run_qc']
    if 'mb_make_sce_object' in config['resources']:
      MB_MAKE_SCE_OBJECT              = config['resources']['mb_make_sce_object']
    if 'mb_run_harmony' in config['resources']:
      MB_RUN_HARMONY                  = config['resources']['mb_run_harmony']
    if 'mb_run_marker_genes' in config['resources']:
      MB_RUN_MARKER_GENES             = config['resources']['mb_run_marker_genes']
    if 'mb_html_marker_genes' in config['resources']:
      MB_HTML_MARKER_GENES            = config['resources']['mb_html_marker_genes']
    if 'mb_lbl_label_celltypes' in config['resources']:
      MB_LBL_LABEL_CELLTYPES          = config['resources']['mb_lbl_label_celltypes']
    if 'mb_lbl_save_subset_sces' in config['resources']:
      MB_LBL_SAVE_SUBSET_SCES         = config['resources']['mb_lbl_save_subset_sces']
    if 'mb_lbl_render_template_rmd' in config['resources']:
      MB_LBL_RENDER_TEMPLATE_RMD      = config['resources']['mb_lbl_render_template_rmd']
    if 'mb_meta_save_metacells' in config['resources']:
      MB_META_SAVE_METACELLS          = config['resources']['mb_meta_save_metacells']
    if 'mb_pb_make_pbs' in config['resources']:
      MB_PB_MAKE_PBS                  = config['resources']['mb_pb_make_pbs']
    if 'mb_pb_calc_empty_genes' in config['resources']:
      MB_PB_CALC_EMPTY_GENES          = config['resources']['mb_pb_calc_empty_genes']
    if 'mb_zoom_run_zoom' in config['resources']:
      MB_ZOOM_RUN_ZOOM                = config['resources']['mb_zoom_run_zoom']
    if 'mb_zoom_render_template_rmd' in config['resources']:
      MB_ZOOM_RENDER_TEMPLATE_RMD     = config['resources']['mb_zoom_render_template_rmd']

  return RETRIES, MB_RUN_ALEVIN_FRY, MB_SAVE_ALEVIN_TO_H5, \
    MB_RUN_AMBIENT, MB_GET_BARCODE_QC_METRICS, \
    MB_RUN_SCDBLFINDER, MB_COMBINE_SCDBLFINDER_OUTPUTS, \
    MB_RUN_QC, \
    MB_MAKE_SCE_OBJECT, \
    MB_RUN_HARMONY, \
    MB_RUN_MARKER_GENES, MB_HTML_MARKER_GENES, \
    MB_LBL_LABEL_CELLTYPES, MB_LBL_SAVE_SUBSET_SCES, MB_LBL_RENDER_TEMPLATE_RMD, \
    MB_META_SAVE_METACELLS, \
    MB_PB_MAKE_PBS, MB_PB_CALC_EMPTY_GENES, \
    MB_ZOOM_RUN_ZOOM, MB_ZOOM_RENDER_TEMPLATE_RMD


