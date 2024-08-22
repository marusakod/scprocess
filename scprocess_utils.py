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


# ambient
def get_ambient_parameters(config):
  # set default values
  AMBIENT_METHOD                = 'cellbender'
  CELLBENDER_VERSION            = 'v0.3.0'
  CELLBENDER_PROP_MAX_KEPT      = 0.9
  CUSTOM_PARAMS_F               = None
  FORCE_EXPECTED_CELLS          = None
  FORCE_TOTAL_DROPLETS_INCLUDED = None
  FORCE_LOW_COUNT_THRESHOLD     = None
  CELLBENDER_LEARNING_RATE      = 1e-4
  CELL_CALLS_METHOD             = 'barcodeRanks'

  # change defaults if specified
  if ('ambient' in config) and (config['ambient'] is not None):
    if 'ambient_method' in config['ambient']:
      AMBIENT_METHOD                 = config['ambient']['ambient_method']
    if 'cellbender_version' in config['ambient']:
      CELLBENDER_VERSION            = config['ambient']['cellbender_version']
    if 'cb_max_prop_kept' in config['ambient']:
      CELLBENDER_PROP_MAX_KEPT      = config['ambient']['cb_max_prop_kept']
    if 'custom_params' in config['ambient']:
      CUSTOM_PARAMS_F    = config['ambient']['custom_params']
    if 'force_expected_cells' in config['ambient']:
      FORCE_EXPECTED_CELLS          = config['ambient']['force_expected_cells']
    if 'force_total_droplets_included' in config['ambient']:
      FORCE_TOTAL_DROPLETS_INCLUDED = config['ambient']['force_total_droplets_included']
    if 'force_low_count_threshold' in config['ambient']:
      FORCE_LOW_COUNT_THRESHOLD     = config['ambient']['force_low_count_threshold']
    if 'learning_rate' in config['ambient']:
      CELLBENDER_LEARNING_RATE      = config['ambient']['learning_rate']

  # get cellbender image (maybe skip this if cellbender is not selected?)
  if CELLBENDER_VERSION == 'v0.3.0':
    CELLBENDER_IMAGE              = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.0'
  elif CELLBENDER_VERSION == 'v0.2.0':
    CELLBENDER_IMAGE              = 'docker://us.gcr.io/broad-dsde-methods/cellbender:0.2.0'
  else:
    raise ValueError(f"selected cellbender version {CELLBENDER_VERSION} not supported")

  # some checks on custom parameters for cellbender
  if CUSTOM_PARAMS_F is not None: 
    assert os.path.exists(CUSTOM_PARAMS_F), \
      f"specified path {CUSTOM_PARAMS_F} does not exist"
    
    # get params  
    params_df   = pd.read_csv(CUSTOM_PARAMS_F)
    meta_df     = pd.read_csv(config["sample_metadata"])

    # check if sample ids in custom params file match sample ids in metadata
    for s in params_df['sample_id'].tolist():
        assert s in meta_df['sample_id'].tolist(), \
         f"sample_id {s} in {CUSTOM_PARAMS_F} not present in metadata"
    
    # depending on ambient method and cell calls method check which columns need to be present in the custom file
    # Check the columns based on AMBIENT_METHOD and CELL_CALLS_METHOD
    if AMBIENT_METHOD == 'cellbender':
        expected_columns = ['sample_id', 'total_droplets_included', 'expected_cells', 'low_count_threshold', 'learning_rate', 'empty_start', 'empty_end']
    else:
        if CELL_CALLS_METHOD == 'emptyDrops':
            expected_columns = ['retain', 'empty_start', 'empty_end']
        elif CELL_CALLS_METHOD == 'barcodeRanks':
            expected_columns = ['expected_cells', 'empty_start', 'empty_end']
        else:
            raise ValueError(f"Unsupported CELL_CALLS_METHOD: {CELL_CALLS_METHOD}")

    # Verify that the columns in params_df match the expected columns
    assert all(params_df.columns.values == expected_columns), \
        f"column names in {CUSTOM_PARAMS_F} are not correct. Expected columns: {expected_columns}"
      
      
  return CELLBENDER_IMAGE, CELLBENDER_PROP_MAX_KEPT, AMBIENT_METHOD, CUSTOM_PARAMS_F, CELL_CALLS_METHOD, \
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
  QC_MAX_SPLICE       = 1
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
def get_label_celltypes_parameters(config, SPECIES): 
  # set some more default values
  LBL_XGB_F       = None
  LBL_XGB_CLS_F   = None
  LBL_TISSUE      = "brain_cns"
  LBL_GENE_VAR    = "gene_id"
  LBL_SEL_RES_CL  = "RNA_snn_res.2"
  LBL_MIN_PRED    = 0.8
  LBL_MIN_CL_PROP = 0.5
  LBL_MIN_CL_SIZE = 100
  LBL_SCE_SUBSETS = None

  # change defaults if specified
  if ('label_celltypes' in config) and (config['label_celltypes'] is not None):
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

    # pick labeller
    xgb_dir         = '/projects/site/pred/neurogenomics/resources/scprocess_data/data/xgboost'
    if SPECIES == 'human':
      if LBL_TISSUE == 'brain_cns':    
        # get cluster levels
        # LBL_XGB_F       = "/projects/site/pred/neurogenomics/resources/scprocess_data/data/xgboost/xgboost_obj_hvgs_Bryois_2022_2023-09-18.rds"
        # LBL_XGB_CLS_F   = "/projects/site/pred/neurogenomics/resources/scprocess_data/data/xgboost/allowed_cls_Bryois_2022_2023-09-18.csv"
        # LBL_XGB_F       = os.path.join(xgb_dir, "Macnair_2022_2023-10-18/xgboost_obj_hvgs_Macnair_2022_2023-10-18.rds")
        # LBL_XGB_CLS_F   = os.path.join(xgb_dir, "Macnair_2022_2023-10-18/allowed_cls_Macnair_2022_2023-10-18.csv")
        # LBL_XGB_F       = os.path.join(xgb_dir, "Siletti_Macnair-2024-03-06/xgboost_obj_hvgs_Siletti_Macnair_2024-03-06.rds")
        # LBL_XGB_CLS_F   = os.path.join(xgb_dir, "Siletti_Macnair-2024-03-06/allowed_cls_Siletti_Macnair_2024-03-06.csv")
        LBL_XGB_F       = os.path.join(xgb_dir, "Siletti_Macnair-2024-03-11/xgboost_obj_hvgs_Siletti_Macnair_2024-03-11.rds")
        LBL_XGB_CLS_F   = os.path.join(xgb_dir, "Siletti_Macnair-2024-03-11/allowed_cls_Siletti_Macnair_2024-03-11.csv")
      else:
        raise ValueError(f"Unknown 'lbl_tissue' value in species {SPECIES}: {LBL_TISSUE}")
    elif SPECIES == 'mouse':
      raise ValueError('sorry no mouse classifiers implemented yet, please bug Will...')
    else:
      raise ValueError(f'sorry species {SPECIES} not recognised')

    # check these are ok
    assert os.path.isfile(LBL_XGB_F)
    assert os.path.isfile(LBL_XGB_CLS_F)
 
  return LBL_XGB_F, LBL_XGB_CLS_F, LBL_GENE_VAR, LBL_SEL_RES_CL, LBL_MIN_PRED, LBL_MIN_CL_PROP, LBL_MIN_CL_SIZE, LBL_SCE_SUBSETS


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
  PB_CUSTOM_EMPTIES_F = None
  PB_SUBSETS          = []
  PB_DO_ALL           = False

  # change defaults if specified
  if ('pb_empties' in config) and (config['pb_empties'] is not None):
    if 'custom_empties_f' in config['pb_empties']:
      PB_CUSTOM_EMPTIES_F = config['pb_empties']['custom_empties_f']
    if 'subsets' in config['pb_empties']:
      PB_SUBSETS          = config['pb_empties']['subsets']


  # parse out empties and all
  if "all" in PB_SUBSETS:
    PB_DO_ALL     = True
    PB_SUBSETS    = [x for x in PB_SUBSETS if x != "all"]

  # check that files exist
  if PB_CUSTOM_EMPTIES_F is not None:
    # check whether custom empties file exists
    if not os.path.isfile(PB_CUSTOM_EMPTIES_F):
      PROJ_DIR            = config["proj_dir"]
      PB_CUSTOM_EMPTIES_F = os.path.join(PROJ_DIR, PB_CUSTOM_EMPTIES_F)

    # does this exist?
    assert os.path.isfile(PB_CUSTOM_EMPTIES_F), f"custom empties file {PB_CUSTOM_EMPTIES_F} does not exist"

  return PB_CUSTOM_EMPTIES_F, PB_SUBSETS, PB_DO_ALL


# define marker_genes parameters
def get_zoom_parameters(config): 
  if ('zoom' not in config) or (config['zoom'] is None):
    ZOOM_NAMES    = []
    ZOOM_SPEC_LS  = []
  else:
    cl_ls         = __get_cl_ls(config)
    ZOOM_NAMES    = list(config['zoom'].keys())
    ZOOM_SPEC_LS  = dict(zip(
      ZOOM_NAMES,
      [ __get_one_zoom_parameters(config, zoom_name, cl_ls) for zoom_name in config['zoom'] ]
      ))

  return ZOOM_NAMES, ZOOM_SPEC_LS


# get rule resource parameters
def get_resource_parameters(config):
  # set default values
  MB_RUN_ALEVIN_FRY               = 8192
  MB_SAVE_ALEVIN_TO_H5            = 8192
  MB_RUN_CELLBENDER               = 32768
  MB_GET_CELLBENDER_QC_METRICS    = 4096
  MB_RUN_SCDBLFINDER              = 4096
  MB_COMBINE_SCDBLFINDER_OUTPUTS  = 8192
  MB_RUN_QC                       = 16384
  MB_MAKE_SCE_OBJECT              = 16384
  MB_RUN_HARMONY                  = 16384
  MB_RUN_MARKER_GENES             = 24576
  MB_HTML_MARKER_GENES            = 8192
  MB_LBL_LABEL_CELLTYPES          = 16384
  MB_LBL_SAVE_SUBSET_SCES         = 16384
  MB_LBL_RENDER_TEMPLATE_RMD      = 4096
  MB_META_SAVE_METACELLS          = 16384
  MB_PB_MAKE_PBS                  = 8192
  MB_PB_CALC_EMPTY_GENES          = 8192
  MB_ZOOM_RUN_ZOOM                = 16384
  MB_ZOOM_RENDER_TEMPLATE_RMD     = 4096

  # change defaults if specified
  if ('resources' in config) and (config['resources'] is not None):
    if 'mb_run_alevin_fry' in config['resources']:
      MB_RUN_ALEVIN_FRY               = config['resources']['mb_run_alevin_fry']
    if 'mb_save_alevin_to_h5' in config['resources']:
      MB_SAVE_ALEVIN_TO_H5            = config['resources']['mb_save_alevin_to_h5']
    if 'mb_run_cellbender' in config['resources']:
      MB_RUN_CELLBENDER               = config['resources']['mb_run_cellbender']
    if 'mb_get_cellbender_qc_metrics' in config['resources']:
      MB_GET_CELLBENDER_QC_METRICS    = config['resources']['mb_get_cellbender_qc_metrics']
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

  return MB_RUN_ALEVIN_FRY, MB_SAVE_ALEVIN_TO_H5, \
    MB_RUN_CELLBENDER, MB_GET_CELLBENDER_QC_METRICS, \
    MB_RUN_SCDBLFINDER, MB_COMBINE_SCDBLFINDER_OUTPUTS, \
    MB_RUN_QC, \
    MB_MAKE_SCE_OBJECT, \
    MB_RUN_HARMONY, \
    MB_RUN_MARKER_GENES, MB_HTML_MARKER_GENES, \
    MB_LBL_LABEL_CELLTYPES, MB_LBL_SAVE_SUBSET_SCES, MB_LBL_RENDER_TEMPLATE_RMD, \
    MB_META_SAVE_METACELLS, \
    MB_PB_MAKE_PBS, MB_PB_CALC_EMPTY_GENES, \
    MB_ZOOM_RUN_ZOOM, MB_ZOOM_RENDER_TEMPLATE_RMD


