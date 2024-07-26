# load modules
import warnings
import yaml
import pandas as pd
import os
import re
import glob
import datetime
import subprocess
from snakemake.utils import validate, min_version


def get_setup_parameters(config_f):
    with open(config_f, "r") as stream:
        config = yaml.safe_load(stream)

    # check if required arguments are in config  
    for v in ["scprocess_data_dir", "index_dir_name"]:
        assert v in config, f"{v} not in config file"

    # set defaults (if rule setup runs and 'setup' is not specified in the config there should be an error)
    SCPROCESS_DATA_DIR = config["scprocess_data_dir"]
    INDEX_DIR_NAME = config["index_dir_name"]
    FASTA = None 
    GTF = None
    USE_DECOYS = False

    # change defaults if specified in config
    if ('setup' in config) and (config['setup'] is not None):
        if 'fasta' in config['setup']:
            FASTA = config['setup']['fasta']
        if 'gtf' in config['setup']:
            GTF = config['setup']['gtf']
        if 'use_decoys' in config['setup']:
            USE_DECOYS = config['setup']['use_decoys']

    return SCPROCESS_DATA_DIR, INDEX_DIR_NAME, FASTA, GTF, USE_DECOYS
      
  

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
def get_alevin_parameters(config):
  # code to install alevin:
  # ml .testing
  # ml Anaconda3/2023.03-libmamba
  # conda create --prefix /projects/site/pred/neurogenomics/resources/scprocess_data/conda/af python=3.11
  # conda activate /projects/site/pred/neurogenomics/resources/scprocess_data/conda/af
  # conda install simpleaf piscem -c bioconda -c conda-forge
  # ALEVIN_FRY_IMAGE = "docker://combinelab/usefulaf:latest"

  # define locations
  # AF_CONDA_ENV  = "/projects/site/pred/neurogenomics/resources/scprocess_data/conda/af"
  AF_CONDA_ENV  = "/projects/site/pred/neurogenomics/resources/scprocess_data/conda/simpleaf_0.15.1"
  AF_PISCEM_BIN = "/projects/site/pred/neurogenomics/resources/scprocess_data/conda/piscem-x86_64-unknown-linux-gnu/piscem"
  scdata_dir    = "/projects/site/pred/neurogenomics/resources/scprocess_data/data"
  AF_HOME_DIR   = os.path.join(scdata_dir, "alevin_fry_home")

  # get species
  SPECIES       = config["alevin"]["species"]
  assert SPECIES in ["human", "mouse"], f"species {SPECIES} not yet supported"

  # define some defaults
  AF_CHEMISTRY  = "10xv3"
  if SPECIES == "human":
    AF_REF_TXOME  = "human"
    AF_MITO_STR   = "^MT-"
  elif SPECIES == "mouse":
    AF_REF_TXOME  = "mouse"
    AF_MITO_STR   = "^mt-"

  # get any non-default values
  if ('alevin' in config) and (config['alevin'] is not None):
    if 'chemistry' in config['alevin']:
      AF_CHEMISTRY  = config['alevin']['chemistry']
    if 'ref_txome' in config['alevin']:
      AF_REF_TXOME  = config['alevin']['ref_txome']

  # check inputs look ok
  assert AF_CHEMISTRY in ["10xv3", "10xv2", "multiome"], "chemistry not recognised"
  if SPECIES == "human":
    assert AF_REF_TXOME in ["human", "human_ebv", "human_Pool_2023"], "reference transcriptome not recognised"
  elif SPECIES == "mouse":
    assert AF_REF_TXOME in ["mouse", "mouse_Pool_2023"], "reference transcriptome not recognised"
  elif SPECIES == "cyno":
    assert AF_REF_TXOME in ["cyno"], "reference transcriptome not recognised"
  else:
    raise ValueError(f"species {SPECIES} does not have ref txome defined")

  # get alevin parameters
  if AF_REF_TXOME == "human":
    AF_INDEX_DIR  = os.path.join(scdata_dir, "alevin_fry_home/human_2020-A_splici_rRNA")
    AF_GTF_DT_F   = os.path.join(scdata_dir, "reference_genome/refdata-gex-GRCh38-2020-A-rRNA/refdata-gex-GRCh38-2020-A-rRNA_gtf_dt.txt.gz")
    # AF_INDEX_DIR  = os.path.join(scdata_dir, "alevin_fry_home/human_2020-A_splici_clean")
    # AF_GTF_DT_F   = os.path.join(scdata_dir, "reference_genome/refdata-gex-GRCh38-2020-A/gex-GRCh38-2020-A_gtf_dt.txt.gz")

  elif AF_REF_TXOME == "human_ebv":
    AF_INDEX_DIR  = os.path.join(scdata_dir, "alevin_fry_home/ebv_splici_clean")
    AF_GTF_DT_F   = os.path.join(scdata_dir, "reference_genome/refdata-ebv/human_ebv_gtf_dt.txt.gz")

  elif AF_REF_TXOME == "human_Pool_2023":
    AF_INDEX_DIR  = os.path.join(scdata_dir, "alevin_fry_home/human_pool_2023_splici")
    AF_GTF_DT_F   = os.path.join(scdata_dir, "reference_genome/Pool_2023_human/human_Pool_2023_gtf_dt.txt.gz")

  elif AF_REF_TXOME == "mouse":
    AF_INDEX_DIR  = os.path.join(scdata_dir, "alevin_fry_home/mouse_2020-A_splici_clean")
    AF_GTF_DT_F   = os.path.join(scdata_dir, "reference_genome/refdata-gex-mm10-2020-A/gex-mm10-2020-A_gtf_dt.txt.gz")

  elif AF_REF_TXOME == "mouse_Pool_2023":
    AF_INDEX_DIR  = os.path.join(scdata_dir, "alevin_fry_home/mouse_pool_2023_splici")
    AF_GTF_DT_F   = os.path.join(scdata_dir, "reference_genome/Pool_2023_mouse/mouse_Pool_2023_gtf_dt.txt.gz")

  elif AF_REF_TXOME == "cyno":
    AF_INDEX_DIR  = os.path.join(scdata_dir, "alevin_fry_home/cyno_macfas5_splici")
    AF_GTF_DT_F   = os.path.join(scdata_dir, "reference_genome/cyno_macfas5/cyno_macfas5_gtf_dt.txt.gz")

  else:
    raise ValueError("something went wrong with AF_REF_TXOME :(")

  # define list of nice barcodes
  if AF_CHEMISTRY == "10xv3":
    AF_WHITELIST_F  = "/projects/site/pred/neurogenomics/resources/scprocess_data/data/cellranger_ref/cellranger_barcode_whitelist_v3.txt"

  elif AF_CHEMISTRY == "10xv2":
    AF_WHITELIST_F  = "/projects/site/pred/neurogenomics/resources/scprocess_data/data/cellranger_ref/cellranger_barcode_whitelist_v2.txt"

  elif AF_CHEMISTRY == "multiome":
    AF_CHEMISTRY    = "10xv3"
    AF_WHITELIST_F  = "/projects/site/pred/neurogenomics/resources/scprocess_data/data/cellranger_ref/cellranger_barcode_whitelist_multiome_GEX.txt"

  # check that the specified files exist
  for p in [AF_INDEX_DIR, AF_GTF_DT_F, AF_WHITELIST_F]:
    assert os.path.exists(p), f"specified path {p} does not exist"

  return SPECIES, AF_MITO_STR, AF_CONDA_ENV, AF_PISCEM_BIN, AF_HOME_DIR, AF_CHEMISTRY, AF_INDEX_DIR, AF_GTF_DT_F, AF_WHITELIST_F


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


# define make_sce parameters
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


# define QC parameters
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
def get_integration_parameters(config): 
  # set default values
  SPECIES       = config["alevin"]["species"]
  if SPECIES == "human":
    INT_EXC_REGEX   = "^MT-"
  elif SPECIES == "mouse":
    INT_EXC_REGEX   = "^mt-"
  elif SPECIES == "human_ebv":
    INT_EXC_REGEX   = "^MT-"
  elif SPECIES == "cyno":
    INT_EXC_REGEX   = "^(ATP6|ATP8|COX1|COX2|COX3|CYTB|ND1|ND2|ND3|ND4|ND4L|ND5|ND6)$"

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
def get_marker_genes_parameters(config, SPECIES): 
  # set some more default values
  MKR_GSEA_DIR    = "/projects/site/pred/neurogenomics/resources/scprocess_data/data/gmt_pathways"
  MKR_MIN_CL_SIZE = 1e2
  MKR_MIN_CELLS   = 10
  MKR_NOT_OK_RE   = "(lincRNA|lncRNA|pseudogene|antisense)"
  MKR_MIN_CPM_MKR = 50
  MKR_MIN_CPM_GO  = 1
  MKR_MAX_ZERO_P  = 0.5
  MKR_GSEA_CUT    = 0.1

  # specify canonical marker file
  if SPECIES == "human":
    MKR_CANON_F     = "/projects/site/pred/neurogenomics/resources/scprocess_data/data/marker_genes/canonical_brain_celltype_markers_human_2023-10-17.txt"
  elif SPECIES == "mouse":
    MKR_CANON_F     = "/projects/site/pred/neurogenomics/resources/scprocess_data/data/marker_genes/canonical_brain_celltype_markers_mouse_2023-10-17.txt"
  elif SPECIES == "human_ebv":
    MKR_CANON_F     = "/projects/site/pred/neurogenomics/resources/scprocess_data/data/marker_genes/canonical_brain_celltype_markers_human_2023-10-17.txt"
  elif SPECIES == "cyno":
    MKR_CANON_F     = "/projects/site/pred/neurogenomics/resources/scprocess_data/data/marker_genes/canonical_brain_celltype_markers_human_2023-10-17.txt"
  else:
    ValueError(f"MKR_CANON_F not defined for species {SPECIES}")

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


def __get_cl_ls(config):
  # get parameters
  PROJ_DIR, SHORT_TAG, FULL_TAG, _, _, _, _, _, _, DATE_STAMP = \
    get_project_parameters(config)
  _, _, _, _, _, _, _, INT_SEL_RES = \
    get_integration_parameters(config)

  # specify harmony outputs
  int_dir     = f"{PROJ_DIR}/output/{SHORT_TAG}06_integration"
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
    if 'mb_zoom_run_zoom' in config['resources']:
      MB_ZOOM_RUN_ZOOM                = config['resources']['mb_zoom_run_zoom']
    if 'mb_zoom_render_template_rmd' in config['resources']:
      MB_ZOOM_RENDER_TEMPLATE_RMD     = config['resources']['mb_zoom_render_template_rmd']

  return MB_RUN_ALEVIN_FRY, MB_SAVE_ALEVIN_TO_H5, MB_RUN_CELLBENDER, MB_GET_CELLBENDER_QC_METRICS, \
    MB_RUN_SCDBLFINDER, MB_COMBINE_SCDBLFINDER_OUTPUTS, MB_RUN_QC, MB_MAKE_SCE_OBJECT, MB_RUN_HARMONY, \
    MB_RUN_MARKER_GENES, MB_HTML_MARKER_GENES, MB_LBL_LABEL_CELLTYPES, MB_LBL_SAVE_SUBSET_SCES, \
    MB_LBL_RENDER_TEMPLATE_RMD, MB_ZOOM_RUN_ZOOM, MB_ZOOM_RENDER_TEMPLATE_RMD


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
def exclude_samples_without_fastq_files(PROJ_DIR, SAMPLES):
  # get parameters
  fastqs_dir    = f"{PROJ_DIR}/data/fastqs"

  # get fastq files for each sample
  chk_samples   = []
  for sample in SAMPLES:
    R1_fs = find_fastq_files(fastqs_dir, sample, "R1")
    if len(R1_fs) > 0:
      chk_samples.append(sample)
    else:
      print(f"WARNING: no fastq files found for sample {sample}; excluded")

  return chk_samples


def get_mem_mb(wildcards, attempt, threads, mem_chunk = 1048):
  return attempt * threads * mem_chunk


def render_html(proj_dir, rlibs_dir, template_f, template_dict, rmd_f):
  if not os.path.isfile(rmd_f):
    # read the contents of the template_f file
    with open(template_f, 'r') as f:
      template_str = f.read()

    # create a string.Template object using the contents of the file
    from string import Template
    template    = Template(template_str)
    filled_str  = template.substitute(template_dict)
    with open(rmd_f, 'w') as f:
      f.write(filled_str)

  # render rmd file via Rscript
  bash_str = f"""
    # set up R
    ml purge
    ml .testing
    ml R-bundle-Bioconductor/3.18-foss-2020a-R-4.3.2
    ml libgit2/1.1.0-GCCcore-9.3.0

    # render rmd files
    export R_LIBS_USER='{rlibs_dir}'
    Rscript -e "source('scripts/render_reports.R'); \
      render_reports('{proj_dir}', rmd_ls_concat = '{rmd_f}')"
    """
  subprocess.run(bash_str, shell = True)
