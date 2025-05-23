# load modules
import yaml
import pandas as pd
import os
import re
import warnings
import glob
from snakemake.utils import validate, min_version

from scprocess_utils import *

SCPROCESS_DATA_DIR = os.getenv('SCPROCESS_DATA_DIR')

# get parameters
PROJ_DIR, FASTQ_DIR, SHORT_TAG, FULL_TAG, YOUR_NAME, AFFILIATION, METADATA_F, METADATA_VARS, \
EXC_SAMPLES, SAMPLES, DATE_STAMP, CUSTOM_SAMPLE_PARAMS_F, SPECIES, \
DEMUX_TYPE, HTO_FASTQ_DIR, FEATURE_REF, DEMUX_F, BATCH_VAR, EXC_POOLS, POOL_IDS, SAMPLE_VAR, SAMPLE_MAPPING = \
  get_project_parameters(config, SCPROCESS_DATA_DIR)
AF_MITO_STR, AF_HOME_DIR, AF_INDEX_DIR, AF_GTF_DT_F, CHEMISTRY = \
  get_alevin_parameters(config, SCPROCESS_DATA_DIR, SPECIES)
CELLBENDER_IMAGE, CELLBENDER_PROP_MAX_KEPT, AMBIENT_METHOD, CELL_CALLS_METHOD, \
FORCE_EXPECTED_CELLS, FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE = \
  get_ambient_parameters(config)
QC_HARD_MIN_COUNTS, QC_HARD_MIN_FEATS, QC_HARD_MAX_MITO, QC_MIN_COUNTS, QC_MIN_FEATS, \
  QC_MIN_MITO, QC_MAX_MITO, QC_MIN_SPLICE, QC_MAX_SPLICE, QC_MIN_CELLS, DBL_MIN_FEATS = \
  get_qc_parameters(config)
HVG_METHOD, HVG_GROUP_VAR, HVG_CHUNK_SIZE, NUM_CHUNKS, GROUP_NAMES, CHUNK_NAMES, N_HVGS, EXCLUDE_AMBIENT_GENES = \
  get_hvg_parameters(config, METADATA_F, AF_GTF_DT_F)
INT_EXC_REGEX, INT_N_DIMS, INT_DBL_RES, INT_DBL_CL_PROP, INT_THETA, INT_RES_LS, INT_SEL_RES = \
  get_integration_parameters(config, AF_MITO_STR)
MKR_GSEA_DIR, MKR_MIN_CL_SIZE, MKR_MIN_CELLS, MKR_NOT_OK_RE, MKR_MIN_CPM_MKR, MKR_MIN_CPM_GO, MKR_MAX_ZERO_P, MKR_GSEA_CUT, CUSTOM_MKR_NAMES, CUSTOM_MKR_PATHS = \
  get_marker_genes_parameters(config, PROJ_DIR, SCPROCESS_DATA_DIR)
LBL_XGB_F, LBL_XGB_CLS_F, LBL_GENE_VAR, LBL_SEL_RES_CL, LBL_MIN_PRED, LBL_MIN_CL_PROP, LBL_MIN_CL_SIZE, LBL_SCE_SUBSETS, LBL_TISSUE, CUSTOM_LABELS_F = \
 get_label_celltypes_parameters(config, SPECIES, SCPROCESS_DATA_DIR)
ZOOM_NAMES, ZOOM_SPEC_LS = get_zoom_parameters(config, AF_MITO_STR, SCPROCESS_DATA_DIR)
META_SUBSETS, META_MAX_CELLS = get_metacells_parameters(config)
AMBIENT_GENES_GRP_NAMES, AMBIENT_GENES_GRP_VAR, AMBIENT_GENES_LOGFC_THR, AMBIENT_GENES_FDR_THR = get_pb_empties_parameters(config, HVG_METHOD, GROUP_NAMES, HVG_GROUP_VAR)
RETRIES, MB_RUN_ALEVIN_FRY, MB_SAVE_ALEVIN_TO_H5, \
  MB_RUN_AMBIENT, MB_GET_BARCODE_QC_METRICS, \
  MB_RUN_SCDBLFINDER, MB_COMBINE_SCDBLFINDER_OUTPUTS, \
  MB_RUN_QC, MB_RUN_HVGS, \
  MB_MAKE_SCE_OBJECT, \
  MB_RUN_HARMONY, \
  MB_RUN_MARKER_GENES, MB_HTML_MARKER_GENES, \
  MB_LBL_LABEL_CELLTYPES, MB_LBL_SAVE_SUBSET_SCES, MB_LBL_RENDER_TEMPLATE_RMD, \
  MB_META_SAVE_METACELLS, \
  MB_PB_MAKE_PBS, MB_PB_CALC_EMPTY_GENES, \
  MB_ZOOM_RUN_ZOOM, MB_ZOOM_RENDER_TEMPLATE_RMD = \
  get_resource_parameters(config)

# specify locations
fastqs_dir    = f"{PROJ_DIR}/data/fastqs"
code_dir      = f"{PROJ_DIR}/code"
af_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_alevin_fry"
af_rna_dir    = 'rna/' if DEMUX_TYPE == 'af' else ''
amb_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_ambient"
demux_dir     = f"{PROJ_DIR}/output/{SHORT_TAG}_demultiplexing"
dbl_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_doublet_id"
qc_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_qc"
hvg_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_hvg"
int_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_integration"
mkr_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_marker_genes"
lbl_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_label_celltypes"
meta_dir      = f"{PROJ_DIR}/output/{SHORT_TAG}_metacells"
pb_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_pseudobulk"
empty_dir     = f"{PROJ_DIR}/output/{SHORT_TAG}_empties"
zoom_dir      = f"{PROJ_DIR}/output/{SHORT_TAG}_zoom"
rmd_dir       = f"{PROJ_DIR}/analysis"
docs_dir      = f"{PROJ_DIR}/public"


# make nice zoom variables
zoom_df       = pd.DataFrame({ \
 'zoom_name': ZOOM_NAMES, \
 'zoom_res': [ ZOOM_SPEC_LS[ zn ][ 'zoom_res' ] for zn in ZOOM_NAMES] \
 })

# exclude all samples without fastq files
if DEMUX_TYPE != "":
 POOL_IDS = exclude_samples_without_fastq_files(FASTQ_DIR, POOL_IDS, HTO=False)
else:
 SAMPLES  = exclude_samples_without_fastq_files(FASTQ_DIR, SAMPLES, HTO=False)

# exclude all samples without hto fastq files
if DEMUX_TYPE == "af":
 POOL_IDS = exclude_samples_without_fastq_files(HTO_FASTQ_DIR, POOL_IDS, HTO=True)

# exclude pools and samples from sample mapping dictionary
SAMPLE_MAPPING = filter_sample_mapping(SAMPLE_MAPPING, POOL_IDS, SAMPLES)

# join all samples into a single string
POOL_STR   = ','.join(POOL_IDS)
SAMPLE_STR = ','.join(SAMPLES)

runs = POOL_IDS if DEMUX_TYPE != "" else SAMPLES
RUNS_STR = ','.join(runs)

# alevin hto index outputs (optional)
hto_index_outs = [
    af_dir + '/hto.tsv',
    af_dir + '/t2g_hto.tsv',
    af_dir + '/hto_index/ref_indexing.log'
  ] if DEMUX_TYPE == "af" else []

# alevin hto quantification outputs (optional)
hto_af_outs = expand(
  [
  af_dir    + '/af_{run}/hto/af_quant/',
  af_dir    + '/af_{run}/hto/af_quant/alevin/quants_mat.mtx',
  af_dir    + '/af_{run}/hto/af_quant/alevin/quants_mat_cols.txt',
  af_dir    + '/af_{run}/hto/af_quant/alevin/quants_mat_rows.txt',
  af_dir    + '/af_{run}/hto/af_hto_counts_mat.h5',
  af_dir    + '/af_{run}/hto/knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz'
  ], run = runs
) if DEMUX_TYPE == 'af' else []

# seurat demultiplexing outputs (optional)
hto_sce_fs = expand(
  demux_dir + '/sce_cells_htos_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
  run = runs
  ) if DEMUX_TYPE == "af" else []

# mutliplxing report (optional)
hto_rmd_f  = (rmd_dir   + '/' + SHORT_TAG + '_demultiplexing.Rmd') if DEMUX_TYPE == 'af' else []
hto_html_f = (docs_dir  + '/' + SHORT_TAG + '_demultiplexing.html') if DEMUX_TYPE == 'af' else []

# fgsea outputs (optional)
fgsea_outs = [
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'go_bp_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'go_cc_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'go_mf_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'paths_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'hlmk_' + DATE_STAMP + '.txt.gz', 
] if SPECIES in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020'] else []


# one rule to rule them all
rule all:
  input:
    # hto outputs
    hto_index_outs, 
    hto_af_outs, 
    expand(
      [
      # alevin_fry
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_quant/',
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_quant/alevin/quants_mat.mtx',
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_quant/alevin/quants_mat_cols.txt',
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_quant/alevin/quants_mat_rows.txt',
      af_dir    + '/af_{run}/' + af_rna_dir + 'af_counts_mat.h5',
      af_dir    + '/af_{run}/' + af_rna_dir + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz',
      af_dir    + '/af_{run}/' + af_rna_dir + 'ambient_params_{run}_' + DATE_STAMP + '.yaml',
      # ambient (cellbender, decontx or nothing)
      amb_dir   + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml',
      # barcode qc metrics
      amb_dir   + '/ambient_{run}/barcodes_qc_metrics_{run}_' + DATE_STAMP + '.txt.gz',
      # qc
      #qc_dir  + '/qc_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
      #qc_dir  + '/coldata_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
      #qc_dir  + '/rowdata_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
      # doublet id
      dbl_dir + '/dbl_{run}/scDblFinder_{run}_outputs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
      dbl_dir + '/dbl_{run}/scDblFinder_{run}_dimreds_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
      ], run =  runs), 
    # ambient sample statistics
    amb_dir + '/ambient_sample_statistics_' + DATE_STAMP + '.txt',  
    # demultiplexing
    hto_sce_fs,  
    # qc
    qc_dir  + '/qc_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    qc_dir  + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    qc_dir  + '/qc_sample_statistics_' + DATE_STAMP + '.txt', 
    qc_dir  + '/sce_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml', 
    # pseudobulks and empties
    pb_dir  + '/af_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    pb_dir  + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
    pb_dir  + '/pb_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    empty_dir + '/edger_empty_genes_' + FULL_TAG + '_all_' + DATE_STAMP + '.txt.gz', 
    # hvgs
    hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    hvg_dir + '/standardized_variance_stats_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    hvg_dir + '/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    hvg_dir + '/top_hvgs_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5'


    
      # integration
      #int_dir   + '/sce_clean_'           + FULL_TAG + '_' + DATE_STAMP + '.rds',
      #int_dir   + '/integrated_dt_'       + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
      #int_dir   + '/harmony_hvgs_'        + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
      # marker genes
      #mkr_dir   + '/pb_'              + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.rds',
      #mkr_dir   + '/pb_marker_genes_' + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.txt.gz',
      #mkr_dir   + '/pb_hvgs_'         + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.txt.gz',
      # fgsea outputs
      #fgsea_outs, 
      # code
      #code_dir  + '/utils.R',
      #code_dir  + '/ambient.R',
      #code_dir  + '/qc.R', 
      #code_dir  + '/make_sce.R',
      #code_dir  + '/doublet_id.R',
      #code_dir  + '/integration.R', 
      #code_dir  + '/marker_genes.R',
      #code_dir  + '/zoom.R', 
      #code_dir  + '/multiplexing.R',
      # markdowns
      #rmd_dir   + '/' + SHORT_TAG + '_alevin_fry.Rmd',
      #rmd_dir   + '/' + SHORT_TAG + '_ambient.Rmd', 
      #rmd_dir   + '/' + SHORT_TAG + '_qc.Rmd', 
      #rmd_dir   + '/' + SHORT_TAG + '_integration.Rmd', 
      #rmd_dir   + '/' + SHORT_TAG + f'_marker_genes_{INT_SEL_RES}.Rmd', 
      #hto_rmd_f, 
      # reports
      #docs_dir  + '/' + SHORT_TAG + '_alevin_fry.html',
      #docs_dir  + '/' + SHORT_TAG + '_ambient.html', 
      #docs_dir  + '/' + SHORT_TAG + '_qc.html', 
      #docs_dir  + '/' + SHORT_TAG + '_integration.html',
      #docs_dir  + '/' + SHORT_TAG + f'_marker_genes_{INT_SEL_RES}.html', 
      #hto_html_f 

rule simpleaf:
  input:
    hto_index_outs, 
    hto_af_outs,
    expand( 
      [
      af_dir    + '/af_{sample}/' + af_rna_dir + 'af_quant/',
      af_dir    + '/af_{sample}/' + af_rna_dir + 'af_quant/alevin/quants_mat.mtx',
      af_dir    + '/af_{sample}/' + af_rna_dir + 'af_quant/alevin/quants_mat_cols.txt',
      af_dir    + '/af_{sample}/' + af_rna_dir + 'af_quant/alevin/quants_mat_rows.txt',
      af_dir    + '/af_{sample}/' + af_rna_dir + 'af_counts_mat.h5',
      af_dir    + '/af_{sample}/' + af_rna_dir + 'knee_plot_data_{sample}_' + DATE_STAMP + '.txt.gz',
      af_dir    + '/af_{sample}/' + af_rna_dir + 'ambient_params_{sample}_' + DATE_STAMP + '.yaml',
      ],
     sample = runs), 
     rmd_dir   + '/' + SHORT_TAG + '_alevin_fry.Rmd',
     docs_dir  + '/' + SHORT_TAG + '_alevin_fry.html'


rule ambient:
  input: 
    expand(amb_dir + '/ambient_{sample}/ambient_{sample}_' + DATE_STAMP + '_output_paths.yaml',
    sample = POOL_IDS if DEMUX_TYPE != "" else SAMPLES),
    amb_dir + '/ambient_{sample}/barcodes_qc_metrics_{sample}_' + DATE_STAMP + '.txt.gz', 
    amb_dir + '/ambient_sample_statistics_' + DATE_STAMP + '.txt', 
    rmd_dir   + '/' + SHORT_TAG + '_ambient.Rmd', 
    docs_dir  + '/' + SHORT_TAG + '_ambient.html'

rule qc:
  input:     
    qc_dir    + '/qc_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    qc_dir    + '/keep_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    rmd_dir   + '/' + SHORT_TAG + '_qc.Rmd',
    docs_dir  + '/' + SHORT_TAG + '_qc.html'


rule integration:
  input:     
    int_dir   + '/sce_clean_'           + FULL_TAG + '_' + DATE_STAMP + '.rds',
    int_dir   + '/integrated_dt_'       + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    int_dir   + '/harmony_hvgs_'        + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    rmd_dir   + '/' + SHORT_TAG + '_integration.Rmd',
    docs_dir  + '/' + SHORT_TAG + '_integration.html'


rule marker_genes:
  input:     
    mkr_dir   + '/pb_'              + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.rds',
    mkr_dir   + '/pb_marker_genes_' + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/pb_hvgs_'         + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'go_bp_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'go_cc_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'go_mf_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'paths_' + DATE_STAMP + '.txt.gz',
    mkr_dir   + '/fgsea_'           + FULL_TAG + f'_{INT_SEL_RES}_' + 'hlmk_' + DATE_STAMP + '.txt.gz',
    rmd_dir   + '/' + SHORT_TAG + f'_marker_genes_{INT_SEL_RES}.Rmd',
    docs_dir  + '/' + SHORT_TAG + f'_marker_genes_{INT_SEL_RES}.html'


rule label_celltypes:
  input:
    lbl_dir + '/hvg_mat_for_labelling_' + LBL_GENE_VAR + '_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    lbl_dir + '/cell_annotations_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    lbl_dir   + '/sce_subset_specifications_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    expand([
      lbl_dir   +'/sce_subset_' + FULL_TAG + '_{s}_' + DATE_STAMP + '.rds'
    ], s = [] if LBL_SCE_SUBSETS is None else [*LBL_SCE_SUBSETS] ), 
    code_dir  + '/label_celltypes.R',
    rmd_dir   + '/' + SHORT_TAG + '_label_celltypes.Rmd', 
    docs_dir  + '/' + SHORT_TAG + '_label_celltypes.html'


rule zoom:
  input:
    # zoom_imputed_dt
    expand('%s/{zoom_name}/zoom_imputed_dt_%s_{zoom_name}_{zoom_res}_%s.txt.gz' % \
           (zoom_dir, FULL_TAG, DATE_STAMP), \
           zip, zoom_name=zoom_df['zoom_name'], zoom_res=zoom_df['zoom_res']),
    # zoom_sce_clean
    expand('%s/{zoom_name}/zoom_sce_clean_%s_{zoom_name}_{zoom_res}_%s.rds' % \
           (zoom_dir, FULL_TAG, DATE_STAMP), \
           zip, zoom_name=zoom_df['zoom_name'], zoom_res=zoom_df['zoom_res']),
    # zoom_integrated_dt
    expand('%s/{zoom_name}/zoom_integrated_dt_%s_{zoom_name}_{zoom_res}_%s.txt.gz' % \
           (zoom_dir, FULL_TAG, DATE_STAMP), \
           zip, zoom_name=zoom_df['zoom_name'], zoom_res=zoom_df['zoom_res']),
    # zoom_pb
    expand('%s/{zoom_name}/zoom_pb_%s_{zoom_name}_{zoom_res}_%s.rds' % \
           (zoom_dir, FULL_TAG, DATE_STAMP), \
           zip, zoom_name=zoom_df['zoom_name'], zoom_res=zoom_df['zoom_res']),
    # zoom_pb_marker_genes
    expand('%s/{zoom_name}/zoom_pb_marker_genes_%s_{zoom_name}_{zoom_res}_%s.txt.gz' % \
           (zoom_dir, FULL_TAG, DATE_STAMP), \
           zip, zoom_name=zoom_df['zoom_name'], zoom_res=zoom_df['zoom_res']),
    # zoom_pb_hvgs
    expand('%s/{zoom_name}/zoom_pb_hvgs_%s_{zoom_name}_{zoom_res}_%s.txt.gz' % \
           (zoom_dir, FULL_TAG, DATE_STAMP), \
           zip, zoom_name=zoom_df['zoom_name'], zoom_res=zoom_df['zoom_res']),
    # zoom_fgsea_go_bp
    expand('%s/{zoom_name}/zoom_fgsea_%s_{zoom_name}_{zoom_res}_go_bp_%s.txt.gz' % \
           (zoom_dir, FULL_TAG, DATE_STAMP), \
           zip, zoom_name=zoom_df['zoom_name'], zoom_res=zoom_df['zoom_res']),
    # zoom_fgsea_go_cc
    expand('%s/{zoom_name}/zoom_fgsea_%s_{zoom_name}_{zoom_res}_go_cc_%s.txt.gz' % \
           (zoom_dir, FULL_TAG, DATE_STAMP), \
           zip, zoom_name=zoom_df['zoom_name'], zoom_res=zoom_df['zoom_res']),
    # zoom_fgsea_go_mf
    expand('%s/{zoom_name}/zoom_fgsea_%s_{zoom_name}_{zoom_res}_go_mf_%s.txt.gz' % \
           (zoom_dir, FULL_TAG, DATE_STAMP), \
           zip, zoom_name=zoom_df['zoom_name'], zoom_res=zoom_df['zoom_res']),
    # zoom_fgsea_paths
    expand('%s/{zoom_name}/zoom_fgsea_%s_{zoom_name}_{zoom_res}_paths_%s.txt.gz' % \
           (zoom_dir, FULL_TAG, DATE_STAMP), \
           zip, zoom_name=zoom_df['zoom_name'], zoom_res=zoom_df['zoom_res']),
    # zoom_fgsea_hlmk
    expand('%s/{zoom_name}/zoom_fgsea_%s_{zoom_name}_{zoom_res}_hlmk_%s.txt.gz' % \
           (zoom_dir, FULL_TAG, DATE_STAMP), \
           zip, zoom_name=zoom_df['zoom_name'], zoom_res=zoom_df['zoom_res']), 
    # Rmd and html files
    expand('%s/%s_zoom_{zoom_name}_{zoom_res}.Rmd' % (rmd_dir, SHORT_TAG), \
           zip, zoom_name=zoom_df['zoom_name'], zoom_res=zoom_df['zoom_res']),
    expand('%s/%s_zoom_{zoom_name}_{zoom_res}.html' % (docs_dir, SHORT_TAG), \
           zip, zoom_name=zoom_df['zoom_name'], zoom_res=zoom_df['zoom_res'])


#rule pb_empties:
# input:
#    pb_dir + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
#    expand([
#     pb_dir + '/pb_subset_' + FULL_TAG + '_{subset}_' + DATE_STAMP + '.rds'
#     ], subset = PB_SUBSETS ),
#    expand([
#     empty_dir + '/edger_empty_genes_' + FULL_TAG + '_{subset}_' + DATE_STAMP + '.txt.gz'
#     ], subset = PB_SUBSETS ),
#     (pb_dir + '/pb_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds') if PB_DO_ALL else [],
#     (empty_dir + '/edger_empty_genes_' + FULL_TAG + '_all_' + DATE_STAMP + '.txt.gz') if PB_DO_ALL else []


rule metacells:
 input:
    expand([
     meta_dir  + '/metacells_sce_' + FULL_TAG + '_{subset}_{max_cells}_' + DATE_STAMP + '.rds'
     ], subset = META_SUBSETS, max_cells = META_MAX_CELLS ),
    expand([
     meta_dir  + '/metacells_map_' + FULL_TAG + '_{subset}_{max_cells}_' + DATE_STAMP + '.txt.gz'
     ], subset = META_SUBSETS, max_cells = META_MAX_CELLS )

      

# define rules that are needed
include: "rules/alevin_fry.smk"
include: "rules/ambient.smk"
include: "rules/make_sce.smk"
include: "rules/doublet_id.smk"
include: "rules/qc.smk"
include: "rules/pb_empties.smk"
include: "rules/hvgs.smk"
#include: "rules/integration.smk"
#include: "rules/marker_genes.smk"
#include: "rules/render_htmls.smk"
#include: "rules/label_and_subset.smk"
#include: "rules/zoom.smk"
#include: "rules/metacells.smk"
#include: "rules/pb_empties.smk"
