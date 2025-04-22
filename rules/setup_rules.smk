import pandas as pd
import os
import re
import glob
import subprocess
import numpy as np
import sys
import warnings

sys.path.append('./scripts')
from setup import get_af_index_parameters
# from setup_utils import check_valid_index

SCPROCESS_DATA_DIR = os.getenv('SCPROCESS_DATA_DIR')

#configfile = '/Users/marusa/Projects/scprocess_test/configs/config-setup-test.yaml'
# get parameters from configfile
IDX_PARAMS_LS = get_af_index_parameters(config) 
GENOMES       = list(IDX_PARAMS_LS.keys())

# define simpleaf index files
AF_INDEX_FS = [
  'index/piscem_idx.ctab',
  'index/piscem_idx.ectab',
  'index/piscem_idx.json',
  'index/piscem_idx.refinfo',
  'index/piscem_idx.sshash',
  'index/piscem_idx_cfish.json',
  'index/simpleaf_index.json',
  'index/t2g_3col.tsv',
  'ref/gene_id_to_name.tsv',
  'ref/roers_make-ref.json',
  'ref/roers_ref.fa',
  'ref/t2g_3col.tsv',
  'simpleaf_index_log.json'
]

# define tiny rule
localrules: save_index_parameters_csv

# define top level rule
rule all:
  input:
    # rule download scprocess repo data
    SCPROCESS_DATA_DIR + '/marker_genes/canonical_brain_celltype_markers_human_2023-10-17.txt',
    SCPROCESS_DATA_DIR + '/marker_genes/canonical_brain_celltype_markers_mouse_2023-10-17.txt',
    SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3LT.txt',
    SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3v2_5v1_5v2.txt',
    SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3v3.txt',
    SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_multiome_gex.txt',
    SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3v4.txt',
    SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_5v3.txt',
    SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_whitelists.csv',
    SCPROCESS_DATA_DIR + '/gmt_pathways/c1.all.v2023.1.Hs.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/c2.all.v2023.1.Hs.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/c2.cp.biocarta.v2023.1.Hs.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/c2.cp.kegg.v2023.1.Hs.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/c2.cp.v2023.1.Hs.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/c3.mir.v2023.1.Hs.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/c3.tft.v2023.1.Hs.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/c5.go.bp.v2023.1.Hs.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/c5.go.cc.v2023.1.Hs.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/c5.go.mf.v2023.1.Hs.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/c5.hpo.v2023.1.Hs.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/c7.immunesigdb.v2023.1.Hs.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/c8.all.v2023.1.Hs.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/h.all.v2023.1.Hs.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/m2.cp.biocarta.v2023.1.Mm.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/m2.cp.v2023.1.Mm.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/m5.go.bp.v2023.1.Mm.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/m5.go.cc.v2023.1.Mm.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/m5.go.mf.v2023.1.Mm.symbols.gmt',
    SCPROCESS_DATA_DIR + '/gmt_pathways/mh.all.v2023.1.Mm.symbols.gmt',
    SCPROCESS_DATA_DIR + '/xgboost/Siletti_Macnair-2024-03-11/allowed_cls_Siletti_Macnair_2024-03-11.csv',
    SCPROCESS_DATA_DIR + '/xgboost/Siletti_Macnair-2024-03-11/xgboost_obj_hvgs_Siletti_Macnair_2024-03-11.rds',
    # rule download_or_build_af_indices
    # expand(SCPROCESS_DATA_DIR + '/alevin_fry_home/{genome}/index/{file}', genome=GENOMES, file=AF_INDEX_FS)
    expand([ SCPROCESS_DATA_DIR + '/alevin_fry_home/{genome}/' + f'{file}' for file in AF_INDEX_FS], genome=GENOMES),
    expand(SCPROCESS_DATA_DIR + '/alevin_fry_home/{genome}/{genome}_index_params.yaml', genome=GENOMES),
    # rule get_reference_genome_data 
    SCPROCESS_DATA_DIR + '/index_parameters.csv'


# rule for getting scprocess data from github repo (maybe not a good idea to have all files as outputs)
rule download_scprocess_files:
  output:
    brain_mkrs_hsa   = SCPROCESS_DATA_DIR + '/marker_genes/canonical_brain_celltype_markers_human_2023-10-17.txt',
    brain_mkrs_mmu   = SCPROCESS_DATA_DIR + '/marker_genes/canonical_brain_celltype_markers_mouse_2023-10-17.txt',
    wl_3lt           = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3LT.txt',
    wl_3v2_5v1_5v2   = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3v2_5v1_5v2.txt',
    wl_3v3           = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3v3.txt',
    wl_3v4	         = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3v4.txt',
    wl_5v3	         = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_5v3.txt',
    wl_multiome      = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_multiome_gex.txt',
    all_wl           = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_whitelists.csv',
    gmt_f_1          = SCPROCESS_DATA_DIR + '/gmt_pathways/c1.all.v2023.1.Hs.symbols.gmt',
    gmt_f_2          = SCPROCESS_DATA_DIR + '/gmt_pathways/c2.all.v2023.1.Hs.symbols.gmt',
    gmt_f_3          = SCPROCESS_DATA_DIR + '/gmt_pathways/c2.cp.biocarta.v2023.1.Hs.symbols.gmt',
    gmt_f_4          = SCPROCESS_DATA_DIR + '/gmt_pathways/c2.cp.kegg.v2023.1.Hs.symbols.gmt',
    gmt_f_5          = SCPROCESS_DATA_DIR + '/gmt_pathways/c2.cp.v2023.1.Hs.symbols.gmt',
    gmt_f_6          = SCPROCESS_DATA_DIR + '/gmt_pathways/c3.mir.v2023.1.Hs.symbols.gmt',
    gmt_f_7          = SCPROCESS_DATA_DIR + '/gmt_pathways/c3.tft.v2023.1.Hs.symbols.gmt',
    gmt_f_8          = SCPROCESS_DATA_DIR + '/gmt_pathways/c5.go.bp.v2023.1.Hs.symbols.gmt',
    gmt_f_9          = SCPROCESS_DATA_DIR + '/gmt_pathways/c5.go.cc.v2023.1.Hs.symbols.gmt',
    gmt_f_10         = SCPROCESS_DATA_DIR + '/gmt_pathways/c5.go.mf.v2023.1.Hs.symbols.gmt',
    gmt_f_11         = SCPROCESS_DATA_DIR + '/gmt_pathways/c5.hpo.v2023.1.Hs.symbols.gmt',
    gmt_f_12         = SCPROCESS_DATA_DIR + '/gmt_pathways/c7.immunesigdb.v2023.1.Hs.symbols.gmt',
    gmt_f_13         = SCPROCESS_DATA_DIR + '/gmt_pathways/c8.all.v2023.1.Hs.symbols.gmt',
    gmt_f_14         = SCPROCESS_DATA_DIR + '/gmt_pathways/h.all.v2023.1.Hs.symbols.gmt',
    gmt_f_15         = SCPROCESS_DATA_DIR + '/gmt_pathways/m2.cp.biocarta.v2023.1.Mm.symbols.gmt',
    gmt_f_16         = SCPROCESS_DATA_DIR + '/gmt_pathways/m2.cp.v2023.1.Mm.symbols.gmt',
    gmt_f_17         = SCPROCESS_DATA_DIR + '/gmt_pathways/m5.go.bp.v2023.1.Mm.symbols.gmt',
    gmt_f_18         = SCPROCESS_DATA_DIR + '/gmt_pathways/m5.go.cc.v2023.1.Mm.symbols.gmt',
    gmt_f_19         = SCPROCESS_DATA_DIR + '/gmt_pathways/m5.go.mf.v2023.1.Mm.symbols.gmt',
    gmt_f_20         = SCPROCESS_DATA_DIR + '/gmt_pathways/mh.all.v2023.1.Mm.symbols.gmt',
    xgb_csv_f        = SCPROCESS_DATA_DIR + '/xgboost/Siletti_Macnair-2024-03-11/allowed_cls_Siletti_Macnair_2024-03-11.csv',
    xgb_rds_f        = SCPROCESS_DATA_DIR + '/xgboost/Siletti_Macnair-2024-03-11/xgboost_obj_hvgs_Siletti_Macnair_2024-03-11.rds'
  conda:
    '../envs/py_env.yml'
  threads: 1
  shell:
    """
    python3 scripts/setup.py get_scprocess_data {SCPROCESS_DATA_DIR}
    """


# rule for downloading reference genome files from 10x and dealing with custom genomes
rule set_up_one_af_index:
  output:
    [ SCPROCESS_DATA_DIR + '/alevin_fry_home/{genome}/' + f'{file}' for file in AF_INDEX_FS],
    SCPROCESS_DATA_DIR + "/alevin_fry_home/{genome}/{genome}_index_params.yaml"
  params:
    fasta       = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ]['fasta'],
    gtf         = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ]['gtf'],
    index_dir   = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ]['index_dir'],
    mito_str    = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ]['mito_str'],
    is_prebuilt = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ]['is_prebuilt'],
    is_tenx     = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ]['is_tenx'],
    has_decoys  = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ]['has_decoys'],
    has_rrna    = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ]['has_rrna']
  conda:
    '../envs/alevin_fry.yml'
  resources:
    mem_mb = 8192
  threads: 8
  shell:
    """  
    python3 scripts/setup.py set_up_af_index {SCPROCESS_DATA_DIR} {wildcards.genome} \
      {params.fasta} {params.gtf} {params.index_dir} {params.mito_str} \
      {params.is_prebuilt} {params.is_tenx} {params.has_decoys} {params.has_rrna} {threads}
    """


rule save_index_parameters_csv:
  input:
    yamls   = expand(SCPROCESS_DATA_DIR + "/alevin_fry_home/{genome}/{genome}_index_params.yaml", genome = GENOMES)
  output:
    csv     = SCPROCESS_DATA_DIR + '/index_parameters.csv'
  conda:
    '../envs/py_env.yml'
  resources:
    mem_mb = 512
  threads: 1
  shell:
    """
    python3 scripts/setup.py save_index_params_csv {output.csv} {input.yamls}
    """
