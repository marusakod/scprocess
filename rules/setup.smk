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

# get parameters
SCDATA_DIR    = os.getenv('SCPROCESS_DATA_DIR')
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
    f'{SCDATA_DIR}/marker_genes/human_brain.csv',
    f'{SCDATA_DIR}/marker_genes/mouse_brain.csv',
    f'{SCDATA_DIR}/cellranger_ref/cellranger_gex_barcode_whitelist_3LT.txt',
    f'{SCDATA_DIR}/cellranger_ref/cellranger_gex_barcode_whitelist_3v2_5v1_5v2.txt',
    f'{SCDATA_DIR}/cellranger_ref/cellranger_gex_barcode_whitelist_3v3.txt',
    f'{SCDATA_DIR}/cellranger_ref/cellranger_gex_barcode_whitelist_3v4.txt',
    f'{SCDATA_DIR}/cellranger_ref/cellranger_gex_barcode_whitelist_5v3.txt',
    f'{SCDATA_DIR}/cellranger_ref/cellranger_gex_barcode_whitelist_multiome_gex.txt',
    f'{SCDATA_DIR}/cellranger_ref/cellranger_hto_barcode_whitelist_3LT.txt',
    f'{SCDATA_DIR}/cellranger_ref/cellranger_hto_barcode_whitelist_3v3.txt',
    f'{SCDATA_DIR}/cellranger_ref/cellranger_hto_barcode_whitelist_3v4.txt',
    f'{SCDATA_DIR}/cellranger_ref/cellranger_whitelist_translation_3LT.txt',
    f'{SCDATA_DIR}/cellranger_ref/cellranger_whitelist_translation_3v3.txt',
    f'{SCDATA_DIR}/cellranger_ref/cellranger_whitelist_translation_3v4.txt',
    f'{SCDATA_DIR}/cellranger_ref/cellranger_whitelists.csv',
    f'{SCDATA_DIR}/gmt_pathways/c1.all.v2023.1.Hs.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/c2.all.v2023.1.Hs.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/c2.cp.biocarta.v2023.1.Hs.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/c2.cp.kegg.v2023.1.Hs.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/c2.cp.v2023.1.Hs.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/c3.mir.v2023.1.Hs.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/c3.tft.v2023.1.Hs.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/c5.go.bp.v2023.1.Hs.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/c5.go.cc.v2023.1.Hs.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/c5.go.mf.v2023.1.Hs.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/c5.hpo.v2023.1.Hs.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/c7.immunesigdb.v2023.1.Hs.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/c8.all.v2023.1.Hs.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/h.all.v2023.1.Hs.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/m2.cp.biocarta.v2023.1.Mm.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/m2.cp.v2023.1.Mm.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/m5.go.bp.v2023.1.Mm.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/m5.go.cc.v2023.1.Mm.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/m5.go.mf.v2023.1.Mm.symbols.gmt',
    f'{SCDATA_DIR}/gmt_pathways/mh.all.v2023.1.Mm.symbols.gmt',
    f'{SCDATA_DIR}/xgboost/Siletti_Macnair-2025-07-23/allowed_cls_Siletti_Macnair_2025-07-23.csv',
    f'{SCDATA_DIR}/xgboost/Siletti_Macnair-2025-07-23/xgboost_obj_hvgs_Siletti_Macnair_2025-07-23.rds',
    # rule download_or_build_af_indices
    expand([ f'{SCDATA_DIR}/alevin_fry_home/{{genome}}/{file}' for file in AF_INDEX_FS], genome=GENOMES),
    expand(f'{SCDATA_DIR}/alevin_fry_home/{{genome}}/{{genome}}_index_params.yaml', genome=GENOMES),
    f'{SCDATA_DIR}/celltypist/celltypist_models.csv', 
    # rule get_reference_genome_data 
    f'{SCDATA_DIR}/index_parameters.csv'


# rule for getting scprocess data from github repo (maybe not a good idea to have all files as outputs)
rule download_scprocess_files:
  output:
    brain_mkrs_hsa   = f'{SCDATA_DIR}/marker_genes/human_brain.csv',
    brain_mkrs_mmu   = f'{SCDATA_DIR}/marker_genes/mouse_brain.csv',
    wl_3lt           = f'{SCDATA_DIR}/cellranger_ref/cellranger_gex_barcode_whitelist_3LT.txt',
    wl_3v2_5v1_5v2   = f'{SCDATA_DIR}/cellranger_ref/cellranger_gex_barcode_whitelist_3v2_5v1_5v2.txt',
    wl_3v3           = f'{SCDATA_DIR}/cellranger_ref/cellranger_gex_barcode_whitelist_3v3.txt',
    wl_3v4	         = f'{SCDATA_DIR}/cellranger_ref/cellranger_gex_barcode_whitelist_3v4.txt',
    wl_5v3	         = f'{SCDATA_DIR}/cellranger_ref/cellranger_gex_barcode_whitelist_5v3.txt',
    wl_multiome      = f'{SCDATA_DIR}/cellranger_ref/cellranger_gex_barcode_whitelist_multiome_gex.txt',
    wl_hto_3lt       = f'{SCDATA_DIR}/cellranger_ref/cellranger_hto_barcode_whitelist_3LT.txt',
    wl_hto_3v3       = f'{SCDATA_DIR}/cellranger_ref/cellranger_hto_barcode_whitelist_3v3.txt',
    wl_hto_3v4	     = f'{SCDATA_DIR}/cellranger_ref/cellranger_hto_barcode_whitelist_3v4.txt',
    wl_translation_3lt =  f'{SCDATA_DIR}/cellranger_ref/cellranger_whitelist_translation_3LT.txt',
    wl_translation_3v3 =  f'{SCDATA_DIR}/cellranger_ref/cellranger_whitelist_translation_3v3.txt',
    wl_translation_3v4 =  f'{SCDATA_DIR}/cellranger_ref/cellranger_whitelist_translation_3v4.txt',
    all_wl           = f'{SCDATA_DIR}/cellranger_ref/cellranger_whitelists.csv',
    gmt_f_1          = f'{SCDATA_DIR}/gmt_pathways/c1.all.v2023.1.Hs.symbols.gmt',
    gmt_f_2          = f'{SCDATA_DIR}/gmt_pathways/c2.all.v2023.1.Hs.symbols.gmt',
    gmt_f_3          = f'{SCDATA_DIR}/gmt_pathways/c2.cp.biocarta.v2023.1.Hs.symbols.gmt',
    gmt_f_4          = f'{SCDATA_DIR}/gmt_pathways/c2.cp.kegg.v2023.1.Hs.symbols.gmt',
    gmt_f_5          = f'{SCDATA_DIR}/gmt_pathways/c2.cp.v2023.1.Hs.symbols.gmt',
    gmt_f_6          = f'{SCDATA_DIR}/gmt_pathways/c3.mir.v2023.1.Hs.symbols.gmt',
    gmt_f_7          = f'{SCDATA_DIR}/gmt_pathways/c3.tft.v2023.1.Hs.symbols.gmt',
    gmt_f_8          = f'{SCDATA_DIR}/gmt_pathways/c5.go.bp.v2023.1.Hs.symbols.gmt',
    gmt_f_9          = f'{SCDATA_DIR}/gmt_pathways/c5.go.cc.v2023.1.Hs.symbols.gmt',
    gmt_f_10         = f'{SCDATA_DIR}/gmt_pathways/c5.go.mf.v2023.1.Hs.symbols.gmt',
    gmt_f_11         = f'{SCDATA_DIR}/gmt_pathways/c5.hpo.v2023.1.Hs.symbols.gmt',
    gmt_f_12         = f'{SCDATA_DIR}/gmt_pathways/c7.immunesigdb.v2023.1.Hs.symbols.gmt',
    gmt_f_13         = f'{SCDATA_DIR}/gmt_pathways/c8.all.v2023.1.Hs.symbols.gmt',
    gmt_f_14         = f'{SCDATA_DIR}/gmt_pathways/h.all.v2023.1.Hs.symbols.gmt',
    gmt_f_15         = f'{SCDATA_DIR}/gmt_pathways/m2.cp.biocarta.v2023.1.Mm.symbols.gmt',
    gmt_f_16         = f'{SCDATA_DIR}/gmt_pathways/m2.cp.v2023.1.Mm.symbols.gmt',
    gmt_f_17         = f'{SCDATA_DIR}/gmt_pathways/m5.go.bp.v2023.1.Mm.symbols.gmt',
    gmt_f_18         = f'{SCDATA_DIR}/gmt_pathways/m5.go.cc.v2023.1.Mm.symbols.gmt',
    gmt_f_19         = f'{SCDATA_DIR}/gmt_pathways/m5.go.mf.v2023.1.Mm.symbols.gmt',
    gmt_f_20         = f'{SCDATA_DIR}/gmt_pathways/mh.all.v2023.1.Mm.symbols.gmt',
    xgb_csv_f        = f'{SCDATA_DIR}/xgboost/Siletti_Macnair-2025-07-23/allowed_cls_Siletti_Macnair_2025-07-23.csv',
    xgb_rds_f        = f'{SCDATA_DIR}/xgboost/Siletti_Macnair-2025-07-23/xgboost_obj_hvgs_Siletti_Macnair_2025-07-23.rds'
  conda:
    '../envs/py_env.yaml'
  threads: 1
  shell: """
    python3 scripts/setup.py get_scprocess_data {SCDATA_DIR} {output.all_wl}
    """


# rule for downloading reference genome files from 10x and dealing with custom genomes
rule set_up_one_af_index:
  output:
    [ f'{SCDATA_DIR}/alevin_fry_home/{{genome}}/{file}' for file in AF_INDEX_FS],
    f'{SCDATA_DIR}/alevin_fry_home/{{genome}}/{{genome}}_index_params.yaml'
  params:
    fasta       = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ].get('fasta', []),
    gtf         = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ].get('gtf', []),
    index_dir   = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ].get('index_dir', None),
    mito_str    = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ].get('mito_str', []),
    is_prebuilt = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ].get('is_prebuilt', False),
    is_tenx     = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ].get('is_tenx', False),
    has_decoys  = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ].get('decoys', True),
    has_rrna    = lambda wildcards: IDX_PARAMS_LS[ wildcards.genome ].get('rrnas', True)
  conda:
    '../envs/alevin_fry.yaml'
  resources:
    mem_mb = 8192
  threads: 8
  shell: """
    python3 scripts/setup.py set_up_af_index {SCDATA_DIR} {wildcards.genome} \
      {params.fasta} {params.gtf} {params.index_dir} {params.mito_str} \
      {params.is_prebuilt} {params.is_tenx} {params.has_decoys} {params.has_rrna} {threads}
    """


rule save_index_parameters_csv:
  input:
    yamls   = expand(f'{SCDATA_DIR}/alevin_fry_home/{{genome}}/{{genome}}_index_params.yaml', genome = GENOMES)
  output:
    csv     = f'{SCDATA_DIR}/index_parameters.csv'
  conda:
    '../envs/py_env.yaml'
  resources:
    mem_mb = 512
  threads: 1
  shell: """
    python3 scripts/setup.py save_index_params_csv {output.csv} {input.yamls}
    """

    
rule download_celltypist_models:
  output:
    models_f  = f'{SCDATA_DIR}/celltypist/celltypist_models.csv'
  conda:
    '../envs/celltypist.yaml'
  threads: 1
  shell:"""
    # download celltypist models
    python3 scripts/label_celltypes.py download_models {output.models_f}
    """
