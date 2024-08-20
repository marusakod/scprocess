import pandas as pd
import os
import re
import glob
import subprocess
import numpy as np
import sys
import warnings


from setup_utils import get_setup_parameters

SCPROCESS_DATA_DIR = os.getenv('SCPROCESS_DATA_DIR')

#configfile = '/Users/marusa/Projects/scprocess_test/configs/config-setup-test.yaml'
# get parameters from configfile
GENOMES_STR, FASTA_FS, GTF_FS, MITO_STRS, DECOYS = get_setup_parameters(config) 
GENOMES = GENOMES_STR.split(',')

rule all:
    input:
        # rule download scprocess repo data
        SCPROCESS_DATA_DIR + '/marker_genes/canonical_brain_celltype_markers_human_2023-10-17.txt',
        SCPROCESS_DATA_DIR + '/marker_genes/canonical_brain_celltype_markers_mouse_2023-10-17.txt',
        SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3LT.txt',
        SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3v1.txt',
        SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3v2_5v1_5v2.txt',
        SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3v3.txt',
        SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_flex.txt',
        SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_multiome_gex.txt',
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
        # rule get_reference_genomes 
        SCPROCESS_DATA_DIR + '/setup_parameters.csv', 
         # rule build_af_indices
        expand(SCPROCESS_DATA_DIR + '/alevin_fry_home/{genome}/simpleaf_index_log.json', genome=GENOMES)




# rule for getting scprocess data from github repo (maybe not a good idea to have all files as outputs)
rule download_scprocess_repo_data:
  output:
    brain_mkrs_hsa    = SCPROCESS_DATA_DIR + '/marker_genes/canonical_brain_celltype_markers_human_2023-10-17.txt',
    brain_mkrs_mmu    = SCPROCESS_DATA_DIR + '/marker_genes/canonical_brain_celltype_markers_mouse_2023-10-17.txt',
    wl_3lt            = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3LT.txt',
    wl_3v1            = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3v1.txt',
    wl_3v2_5v1_5v2    = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3v2_5v1_5v2.txt',
    wl_3v3            = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_3v3.txt',
    wl_flex           = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_flex.txt',
    wl_multiome       = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_barcode_whitelist_multiome_gex.txt',
    all_wl            = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_whitelists.csv',
    gmt_f_1           = SCPROCESS_DATA_DIR + '/gmt_pathways/c1.all.v2023.1.Hs.symbols.gmt',
    gmt_f_2           = SCPROCESS_DATA_DIR + '/gmt_pathways/c2.all.v2023.1.Hs.symbols.gmt',
    gmt_f_3           = SCPROCESS_DATA_DIR + '/gmt_pathways/c2.cp.biocarta.v2023.1.Hs.symbols.gmt',
    gmt_f_4           = SCPROCESS_DATA_DIR + '/gmt_pathways/c2.cp.kegg.v2023.1.Hs.symbols.gmt',
    gmt_f_5           = SCPROCESS_DATA_DIR + '/gmt_pathways/c2.cp.v2023.1.Hs.symbols.gmt',
    gmt_f_6           = SCPROCESS_DATA_DIR + '/gmt_pathways/c3.mir.v2023.1.Hs.symbols.gmt',
    gmt_f_7           = SCPROCESS_DATA_DIR + '/gmt_pathways/c3.tft.v2023.1.Hs.symbols.gmt',
    gmt_f_8           = SCPROCESS_DATA_DIR + '/gmt_pathways/c5.go.bp.v2023.1.Hs.symbols.gmt',
    gmt_f_9           = SCPROCESS_DATA_DIR + '/gmt_pathways/c5.go.cc.v2023.1.Hs.symbols.gmt',
    gmt_f_10          = SCPROCESS_DATA_DIR + '/gmt_pathways/c5.go.mf.v2023.1.Hs.symbols.gmt',
    gmt_f_11          = SCPROCESS_DATA_DIR + '/gmt_pathways/c5.hpo.v2023.1.Hs.symbols.gmt',
    gmt_f_12          = SCPROCESS_DATA_DIR + '/gmt_pathways/c7.immunesigdb.v2023.1.Hs.symbols.gmt',
    gmt_f_13          = SCPROCESS_DATA_DIR + '/gmt_pathways/c8.all.v2023.1.Hs.symbols.gmt',
    gmt_f_14          = SCPROCESS_DATA_DIR + '/gmt_pathways/h.all.v2023.1.Hs.symbols.gmt',
    gmt_f_15          = SCPROCESS_DATA_DIR + '/gmt_pathways/m2.cp.biocarta.v2023.1.Mm.symbols.gmt',
    gmt_f_16          = SCPROCESS_DATA_DIR + '/gmt_pathways/m2.cp.v2023.1.Mm.symbols.gmt',
    gmt_f_17          = SCPROCESS_DATA_DIR + '/gmt_pathways/m5.go.bp.v2023.1.Mm.symbols.gmt',
    gmt_f_18          = SCPROCESS_DATA_DIR + '/gmt_pathways/m5.go.cc.v2023.1.Mm.symbols.gmt',
    gmt_f_19          = SCPROCESS_DATA_DIR + '/gmt_pathways/m5.go.mf.v2023.1.Mm.symbols.gmt',
    gmt_f_20          = SCPROCESS_DATA_DIR + '/gmt_pathways/mh.all.v2023.1.Mm.symbols.gmt',
    xgb_csv_f         = SCPROCESS_DATA_DIR + '/xgboost/Siletti_Macnair-2024-03-11/allowed_cls_Siletti_Macnair_2024-03-11.csv',
    xgb_rds_f         = SCPROCESS_DATA_DIR + '/xgboost/Siletti_Macnair-2024-03-11/xgboost_obj_hvgs_Siletti_Macnair_2024-03-11.rds'
  conda:
    'envs/py_env.yml'
  threads: 1
  shell:
    """
    python3 scripts/setup.py get_scprocess_data {SCPROCESS_DATA_DIR}

    """


# rule for downloading reference genome files from 10x and dealing with custom genomes
rule get_reference_genomes:
  output:
    ref_params_f = SCPROCESS_DATA_DIR + '/setup_parameters.csv'
  conda:
    'envs/py_env.yml'
  resources:
    mem_mb = 8192
  threads: 1
  shell:
    """  
    python3 scripts/setup.py get_genome_params {GENOMES_STR} {FASTA_FS} {GTF_FS} {MITO_STRS} {DECOYS} {SCPROCESS_DATA_DIR}
    
    """

rule build_af_indices:
  input: 
    ref_params_f = SCPROCESS_DATA_DIR + '/setup_parameters.csv'
  output:
    simpleaf_log_f = SCPROCESS_DATA_DIR + '/alevin_fry_home/{genome}/simpleaf_index_log.json'
  conda:
    'envs/alevin_fry.yml'
  resources:
    mem_mb = 16384
  threads: 16
  shell:
    """
    python3 scripts/build_af_index.py {wildcards.genome} {input.ref_params_f} {SCPROCESS_DATA_DIR} {threads}

    """
