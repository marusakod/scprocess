import pandas as pd
import os
import re
import glob
import subprocess
import numpy as np
import sys
import warnings
import pathlib
import json
import shlex

scprocess_dir = pathlib.Path(config.pop('scprocess_dir'))
sys.path.append(str(scprocess_dir / 'scripts'))
from setup import get_txome_index_parameters, get_probe_set_parameters, get_xgboost_parameters, get_xgboost_source_files, CELLRANGER_RESOURCE_FILES

# get parameters
SCDATA_DIR       = os.getenv('SCPROCESS_DATA_DIR')
IDX_PARAMS_LS    = get_txome_index_parameters(config) 
PROBE_SET_LS     = get_probe_set_parameters(config)
RANGER_URL       = config['cellranger_url']
REF_TXOMES       = list(IDX_PARAMS_LS.keys())
PROBE_SET_NAMES  = list(PROBE_SET_LS.keys())
XGB_PARAMS_LS    = get_xgboost_parameters(config)
XGB_NAMES        = list(XGB_PARAMS_LS.keys())
logs_dir         = f'{SCDATA_DIR}/.log/setup'

# The manifest is the single list of marker/gene-set files supplied by Zenodo.
SUPPORT_MANIFEST_F = str(scprocess_dir / 'resources/scprocess_data.json')
SUPPORT_MANIFEST = json.loads(pathlib.Path(SUPPORT_MANIFEST_F).read_text())
SUPPORT_FILES = [str(pathlib.Path(SCDATA_DIR) / filename) for filename in SUPPORT_MANIFEST['files']]
CELLRANGER_FILES = [str(pathlib.Path(SCDATA_DIR) / 'cellranger_ref' / name) for name in CELLRANGER_RESOURCE_FILES]

# define simpleaf index files for standard txome indices
TXOME_INDEX_FS = [
  'index/piscem_idx.ctab',
  'index/piscem_idx.ectab',
  'index/piscem_idx.refinfo',
  'index/piscem_idx.ssi',
  'index/piscem_idx_cfish.json',
  'index/piscem_idx_ver.json',
  'index/simpleaf_index.json',
  'index/t2g_3col.tsv',
  'ref/gene_id_to_name.tsv',
  'ref/roers_make-ref.json',
  'ref/roers_ref.fa',
  'ref/t2g_3col.tsv',
  'simpleaf_index_log.json'
]

# define simpleaf index files for probe set (direct-ref) indices
PROBE_SET_INDEX_FS = [
  'index/piscem_idx.ctab',
  'index/piscem_idx.ectab',
  'index/piscem_idx.sigs.json',
  'index/piscem_idx.refinfo',
  'index/piscem_idx_cfish.json',
  'index/simpleaf_index.json',
  'index/t2g_3col.tsv',
  'simpleaf_index_log.json'
]

# define tiny rule
localrules: save_index_parameters_csv, add_xgboost, save_available_xgboost_csv

# define top level rule
rule all:
  input:
    SUPPORT_FILES,
    f'{SCDATA_DIR}/scprocess_support_manifest.json',
    CELLRANGER_FILES,
    # xgboost classifiers (prebuilt + user-trained)
    expand(f'{SCDATA_DIR}/xgboost/{{xgb_name}}/{{xgb_name}}_manifest.json', xgb_name=XGB_NAMES),
    expand(f'{SCDATA_DIR}/xgboost/{{xgb_name}}/{{xgb_name}}_allowed_cls.csv', xgb_name=XGB_NAMES),
    expand(f'{SCDATA_DIR}/xgboost/{{xgb_name}}/{{xgb_name}}_selected_genes.txt', xgb_name=XGB_NAMES),
    f'{SCDATA_DIR}/xgboost/xgboost_models.csv',
    f'{SCDATA_DIR}/alevin_fry_home/chemistries.json',
    # rule download_or_build_txome_indices
    expand([ f'{SCDATA_DIR}/alevin_fry_home/ref_txomes/{{ref_txome}}/{file}' for file in TXOME_INDEX_FS], ref_txome=REF_TXOMES),
    expand(f'{SCDATA_DIR}/alevin_fry_home/ref_txomes/{{ref_txome}}/{{ref_txome}}_index_params.yaml', ref_txome=REF_TXOMES),
    # rule set_up_one_probe_set_index
    expand([ f'{SCDATA_DIR}/alevin_fry_home/probe_sets/{{probe_set}}/{file}' for file in PROBE_SET_INDEX_FS], probe_set=PROBE_SET_NAMES),
    expand(f'{SCDATA_DIR}/alevin_fry_home/probe_sets/{{probe_set}}/{{probe_set}}_index_params.yaml', probe_set=PROBE_SET_NAMES),
    f'{SCDATA_DIR}/celltypist/celltypist_models.csv',
    # rule get_reference_genome_data
    f'{SCDATA_DIR}/index_parameters.csv'


# Support data and vendor files have independent downloads and dependencies.
rule download_scprocess_files:
  input:
    support_manifest = SUPPORT_MANIFEST_F,
    setup_script = str(scprocess_dir / 'scripts/setup.py')
  output:
    support_files = update(SUPPORT_FILES),
    installed_manifest = update(f'{SCDATA_DIR}/scprocess_support_manifest.json')
  conda:
    '../envs/py_env.yaml'
  threads: 1
  log:
    f'{logs_dir}/download_scprocess_files.log'
  shell: """
    exec &>> {log}
    python3 {input.setup_script:q} get_scprocess_data {SCDATA_DIR:q}
    """


rule download_cellranger_files:
  output:
    resource_files = update(CELLRANGER_FILES)
  params:
    ranger_version = config['cellranger_version'],
    ranger_url_arg = shlex.quote(RANGER_URL),
    whitelists_f = f'{SCDATA_DIR}/cellranger_ref/cellranger_whitelists.csv',
    version_f = f'{SCDATA_DIR}/cellranger_ref/.cellranger_version'
  conda:
    '../envs/py_env.yaml'
  threads: 1
  log:
    f'{logs_dir}/download_cellranger_files.log'
  shell: """
    exec &>> {log}
    python3 {scprocess_dir}/scripts/setup.py get_cellranger_data \
      {SCDATA_DIR:q} \
      {params.ranger_url_arg} \
      {params.whitelists_f:q} \
      {params.version_f:q}
    """


rule refresh_simpleaf_chems:
  output:
    chem_json = f'{SCDATA_DIR}/alevin_fry_home/chemistries.json',
  conda:
    '../envs/alevin_fry.yaml'
  threads: 1
  log:
    f'{logs_dir}/refresh_simpleaf_chems.log'
  shell: """
    exec &>> {log}

    ALEVIN_FRY_HOME="{SCDATA_DIR}/alevin_fry_home" simpleaf set-paths
    ALEVIN_FRY_HOME="{SCDATA_DIR}/alevin_fry_home" simpleaf chemistry refresh
    """


# rule for downloading reference genome files from 10x and dealing with custom genomes
rule set_up_one_txome_index:
  output:
    [ f'{SCDATA_DIR}/alevin_fry_home/ref_txomes/{{ref_txome}}/{file}' for file in TXOME_INDEX_FS],
    f'{SCDATA_DIR}/alevin_fry_home/ref_txomes/{{ref_txome}}/{{ref_txome}}_index_params.yaml'
  params:
    fasta       = lambda wildcards: IDX_PARAMS_LS[ wildcards.ref_txome ].get('fasta', []),
    gtf         = lambda wildcards: IDX_PARAMS_LS[ wildcards.ref_txome ].get('gtf', []),
    index_dir   = lambda wildcards: IDX_PARAMS_LS[ wildcards.ref_txome ].get('index_dir', None),
    mito_str    = lambda wildcards: IDX_PARAMS_LS[ wildcards.ref_txome ].get('mito_str', []),
    is_prebuilt = lambda wildcards: IDX_PARAMS_LS[ wildcards.ref_txome ].get('is_prebuilt', False),
    is_tenx     = lambda wildcards: IDX_PARAMS_LS[ wildcards.ref_txome ].get('is_tenx', False),
    has_decoys  = lambda wildcards: IDX_PARAMS_LS[ wildcards.ref_txome ].get('decoys', True),
    has_rrna    = lambda wildcards: IDX_PARAMS_LS[ wildcards.ref_txome ].get('rrnas', True)
  conda:
    '../envs/alevin_fry.yaml'
  resources:
    mem_mb = 8192
  threads: 8
  log: 
    f'{logs_dir}/set_up_one_txome_index_{{ref_txome}}.log'
  shell: """
    exec &>> {log}

    python3 {scprocess_dir}/scripts/setup.py set_up_txome_index {SCDATA_DIR} {wildcards.ref_txome} \
      {params.fasta} {params.gtf} {params.index_dir} {params.mito_str} \
      {params.is_prebuilt} {params.is_tenx} {params.has_decoys} {params.has_rrna} {threads}
    """


rule set_up_one_probe_set_index:
  output:
    [ f'{SCDATA_DIR}/alevin_fry_home/probe_sets/{{probe_set}}/{file}' for file in PROBE_SET_INDEX_FS],
    f'{SCDATA_DIR}/alevin_fry_home/probe_sets/{{probe_set}}/{{probe_set}}_index_params.yaml'
  conda:
    '../envs/alevin_fry.yaml'
  resources:
    mem_mb = 4096
  threads: 8
  log:
    f'{logs_dir}/set_up_one_probe_set_index_{{probe_set}}.log'
  shell: """
    exec &>> {log}

    python3 {scprocess_dir}/scripts/setup.py set_up_probe_set_index {SCDATA_DIR} {wildcards.probe_set} {threads}
    """


rule save_index_parameters_csv:
  input:
    yamls   = expand(f'{SCDATA_DIR}/alevin_fry_home/ref_txomes/{{ref_txome}}/{{ref_txome}}_index_params.yaml', ref_txome = REF_TXOMES) +
     expand(f'{SCDATA_DIR}/alevin_fry_home/probe_sets/{{probe_set}}/{{probe_set}}_index_params.yaml', probe_set = PROBE_SET_NAMES)
  output:
    csv     = f'{SCDATA_DIR}/index_parameters.csv'
  conda:
    '../envs/py_env.yaml'
  resources:
    mem_mb = 512
  threads: 1
  log:
    f'{logs_dir}/save_index_parameters_csv.log'
  shell: """
    exec &>> {log}

    python3 {scprocess_dir}/scripts/setup.py save_index_params_csv {output.csv} {input.yamls}
    """


rule download_celltypist_models:
  output:
    models_f  = f'{SCDATA_DIR}/celltypist/celltypist_models.csv'
  conda:
    '../envs/label_celltypes.yaml'
  threads: 1
  log:
    f'{logs_dir}/download_celltypist_models.log'
  shell:"""
    exec &>> {log}

    # download celltypist models
    python3 {scprocess_dir}/scripts/label_celltypes.py download_models {output.models_f}
    """


rule add_xgboost:
  input:
    setup_script = str(scprocess_dir / 'scripts/setup.py'),
    source_files = lambda wildcards: get_xgboost_source_files(XGB_PARAMS_LS[wildcards.xgb_name]),
    manifest = lambda wildcards: ([str(scprocess_dir / 'resources/xgboost' / f'{wildcards.xgb_name}.json')]
      if (scprocess_dir / 'resources/xgboost' / f'{wildcards.xgb_name}.json').is_file() else []),
  output:
    model_f = update(f'{SCDATA_DIR}/xgboost/{{xgb_name}}/{{xgb_name}}_xgboost_model.json'),
    cls_f   = update(f'{SCDATA_DIR}/xgboost/{{xgb_name}}/{{xgb_name}}_allowed_cls.csv'),
    genes_f = update(f'{SCDATA_DIR}/xgboost/{{xgb_name}}/{{xgb_name}}_selected_genes.txt'),
    identity_f = update(f'{SCDATA_DIR}/xgboost/{{xgb_name}}/{{xgb_name}}_manifest.json'),
  params:
    ref_tag     = lambda wildcards: XGB_PARAMS_LS[wildcards.xgb_name].get('ref_tag', wildcards.xgb_name),
    src_dir     = lambda wildcards: XGB_PARAMS_LS[wildcards.xgb_name].get('src_dir', None),
    is_prebuilt = lambda wildcards: XGB_PARAMS_LS[wildcards.xgb_name].get('is_prebuilt', False),
  threads: 1
  log:
    f'{logs_dir}/add_xgboost_{{xgb_name}}.log'
  shell: """
    exec &>> {log}

    python3 {scprocess_dir}/scripts/setup.py add_xgboost_classifier {SCDATA_DIR} {wildcards.xgb_name} {params.ref_tag} {params.src_dir} {params.is_prebuilt}
    """


rule save_available_xgboost_csv:
  input:
    setup_script = str(scprocess_dir / 'scripts/setup.py'),
    model_manifests = expand(f'{SCDATA_DIR}/xgboost/{{xgb_name}}/{{xgb_name}}_manifest.json', xgb_name=XGB_NAMES),
  output:
    csv_f = update(f'{SCDATA_DIR}/xgboost/xgboost_models.csv'),
  params:
    xgb_names = ' '.join(XGB_NAMES),
  threads: 1
  log:
    f'{logs_dir}/save_available_xgboost_csv.log'
  shell: """
    exec &>> {log}

    python3 {scprocess_dir}/scripts/setup.py save_available_xgboost_csv {output.csv_f} {SCDATA_DIR} {params.xgb_names}
    """
