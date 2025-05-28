# snakemake rule for running alevin-fry

import glob
import re
import yaml
import pandas as pd
import os


# get alevin params


# check if custom chemistry and knees are defined for a sample
def parse_alevin_params(CUSTOM_SAMPLE_PARAMS_F, CHEMISTRY, SCPROCESS_DATA_DIR, sample):
  # set defaults
  SAMPLE_CHEMISTRY = CHEMISTRY
  
  if CUSTOM_SAMPLE_PARAMS_F is not None: 
    with open(CUSTOM_SAMPLE_PARAMS_F) as f:
      custom_smpl_params = yaml.load(f, Loader=yaml.FullLoader)
    # get all samples with custom params
      custom_smpls = list(custom_smpl_params.keys())

      valid_chems = ['3LT', '3v2', '5v1', '5v2', '3v3', 'multiome', '3v4', '5v3']

      if sample in custom_smpls:
        # check if chemistry is defined
        if 'chemistry' in custom_smpl_params[sample] and (custom_smpl_params[sample]['chemistry'] is not None):
          SAMPLE_CHEMISTRY = custom_smpl_params[sample]['chemistry']
          # check if valid
          assert SAMPLE_CHEMISTRY in valid_chems, \
            f"chemistry not valid for sample {sample}"

  # get expected ori, af chemistry and whitelist f
  if SAMPLE_CHEMISTRY in ['3v2', '5v1', '5v2']:
    AF_CHEMISTRY = '10xv2' 
  else: 
    AF_CHEMISTRY = '10xv3'

  if SAMPLE_CHEMISTRY in ['5v1', '5v2', '5v3']:
    EXPECTED_ORI = 'rc'
  else:
    EXPECTED_ORI = 'fw'

  wl_df_f = os.path.join(SCPROCESS_DATA_DIR, 'cellranger_ref/cellranger_whitelists.csv')
  wl_df = pd.read_csv(wl_df_f)
  wl_f = wl_df.loc[wl_df['chemistry'] == SAMPLE_CHEMISTRY, 'barcodes_f'].values[0]
  WHITELIST_F = os.path.join(SCPROCESS_DATA_DIR, 'cellranger_ref', wl_f)

  return AF_CHEMISTRY, EXPECTED_ORI, WHITELIST_F


def parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, sample, 
  FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD):

  # set defaults
  KNEE_1 = ""
  INF_1 = ""
  KNEE_2 = ""
  INF_2 = ""
  TOTAL_DROPLETS_INCLUDED = ""
  EXPECTED_CELLS = ""
  LOW_COUNT_THRESHOLD = ""


  # check for global cellbender params defined in config
  if AMBIENT_METHOD == 'cellbender':
    if FORCE_EXPECTED_CELLS is not None:
       EXPECTED_CELLS = FORCE_EXPECTED_CELLS
    if FORCE_TOTAL_DROPLETS_INCLUDED is not None:
       TOTAL_DROPLETS_INCLUDED = FORCE_TOTAL_DROPLETS_INCLUDED
    if FORCE_LOW_COUNT_THRESHOLD is not None:
       LOW_COUNT_THRESHOLD = FORCE_LOW_COUNT_THRESHOLD

  # check custom_sample_params_f for sample specific params
  if CUSTOM_SAMPLE_PARAMS_F is not None: 
    with open(CUSTOM_SAMPLE_PARAMS_F) as f:
      custom_smpl_params = yaml.load(f, Loader=yaml.FullLoader)
      # get all samples with custom params
      custom_smpls = list(custom_smpl_params.keys())

      amb_param_names = ['knee1', 'shin1', 'knee2', 'shin2']
      bender_param_names = ['expected_cells', 'total_droplets_included', 'low_count_threshold']

      if sample in custom_smpls:
        # check if knees and shins are defined for a sample
        if 'ambient' in custom_smpl_params[sample] and (custom_smpl_params[sample]['ambient'] is not None):
          amb_smpl_params = custom_smpl_params[sample]['ambient']
          for p_names in amb_param_names:
            assert p_names in list(amb_smpl_params.keys()), \
              f"{p_names} value missing from 'ambient' for sample {sample}"
          
          KNEE_1 = amb_smpl_params['knee1']
          INF_1 = amb_smpl_params['shin1']
          KNEE_2 = amb_smpl_params['knee2']
          INF_2 = amb_smpl_params['shin2']

        # check if specific cellbender params are defined for a sample
        if AMBIENT_METHOD == 'cellbender' and 'cellbender' in custom_smpl_params[sample] and (custom_smpl_params[sample]['cellbender'] is not None):
          if 'expected_cells' in custom_smpl_params[sample]['cellbender']:
            EXPECTED_CELLS = custom_smpl_params[sample]['cellbender']['expected_cells']
          if 'total_droplets_included' in custom_smpl_params[sample]['cellbender']:
            TOTAL_DROPLETS_INCLUDED = custom_smpl_params[sample]['cellbender']['total_droplets_included']
          if 'low_count_threshould' in custom_smpl_params[sample]['cellbender']:
            LOW_COUNT_THRESHOLD =  custom_smpl_params[sample]['cellbender']['low_count_threshold']

  return KNEE_1, INF_1, KNEE_2, INF_2, EXPECTED_CELLS, TOTAL_DROPLETS_INCLUDED, LOW_COUNT_THRESHOLD


rule run_alevin_fry:
  input:
    R1_fs      = lambda wildcards: find_fastq_files(FASTQ_DIR, wildcards.sample, "R1"),
    R2_fs      = lambda wildcards: find_fastq_files(FASTQ_DIR, wildcards.sample, "R2")
  threads: 32
  retries: RETRIES
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * MB_RUN_ALEVIN_FRY
  params:
    af_chemistry  = lambda wildcards: parse_alevin_params(CUSTOM_SAMPLE_PARAMS_F, CHEMISTRY, SCPROCESS_DATA_DIR, wildcards.sample)[0],
    exp_ori       = lambda wildcards: parse_alevin_params(CUSTOM_SAMPLE_PARAMS_F, CHEMISTRY, SCPROCESS_DATA_DIR, wildcards.sample)[1],
    whitelist_f   = lambda wildcards: parse_alevin_params(CUSTOM_SAMPLE_PARAMS_F, CHEMISTRY, SCPROCESS_DATA_DIR, wildcards.sample)[2]
  output:
    fry_dir     = directory(af_dir + '/af_{sample}/af_quant/'),
    rad_f       = temp(af_dir + '/af_{sample}/af_map/map.rad'),
    mtx_f       = af_dir + '/af_{sample}/af_quant/alevin/quants_mat.mtx',
    cols_f      = af_dir + '/af_{sample}/af_quant/alevin/quants_mat_cols.txt',
    rows_f      = af_dir + '/af_{sample}/af_quant/alevin/quants_mat_rows.txt'
  conda:
    '../envs/alevin_fry.yml'
  shell:
    """
    # mess about with input strings
    R1_fs=$(echo {input.R1_fs} | sed "s/ /,/g")
    R2_fs=$(echo {input.R2_fs} | sed "s/ /,/g")

    export ALEVIN_FRY_HOME="{AF_HOME_DIR}"
    simpleaf set-paths
  
    # simpleaf quantfication
    simpleaf quant \
      --reads1 $R1_fs  \
      --reads2 $R2_fs  \
      --threads {threads} \
      --index {AF_INDEX_DIR}/index \
      --chemistry {params.af_chemistry} --resolution cr-like \
      --expected-ori {params.exp_ori} \
      --t2g-map {AF_INDEX_DIR}/index/t2g_3col.tsv \
      --unfiltered-pl {params.whitelist_f} --min-reads 1 \
      --output {af_dir}/af_{wildcards.sample}

    """


rule save_alevin_to_h5:
  input: 
    fry_dir     = af_dir + '/af_{sample}/af_quant/'
  output: 
    h5_f        = af_dir + '/af_{sample}/af_counts_mat.h5',
    amb_yaml_f   = af_dir + '/af_{sample}/ambient_params_{sample}_' + DATE_STAMP + '.yaml',
    knee_data_f = af_dir + '/af_{sample}/knee_plot_data_{sample}_' + DATE_STAMP + '.txt.gz'
  params:
    knee1         = lambda wildcards: parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, wildcards.sample, 
      FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[0],
    inf1          = lambda wildcards: parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, wildcards.sample, 
      FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[1],
    knee2         = lambda wildcards: parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, wildcards.sample, 
      FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[2], 
    inf2          = lambda wildcards: parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, wildcards.sample, 
      FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[3],
    exp_cells     = lambda wildcards: parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, wildcards.sample, 
      FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[4], 
    total_inc     = lambda wildcards: parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, wildcards.sample, 
      FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[5], 
    low_count_thr = lambda wildcards: parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, wildcards.sample, 
      FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[6]
  threads: 1
  retries: RETRIES
  resources:
    mem_mb = lambda wildcards, attempt: attempt * MB_SAVE_ALEVIN_TO_H5
  conda: 
   '../envs/rlibs.yml'
  shell:
    """
    Rscript -e "source('scripts/alevin_fry.R'); \
      save_alevin_h5_calculate_amb_params('{wildcards.sample}', '{input.fry_dir}', \
        '{output.h5_f}', '{output.amb_yaml_f}', '{output.knee_data_f}', \
        '{params.knee1}', '{params.inf1}', '{params.knee2}', '{params.inf2}',
        '{params.exp_cells}', '{params.total_inc}', '{params.low_count_thr}')"
    """

