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

  # sort out custom sample parameters
  if CUSTOM_SAMPLE_PARAMS_F is not None: 
    with open(CUSTOM_SAMPLE_PARAMS_F) as f:
      # get all samples with custom params
      custom_smpl_params = yaml.load(f, Loader=yaml.FullLoader)
      custom_smpls = list(custom_smpl_params.keys())

      # check each sample
      valid_chems = ['3LT', '3v2', '3v3', '3v4', '5v1', '5v2', '5v3', 'multiome']
      if sample in custom_smpls:
        # check if chemistry is defined
        if 'chemistry' in custom_smpl_params[sample] and (custom_smpl_params[sample]['chemistry'] is not None):
          SAMPLE_CHEMISTRY = custom_smpl_params[sample]['chemistry']
          # check if valid
          assert SAMPLE_CHEMISTRY in valid_chems, \
            f"chemistry not valid for sample {sample}"

  # get af chemistry and expected orientation
  if SAMPLE_CHEMISTRY in ['3v2', '5v1', '5v2']:
    AF_CHEMISTRY = '10xv2' 
  else: 
    AF_CHEMISTRY = '10xv3'
  if SAMPLE_CHEMISTRY in ['5v1', '5v2', '5v3']:
    EXPECTED_ORI = 'rc'
  else:
    EXPECTED_ORI = 'fw'

  # sort out whitelist file
  wl_df_f = os.path.join(SCPROCESS_DATA_DIR, 'cellranger_ref/cellranger_whitelists.csv')
  wl_df = pd.read_csv(wl_df_f)
  wl_f = wl_df.loc[wl_df['chemistry'] == SAMPLE_CHEMISTRY, 'barcodes_f'].values[0]
  wl_trans_f = wl_df.loc[wl_df['chemistry'] == SAMPLE_CHEMISTRY, 'translation_f'].values[0]
  if type(wl_trans_f) == str:
    WHITELIST_TRANS_F = os.path.join(SCPROCESS_DATA_DIR, 'cellranger_ref', wl_trans_f)
  else:
    WHITELIST_TRANS_F = None
  WHITELIST_F = os.path.join(SCPROCESS_DATA_DIR, 'cellranger_ref', wl_f)

  return AF_CHEMISTRY, EXPECTED_ORI, WHITELIST_F, WHITELIST_TRANS_F


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


# build hto index if at least one multiplexed sample with htos
if DEMUX_TYPE == "hto":
  rule build_hto_index:
    input:
      feature_ref_f = FEATURE_REF
    threads: 8
    retries: RETRIES
    resources:
      mem_mb = 4096
    output: 
      hto_f       = af_dir + '/hto.tsv',
      t2g_f       = af_dir + '/t2g_hto.tsv',
      idx_log_f   = af_dir + '/hto_index/ref_indexing.log'
    conda:
      '../envs/alevin_fry.yaml'
    shell:
      """
      cd {af_dir}
      # create a tsv file from feature ref file
      awk -F "," '
        NR == 1 {{
          for (i = 1; i <= NF; i++) {{
            if ($i == "hto_id") hto_id_col = i
            if ($i == "sequence") seq_col = i
          }}
        }}
        NR > 1 {{
          print $hto_id_col "\t" $seq_col
        }}' {input.feature_ref_f} > {output.hto_f}

      # 
      salmon index -t {output.hto_f} -i hto_index --features -k7
      
      # Make gene-to-transcript mapping file
      awk '{{print $1"\t"$1;}}' {output.hto_f} > {output.t2g_f}
      """


rule run_mapping:
  threads: 8
  retries: RETRIES
  resources:
    mem_mb        = lambda wildcards, attempt: attempt * MB_RUN_MAPPING
  params:
    af_chemistry  = lambda wildcards: parse_alevin_params(CUSTOM_SAMPLE_PARAMS_F, CHEMISTRY, SCPROCESS_DATA_DIR, wildcards.run)[0],
    exp_ori       = lambda wildcards: parse_alevin_params(CUSTOM_SAMPLE_PARAMS_F, CHEMISTRY, SCPROCESS_DATA_DIR, wildcards.run)[1],
    whitelist_f   = lambda wildcards: parse_alevin_params(CUSTOM_SAMPLE_PARAMS_F, CHEMISTRY, SCPROCESS_DATA_DIR, wildcards.run)[2],
    where         = lambda wildcards: SAMPLE_FQS[wildcards.run]["where"],
    R1_fs         = lambda wildcards: SAMPLE_FQS[wildcards.run]["R1_fs"],
    R2_fs         = lambda wildcards: SAMPLE_FQS[wildcards.run]["R2_fs"]
  output:
    rad_f         = temp(af_dir + '/af_{run}/' + af_rna_dir + 'af_map/map.rad'),
    collate_rad_f = temp(af_dir + '/af_{run}/' + af_rna_dir + 'af_quant/map.collated.rad'), 
    fry_dir       = directory(af_dir + '/af_{run}/' +  af_rna_dir + 'af_quant/'),
    mtx_f         = af_dir + '/af_{run}/'  + af_rna_dir + 'af_quant/alevin/quants_mat.mtx',
    cols_f        = af_dir + '/af_{run}/' + af_rna_dir +'af_quant/alevin/quants_mat_cols.txt',
    rows_f        = af_dir + '/af_{run}/' + af_rna_dir +'af_quant/alevin/quants_mat_rows.txt'
  conda:
    '../envs/alevin_fry.yaml'
  shell:"""
    python3 scripts/mapping.py {wildcards.run} \
      --af_dir {af_dir} \
      --demux_type {DEMUX_TYPE} \
      --what rna \
      --af_home_dir {AF_HOME_DIR} \
      --where {params.where} \
      --R1_fs {params.R1_fs} \
      --R2_fs {params.R2_fs} \
      --threads {threads} \
      --af_index_dir {AF_INDEX_DIR} \
      --tenx_chemistry {params.af_chemistry} \
      --exp_ori {params.exp_ori} \
      --whitelist_f {params.whitelist_f}
    """


if DEMUX_TYPE == "hto":
  rule run_mapping_hto:
    input:
      hto_R1_fs     = lambda wildcards: find_fastq_files(HTO_FASTQ_DIR, wildcards.run, "R1"), 
      hto_R2_fs     = lambda wildcards: find_fastq_files(HTO_FASTQ_DIR, wildcards.run, "R2"),
      t2g_f         = af_dir + '/t2g_hto.tsv'
    threads: 8
    retries: RETRIES
    resources:
      mem_mb      = lambda wildcards, attempt: attempt * MB_RUN_MAPPING
    params:
      af_chemistry  = lambda wildcards: parse_alevin_params(CUSTOM_SAMPLE_PARAMS_F, CHEMISTRY, SCPROCESS_DATA_DIR, wildcards.run)[0],
      whitelist_f   = lambda wildcards: parse_alevin_params(CUSTOM_SAMPLE_PARAMS_F, CHEMISTRY, SCPROCESS_DATA_DIR, wildcards.run)[2]
    output:
      fry_dir     = directory(af_dir + '/af_{run}/hto/af_quant/'),
      rad_f       = temp(af_dir + '/af_{run}/hto/af_map/map.rad'),
      mtx_f       = af_dir + '/af_{run}/hto/af_quant/alevin/quants_mat.mtx',
      cols_f      = af_dir + '/af_{run}/hto/af_quant/alevin/quants_mat_cols.txt',
      rows_f      = af_dir + '/af_{run}/hto/af_quant/alevin/quants_mat_rows.txt'
    conda:
      '../envs/alevin_fry.yaml'
    shell:"""
      python3 scripts/mapping.py {wildcards.run} \
        --af_dir {af_dir}\
        --demux_type {DEMUX_TYPE} \
        --what hto \
        --af_home_dir {AF_HOME_DIR} \
        --where {params.fastqs["where"]} \
        --R1_fs {params.fastqs["R1_fs"]} \
        --R2_fs {params.fastqs["R2_fs"]} \
        --threads {threads} \
        --af_index_dir {AF_INDEX_DIR} \
        --tenx_chemistry {params.af_chemistry} \
        --exp_ori fw \
        --whitelist_f {params.whitelist_f}
      """


rule save_alevin_to_h5:
  input: 
    fry_dir     = af_dir + '/af_{run}/' + af_rna_dir + 'af_quant/'
  output: 
    h5_f        = af_dir + '/af_{run}/' + af_rna_dir + 'af_counts_mat.h5',
    amb_yaml_f  = af_dir + '/af_{run}/' + af_rna_dir + 'ambient_params_{run}_' + DATE_STAMP + '.yaml',
    knee_data_f = af_dir + '/af_{run}/' + af_rna_dir + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz'
  params:
    knee1         = lambda wildcards: parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, wildcards.run, 
      FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[0],
    inf1          = lambda wildcards: parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, wildcards.run, 
      FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[1],
    knee2         = lambda wildcards: parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, wildcards.run, 
      FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[2], 
    inf2          = lambda wildcards: parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, wildcards.run, 
      FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[3],
    exp_cells     = lambda wildcards: parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, wildcards.run, 
      FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[4], 
    total_inc     = lambda wildcards: parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, wildcards.run, 
      FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[5], 
    low_count_thr = lambda wildcards: parse_knee_finder_params(CUSTOM_SAMPLE_PARAMS_F, AMBIENT_METHOD, wildcards.run, 
      FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[6]
  threads: 1
  retries: RETRIES
  resources:
    mem_mb = lambda wildcards, attempt: attempt * MB_SAVE_ALEVIN_TO_H5
  conda: 
   '../envs/rlibs.yaml'
  shell:
    """
    Rscript -e "source('scripts/mapping.R');
      save_alevin_h5_ambient_params(
        sample        = '{wildcards.run}',
        fry_dir       = '{input.fry_dir}',
        h5_f          = '{output.h5_f}',
        cb_yaml_f     = '{output.amb_yaml_f}',
        knee_data_f   = '{output.knee_data_f}',
        sample_var    = '{SAMPLE_VAR}',
        knee1         = '{params.knee1}',
        inf1          = '{params.inf1}',
        knee2         = '{params.knee2}',
        inf2          = '{params.inf2}',
        exp_cells     = '{params.exp_cells}',
        total_included= '{params.total_inc}',
        low_count_thr = '{params.low_count_thr}')"
    """


if DEMUX_TYPE == "hto":
  rule save_alevin_hto_to_h5:
    input: 
      fry_dir     = af_dir + '/af_{run}/hto/af_quant/'
    output: 
      h5_f        = af_dir + '/af_{run}/hto/af_hto_counts_mat.h5',
      knee_data_f = af_dir + '/af_{run}/hto/knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz'
    threads: 1
    retries: RETRIES
    resources:
      mem_mb = lambda wildcards, attempt: attempt * MB_SAVE_ALEVIN_TO_H5
    conda: 
      '../envs/rlibs.yaml'
    shell:
      """
      Rscript -e "source('scripts/mapping.R');
        save_alevin_h5_knee_params_df(
          sample      = '{wildcards.run}', 
          fry_dir     = '{input.fry_dir}', 
          hto_mat     = 1, 
          sample_var  = '{SAMPLE_VAR}',
          h5_f        = '{output.h5_f}',
          knee_data_f = '{output.knee_data_f}')"
      """
