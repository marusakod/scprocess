# snakemake rule for running alevin-fry

import glob
import re
import yaml
import polars as pl
import os
from math import ceil

# mem_mb        = lambda wildcards, attempt: attempt * config['resources']['gb_run_mapping'] * MB_PER_GB
# Attempt dynamic memory based on size of R1 fastq file, but at least 32GB. Currently set to 4x size of R1 file, usually in the range of 10-15 GB.
rule run_mapping:
  params:
    demux_type    = config['multiplexing']['demux_type'],
    af_home_dir   = config['mapping']['alevin_fry_home'],
    af_index_dir  = config['mapping']['af_index_dir'],
    af_chemistry  = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["af_chemistry"],
    exp_ori       = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["expected_ori"],
    whitelist_f   = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["whitelist_f"],
    where         = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["where"],
    R1_fs         = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["R1_fs"],
    R2_fs         = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["R2_fs"]
  output:
    rad_f         = temp(af_dir + '/af_{run}/' + af_rna_dir + 'af_map/map.rad'),
    collate_rad_f = temp(af_dir + '/af_{run}/' + af_rna_dir + 'af_quant/map.collated.rad'), 
    fry_dir       = directory(af_dir + '/af_{run}/' +  af_rna_dir + 'af_quant/'),
    mtx_f         = af_dir + '/af_{run}/'  + af_rna_dir + 'af_quant/alevin/quants_mat.mtx',
    cols_f        = af_dir + '/af_{run}/' + af_rna_dir +'af_quant/alevin/quants_mat_cols.txt',
    rows_f        = af_dir + '/af_{run}/' + af_rna_dir +'af_quant/alevin/quants_mat_rows.txt'
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_mapping/run_mapping_{run}_' + DATE_STAMP + '.benchmark.txt'
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb        = lambda wildcards: max(ceil(config['resources']['gb_run_mapping_per_gb_fq'] * RUN_PARAMS[wildcards.run]["mapping"]["R1_fs_size_gb"] * MB_PER_GB), 32 * MB_PER_GB)
  conda:
    '../envs/alevin_fry.yaml'
  shell:"""
    # check whether doing arvados
    ARV_REGEX="^arkau-[0-9a-z]{{5}}-[0-9a-z]{{15}}$"
    if [[ "{params.where}" =~ $ARV_REGEX ]]; then
      ml arvados
      arv-env arkau
    fi
    # run mapping
    python3 scripts/mapping.py {wildcards.run} \
      --af_dir          {af_dir} \
      --demux_type      {params.demux_type} \
      --what            "rna" \
      --af_home_dir     {params.af_home_dir} \
      --where           {params.where} \
      --R1_fs           {params.R1_fs} \
      --R2_fs           {params.R2_fs} \
      --threads         {threads} \
      --af_index_dir    {params.af_index_dir} \
      --tenx_chemistry  {params.af_chemistry} \
      --exp_ori         {params.exp_ori} \
      --whitelist_f     {params.whitelist_f}
    """


rule save_alevin_to_h5:
  input: 
    fry_dir     = af_dir + '/af_{run}/' + af_rna_dir + 'af_quant/'
  output: 
    h5_f        = af_dir + '/af_{run}/' + af_rna_dir + 'af_counts_mat.h5',
    amb_yaml_f  = af_dir + '/af_{run}/' + af_rna_dir + 'ambient_params_{run}_' + DATE_STAMP + '.yaml',
    knee_data_f = af_dir + '/af_{run}/' + af_rna_dir + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz'
  params:
    knee1         = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["knee1"],
    shin1         = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["shin1"],
    knee2         = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["knee2"],
    shin2         = lambda wildcards: RUN_PARAMS[wildcards.run]["mapping"]["shin2"],
    exp_cells     = lambda wildcards: RUN_PARAMS[wildcards.run]["ambient"]["cb_expected_cells"],
    total_inc     = lambda wildcards: RUN_PARAMS[wildcards.run]["ambient"]["cb_total_droplets_included"],
    low_count_thr = lambda wildcards: RUN_PARAMS[wildcards.run]["ambient"]["cb_low_count_threshold"]
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb = lambda wildcards, attempt: attempt * config['resources']['gb_save_alevin_to_h5'] * MB_PER_GB
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_mapping/save_alevin_to_h5_{run}_' + DATE_STAMP + '.benchmark.txt'
  conda: 
   '../envs/rlibs.yaml'
  shell: """
    Rscript -e "source('scripts/mapping.R');
      save_alevin_h5_ambient_params(
        run           = '{wildcards.run}',
        fry_dir       = '{input.fry_dir}',
        h5_f          = '{output.h5_f}',
        cb_yaml_f     = '{output.amb_yaml_f}',
        knee_data_f   = '{output.knee_data_f}',
        run_var       = '{RUN_VAR}',
        knee1         = '{params.knee1}',
        shin1         = '{params.shin1}',
        knee2         = '{params.knee2}',
        shin2         = '{params.shin2}',
        exp_cells     = '{params.exp_cells}',
        total_included= '{params.total_inc}',
        low_count_thr = '{params.low_count_thr}'
      )"
    """

