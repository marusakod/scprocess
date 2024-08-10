# snakemake rule for running alevin-fry

import glob
import re

rule detect_chemistry:
  input: 
    wl_f  = SCPROCESS_DATA_DIR + 'cellranger_ref/cellranger_whitelists.csv'
  threads: 8
  retries: 5
  resources:    
    mem_mb = 8192
  conda: 
    'envs/chem_detection.yml'
  output:
    chem_f = af_dir + '/chemistry_stats.csv'
  shell:
    """  
    python3 scripts/detect_chemistry.py {FASTQ_DIR} {SAMPLE_STR} {input.wl_f} {output.chem_f} {threads}
    
    """
    

rule run_alevin_fry:
  input:
    chem_stats  = af_dir + 'chemistry_stats.csv'
    R1_fs      = lambda wildcards: ",".join(find_fastq_files(fastqs_dir, wildcards.sample, "R1")),
    R2_fs      = lambda wildcards: ",".join(find_fastq_files(fastqs_dir, wildcards.sample, "R2"))
  threads: 16
  retries: 5
  resources:
    mem_mb      = 16384
  output:
    fry_dir     = directory(af_dir + '/af_{sample}/af_quant/'),
    rad_f       = temp(af_dir + '/af_{sample}/af_map/map.rad'),
    mtx_f       = af_dir + '/af_{sample}/af_quant/alevin/quants_mat.mtx',
    cols_f      = af_dir + '/af_{sample}/af_quant/alevin/quants_mat_cols.txt',
    rows_f      = af_dir + '/af_{sample}/af_quant/alevin/quants_mat_rows.txt'
  conda:
    'envs/alevin_fry.yml'
  shell:
    """
    python3 scripts/simpleaf_quant.py \
    {wildcards.sample} \
    {input.chem_stats} \
    {ALEVIN_FRY_HOME} \
    {af_dir} \
    {AF_INDEX_DIR} \
    {input.R1_fs} \
    {input.R2_fs} \{threads}

    """


rule save_alevin_to_h5:
  input: 
    fry_dir     = af_dir + '/af_{sample}/af_quant/'
  output: 
    h5_f        = af_dir + '/af_{sample}/af_counts_mat.h5',
    cb_yaml_f   = af_dir + '/af_{sample}/bender_params_{sample}_' + DATE_STAMP + '.yaml',
    knee_data_f = af_dir + '/af_{sample}/knee_plot_data_{sample}_' + DATE_STAMP + '.txt.gz'
  threads: 1
  retries: 5
  resources:
    mem_mb      = 8192
  conda: 
   'envs/rlibs.yaml'
  shell:
    """
    Rscript -e "source('scripts/alevin_fry.R'); \
      save_alevin_h5_calculate_cb_params('{wildcards.sample}', '{input.fry_dir}', \
        '{output.h5_f}', '{output.cb_yaml_f}', '{output.knee_data_f}')"
    """

