# snakemake rule for running alevin-fry

import glob
import re

rule detect_chemistry:
  input: 
    wl_f  = SCPROCESS_DATA_DIR + '/cellranger_ref/cellranger_whitelists.csv'
  threads: 8
  resources:    
    mem_mb = 8192
  conda: 
    '../envs/chem_detection.yml'
  output:
    chem_f = af_dir + '/chemistry_stats.csv'
  shell:
    """ 
    mkdir -p {af_dir}  
    python3 scripts/detect_chemistry.py -f {FASTQ_DIR} -s {SAMPLE_STR} -b {input.wl_f} -o {output.chem_f} -c {threads}
    
    """
    

rule run_alevin_fry:
  input:
    chem_stats  = af_dir + '/chemistry_stats.csv',
    R1_fs      = lambda wildcards: find_fastq_files(FASTQ_DIR, wildcards.sample, "R1"),
    R2_fs      = lambda wildcards: find_fastq_files(FASTQ_DIR, wildcards.sample, "R2")
  threads: 8
  resources:
    mem_mb      = 16384
  params:
    R1_fs_str= lambda wildcards, input: ','.join(input.R1_fs),
    R2_fs_str= lambda wildcards, input: ','.join(input.R2_fs)
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

    python3 scripts/simpleaf_quant.py {wildcards.sample} {input.chem_stats} {AF_HOME_DIR} {af_dir} {AF_INDEX_DIR} {params.R1_fs_str} {params.R2_fs_str} {threads}

    """


rule save_alevin_to_h5:
  input: 
    fry_dir     = af_dir + '/af_{sample}/af_quant/'
  output: 
    h5_f        = af_dir + '/af_{sample}/af_counts_mat.h5',
    amb_yaml_f   = af_dir + '/af_{sample}/ambient_params_{sample}_' + DATE_STAMP + '.yaml',
    knee_data_f = af_dir + '/af_{sample}/knee_plot_data_{sample}_' + DATE_STAMP + '.txt.gz'
  threads: 1
  resources:
    mem_mb      = 8192
  conda: 
   '../envs/rlibs.yml'
  shell:
    """
    Rscript -e "source('scripts/alevin_fry.R'); \
      save_alevin_h5_calculate_cb_params('{wildcards.sample}', '{input.fry_dir}', \
        '{output.h5_f}', '{output.amb_yaml_f}', '{output.knee_data_f}')"
    """

