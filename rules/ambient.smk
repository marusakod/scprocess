# snakemake rule for running cellbender

import numpy as np
import yaml
import pandas as pd

localrules: get_ambient_sample_statistics

def parse_ambient_params(AMBIENT_METHOD, CUSTOM_SAMPLE_PARAMS_F, sample, amb_yaml_f, CELLBENDER_LEARNING_RATE):
  # get cellbender parameters from yaml file
  with open(amb_yaml_f) as f:
    amb_params = yaml.load(f, Loader=yaml.FullLoader)
      
  # get parameters for cellbender
  EXPECTED_CELLS = amb_params['expected_cells']
  TOTAL_DROPLETS_INCLUDED = amb_params['total_droplets_included']
  LOW_COUNT_THRESHOLD = amb_params['low_count_threshold']
  KNEE_1 = amb_params['knee_1']
  INFLECTION_1 = amb_params['inflection_1']
  KNEE_2 = amb_params['knee_2']
  INFLECTION_2 = amb_params['inflection_2']
  LEARNING_RATE = CELLBENDER_LEARNING_RATE

  # check if cellbender specific parameters are defined for individual samples 
  if AMBIENT_METHOD == 'cellbender':
    if CUSTOM_SAMPLE_PARAMS_F is not None:
      with open(CUSTOM_SAMPLE_PARAMS_F) as f:
        custom_smpl_params = yaml.load(f, Loader=yaml.FullLoader)
      # get all samples with custom params
      custom_smpls = list(custom_smpl_params.keys())

      if sample in custom_smpls:
        # check if cellbender is defined
        if 'cellbender' in custom_smpl_params[sample] and (custom_smpl_params[sample]['cellbender'] is not None):
          if 'learning_rate' in custom_smpl_params[sample]['cellbender']: 
            LEARNING_RATE = custom_smpl_params[sample]['cellbender']['learning_rate']

  return EXPECTED_CELLS, TOTAL_DROPLETS_INCLUDED, LOW_COUNT_THRESHOLD, LEARNING_RATE, \
    KNEE_1, INFLECTION_1, KNEE_2, INFLECTION_2




# metrics_fs_ls is the list of knee files (this should exist for all ambient methods)
# this function should run in the last rule
def extract_ambient_sample_statistics(AMBIENT_METHOD, SAMPLE_VAR, samples_ls, metrics_fs_ls, ambient_outs_yamls, custom_f, max_kept=0.9):

  # loop through samples
  kept_arr    = []
  totals_arr  = []
  for sample, metrics_f, ambient_outs_yaml in zip(samples_ls, metrics_fs_ls, ambient_outs_yamls):
    # Load ambient outs yaml file
    with open(ambient_outs_yaml) as f:
      amb_outs = yaml.load(f, Loader=yaml.FullLoader)

    bc_f = amb_outs['bcs_f']
   
    # count the number of barcodes
    barcode_count = pd.read_csv(bc_f, header=None).shape[0]
    kept_arr.append(barcode_count)

    if AMBIENT_METHOD == 'cellbender':
      # get the number of total droplets included
      total_droplets = pd.read_csv(metrics_f)['total_droplets_included'][0]
      totals_arr.append(total_droplets)

  kept_arr = np.array(kept_arr)

  if AMBIENT_METHOD != 'cellbender':
    sample_df = pd.DataFrame({
      SAMPLE_VAR : samples_ls,
      'kept_droplets': kept_arr
    })
  else:
    totals_arr = np.array(totals_arr)

    # replace dodgy totals values with custom if need be
    if custom_f is not None and os.path.isfile(custom_f):
      # open custom yaml file 
      with open(custom_f) as f:
        all_params  = yaml.load(f, Loader=yaml.FullLoader)

      # get all samples with custom params
      custom_ls   = [ s for s in all_params.keys() if 'ambient' in all_params[s] ]
      amb_params  = [ all_params[s]['ambient'] for s in custom_ls ]
      custom_df   = pd.DataFrame(columns = [SAMPLE_VAR, 'total_droplets_included'])
      for i, s in enumerate(custom_ls):
        this_params = amb_params[ i ]
        if 'total_droplets_included' in this_params:
          tmp_df      = pd.DataFrame({SAMPLE_VAR: s, 'total_droplets_included': this_params['total_droplets_included']})
          custom_df   = custom_df.append( tmp_df )
      
      # use these to edit totals_arr
      samples_arr = np.array(samples_ls)
      for idx, row in custom_df.iterrows():
        match_idx = np.where(samples_arr == row[SAMPLE_VAR])
        totals_arr[match_idx] = row['total_droplets_included']

    # do some calculations
    prop_kept = kept_arr / totals_arr
    bad_idx = prop_kept > max_kept

    # assemble into dataframe
    sample_df = pd.DataFrame({
      SAMPLE_VAR : samples_ls,
      'total_droplets': totals_arr,
      'kept_droplets': kept_arr,
      'prop_kept_by_cb': prop_kept,
      'bad_sample': bad_idx
    })

  return sample_df


if AMBIENT_METHOD == 'cellbender':
  rule run_cellbender:
    input:
      h5_f      = af_dir + '/af_{run}/' + af_rna_dir  + 'af_counts_mat.h5',
      amb_yaml_f = af_dir + '/af_{run}/' + af_rna_dir + 'ambient_params_{run}_' + DATE_STAMP + '.yaml'
    params:
      expected_cells          = lambda wildcards: parse_ambient_params(AMBIENT_METHOD, CUSTOM_SAMPLE_PARAMS_F, wildcards.run,
        af_dir + f'/af_{wildcards.run}/' + af_rna_dir + f'ambient_params_{wildcards.run}_{DATE_STAMP}.yaml', CELLBENDER_LEARNING_RATE)[0],
      total_droplets_included = lambda wildcards: parse_ambient_params(AMBIENT_METHOD, CUSTOM_SAMPLE_PARAMS_F, wildcards.run,
        af_dir + f'/af_{wildcards.run}/' + af_rna_dir + f'ambient_params_{wildcards.run}_{DATE_STAMP}.yaml', CELLBENDER_LEARNING_RATE)[1],
      low_count_threshold     = lambda wildcards: parse_ambient_params(AMBIENT_METHOD, CUSTOM_SAMPLE_PARAMS_F, wildcards.run,
        af_dir + f'/af_{wildcards.run}/' + af_rna_dir + f'ambient_params_{wildcards.run}_{DATE_STAMP}.yaml', CELLBENDER_LEARNING_RATE)[2],
      learning_rate           = lambda wildcards: parse_ambient_params(AMBIENT_METHOD, CUSTOM_SAMPLE_PARAMS_F, wildcards.run,
        af_dir + f'/af_{wildcards.run}/' + af_rna_dir + f'ambient_params_{wildcards.run}_{DATE_STAMP}.yaml', CELLBENDER_LEARNING_RATE)[3]
    output:
        ambient_yaml_out = amb_dir + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml',
        tmp_f            = temp(amb_dir + '/ambient_{run}/ckpt.tar.gz')
    threads: 1
    retries: RETRIES
    resources:
      mem_mb      = lambda wildcards, attempt: attempt * MB_RUN_AMBIENT
    singularity:
      CELLBENDER_IMAGE
    shell:
      """
      # get parameters for cellbender
      EXPECTED_CELLS={params.expected_cells}
      TOTAL_DROPLETS_INCLUDED={params.total_droplets_included}
      LOW_COUNT_THRESHOLD={params.low_count_threshold}
      LEARNING_RATE={params.learning_rate}

      # create main ambient directory
      mkdir -p {amb_dir}

      # change to sample ambient directory
      amb_dir=$(dirname {output.ambient_yaml_out})
      mkdir -p $amb_dir
      cd $amb_dir

      # define output files
      raw_counts_f="{amb_dir}/ambient_{wildcards.run}/bender_{wildcards.run}_{DATE_STAMP}.h5"
      filt_counts_f="{amb_dir}/ambient_{wildcards.run}/bender_{wildcards.run}_{DATE_STAMP}_filtered.h5"
      bcs_f="{amb_dir}/ambient_{wildcards.run}/bender_{wildcards.run}_{DATE_STAMP}_cell_barcodes.csv"
      tmp_f="{output.tmp_f}"

      
      if [ $EXPECTED_CELLS -eq 0 ]; then
        echo "running cellbender without expected_cells"
        # run cellbender
        cellbender remove-background \
          --input {input.h5_f} \
          --output $cb_full_f \
          --total-droplets-included $TOTAL_DROPLETS_INCLUDED \
          --low-count-threshold $LOW_COUNT_THRESHOLD \
          --learning-rate $LEARNING_RATE \
          --checkpoint-mins 30.0 \
          --cuda
      else
        echo "running cellbender with expected_cells"
        # run cellbender
        cellbender remove-background \
          --input {input.h5_f} \
          --output $raw_counts_f \
          --expected-cells $EXPECTED_CELLS \
          --total-droplets-included $TOTAL_DROPLETS_INCLUDED \
          --low-count-threshold $LOW_COUNT_THRESHOLD \
          --learning-rate $LEARNING_RATE \
          --checkpoint-mins 30.0 \
          --cuda
      fi

      # Create the output yaml file
      echo "raw_counts_f: $raw_counts_f" >> {output.ambient_yaml_out}
      echo "filt_counts_f: $filt_counts_f" >> {output.ambient_yaml_out}
      echo "bcs_f: $bcs_f" >> {output.ambient_yaml_out}
      #echo "tmp_f: $tmp_f" >> {output.ambient_yaml_out}

      # check whether temp file was actually made; if not, make an empty one
      if [ ! -f $tmp_f ]; then
        touch $tmp_f
      fi
      """

if AMBIENT_METHOD == 'decontx':
  rule run_decontx:
    input:
      h5_f        = af_dir + '/af_{run}/' + af_rna_dir + 'af_counts_mat.h5',
      amb_yaml_f  = af_dir + '/af_{run}/' + af_rna_dir + 'ambient_params_{run}_' + DATE_STAMP + '.yaml', 
      knee_data_f = af_dir + '/af_{run}/' + af_rna_dir + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz'
    output:
      ambient_yaml_out = amb_dir + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml'
    threads: 4
    retries: RETRIES
    resources:
      mem_mb    = lambda wildcards, attempt: attempt * MB_RUN_AMBIENT
    conda: 
      '../envs/rlibs.yml'
    shell:
      """
      # create main ambient directory
      mkdir -p {amb_dir}

      # change to sample ambient directory
      amb_dir=$(dirname {output.ambient_yaml_out})
      mkdir -p $amb_dir
  
      # define output file names
      filt_counts_f="{amb_dir}/ambient_{wildcards.run}/decontx_{wildcards.run}_{DATE_STAMP}_filtered.h5"
      bcs_f="{amb_dir}/ambient_{wildcards.run}/decontx_{wildcards.run}_{DATE_STAMP}_cell_barcodes.csv"
      dcx_params_f="{amb_dir}/ambient_{wildcards.run}/decontx_{wildcards.run}_{DATE_STAMP}_params.txt.gz"

      # run cell calling and decontamination
   
      Rscript -e "source('scripts/ambient.R'); 
      get_cell_mat_and_barcodes(
      out_mat_f = '$filt_counts_f',
      out_bcs_f = '$bcs_f',
      out_dcx_f = '$dcx_params_f',
      sel_s     = '{wildcards.run}',
      af_mat_f  = '{input.h5_f}',
      knee_f    = '{input.knee_data_f}',
      ncores    = {threads},
      cell_calls_method = '{CELL_CALLS_METHOD}',
      ambient_method    = '{AMBIENT_METHOD}')"

      # Create the output yaml file
      echo "filt_counts_f: $filt_counts_f" >> {output.ambient_yaml_out}
      echo "bcs_f: $bcs_f" >> {output.ambient_yaml_out}
      echo "dcx_params_f: $dcx_params_f" >> {output.ambient_yaml_out}

      """


if AMBIENT_METHOD == 'none':
  rule run_cell_calling:
    input:
      h5_f        = af_dir + '/af_{run}/' + af_rna_dir + 'af_counts_mat.h5',
      amb_yaml_f  = af_dir + '/af_{run}/' + af_rna_dir + 'ambient_params_{run}_' + DATE_STAMP + '.yaml',
      knee_data_f = af_dir + '/af_{run}/' + af_rna_dir + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz'
    output:
      ambient_yaml_out = amb_dir + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml'
    threads: 4
    retries: RETRIES
    conda:
      '../envs/rlibs.yml'
    resources:
      mem_mb    = lambda wildcards, attempt: attempt * MB_RUN_AMBIENT
    shell:
      """
      # create main ambient directory
      mkdir -p {amb_dir}

      # change to sample ambient directory
      amb_dir=$(dirname {output.ambient_yaml_out})
      mkdir -p $amb_dir
      
      # define output file names
      filt_counts_f="{amb_dir}/ambient_{wildcards.run}/uncorrected_{wildcards.run}_{DATE_STAMP}_filtered.h5"
      bcs_f="{amb_dir}/ambient_{wildcards.run}/uncorrected_{wildcards.run}_{DATE_STAMP}_cell_barcodes.csv"

      # run cell calling and decontamination
      Rscript -e "source('scripts/ambient.R');
      get_cell_mat_and_barcodes(
      out_mat_f = '$filt_counts_f',
      out_bcs_f = '$bcs_f', 
      sel_s     = '{wildcards.run}', 
      af_mat_f  = '{input.h5_f}', 
      knee_f    = '{input.knee_data_f}', 
      ncores    = {threads}, 
      cell_calls_method = '{CELL_CALLS_METHOD}', 
      ambient_method    = '{AMBIENT_METHOD}')"

      # Create the output yaml file
      echo "filt_counts_f: $filt_counts_f" >> {output.ambient_yaml_out}
      echo "bcs_f: $bcs_f" >> {output.ambient_yaml_out}

      """


# rule get_barcode_qc_metrics:
#   input:
#     af_h5_f     = af_dir  + '/af_{run}/' + af_rna_dir + 'af_counts_mat.h5',
#     amb_yaml_f  = amb_dir + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml',
#     knee_yaml_f = af_dir  + '/af_{run}/' + af_rna_dir + 'ambient_params_{run}_' + DATE_STAMP + '.yaml'
#   params:
#     expected_cells = lambda wildcards: parse_ambient_params(AMBIENT_METHOD, CUSTOM_SAMPLE_PARAMS_F, wildcards.run,
#         af_dir + f'/af_{wildcards.run}/' + af_rna_dir + f'ambient_params_{wildcards.run}_{DATE_STAMP}.yaml', CELLBENDER_LEARNING_RATE)[0]
#   output:
#     bc_qc_f     = amb_dir + '/ambient_{run}/barcodes_qc_metrics_{run}_' + DATE_STAMP + '.txt.gz'
#   threads: 1
#   retries: RETRIES
#   conda:
#     '../envs/rlibs.yml'
#   resources:
#     mem_mb      = lambda wildcards, attempt: attempt * MB_GET_BARCODE_QC_METRICS
#   shell:
#     """
#     # save barcode stats
#     Rscript -e "source('scripts/ambient.R'); \
#       save_barcode_qc_metrics('{input.af_h5_f}', '{input.amb_yaml_f}', \
#         '{output.bc_qc_f}', {params.expected_cells}, '{AMBIENT_METHOD}')"
#     """


rule get_ambient_sample_statistics:
  input:
    metrics_fs  = expand(af_dir + '/af_{run}/' + af_rna_dir + 'knee_plot_data_{run}_' + DATE_STAMP + '.txt.gz', run=runs),
    amb_yaml_fs = expand(amb_dir + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml', run=runs)
  output:
    smpl_stats_f  = amb_dir + '/ambient_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  run:
    sample_stats_df   = extract_ambient_sample_statistics(AMBIENT_METHOD, SAMPLE_VAR, runs, input.metrics_fs, input.amb_yaml_fs,
      CUSTOM_SAMPLE_PARAMS_F, CELLBENDER_PROP_MAX_KEPT)
    sample_stats_df.to_csv(output.smpl_stats_f, index = False)
