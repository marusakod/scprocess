# snakemake rule for running cellbender

import numpy as np
import yaml
import pandas as pd

localrules: exclude_bad_cellbender_samples

def parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, sample, cb_yaml_f,
                         FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE):

    # get cellbender parameters from yaml file
    with open(cb_yaml_f) as f:
        cb_params = yaml.load(f, Loader=yaml.FullLoader)
        
    # get parameters for cellbender
    EXPECTED_CELLS = cb_params['expected_cells']
    TOTAL_DROPLETS_INCLUDED = cb_params['total_droplets_included']
    LOW_COUNT_THRESHOLD = cb_params['low_count_threshold']
    KNEE_1 = cb_params['knee_1']
    INFLECTION_1 = cb_params['inflection_1']
    KNEE_2 = cb_params['knee_2']
    INFLECTION_2 = cb_params['inflection_2']
    LEARNING_RATE = CELLBENDER_LEARNING_RATE
    
    # additional custom parameters for decontx and none
    EMPTY_START = None
    EMPTY_END = None

    # force parameters if needed (this is ignored if method is not cellbender)
    if FORCE_EXPECTED_CELLS is not None:
        EXPECTED_CELLS = FORCE_EXPECTED_CELLS
    if FORCE_TOTAL_DROPLETS_INCLUDED is not None:
        TOTAL_DROPLETS_INCLUDED = FORCE_TOTAL_DROPLETS_INCLUDED
    if FORCE_LOW_COUNT_THRESHOLD is not None:
        LOW_COUNT_THRESHOLD = FORCE_LOW_COUNT_THRESHOLD

    # load up custom parameters if exists
    if CUSTOM_PARAMS_F is not None:
        params_df = pd.read_csv(CUSTOM_PARAMS_F)
        
        # if sample is in custom parameters, then use parameters from there
        if params_df['sample_id'].str.contains(sample).any():
            # subset custom params df 
            filt_params_df = params_df[params_df['sample_id'] == sample]
            filt_params_df = filt_params_df.reset_index(drop=True)
            fasta_f = filt_params_df.loc[0, 'fasta_f']

            print(f'using custom parameters for {sample}')
            
            if AMBIENT_METHOD == 'cellbender':
                EXPECTED_CELLS = int(filt_params_df.loc[0, 'expected_cells'])
                TOTAL_DROPLETS_INCLUDED = int(filt_params_df.loc[0, 'total_droplets_included'])
                LOW_COUNT_THRESHOLD = int(filt_params_df.loc[0, 'low_count_threshold'])
                LEARNING_RATE = filt_params_df.loc[0, 'learning_rate']
                
            else:
                if CELL_CALLS_METHOD == 'emptyDrops':
                    KNEE_1 = int(filt_params_df.loc[0, 'retain'])
                    EMPTY_START = int(filt_params_df.loc[0, 'empty_start'])
                    EMPTY_END = int(filt_params_df.loc[0, 'empty_end'])
                else:
                    EXPECTED_CELLS = int(filt_params_df.loc[0, 'expected_cells'])
                    EMPTY_START = int(filt_params_df.loc[0, 'empty_start'])
                    EMPTY_END = int(filt_params_df.loc[0, 'empty_end'])

    return EXPECTED_CELLS, TOTAL_DROPLETS_INCLUDED, LOW_COUNT_THRESHOLD, LEARNING_RATE, \
           KNEE_1, INFLECTION_1, KNEE_2, INFLECTION_2, EMPTY_START, EMPTY_END




# metrics_fs_ls is the list of knee files (this should exist for all ambient methods)
# this function should run in the last rule
def extract_sample_statistics(AMBIENT_METHOD, samples_ls, metrics_fs_ls, ambient_outs_yamls, custom_f, max_kept=0.9):
    kept_arr = []
    totals_arr = []

    for sample, metrics_f, ambient_outs_yaml in zip(samples_ls, metrics_fs_ls, ambient_outs_yamls):
        # Load ambient outs YAML file
        with open(ambient_outs_yaml) as f:
            amb_outs = yaml.load(f, Loader=yaml.FullLoader)

        # Determine the correct barcode file path based on AMBIENT_METHOD
        if AMBIENT_METHOD == 'cellbender':
            bc_f = amb_outs['cb_bcs_f']
        elif AMBIENT_METHOD == 'decontx':
            bc_f = amb_outs['dcx_bcs_f']
        else:
            bc_f = amb_outs['cell_bcs_f']

        # Read the CSV file and count the number of barcodes
        barcode_count = pd.read_csv(bc_f, header=None).shape[0]
        kept_arr.append(barcode_count)

        if AMBIENT_METHOD == 'cellbender':
            # Read the metrics file and get the number of cells called as barcodes
            total_droplets = pd.read_csv(metrics_f)['total_droplets_included'][0]
            totals_arr.append(total_droplets)

    kept_arr = np.array(kept_arr)

    if AMBIENT_METHOD != 'cellbender':
        sample_df = pd.DataFrame({
            'sample_id': samples_ls,
            'kept_droplets': kept_arr
        })
    else:
        totals_arr = np.array(totals_arr)

        # Replace dodgy totals values with custom if need be
        if custom_f is not None and os.path.isfile(custom_f):
            # Load up custom parameters
            custom_df = pd.read_csv(custom_f)[['sample_id', 'total_droplets_included']]

            # Iterate through rows
            samples_arr = np.array(samples_ls)
            for idx, row in custom_df.iterrows():
                match_idx = np.where(samples_arr == row['sample_id'])
                totals_arr[match_idx] = row['total_droplets_included']

        # Do some calculations
        prop_kept = kept_arr / totals_arr
        bad_idx = prop_kept > max_kept

        # Assemble into dataframe
        sample_df = pd.DataFrame({
            'sample_id': samples_ls,
            'total_droplets': totals_arr,
            'kept_droplets': kept_arr,
            'prop_kept_by_cb': prop_kept,
            'bad_sample': bad_idx
        })

    return sample_df




if AMBIENT_METHOD == 'cellbender':
  rule run_ambient:
    input:
      h5_f      = af_dir + '/af_{sample}/af_counts_mat.h5',
      cb_yaml_f = af_dir + '/af_{sample}/ambient_params_{sample}_' + DATE_STAMP + '.yaml'
    params:
      expected_cells          = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[0],
      total_droplets_included = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[1],
      low_count_threshold     = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[2],
      learning_rate           = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[3]
    output:
        ambient_yaml_out = amb_dir + '/ambient_{sample}/ambient_{sample}_' + DATE_STAMP + '_output_paths.yaml'
    threads: 1
    resources:
      mem_mb      = 8192,
      nvidia_gpu  = 1
    singularity:
      CELLBENDER_IMAGE
    shell:
      """
      # get parameters for cellbender
      EXPECTED_CELLS={params.expected_cells}
      TOTAL_DROPLETS_INCLUDED={params.total_droplets_included}
      LOW_COUNT_THRESHOLD={params.low_count_threshold}
      LEARNING_RATE={params.learning_rate}

      # change to cellbender directory
      amb_dir=$(dirname {output.cb_full_f})
      mkdir -p $amb_dir
      cd $amb_dir

      # define output files
      cb_full_f = amb_dir + '/ambient_{sample}/bender_{sample}_' + DATE_STAMP + '.h5',
      cb_filt_f = amb_dir + '/ambient_{sample}/bender_{sample}_' + DATE_STAMP + '_filtered.h5',
      cb_bcs_f  = amb_dir + '/ambient_{sample}/bender_{sample}_' + DATE_STAMP + '_cell_barcodes.csv'
      tmp_f     = temp(amb_dir + '/ambient_{sample}/ckpt.tar.gz')

      
      if [ $EXPECTED_CELLS -eq 0 ]; then
        echo "running cellbender without expected_cells"
        # run cellbender
        cellbender remove-background \
          --input {input.h5_f} \
          --output {output.cb_full_f} \
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
          --output {output.cb_full_f} \
          --expected-cells $EXPECTED_CELLS \
          --total-droplets-included $TOTAL_DROPLETS_INCLUDED \
          --low-count-threshold $LOW_COUNT_THRESHOLD \
          --learning-rate $LEARNING_RATE \
          --checkpoint-mins 30.0 \
          --cuda
      fi

      # Create the output yaml file
      echo "  cb_full_f: $cb_full_f" >> {output.ambient_yaml_out}
      echo "  cb_filt_f: $cb_filt_f" >> {output.ambient_yaml_out}
      echo "  cb_bcs_f: $cb_bcs_f" >> {output.ambient_yaml_out}
      echo "  tmp_f: $tmp_f" >> {output.ambient_yaml_out}

      # check whether temp file was actually made; if not, make an empty one
      if [ ! -f $tmp_f ]; then
        touch $tmp_f
      fi
      """
elif AMBIENT_METHOD == 'decontx':
  localrules: run_ambient
  rule run_ambient:
    input:
      h5_f      = af_dir + '/af_{sample}/af_counts_mat.h5',
      cb_yaml_f = af_dir + '/af_{sample}/ambient_params_{sample}_' + DATE_STAMP + '.yaml'
    params:
      expected_cells          = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[0],
      total_droplets_included = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[1],
      knee_1                  = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[4],
      inflection_1            = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[5],
      knee_2                  = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[6],
      empty_start             = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[8],
      empty_end               = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[9]
    output:
      ambient_yaml_out = amb_dir + '/ambient_{sample}/ambient_{sample}_' + DATE_STAMP + '_output_paths.yaml'
    threads: 4
    resources:
      mem_mb    = 8192
    conda: 
      '../envs/rlibs.yml'
    shell:
      """
      # define output file names
      dcx_filt_f = amb_dir + '/ambient_{sample}/decontx_{sample}_' + DATE_STAMP + '_filtered.h5',
      dcx_bcs_f  = amb_dir + '/ambient_{sample}/decontx_{sample}_' + DATE_STAMP + '_cell_barcodes.csv'
      dcx_params_f = amb_dir + '/ambient_{sample}/decontx_{sample}_' + DATE_STAMP + '_params.txt.gz'

      # run cell calling and decontamination
   
      Rscript -e "source('scripts/cellbender.R'); \
      get_cell_mat_and_barcodes(
      out_mat_f = $dcx_filt_f, \
      out_bcs_f = $dcx_bcs_f, \
      out_dcx_f = $dcx_params_f, \
      sel_s = '{sample}', \
      af_mat_f = '{input.h5_f}', \
      knee_1 = '{params.knee_1}', \
      knee_2 = '{params.knee_2}', \
      inf_1 = '{params.inflection_1}', \
      n_cores = {threads}, \
      total_included = '{params.total_droplets_included}', \
      exp_cells = '{params.expected_cells}', \
      empty_start = '{params.empty_start}', \
      empty_end = '{params.empty_end}', \
      cell_calls_method = '{CELL_CALLS_METHOD}', \
      ambient_method = '{AMBIENT_METHOD}')

      # Create the output yaml file
      echo "  dcx_filt_f: $dcx_filt_f" >> {output.ambient_yaml_out}
      echo "  dcx_bcs_f: $dcx_bcs_f" >> {output.ambient_yaml_out}
      echo "  dcx_params_f: $dcx_params_f" >> {output.ambient_yaml_out}

      """
else:
  localrules: run_ambient
  rule run_ambient:
    input:
      h5_f      = af_dir + '/af_{sample}/af_counts_mat.h5',
      cb_yaml_f = af_dir + '/af_{sample}/ambient_params_{sample}_' + DATE_STAMP + '.yaml'
    params:
      expected_cells          = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[0],
      total_droplets_included = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[1],
      knee_1                  = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[4],
      inflection_1            = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[5],
      knee_2                  = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[6],
      empty_start             = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[8],
      empty_end               = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[9]
    output:
      ambient_yaml_out = amb_dir + '/ambient_{sample}/ambient_{sample}_' + DATE_STAMP + '_output_paths.yaml'
    threads: 4
    conda:
      '../envs/rlibs.yml'
    resources:
      mem_mb    = 8192
    shell:
      """
      # define output file names
      cell_filt_f = amb_dir + '/ambient_{sample}/uncorrected_{sample}_' + DATE_STAMP + '_filtered.h5',
      cell_bcs_f  = amb_dir + '/ambient_{sample}/uncorrected_{sample}_' + DATE_STAMP + '_cell_barcodes.csv'

      # run cell calling and decontamination
      Rscript -e "source('scripts/cellbender.R'); \
      get_cell_mat_and_barcodes(
      out_mat_f = $cell_filt_f, \
      out_bcs_f = $cell_bcs_f, \
      sel_s = '{sample}', \
      af_mat_f = '{input.h5_f}', \
      knee_1 = '{params.knee_1}', \
      knee_2 = '{params.knee_2}', \
      inf_1 = '{params.inflection_1}', \
      total_included = '{params.total_droplets_included}', \
      exp_cells = '{params.expected_cells}', \
      empty_start = '{params.empty_start}', \
      empty_end = '{params.empty_end}', \
      n_cores = {threads}, \
      cell_calls_method = '{CELL_CALLS_METHOD}', \
      ambient_method = '{AMBIENT_METHOD}')

      # Create the output yaml file
      echo "  cell_filt_f: $cell_filt_f" >> {output.ambient_yaml_out}
      echo "  cell_bcs_f: $cell_bcs_f" >> {output.ambient_yaml_out}

      """



rule get_barcode_qc_metrics:
  input:
    af_h5_f     = af_dir + '/af_{sample}/af_counts_mat.h5',
    amb_yaml_f  = amb_dir + '/ambient_{sample}/ambient_{sample}_' + DATE_STAMP + '_output_paths.yaml'
  params:
    expected_cells = lambda wildcards: parse_ambient_params(CUSTOM_PARAMS_F, AMBIENT_METHOD, CELL_CALLS_METHOD, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/ambient_params_{wildcards.sample}_{DATE_STAMP}.yaml',
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD, CELLBENDER_LEARNING_RATE)[0]
  output:
    bc_qc_f     = amb_dir + '/ambient_{sample}/barcodes_qc_metrics_{sample}_' + DATE_STAMP + '.txt.gz'
  threads: 1
  conda:
    '../envs/rlibs.yml'
  resources:
    mem_mb      = 8192
  shell:
    """
    # save barcode stats
    export R_LIBS_USER='{RLIBS_DIR}'
    Rscript -e "source('scripts/cellbender.R'); \
      save_barcode_qc_metrics('{input.af_h5_f}', '{input.amb_yaml_f}', \
        '{output.bc_qc_f}', '{params.expected_cells}', '{AMBIENT_METHOD}')"
    """


rule get_ambient_sample_statistics:
  input:
    metrics_fs  = expand(af_dir + '/af_{sample}/knee_plot_data_{sample}_' + DATE_STAMP + '.txt.gz', sample=SAMPLES),
    amb_yaml_fs = expand(amb_dir + '/ambient_{sample}/ambient_{sample}_' + DATE_STAMP + '_output_paths.yaml', sample=SAMPLES)
  output:
    smpl_stats_f    = amb_dir + '/ambient_sample_statistics_' + DATE_STAMP + '.txt'
  run:
    sample_stats_df   = extract_sample_statistics(AMBIENT_METHOD, SAMPLES, input.metrics_fs, input.amb_yaml_fs,
      CUSTOM_PARAMS_F, CELLBENDER_PROP_MAX_KEPT)
    sample_stats_df.to_csv(output.smpl_stats_f, sep = '\t', index = False)