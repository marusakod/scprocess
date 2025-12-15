# snakemake rule for running cellbender

import os
import numpy as np
import yaml
import pandas as pd
import polars as pl

localrules: get_ambient_run_statistics, make_paths_h5_csv

def parse_ambient_params(run, config, RUN_PARAMS, af_dir, af_rna_dir):
  # get cellbender parameters from yaml file
  amb_yaml_f  = f'{af_dir}/af_{run}/{af_rna_dir}ambient_params_{run}_{config['project']['date_stamp']}.yaml'
  with open(amb_yaml_f) as f:
    amb_params  = yaml.load(f, Loader=yaml.FullLoader)

  # make dictionary for output
  params_dc = {
    "knee1":                      amb_params['knee1'], 
    "shin1":                      amb_params['shin1'], 
    "knee2":                      amb_params['knee2'], 
    "shin2":                      amb_params['shin2'], 
    "cb_expected_cells":          amb_params['cb_expected_cells'], 
    "cb_total_droplets_included": amb_params['cb_total_droplets_included'], 
    "cb_low_count_threshold":     amb_params['cb_low_count_threshold'], 
    "cb_learning_rate":           config['ambient']['cb_learning_rate'], 
    "cb_empty_training_fraction": config['ambient']['cb_empty_training_fraction'], 
    "cb_posterior_batch_size":    config['ambient']['cb_posterior_batch_size']
  }

  # check if cellbender specific parameters are defined for individual samples 
  if run in RUN_PARAMS:
    run_dc    = RUN_PARAMS[run]
    if "mapping" in run_dc:
      for v in ["knee1", "shin1", "knee2", "shin2"]:
        if v in run_dc["mapping"]:
          params_dc[v] = run_dc["mapping"][v]
    if "ambient" in run_dc:
      for v in ['cb_expected_cells', 'cb_total_droplets_included', 'cb_low_count_threshold', 
        'cb_learning_rate', 'cb_empty_training_fraction', 'cb_posterior_batch_size']:
        if v in run_dc["ambient"]:
          val = run_dc["ambient"][v]
          if val != "":
            params_dc[v] = val

  return params_dc


# metrics_fs_ls is the list of knee files (this should exist for all ambient methods)
# this function should run in the last rule
def extract_ambient_run_statistics(config, RUNS, RUN_PARAMS, metrics_fs_ls, ambient_outs_yamls, RUN_VAR):
  # unpack
  do_cellbender = config['ambient']['ambient_method'] == "cellbender"

  # loop through samples
  kept_arr    = []
  totals_arr  = []
  for sample, metrics_f, ambient_outs_yaml in zip(RUNS, metrics_fs_ls, ambient_outs_yamls):
    # get barcodes file from ambient outs yaml file
    with open(ambient_outs_yaml) as f:
      amb_outs = yaml.load(f, Loader=yaml.FullLoader)
   
    # count the number of barcodes
    bc_f          = amb_outs['bcs_f']
    barcode_count = pl.read_csv(bc_f, has_header = False).shape[0]
    kept_arr.append(barcode_count)

    # if cellbender, get the number of total droplets included
    if do_cellbender:
      total_droplets = pl.read_csv(metrics_f, schema_overrides={"rank": pl.Float32})['total_droplets_included'][0]
      totals_arr.append(total_droplets)

  # tidy up into array
  kept_arr = np.array(kept_arr)

  # if no cellbender it's easy
  if not do_cellbender:
    amb_stats_df = pl.DataFrame({
      RUN_VAR:          RUNS,
      'kept_droplets':  kept_arr
    })
  # otherwise add details of what cellbender used
  else:
    # convert to array
    totals_arr  = np.array(totals_arr)

    # get all samples with custom params
    custom_ls   = [ run_name for run_name in RUN_PARAMS if 'ambient' in RUN_PARAMS[run_name] ]
    amb_params  = [ RUN_PARAMS[run_name]['ambient'] for run_name in custom_ls ]

    # check whether to add total_droplet_included
    custom_df   = pl.DataFrame( schema = {RUN_VAR: pl.String, 'cb_total_droplets_included': pl.Int32} )
    for i, run_name in enumerate(custom_ls):
      # get this set of parameters
      this_params = amb_params[ i ]

      # if total droplets specified, use this instead
      if not this_params['cb_total_droplets_included'] == '':
        tmp_df      = pl.DataFrame({
          RUN_VAR: run_name, 
          'cb_total_droplets_included': this_params['cb_total_droplets_included']
        }).with_columns(
          # Cast the column 'cb_total_droplets_included' to Int32
          pl.col('cb_total_droplets_included').cast(pl.Int32)
        )
        custom_df   = custom_df.vstack( tmp_df )
    
    # use these to edit totals_arr
    samples_arr = np.array(RUNS)
    for row in custom_df.rows(named = True):
      match_idx = np.where(samples_arr == row[RUN_VAR])
      totals_arr[match_idx] = row['cb_total_droplets_included']

    # do some calculations
    prop_kept = kept_arr / totals_arr
    bad_idx   = prop_kept > config['ambient']['cb_max_prop_kept']

    # assemble into dataframe
    amb_stats_df = pl.DataFrame({
      RUN_VAR:            RUNS,
      'total_droplets':   totals_arr,
      'kept_droplets':    kept_arr,
      'prop_kept_by_cb':  prop_kept,
      'bad_run':       bad_idx
    })

  return amb_stats_df


if config['ambient']['ambient_method'] == 'cellbender':
  rule run_cellbender:
    input:
      af_h5_f     = f'{af_dir}/af_{{run}}/{af_rna_dir}af_counts_mat.h5',
      amb_yaml_f  = f'{af_dir}/af_{{run}}/{af_rna_dir}ambient_params_{{run}}_{DATE_STAMP}.yaml'
    params:
      cb_version                  = config['ambient']['cb_version'],
      cb_expected_cells           = lambda wildcards: parse_ambient_params(wildcards.run, config, RUN_PARAMS, af_dir, af_rna_dir)['cb_expected_cells'],
      cb_total_droplets_included  = lambda wildcards: parse_ambient_params(wildcards.run, config, RUN_PARAMS, af_dir, af_rna_dir)['cb_total_droplets_included'],
      cb_low_count_threshold      = lambda wildcards: parse_ambient_params(wildcards.run, config, RUN_PARAMS, af_dir, af_rna_dir)['cb_low_count_threshold'],
      cb_learning_rate            = lambda wildcards: parse_ambient_params(wildcards.run, config, RUN_PARAMS, af_dir, af_rna_dir)['cb_learning_rate'],
      cb_empty_training_fraction  = lambda wildcards: parse_ambient_params(wildcards.run, config, RUN_PARAMS, af_dir, af_rna_dir)['cb_empty_training_fraction'],
      cb_posterior_batch_size     = lambda wildcards: parse_ambient_params(wildcards.run, config, RUN_PARAMS, af_dir, af_rna_dir)['cb_posterior_batch_size']
    output:
      ambient_yaml_out = f'{amb_dir}/ambient_{{run}}/ambient_{{run}}_{DATE_STAMP}_output_paths.yaml',
      tmp_f            = temp(f'{amb_dir}/ambient_{{run}}/ckpt.tar.gz')
    threads: 4
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('run_cellbender', rules, 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run),
      runtime = lambda wildcards, attempt, input: get_resources('run_cellbender', rules, 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run)*(1.5**(attempt-1))
    benchmark:
      f'{benchmark_dir}/{SHORT_TAG}_ambient/run_cellbender_{{run}}_{DATE_STAMP}.benchmark.txt'
    container:
      config['ambient']['cellbender_image']
    shell: """
      # get parameters for cellbender
      CB_VERSION={params.cb_version}
      EXPECTED_CELLS={params.cb_expected_cells}
      TOTAL_DROPLETS_INCLUDED={params.cb_total_droplets_included}
      LOW_COUNT_THRESHOLD={params.cb_low_count_threshold}
      LEARNING_RATE={params.cb_learning_rate}
      EMPTY_TRAINING_FRACTION={params.cb_empty_training_fraction}
      POSTERIOR_BATCH_SIZE={params.cb_posterior_batch_size}
      
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

      # get --posterior-batch-size flag
      POSTERIOR_BATCH_SIZE_FLAG=""
      if [ "$CB_VERSION" == "v0.3.2" ]; then
        POSTERIOR_BATCH_SIZE_FLAG="--posterior-batch-size $POSTERIOR_BATCH_SIZE"
      fi

      # switch between cellbender running options
      if [ $EXPECTED_CELLS -eq 0 ]; then
        echo "running cellbender without expected_cells"
        # run cellbender
        cellbender remove-background \
          --input {input.af_h5_f} \
          --output $raw_counts_f \
          --total-droplets-included $TOTAL_DROPLETS_INCLUDED \
          --low-count-threshold $LOW_COUNT_THRESHOLD \
          --empty-drop-training-fraction $EMPTY_TRAINING_FRACTION \
          --learning-rate $LEARNING_RATE \
          --checkpoint-mins 30.0 \
          --cuda \
          $POSTERIOR_BATCH_SIZE_FLAG
      else
        echo "running cellbender with expected_cells"
        # run cellbender
        cellbender remove-background \
          --input {input.af_h5_f} \
          --output $raw_counts_f \
          --expected-cells $EXPECTED_CELLS \
          --total-droplets-included $TOTAL_DROPLETS_INCLUDED \
          --low-count-threshold $LOW_COUNT_THRESHOLD \
          --empty-drop-training-fraction $EMPTY_TRAINING_FRACTION \
          --learning-rate $LEARNING_RATE \
          --checkpoint-mins 30.0 \
          --cuda \
          $POSTERIOR_BATCH_SIZE_FLAG
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


if config['ambient']['ambient_method'] == 'decontx':
  rule run_decontx:
    input:
      af_h5_f         = f'{af_dir}/af_{{run}}/{af_rna_dir}af_counts_mat.h5',
      amb_yaml_f      = f'{af_dir}/af_{{run}}/{af_rna_dir}ambient_params_{{run}}_{DATE_STAMP}.yaml', 
      knee_data_f     = f'{af_dir}/af_{{run}}/{af_rna_dir}knee_plot_data_{{run}}_{DATE_STAMP}.csv.gz'
    output:
      amb_yaml_out    = f'{amb_dir}/ambient_{{run}}/ambient_{{run}}_{DATE_STAMP}_output_paths.yaml'
    params:
      ambient_method  = config['ambient']['ambient_method'],
      cell_calling    = config['ambient']['cell_calling']
    threads: 4
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('run_decontx', rules, 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run),
      runtime = lambda wildcards, input: get_resources('run_decontx', rules, 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run)
    benchmark:
      f'{benchmark_dir}/{SHORT_TAG}_ambient/run_decontx_{{run}}_{DATE_STAMP}.benchmark.txt'
    conda: 
      '../envs/ambientr.yaml'
    shell:"""
      # create main ambient directory
      mkdir -p {amb_dir}

      # change to sample ambient directory
      amb_dir=$(dirname {output.amb_yaml_out})
      mkdir -p $amb_dir
  
      # define output file names
      filt_counts_f="{amb_dir}/ambient_{wildcards.run}/decontx_{wildcards.run}_{DATE_STAMP}_filtered.h5"
      bcs_f="{amb_dir}/ambient_{wildcards.run}/decontx_{wildcards.run}_{DATE_STAMP}_cell_barcodes.csv"
      dcx_params_f="{amb_dir}/ambient_{wildcards.run}/decontx_{wildcards.run}_{DATE_STAMP}_params.csv.gz"

      # run cell calling and decontamination
      Rscript -e "source('scripts/ambient.R'); source('scripts/utils.R'); \
        get_cell_mat_and_barcodes(
          out_mat_f         = '$filt_counts_f',
          out_bcs_f         = '$bcs_f',
          out_dcx_f         = '$dcx_params_f',
          sel_s             = '{wildcards.run}',
          af_mat_f          = '{input.af_h5_f}',
          knee_f            = '{input.knee_data_f}',
          ncores            =  {threads},
          cell_calls_method = '{params.cell_calling}',
          ambient_method    = '{params.ambient_method}'
        )"

      # Create the output yaml file
      echo "filt_counts_f: $filt_counts_f" >> {output.amb_yaml_out}
      echo "bcs_f: $bcs_f" >> {output.amb_yaml_out}
      echo "dcx_params_f: $dcx_params_f" >> {output.amb_yaml_out}

      """


if config['ambient']['ambient_method'] == 'none':
  rule run_cell_calling:
    input:
      h5_f        = f'{af_dir}/af_{{run}}/{af_rna_dir}af_counts_mat.h5',
      amb_yaml_f  = f'{af_dir}/af_{{run}}/{af_rna_dir}ambient_params_{{run}}_{DATE_STAMP}.yaml',
      knee_data_f = f'{af_dir}/af_{{run}}/{af_rna_dir}knee_plot_data_{{run}}_{DATE_STAMP}.csv.gz'
    output:
      ambient_yaml_out = f'{amb_dir}/ambient_{{run}}/ambient_{{run}}_{DATE_STAMP}_output_paths.yaml'
    params:
      ambient_method  = config['ambient']['ambient_method'],
      cell_calling    = config['ambient']['cell_calling']
    threads: 4
    retries: config['resources']['retries']
    conda:
      '../envs/ambientr.yaml'
    resources:
      mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('run_cell_calling', rules, 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run),
      runtime = lambda wildcards, input: get_resources('run_cell_calling', rules, 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run)
    benchmark:
      f'{benchmark_dir}/{SHORT_TAG}_ambient/run_cell_calling_{{run}}_{DATE_STAMP}.benchmark.txt'
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
      Rscript -e "source('scripts/ambient.R'); source('scripts/utils.R');
        get_cell_mat_and_barcodes(
          out_mat_f         = '$filt_counts_f',
          out_bcs_f         = '$bcs_f', 
          sel_s             = '{wildcards.run}', 
          af_mat_f          = '{input.af_h5_f}', 
          knee_f            = '{input.knee_data_f}', 
          ncores            =  {threads}, 
          cell_calls_method = '{params.cell_calling}',
          ambient_method    = '{params.ambient_method}'
        )"

      # Create the output yaml file
      echo "filt_counts_f: $filt_counts_f" >> {output.ambient_yaml_out}
      echo "bcs_f: $bcs_f" >> {output.ambient_yaml_out}
      """


rule get_barcode_qc_metrics:
  input:
    af_h5_f     = f'{af_dir}/af_{{run}}/{af_rna_dir}af_counts_mat.h5',
    amb_yaml_f  = f'{amb_dir}/ambient_{{run}}/ambient_{{run}}_{DATE_STAMP}_output_paths.yaml',
    knee_yaml_f = f'{af_dir}/af_{{run}}/{af_rna_dir}ambient_params_{{run}}_{DATE_STAMP}.yaml'
  params:
    ambient_method  = config['ambient']['ambient_method']
  output:
    bc_qc_f     = f'{amb_dir}/ambient_{{run}}/barcodes_qc_metrics_{{run}}_{DATE_STAMP}.csv.gz'
  threads: 1
  retries: config['resources']['retries']
  conda:
    '../envs/ambientr.yaml'
  resources:
    mem_mb  = lambda wildcards, attempt, input: attempt * get_resources('get_barcode_qc_metrics', rules, 'memory', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run),
    runtime = lambda wildcards, input: get_resources('get_barcode_qc_metrics', rules, 'time', lm_f, config, schema_f, input, BATCHES, RUN_PARAMS, wildcards.run)
  benchmark:
    f'{benchmark_dir}/{SHORT_TAG}_ambient/get_barcode_qc_metrics_{{run}}_{DATE_STAMP}.benchmark.txt'
  shell: """
    # save barcode stats
    Rscript -e "source('scripts/ambient.R'); source('scripts/utils.R'); \
      save_barcode_qc_metrics('{input.af_h5_f}', '{input.amb_yaml_f}', \
        '{output.bc_qc_f}', '{params.ambient_method}')"
    """


rule get_ambient_run_statistics:
  input:
    metrics_fs  = expand(f'{af_dir}/af_{{run}}/{af_rna_dir}knee_plot_data_{{run}}_{DATE_STAMP}.csv.gz', run=RUNS),
    bc_qc_fs    = expand(f'{amb_dir}/ambient_{{run}}/barcodes_qc_metrics_{{run}}_{DATE_STAMP}.csv.gz', run=RUNS),
    amb_yaml_fs = expand(f'{amb_dir}/ambient_{{run}}/ambient_{{run}}_{DATE_STAMP}_output_paths.yaml', run=RUNS)
  output:
    amb_stats_f = f'{amb_dir}/ambient_run_statistics_{FULL_TAG}_{DATE_STAMP}.csv'
  run:
    amb_stats_df = extract_ambient_run_statistics(config, RUNS, 
      RUN_PARAMS, input.metrics_fs, input.amb_yaml_fs, RUN_VAR)
    amb_stats_df.write_csv(output.amb_stats_f)
