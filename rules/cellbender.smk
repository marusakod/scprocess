# snakemake rule for running cellbender
import numpy as np

localrules: exclude_bad_cellbender_samples

# setting up conda environment
# ml .testing Anaconda3/2023.03-libmamba
# conda create -n CellBender python=3.7 --solver=libmamba

def parse_cellbender_params(CUSTOM_CELLBENDER_PARAMS_F, sample, cb_yaml_f, 
  FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD):
  """
  Parse cellbender parameters from yaml file.
  """
  import yaml
  with open(cb_yaml_f) as f:
    cb_params   = yaml.load(f, Loader=yaml.FullLoader)
  # get parameters for cellbender
  EXPECTED_CELLS          = cb_params['expected_cells']
  TOTAL_DROPLETS_INCLUDED = cb_params['total_droplets_included']
  LOW_COUNT_THRESHOLD     = cb_params['low_count_threshold']
  LEARNING_RATE           = CELLBENDER_LEARNING_RATE

  # force parameters if needed
  if FORCE_EXPECTED_CELLS is not None:
    EXPECTED_CELLS          = FORCE_EXPECTED_CELLS
  if FORCE_TOTAL_DROPLETS_INCLUDED is not None:
    TOTAL_DROPLETS_INCLUDED = FORCE_TOTAL_DROPLETS_INCLUDED
  if FORCE_LOW_COUNT_THRESHOLD is not None:
    LOW_COUNT_THRESHOLD     = FORCE_LOW_COUNT_THRESHOLD

  # load up custom parameters if exists
  if CUSTOM_CELLBENDER_PARAMS_F is not None:
    params_df   = pd.read_csv(CUSTOM_CELLBENDER_PARAMS_F)
    # if sample is in custom parameters, then use parameters from there
    if params_df['sample_id'].str.contains(sample).any():
      print(f'using custom cellbender parameters for {sample}')
      row_idx                 = params_df['sample_id'] == sample
      EXPECTED_CELLS          = int(params_df[ row_idx ]['expected_cells'].values[0])
      TOTAL_DROPLETS_INCLUDED = int(params_df[ row_idx ]['total_droplets_included'].values[0])
      LOW_COUNT_THRESHOLD     = int(params_df[ row_idx ]['low_count_threshold'].values[0])
      LEARNING_RATE           = params_df[ row_idx ]['learning_rate'].values[0]

  return EXPECTED_CELLS, TOTAL_DROPLETS_INCLUDED, LOW_COUNT_THRESHOLD, LEARNING_RATE


def extract_bad_cellbender_proportions(do_cellbender, samples_ls, metrics_fs, bcs_fs, custom_f, max_kept = 0.9):
  if not do_cellbender:
    # assemble into dataframe
    cb_bad_df     = pd.DataFrame({
      'sample_id':        samples_ls,
      'total_droplets':   np.array([pd.read_csv(f)['total_droplets_included'][0] for f in metrics_fs]),
      'kept_droplets':    np.array([np.nan for f in metrics_fs]),
      'prop_kept_by_cb':  np.array([np.nan for f in metrics_fs]),
      'bad_sample':       np.array([False for f in metrics_fs])
      })

    return cb_bad_df

  # from each metrics file, get the number of cells called as barcodes
  totals_arr    = np.array([pd.read_csv(f)['total_droplets_included'][0] for f in metrics_fs])

  # replace dodgy totals values with custom if need be
  if (custom_f is not None) and (os.path.isfile(custom_f)):
    # load up custom parameters
    custom_df     = pd.read_csv(custom_f)[[ 'sample_id', 'total_droplets_included' ]]

    # iterate through rows
    samples_arr   = np.array(samples_ls)
    for idx, row in custom_df.iterrows():
      match_idx     = np.where(samples_arr == row['sample_id'])
      totals_arr[ match_idx ] = row['total_droplets_included']

  # from each barcodes file, get the number of cells called as barcodes
  kept_arr      = np.array([pd.read_csv(f, header = None).shape[0] for f in bcs_fs])

  # do some calculations
  prop_kept     = kept_arr / totals_arr
  bad_idx       = prop_kept > max_kept

  # assemble into dataframe
  cb_bad_df     = pd.DataFrame({
    'sample_id':        samples_ls,
    'total_droplets':   totals_arr,
    'kept_droplets':    kept_arr,
    'prop_kept_by_cb':  prop_kept,
    'bad_sample':       bad_idx
    })

  return cb_bad_df


if DO_CELLBENDER:
  rule run_cellbender:
    input:
      h5_f      = af_dir + '/af_{sample}/af_counts_mat.h5',
      cb_yaml_f = af_dir + '/af_{sample}/bender_params_{sample}_' + DATE_STAMP + '.yaml'
    params:
      expected_cells          = lambda wildcards: parse_cellbender_params(CUSTOM_CELLBENDER_PARAMS_F, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/bender_params_{wildcards.sample}_{DATE_STAMP}.yaml', 
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[0],
      total_droplets_included = lambda wildcards: parse_cellbender_params(CUSTOM_CELLBENDER_PARAMS_F, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/bender_params_{wildcards.sample}_{DATE_STAMP}.yaml', 
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[1],
      low_count_threshold     = lambda wildcards: parse_cellbender_params(CUSTOM_CELLBENDER_PARAMS_F, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/bender_params_{wildcards.sample}_{DATE_STAMP}.yaml', 
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[2],
      learning_rate           = lambda wildcards: parse_cellbender_params(CUSTOM_CELLBENDER_PARAMS_F, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/bender_params_{wildcards.sample}_{DATE_STAMP}.yaml', 
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[3]
    output:
      cb_full_f = cb_dir + '/bender_{sample}/bender_{sample}_' + DATE_STAMP + '.h5',
      cb_filt_f = cb_dir + '/bender_{sample}/bender_{sample}_' + DATE_STAMP + '_filtered.h5',
      cb_bcs_f  = cb_dir + '/bender_{sample}/bender_{sample}_' + DATE_STAMP + '_cell_barcodes.csv',
      tmp_f     = temp(cb_dir + '/bender_{sample}/ckpt.tar.gz')
    threads: 1
    resources:
      mem_mb      = lambda wildcards, attempt: attempt * MB_RUN_CELLBENDER,
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
      cb_dir=$(dirname {output.cb_full_f})
      mkdir -p $cb_dir
      cd $cb_dir

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

      # check whether temp file was actually made; if not, make an empty one
      if [ ! -f {output.tmp_f} ]; then
        touch {output.tmp_f}
      fi
      """
else:
  localrules: run_cellbender
  rule run_cellbender:
    input:
      h5_f      = af_dir + '/af_{sample}/af_counts_mat.h5',
      cb_yaml_f = af_dir + '/af_{sample}/bender_params_{sample}_' + DATE_STAMP + '.yaml'
    params:
      expected_cells          = lambda wildcards: parse_cellbender_params(CUSTOM_CELLBENDER_PARAMS_F, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/bender_params_{wildcards.sample}_{DATE_STAMP}.yaml', 
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[0],
      total_droplets_included = lambda wildcards: parse_cellbender_params(CUSTOM_CELLBENDER_PARAMS_F, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/bender_params_{wildcards.sample}_{DATE_STAMP}.yaml', 
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[1],
      low_count_threshold     = lambda wildcards: parse_cellbender_params(CUSTOM_CELLBENDER_PARAMS_F, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/bender_params_{wildcards.sample}_{DATE_STAMP}.yaml', 
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[2],
      learning_rate           = lambda wildcards: parse_cellbender_params(CUSTOM_CELLBENDER_PARAMS_F, wildcards.sample,
        af_dir + f'/af_{wildcards.sample}/bender_params_{wildcards.sample}_{DATE_STAMP}.yaml', 
        FORCE_TOTAL_DROPLETS_INCLUDED, FORCE_EXPECTED_CELLS, FORCE_LOW_COUNT_THRESHOLD)[3]
    output:
      cb_full_f = cb_dir + '/bender_{sample}/bender_{sample}_' + DATE_STAMP + '.h5',
      cb_filt_f = cb_dir + '/bender_{sample}/bender_{sample}_' + DATE_STAMP + '_filtered.h5',
      cb_bcs_f  = cb_dir + '/bender_{sample}/bender_{sample}_' + DATE_STAMP + '_cell_barcodes.csv'
    threads: 1
    resources:
      mem_mb    = 100
    shell:
      """
      touch {output.cb_full_f}
      touch {output.cb_filt_f}
      touch {output.cb_bcs_f}
      """


rule get_cellbender_qc_metrics:
  input:
    knee_data_f = af_dir + '/af_{sample}/knee_plot_data_{sample}_' + DATE_STAMP + '.txt.gz',
    af_h5_f     = af_dir + '/af_{sample}/af_counts_mat.h5',
    cb_full_f   = cb_dir + '/bender_{sample}/bender_{sample}_' + DATE_STAMP + '.h5'
  output:
    cb_qc_f     = cb_dir + '/bender_{sample}/bender_qc_metrics_{sample}_' + DATE_STAMP + '.txt.gz'
  threads: 1
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * MB_GET_CELLBENDER_QC_METRICS
  conda: 
   '../envs/rlibs.yml'
  shell:
    """
  
    Rscript -e "source('scripts/cellbender.R'); \
      save_cellbender_qc_metrics('{input.knee_data_f}', '{input.af_h5_f}', '{input.cb_full_f}', \
        '{output.cb_qc_f}', '{DO_CELLBENDER}')"
    """


rule exclude_bad_cellbender_samples:
  input:
    metrics_fs  = expand( af_dir + '/af_{sample}/knee_plot_data_{sample}_' + DATE_STAMP + '.txt.gz', sample = SAMPLES ),
    bcs_fs      = expand( cb_dir + '/bender_{sample}/bender_{sample}_' + DATE_STAMP + '_cell_barcodes.csv', sample = SAMPLES )
  output:
    cb_bad_f    = cb_dir + '/bender_bad_samples_' + DATE_STAMP + '.txt'
  run:
    cb_bad_df   = extract_bad_cellbender_proportions(DO_CELLBENDER, SAMPLES, input.metrics_fs, input.bcs_fs, 
      CUSTOM_CELLBENDER_PARAMS_F, CELLBENDER_PROP_MAX_KEPT)
    cb_bad_df.to_csv(output.cb_bad_f, sep = '\t', index = False)
