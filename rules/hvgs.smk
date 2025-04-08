
import yaml  
import csv

localrules: make_hvg_df
localrules: merge_sample_stats

def make_hvgs_input_df(DEMUX_TYPE, SAMPLE_VAR, runs, ambient_outs_yamls, SAMPLE_MAPPING, FULL_TAG, DATE_STAMP, hvg_dir):

    df_list = []

    for r, yaml_file in zip(runs, ambient_outs_yamls):
        # get filtered ambient outputs
        with open(yaml_file) as f:
            amb_outs = yaml.load(f, Loader=yaml.FullLoader)

        amb_filt_f = amb_outs['filt_counts_f']

        if DEMUX_TYPE != "":
            # get sample ids for pool
            sample_ids = SAMPLE_MAPPING.get(r, [])

            for sample_id in sample_ids:
                hvg_df = pd.DataFrame({
                    SAMPLE_VAR: [r],
                    'amb_filt_f': [amb_filt_f],
                    'sample_id': [sample_id]
                })

                df_list.append(hvg_df)
        else:
            hvg_df = pd.DataFrame({
                SAMPLE_VAR: [r],
                'amb_filt_f': [amb_filt_f]
            })
            df_list.append(hvg_df)

    # merge dfs for all runs
    hvg_df_full = pd.concat(df_list, ignore_index=True)

    # add path to chunked file
    hvg_df_full['chunked_f'] = hvg_df_full['sample_id'].apply(lambda s: f"{hvg_dir}/chunked_counts_{s}_{FULL_TAG}_{DATE_STAMP}.h5")

    return hvg_df_full



# rule to create df with hvg input files and temporary chunked files

rule make_hvg_df: 
  input:
    ambient_yaml_out   = expand([amb_dir + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml'], run = runs)
  output:
    hvg_paths_f = hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  run:
    hvg_df = make_hvgs_input_df(DEMUX_TYPE, SAMPLE_VAR, runs, input.ambient_yaml_out, SAMPLE_MAPPING, FULL_TAG, DATE_STAMP, hvg_dir)
      # save dataframe
    hvg_df.to_csv(output.hvg_paths_f, index = False)



# create temporary csr h5 files
rule make_temp_counts:
  input:
    hvg_paths_f = hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    qc_f        = qc_dir  + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    rowdata_f   = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    qc_stats_f  = qc_dir + '/qc_sample_statistics_' + DATE_STAMP + '.txt'
  output:
    clean_h5_f  = expand(hvg_dir + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5', sample = SAMPLES)
  threads: 8
  retries: RETRIES
  resources:
    mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
  conda:
    '../envs/hvgs.yml'
  shell:
    """
    python3 scripts/hvgs.py get_csr_counts {input.hvg_paths_f} {input.qc_f} {input.qc_stats_f} {input.rowdata_f} {SAMPLE_VAR} {DEMUX_TYPE}

    """


if HVG_METHOD == 'sample': 
  # calculate stats for each sample separatelly  
  rule calc_hvg_stats_per_sample:
    input: 
      clean_h5_f = hvg_dir + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5'
    output:
      smpl_calcs_f = temp(hvg_dir + '/tmp_calcs_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
    threads: 1
    retries: RETRIES
    resources:
      mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
    conda:
      '../envs/hvgs.yml'
    shell:
      """
      python3 scripts/hvgs.py calculate_stats_per_sample {wildcards.sample} {input.clean_h5_f} {output.smpl_calcs_f}
      """
else:
  rule calc_hvg_stats_per_group:
    input:
      clean_h5_f  = expand(hvg_dir + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5', sample = SAMPLES),
      hvg_paths_f = hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
      rowdata_f   = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
      qc_smpl_stats_f = qc_dir + '/qc_sample_statistics_' + DATE_STAMP + '.txt'
    output: 
      chunk_calcs_f = temp(hvg_dir + '/tmp_calcs_{group}_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
    threads: 8
    resources:
      mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
    conda:
      '../envs/hvgs.yml'
    shell:
      """
      python3 scripts/hvgs.py calculate_stats_per_chunk \
        {input.hvg_paths_f} \
        {input.rowdata_f} \
        {METADATA_F} \
        {input.qc_smpl_stats_f} \
        {output.chunk_calcs_f} \
        {wildcards.chunk} \
        {HVG_METHOD} \
        {HVG_CHUNK_SIZE} \
        --group {wildcards.group} \
        --groupvar {HVG_SPLIT_VAR} \
        --ncores {threads} \
    
      """


if HVG_METHOD == 'sample':
  # combine sample stats
  rule merge_sample_stats:
    input:
      smpl_calcs_f      = expand(hvg_dir + '/tmp_calcs_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', sample = SAMPLES)
    output:
      smpl_calcs_merged = hvg_dir + '/sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
    threads: 1
    retries: RETRIES
    resources:
      mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
    run:
      # read all nonempty input files and concatenate them
      stats_df_ls = [
        pd.read_csv(f, compression='gzip', sep = '\t') 
        for f in input.smpl_calcs_f 
        if gzip.open(f, 'rb').read(1)
      ]
      stats_df_merged = pd.concat(stats_df_ls, ignore_index= True)

      stats_df_merged.to_csv(output.smpl_calcs_merged, sep='\t', index=False, compression='gzip')
else:
  rule merge_chunk_stats:
    input:                 
      chunk_calcs_f = expand(hvg_dir + '/tmp_calcs_{group}_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
      group=GROUP_NAMES, chunk=range(NUM_CHUNKS)),
    output:
      chunk_calcs_merged = hvg_dir + '/group_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
    threads: 1
    retries: RETRIES
    resources:
      mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
    run:
      # read all nonempty input files and concatenate them
      stats_df_ls = [
        pd.read_csv(f, compression='gzip', sep='\t') 
        for f in input.chunk_calcs_f 
        if gzip.open(f, 'rb').read(1)
      ]
      chunk_df_merged = pd.concat(stats_df_ls, ignore_index= True)

      chunk_df_merged.to_csv(output.chunk_calcs_merged, sep='\t', index=False, compression='gzip', quoting=csv.QUOTE_NONE)
    
