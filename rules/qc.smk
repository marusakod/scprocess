# snakemake rule for doing QC on sce object
import pandas as pd
import gzip
import os
import numpy as np

def extract_qc_sample_statistics(ambient_stats_f, qc_merged_f, SAMPLES, SAMPLE_VAR, AMBIENT_METHOD, DEMUX_TYPE, SAMPLE_MAPPING, QC_MIN_CELLS):
    # load the merged qc file
    qc_dt = pd.read_csv(qc_merged_f, sep = '\t', compression='gzip')

    # filter for cells that passed qc
    qc_dt = qc_dt[qc_dt["keep"] == True]

    # count the number of cells per sample
    sample_df = (
        qc_dt.groupby('sample_id')
        .size()
        .reset_index(name='n_cells')
    )

    # identify samples that do not meet the minimum cell threshold
    sample_df['bad_qc'] = sample_df['n_cells'] < QC_MIN_CELLS

    # handle samples excluded after cellbender
    if AMBIENT_METHOD == 'cellbender':
        # load ambient sample stats
        amb_stats = pd.read_csv(ambient_stats_f)
        # get bad pools or samples
        bad_bender = amb_stats.loc[amb_stats['bad_sample'], SAMPLE_VAR].tolist()

        if DEMUX_TYPE != "none":
            bad_bender_samples = []
            for p in bad_bender:
                if p in SAMPLE_MAPPING:
                    bad_bender_samples.extend(SAMPLE_MAPPING[p])
          
            assert all(s in SAMPLES for s in bad_bender_samples), \
                "Some bad bender samples are not in the SAMPLES list."
        else:
            bad_bender_samples = bad_bender

        sample_df['bad_bender'] = False
        if len(bad_bender_samples) != 0: 
           # add bad_bender column to sample_df
          bad_bender_df = pd.DataFrame({
            'sample_id': bad_bender_samples, 
            'n_cells': np.nan, 
            'bad_qc': False, 
            'bad_bender': True
          })

          sample_df = pd.concat([sample_df, bad_bender_df], ignore_index=True)

        assert all([s in sample_df['sample_id'].tolist() for s in SAMPLES])
        # label as bad if bad_bender or bad_qc
        sample_df['bad_sample'] = sample_df['bad_bender'] | sample_df['bad_qc']
    else:
        sample_df['bad_sample'] = sample_df['bad_qc']

    return sample_df


# get output file paths as string
def get_qc_files_str(run, SAMPLE_MAPPING, qc_dir, FULL_TAG, DATE_STAMP):
  if SAMPLE_MAPPING is None:
    sce_str = f"{qc_dir}/sce_cells_tmp_{run}_{FULL_TAG}_{DATE_STAMP}.rds"
    smpl_str= run
  else:
    sce_fs_ls = []
    smpls_ls = []
    for s in SAMPLE_MAPPING[run]:
      sce_fs_ls.append(f"{qc_dir}/sce_cells_tmp_{s}_{FULL_TAG}_{DATE_STAMP}.rds")
      smpls_ls.append(s)

    sce_str  = ','.join(sce_fs_ls)
    smpl_str = ','.join(smpls_ls)
  
  return smpl_str, sce_str
    

rule run_qc:
  input:
    ambient_stats_f = amb_dir + '/ambient_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    amb_yaml_f   = amb_dir + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml',
    demux_f      = (demux_dir + '/sce_cells_htos_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.rds') if DEMUX_TYPE == 'af' else ([DEMUX_F] if DEMUX_TYPE == 'custom' else [])
  output:
    qc_f         = temp(qc_dir  + '/qc_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'), 
    coldata_f    = temp(qc_dir + '/coldata_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'),
    rowdata_f    = temp(qc_dir + '/rowdata_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'), 
    dimred_f     = dbl_dir + '/dbl_{run}/scDblFinder_{run}_dimreds_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    dbl_f        = dbl_dir + '/dbl_{run}/scDblFinder_{run}_outputs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'    
  params:
    all_smpls_str   = lambda wildcards: get_qc_files_str(wildcards.run, SAMPLE_MAPPING, qc_dir, FULL_TAG, DATE_STAMP)[0],
    sce_fs_str      = lambda wildcards: get_qc_files_str(wildcards.run, SAMPLE_MAPPING, qc_dir, FULL_TAG, DATE_STAMP)[1],
    mito_str        = AF_MITO_STR,
    exclude_mito    = EXCLUDE_MITO,
    hard_min_counts = QC_HARD_MIN_COUNTS,
    hard_min_feats  = QC_HARD_MIN_FEATS,
    hard_max_mito   = QC_HARD_MAX_MITO,
    min_counts      = QC_MIN_COUNTS,
    min_feats       = QC_MIN_FEATS,
    min_mito        = QC_MIN_MITO,
    max_mito        = QC_MAX_MITO,
    min_splice      = QC_MIN_SPLICE,
    max_splice      = QC_MAX_SPLICE,
    sample_var      = SAMPLE_VAR,
    demux_type      = DEMUX_TYPE,
    dbl_min_feats   = DBL_MIN_FEATS
  threads: 4
  retries: RETRIES
  resources:
    mem_mb = lambda wildcards, attempt: attempt * MB_RUN_QC
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_qc/run_qc_{run}_' + DATE_STAMP + '.benchmark.txt'
  conda:
    '../envs/rlibs.yaml'
  shell:
    """
    Rscript -e "source('scripts/SampleQC.R'); source('scripts/ambient.R'); \
        main_qc( \
          sel_sample     = '{wildcards.run}', \
          meta_f         = '{METADATA_F}', \
          amb_yaml_f     = '{input.amb_yaml_f}', \
          sample_stats_f = '{input.ambient_stats_f}', \
          demux_f        = '{input.demux_f}', \
          gtf_dt_f       = '{AF_GTF_DT_F}', \
          ambient_method = '{AMBIENT_METHOD}', \
          sce_fs_str     = '{params.sce_fs_str}', \
          all_samples_str = '{params.all_smpls_str}', \
          rowdata_f       = '{output.rowdata_f}', \
          qc_f           = '{output.qc_f}', \
          coldata_f      = '{output.coldata_f}', \
          dimred_f       = '{output.dimred_f}', \
          dbl_f		       = '{output.dbl_f}', \
          mito_str       = '{params.mito_str}', \
          exclude_mito   = '{params.exclude_mito}', \
          hard_min_counts= {params.hard_min_counts}, \
          hard_min_feats = {params.hard_min_feats}, \
          hard_max_mito  = {params.hard_max_mito}, \
          min_counts     = {params.min_counts}, \
          min_feats      = {params.min_feats}, \
          min_mito       = {params.min_mito}, \
          max_mito       = {params.max_mito}, \
          min_splice     = {params.min_splice}, \
          max_splice     = {params.max_splice}, \
          sample_var     = '{params.sample_var}', \
          demux_type     = '{params.demux_type}', \
          dbl_min_feats  = {params.dbl_min_feats} \
        )"

    """


rule merge_qc:
  input:
    qc_fs      = expand(qc_dir  + '/qc_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', run = runs),
    coldata_fs = expand(qc_dir  + '/coldata_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', run = runs)
  output:
    qc_merged_f      = qc_dir  + '/qc_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    coldata_merged_f = qc_dir  + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  threads: 1
  retries: RETRIES
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_qc/merge_qc_' + DATE_STAMP + '.benchmark.txt'
  resources:
    mem_mb = lambda wildcards, attempt: attempt * MB_RUN_QC
  run:
    # read all nonempty input files and concatenate them
    qc_df_ls = [
      pd.read_csv(f, compression='gzip') 
      for f in input.qc_fs 
      if gzip.open(f, 'rb').read(1)  
    ]
    qc_df_merged = pd.concat(qc_df_ls, ignore_index=True)

    coldata_df_ls = [
      pd.read_csv(f, compression='gzip') 
      for f in input.coldata_fs 
      if gzip.open(f, 'rb').read(1)  
    ]
    coldata_df_merged = pd.concat(coldata_df_ls, ignore_index=True)

    # save merged dataframes to output files
    qc_df_merged.to_csv(output.qc_merged_f, sep='\t', index=False, compression='gzip')
    coldata_df_merged.to_csv(output.coldata_merged_f, sep='\t', index=False, compression='gzip')


rule merge_rowdata:
  input:
    rowdata_fs = expand(qc_dir  + '/rowdata_dt_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', run = runs)
  output:
    rowdata_merged_f = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  threads: 1
  retries: RETRIES
  resources:
    mem_mb = lambda wildcards, attempt: attempt * MB_RUN_QC
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_qc/merge_rowdata_' + DATE_STAMP + '.benchmark.txt'
  run:
    # read all nonempty rowdata files 
    rd_df_ls = [
      pd.read_csv(f,compression='gzip') 
      for f in input.rowdata_fs 
      if gzip.open(f, 'rb').read(1)  # Check if file is not empty 
    ]
    
    # check if identical
    first_rd = rd_df_ls[0]
    all_ident= all(first_rd.equals(df) for df in rd_df_ls[1:])
    
    # save only one df
    if all_ident:
      first_rd.to_csv(output.rowdata_merged_f, sep='\t', index=False, compression='gzip')
    else:
      raise ValueError("Error: rowdata for all sce objects not identical.")


rule get_qc_sample_statistics:
  input:
    ambient_stats_f = amb_dir + '/ambient_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    qc_merged_f     = qc_dir  + '/qc_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz' 
  output:
    qc_stats_f      = qc_dir + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  threads: 1
  retries: RETRIES
  resources:
    mem_mb = lambda wildcards, attempt: attempt * MB_RUN_QC
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_qc/get_qc_sample_statistics_' + DATE_STAMP + '.benchmark.txt'
  run:
    sample_stats_df = extract_qc_sample_statistics(input.ambient_stats_f, input.qc_merged_f, SAMPLES, SAMPLE_VAR, AMBIENT_METHOD, DEMUX_TYPE, SAMPLE_MAPPING, QC_MIN_CELLS)
    sample_stats_df.to_csv(output.qc_stats_f, index = False)


# write sce objects paths to a yaml file
rule make_tmp_sce_paths_yaml:
   input:
    qc_stats_f  = qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv' # so that this runs after get_qc_sample_statistics
   output:
    sces_yaml_f = temp(qc_dir  + '/sce_tmp_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml')
   threads: 1
   retries: RETRIES
   run:
    # split paths and sample names
    fs = [f"{qc_dir}/sce_cells_tmp_{s}_{FULL_TAG}_{DATE_STAMP}.rds" for s in SAMPLES]
    
    # check that all files exist
    for f in fs:
     assert os.path.isfile(f), \
      f"File {f} doesn't exist"

    # create a dictionary
    fs_dict = dict(zip(SAMPLES, fs))

    # write to yaml
    with open(output.sces_yaml_f, 'w') as f:
     yaml.dump(fs_dict, f, default_flow_style=False)

