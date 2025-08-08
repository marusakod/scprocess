
def extract_zoom_sample_statistics(qc_stats_f, LABELS_F, LABELS_VAR, LABELS, MIN_N_SAMPLE, AMBIENT_METHOD):
    
    # load inputs
    qc_df   = pd.read_csv(qc_stats_f)
    qc_df   = qc_df.drop('n_cells', axis=1)
    lbls_dt = pd.read_csv(LABELS_F, compression='gzip')
 
    # keep selected labels
    lbls_dt = lbls_dt[lbls_dt[LABELS_VAR].isin(LABELS)]
    
    # count the number of cells per sample
    zoom_sample_stats = (
      
        lbls_dt.groupby('sample_id')
        .size()
        .reset_index(name='n_cells')
    )
    
    # identify samples that do not meet the minimum cell threshold
    zoom_sample_stats['bad_zoom_qc'] = zoom_sample_stats['n_cells'] < MIN_N_SAMPLE
    
    # merge new and existing sample stats
    sample_df = qc_df.merge(zoom_sample_stats, on='sample_id',how='left')
    
    # update 'bad_sample' column
    if AMBIENT_METHOD == 'cellbender':
        sample_df['bad_sample'] = (
            sample_df['bad_bender'] | sample_df['bad_qc'] | sample_df['bad_zoom_qc']
        )
    else:
        sample_df['bad_sample'] = (
            sample_df['bad_qc'] | sample_df['bad_zoom_qc']
        )

    # check that at least 2 good samples remain
    good_smpls_count = (sample_df['bad_sample'] == False).sum()
    assert good_smpls_count >= 2, \
      "Fewer than 2 samples available for this zoom."
    
    return sample_df



rule get_zoom_sample_statistics:
  input:
    qc_stats_f      = qc_dir + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    zoom_stats_f    = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  params:
    zoom_lbls_f     = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['LABELS_F'],
    zoom_lbls_var   = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['LABELS_VAR'], 
    zoom_lbls       = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['LABELS'],
    zoom_min_n_smpl = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['MIN_N_SAMPLE']
  run:
    zoom_stats_df = extract_zoom_sample_statistics(input.qc_stats_f, params.zoom_lbls_f, params.zoom_lbls_var, params.zoom_lbls, params.zoom_min_n_smpl, AMBIENT_METHOD)
    zoom_stats_df.to_csv(output.zoom_stats_f, index = False)



# pseudobulks and empties
rule zoom_make_pb_subset:
  input:
    sces_yaml_f  = int_dir  + '/sce_clean_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml', 
    zoom_stats_f = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    zoom_pb_subset_f  = zoom_dir + '/{zoom_name}/pb_{zoom_name}_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  params:
    zoom_lbls_f     = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['LABELS_F'],
    zoom_lbls_var   = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['LABELS_VAR'],
    zoom_lbls       = lambda wildcards: ','.join(ZOOM_PARAMS_DICT[wildcards.zoom_name]['LABELS'])
  threads: 4
  retries: RETRIES 
  resources:
    mem_mb  = lambda wildcards, attempt: attempt * MB_PB_MAKE_PBS
  conda: 
    '../envs/rlibs.yaml'
  shell:
    """
    Rscript -e "source('scripts/utils.R'); source('scripts/ambient.R'); source('scripts/pseudobulk_and_empties.R'); \
    make_pb_cells( \
      sce_fs_yaml = '{input.sces_yaml_f}',
      qc_stats_f  = '{input.zoom_stats_f}',
      subset_f    = '{params.zoom_lbls_f}',
      subset_col  = '{params.zoom_lbls_var}', 
      subset_str  = '{params.zoom_lbls}', 
      pb_f        = '{output.zoom_pb_subset_f}',
      n_cores     = {threads})"
    """


rule zoom_calculate_ambient_genes:
  input:
    pb_empty_f       = pb_dir + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
    zoom_pb_subset_f = zoom_dir + '/{zoom_name}/pb_{zoom_name}_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  output:
    zoom_empty_gs_f  = zoom_dir + '/{zoom_name}/edger_empty_genes_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  params:
    zoom_fdr_thr     = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['AMBIENT_GENES_FDR_THR'],
    zoom_logfc_thr   = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['AMBIENT_GENES_LOGFC_THR']
  threads: 4
  retries: RETRIES 
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * MB_PB_MAKE_PBS
  conda: 
    '../envs/rlibs.yaml'
  shell:
    """
    Rscript -e "source('scripts/utils.R'); source('scripts/pseudobulk_and_empties.R'); \
    calc_empty_genes(
      pb_cells_f = '{input.zoom_pb_subset_f}',
      pb_empty_f = '{input.pb_empty_f}',
      fdr_thr    = {params.zoom_fdr_thr}, 
      logfc_thr  = {params.zoom_logfc_thr},
      empty_gs_f = '{output.zoom_empty_gs_f}')"
    """

# highly variable genes

rule zoom_make_hvg_df:
    input:
        ambient_yaml_out=expand([amb_dir + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml'], run=runs)
    output:
        hvg_paths_f = zoom_dir + '/{zoom_name}/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv' 
    run:
        hvg_df = make_hvgs_input_df(
            DEMUX_TYPE, SAMPLE_VAR, runs, input.ambient_yaml_out,
            SAMPLE_MAPPING, FULL_TAG, DATE_STAMP, f"{zoom_dir}/{wildcards.zoom_name}"
            )
        hvg_df.to_csv(output.hvg_paths_f, index=False)



rule zoom_make_tmp_csr_matrix:
  input:
    hvg_paths_f        = zoom_dir + '/{zoom_name}/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    smpl_stats_f       = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    rowdata_f          = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    clean_h5_f  = zoom_dir + '/{zoom_name}' + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5'
  params: 
    zoom_lbls_f     = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['LABELS_F'],
    zoom_lbls_var   = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['LABELS_VAR'],
    zoom_lbls       = lambda wildcards: ','.join(ZOOM_PARAMS_DICT[wildcards.zoom_name]['LABELS']), 
    zoom_chunk_size = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['HVG_CHUNK_SIZE']
  threads: 8
  retries: RETRIES
  resources:
    mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
  conda:
    '../envs/hvgs.yaml'
  shell:
    """
    python3 scripts/hvgs.py get_csr_counts \
      {input.hvg_paths_f} \
      {params.zoom_lbls_f} \
      {params.zoom_lbls_var} \
      {params.zoom_lbls} \
      {input.smpl_stats_f} \
      {input.rowdata_f} \
      {SAMPLE_VAR} \
      {DEMUX_TYPE} \
      --size {params.zoom_chunk_size} \
      --ncores {threads}
    """



# if ZOOM_PARAMS_DICT[wildcards.zoom_name]["HVG_METHOD"] == "sample": 
#   localrules: zoom_merge_sample_std_var_stats
#   # calculate stats for each sample separatelly  
#   rule zoom_get_stats_for_std_variance_for_sample:
#     input: 
#       clean_h5_f      = zoom_dir + '/{zoom_name}' + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
#       smpl_stats_f    = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
#       rowdata_f       = qc_dir   + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
#     output:
#       std_var_stats_f = zoom_dir + '/{zoom_name}' + '/tmp_std_var_stats_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
#     threads: 1
#     retries: RETRIES
#     resources:
#       mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
#     conda:
#       '../envs/hvgs.yaml'
#     shell:
#       """
#       python3 scripts/hvgs.py calculate_std_var_stats_for_sample \
#         {wildcards.sample} \
#         {input.smpl_stats_f} \
#         {input.clean_h5_f} \
#         {input.rowdata_f} \
#         {output.std_var_stats_f}
#       """

#   rule zoom_merge_sample_std_var_stats:
#     input:                 
#       std_var_stats_f = expand(zoom_dir + '/{zoom_name}' + '/tmp_std_var_stats_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', sample = SAMPLES)
#     output:
#       std_var_stats_merged_f = zoom_dir + '/{zoom_name}' + '/standardized_variance_stats_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
#     threads: 1
#     retries: RETRIES
#     resources:
#       mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
#     run:
#       merge_tmp_files(input.std_var_stats_f, output.std_var_stats_merged_f)


# else:
#   localrules: merge_group_mean_var, merge_group_std_var_stats

#   rule get_mean_var_for_group:
#     input:
#       clean_h5_f      = expand(hvg_dir + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5', sample = SAMPLES),
#       hvg_paths_f     = hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
#       rowdata_f       = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
#       qc_smpl_stats_f = qc_dir + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
#     output: 
#       mean_var_f      = temp(hvg_dir + '/tmp_mean_var_{group}_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
#     threads: 8
#     resources:
#       mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
#     conda:
#       '../envs/hvgs.yaml'
#     shell:
#       """
#       python3 scripts/hvgs.py calculate_mean_var_for_chunk \
#         {input.hvg_paths_f} \
#         {input.rowdata_f} \
#         {METADATA_F} \
#         {input.qc_smpl_stats_f} \
#         {output.mean_var_f} \
#         {wildcards.chunk} \
#         {HVG_METHOD} \
#         {HVG_CHUNK_SIZE} \
#         --group {wildcards.group} \
#         --groupvar {HVG_GROUP_VAR} \
#         --ncores {threads} 
    
#       """
#   rule merge_group_mean_var:
#     input:                 
#       mean_var_f  = expand(hvg_dir + '/tmp_mean_var_{group}_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
#       group=GROUP_NAMES, chunk=range(NUM_CHUNKS))
#     output:
#       mean_var_merged_f = temp(hvg_dir + '/means_variances_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
#     threads: 1
#     retries: RETRIES
#     resources:
#       mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
#     run:
#       merge_tmp_files(input.mean_var_f, output.mean_var_merged_f)


#   rule get_estimated_variances:
#     input:
#       mean_var_merged_f  = hvg_dir + '/means_variances_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
#     output:
#       estim_vars_f       = temp(hvg_dir + '/estimated_variances_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
#     threads: 1
#     retries: RETRIES
#     conda:
#       '../envs/hvgs.yaml'
#     resources:
#       mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
#     shell:
#       """
#       python3 scripts/hvgs.py calculate_estimated_vars \
#         {output.estim_vars_f} \
#         {HVG_METHOD} \
#         {input.mean_var_merged_f}
#       """

#   rule get_stats_for_std_variance_for_group:
#     input: 
#       clean_h5_fs     = expand(hvg_dir + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5', sample = SAMPLES),
#       estim_vars_f    = hvg_dir + '/estimated_variances_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
#       hvg_paths_f     = hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
#       rowdata_f       = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
#       qc_smpl_stats_f = qc_dir + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
#     output:
#       std_var_stats_f = temp(hvg_dir + '/tmp_std_var_stats_{group}_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
#     threads: 8
#     resources:
#       mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
#     conda:
#       '../envs/hvgs.yaml'
#     shell:
#       """
#       python3 scripts/hvgs.py calculate_std_var_stats_for_chunk \
#         {input.hvg_paths_f} \
#         {input.rowdata_f} \
#         {METADATA_F} \
#         {input.qc_smpl_stats_f} \
#         {output.std_var_stats_f} \
#         {input.estim_vars_f} \
#         {wildcards.chunk} \
#         {HVG_METHOD} \
#         --size {HVG_CHUNK_SIZE} \
#         --group {wildcards.group} \
#         --groupvar {HVG_GROUP_VAR} \
#         --ncores {threads} \
    
#       """

#   rule merge_group_std_var_stats:
#     input:                 
#       std_var_stats_f = expand(hvg_dir + '/tmp_std_var_stats_{group}_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
#       group=GROUP_NAMES, chunk=range(NUM_CHUNKS)),
#     output:
#       std_var_stats_merged_f = temp(hvg_dir + '/standardized_variance_stats_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
#     threads: 1
#     retries: RETRIES
#     resources:
#       mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
#     run:
#       merge_tmp_files(input.std_var_stats_f, output.std_var_stats_merged_f)


# rule get_highly_variable_genes:
#   input:
#     std_var_stats_f = hvg_dir + '/standardized_variance_stats_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
#     empty_gs_fs     = empty_dir + '/edger_empty_genes_' + FULL_TAG + '_all_' + DATE_STAMP + '.txt.gz' 
#   output:
#     hvg_f = hvg_dir + '/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
#   threads: 1
#   retries: RETRIES
#   params:
#     hvg_method = HVG_METHOD,
#     n_hvgs = N_HVGS,
#     exclude_ambient_genes = EXCLUDE_AMBIENT_GENES
#   resources:
#      mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
#   conda:
#     '../envs/hvgs.yaml'
#   shell:
#      """
#      NOAMBIENT_FLAG=""
#      if [ "{params.exclude_ambient_genes}" = "True" ]; then
#        NOAMBIENT_FLAG="--noambient"
#      fi

#      python3 scripts/hvgs.py calculate_hvgs \
#       {input.std_var_stats_f} \
#       {output.hvg_f} \
#       {input.empty_gs_fs} \
#       {params.hvg_method} \
#       {params.n_hvgs} \
#       $NOAMBIENT_FLAG
#      """


# rule create_hvg_matrix:
#   input: 
#     clean_h5_f      = expand(hvg_dir + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5', sample = SAMPLES),
#     qc_smpl_stats_f = qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
#     hvg_paths_f     = hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
#     hvg_f           = hvg_dir + '/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
#   output:
#     hvg_mat_f       = hvg_dir + '/top_hvgs_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5'
#   threads: 1
#   retries: RETRIES
#   resources:
#     mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
#   conda:
#     '../envs/hvgs.yaml'
#   shell:
#     """
#     python3 scripts/hvgs.py read_top_genes \
#       {input.qc_smpl_stats_f} \
#       {input.hvg_paths_f} \
#       {input.hvg_f} \
#       {output.hvg_mat_f} \
#       {SAMPLE_VAR}

#     """


# rule create_doublets_hvg_matrix:
#     input: 
#       hvg_paths_f        = hvg_dir + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
#       hvg_f              = hvg_dir + '/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
#       qc_f               = qc_dir  + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
#       qc_sample_stats_f  = qc_dir  + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
#     output: 
#       dbl_hvg_mat_f      = hvg_dir + '/top_hvgs_doublet_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5'
#     threads: 1
#     retries: RETRIES
#     resources:
#       mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
#     conda: 
#       '../envs/hvgs.yaml'
#     shell: 
#       """
#       python3 scripts/hvgs.py create_doublets_matrix \
#       {input.hvg_paths_f} \
#       {input.hvg_f} \
#       {input.qc_f} \
#       {input.qc_sample_stats_f} \
#       {output.dbl_hvg_mat_f} \
#       {SAMPLE_VAR}

#       """
