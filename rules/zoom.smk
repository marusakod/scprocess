
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


def get_mean_var_input(zoom_name):
    group_names = ZOOM_PARAMS_DICT[zoom_name]['HVG_GROUP_NAMES']
    num_chunks = ZOOM_PARAMS_DICT[zoom_name]['HVG_NUM_CHUNKS']

    return [
        zoom_dir + f'/{zoom_name}/tmp_mean_var_{group}_group_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
        for group in group_names
        for chunk in range(num_chunks)
    ]


def get_tmp_std_var_stats_input(zoom_name):
    hvg_method = ZOOM_PARAMS_DICT[zoom_name]['HVG_METHOD']

    if hvg_method == "sample":
        return [
            zoom_dir + f'/{zoom_name}/tmp_std_var_stats_{sample}_sample_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
            for sample in SAMPLES
        ]
    else:
        group_names = ZOOM_PARAMS_DICT[zoom_name]['HVG_GROUP_NAMES']
        num_chunks = ZOOM_PARAMS_DICT[zoom_name]['HVG_NUM_CHUNKS']

        return [
            zoom_dir + f'/{zoom_name}/tmp_std_var_stats_{group}_group_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
            for group in group_names
            for chunk in range(num_chunks)
        ]


def get_zoom_conditional_outputs(species):
  if species in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020']:
    return {
      'fgsea_go_bp_f': zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_go_bp_' + DATE_STAMP + '.txt.gz', 
      'fgsea_go_cc_f': zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_go_cc_' + DATE_STAMP + '.txt.gz',
      'fgsea_go_mf_f': zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_go_mf_' + DATE_STAMP + '.txt.gz',
      'fgsea_paths_f': zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_paths_' + DATE_STAMP + '.txt.gz',
      'fgsea_hlmk_f':  zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_hlmk_' + DATE_STAMP + '.txt.gz'
    }
  else:
    return {}


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
    clean_h5_f  = temp(zoom_dir + '/{zoom_name}' + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5')
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


rule zoom_get_stats_for_std_variance_for_sample:
    input: 
        clean_h5_f      = zoom_dir + '/{zoom_name}' + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
        smpl_stats_f    = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
        rowdata_f       = qc_dir   + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
    output:
        std_var_stats_f = temp(zoom_dir + '/{zoom_name}' + '/tmp_std_var_stats_{sample}_sample_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
    threads: 1
    retries: RETRIES
    resources:
        mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
    conda:
        '../envs/hvgs.yaml'
    shell:
        """
        python3 scripts/hvgs.py calculate_std_var_stats_for_sample \
          {wildcards.sample} \
          {input.smpl_stats_f} \
          {input.clean_h5_f} \
          {input.rowdata_f} \
          {output.std_var_stats_f}
        """


rule zoom_get_mean_var_for_group:
    input:
        clean_h5_f= expand(
            zoom_dir + '/{zoom_name}' + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5',
            zoom_name = ZOOM_NAMES, sample=SAMPLES
        ),
        hvg_paths_f     = zoom_dir + '/{zoom_name}' + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
        rowdata_f       = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
        smpl_stats_f    = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
    output: 
        mean_var_f      = temp(zoom_dir + '/{zoom_name}' + '/tmp_mean_var_{group}_group_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
    params:
        zoom_hvg_method = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['HVG_METHOD'], 
        zoom_hvg_chunk_size = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['HVG_CHUNK_SIZE'], 
        zoom_group_var = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['HVG_SPLIT_VAR']
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
    conda:
        '../envs/hvgs.yaml'
    shell:
        """
        python3 scripts/hvgs.py calculate_mean_var_for_chunk \
          {input.hvg_paths_f} \
          {input.rowdata_f} \
          {METADATA_F} \
          {input.smpl_stats_f} \
          {output.mean_var_f} \
          {wildcards.chunk} \
          {params.zoom_hvg_method} \
          {params.zoom_hvg_chunk_size} \
          --group {wildcards.group} \
          --groupvar {params.zoom_group_var} \
          --ncores {threads} 
        """


rule zoom_merge_group_mean_var:
    input:                 
        mean_var_f= lambda wildcards: get_mean_var_input(wildcards.zoom_name)
    output:
        mean_var_merged_f = temp(zoom_dir + '/{zoom_name}'  + '/means_variances_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
    threads: 1
    retries: RETRIES
    resources:
        mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
    run:
        merge_tmp_files(input.mean_var_f, output.mean_var_merged_f)


rule zoom_get_estimated_variances:
    input:
        mean_var_merged_f  = zoom_dir + '/{zoom_name}' + '/means_variances_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
    output:
        estim_vars_f       = temp(zoom_dir + '/{zoom_name}' + '/estimated_variances_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
    params: 
        zoom_hvg_method = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['HVG_METHOD']
    threads: 1
    retries: RETRIES
    conda:
        '../envs/hvgs.yaml'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
    shell:
        """
        python3 scripts/hvgs.py calculate_estimated_vars \
          {output.estim_vars_f} \
          {params.zoom_hvg_method} \
          {input.mean_var_merged_f}
        """


rule zoom_get_stats_for_std_variance_for_group:
    input: 
        clean_h5_f= expand(
            zoom_dir + '/{zoom_name}' + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5',
            zoom_name = ZOOM_NAMES, sample=SAMPLES
        ),
        estim_vars_f    = zoom_dir + '/{zoom_name}'  + '/estimated_variances_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
        hvg_paths_f     = zoom_dir + '/{zoom_name}'  + '/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
        rowdata_f       = qc_dir  + '/rowdata_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
        smpl_stats_f    = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
    output:
        std_var_stats_f = temp(zoom_dir + '/{zoom_name}' + '/tmp_std_var_stats_{group}_group_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz')
    params:
        zoom_hvg_method = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['HVG_METHOD'], 
        zoom_hvg_chunk_size = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['HVG_CHUNK_SIZE'], 
        zoom_hvg_group_var = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['HVG_SPLIT_VAR']
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
    conda:
        '../envs/hvgs.yaml'
    shell:
        """
        python3 scripts/hvgs.py calculate_std_var_stats_for_chunk \
          {input.hvg_paths_f} \
          {input.rowdata_f} \
          {METADATA_F} \
          {input.smpl_stats_f} \
          {output.std_var_stats_f} \
          {input.estim_vars_f} \
          {wildcards.chunk} \
          {params.zoom_hvg_method} \
          --size {params.zoom_hvg_chunk_size} \
          --group {wildcards.group} \
          --groupvar {params.zoom_hvg_group_var} \
          --ncores {threads} 
        """


rule zoom_merge_stats_for_std_variance:
    input:
        tmp_std_var_stats_fs=lambda wildcards: get_tmp_std_var_stats_input(wildcards.zoom_name)
    output:
        std_var_stats_merged_f= zoom_dir + '/{zoom_name}/standardized_variance_stats_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * MB_RUN_HVGS
    run:
        merge_tmp_files(input.tmp_std_var_stats_fs, output.std_var_stats_merged_f)

        
rule zoom_get_highly_variable_genes:
  input:
    std_var_stats_f = zoom_dir + '/{zoom_name}/standardized_variance_stats_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    empty_gs_fs     = zoom_dir + '/{zoom_name}/edger_empty_genes_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz' 
  output:
    hvg_f =  zoom_dir + '/{zoom_name}/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  threads: 1
  retries: RETRIES
  params:
    zoom_hvg_method = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['HVG_METHOD'],
    zoom_n_hvgs     = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['N_HVGS'],
    zoom_exclude_ambient_genes = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['EXCLUDE_AMBIENT_GENES']
  resources:
     mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
  conda:
    '../envs/hvgs.yaml'
  shell:
     """
     NOAMBIENT_FLAG=""
     if [ "{params.zoom_exclude_ambient_genes}" = "True" ]; then
       NOAMBIENT_FLAG="--noambient"
     fi

     python3 scripts/hvgs.py calculate_hvgs \
      {input.std_var_stats_f} \
      {output.hvg_f} \
      {input.empty_gs_fs} \
      {params.zoom_hvg_method} \
      {params.zoom_n_hvgs} \
      $NOAMBIENT_FLAG
     """



rule zoom_create_hvg_matrix:
  input: 
    clean_h5_f   = expand(
            zoom_dir + '/{zoom_name}' + '/chunked_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.h5',
            zoom_name = ZOOM_NAMES, sample=SAMPLES
    ),
    smpl_stats_f = zoom_dir  + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    hvg_paths_f  = zoom_dir  + '/{zoom_name}/hvg_paths_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    hvg_f        = zoom_dir  + '/{zoom_name}/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    hvg_mat_f    = zoom_dir  + '/{zoom_name}/top_hvgs_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5'
  threads: 1
  retries: RETRIES
  resources:
    mem_mb = lambda wildcards, attempt: attempt * MB_RUN_HVGS
  conda:
    '../envs/hvgs.yaml'
  shell:
    """
    python3 scripts/hvgs.py read_top_genes \
      {input.smpl_stats_f} \
      {input.hvg_paths_f} \
      {input.hvg_f} \
      {output.hvg_mat_f} \
      {SAMPLE_VAR}

    """


rule zoom_run_integration:
  input:
    hvg_mat_f     = zoom_dir + '/{zoom_name}/top_hvgs_counts_' + FULL_TAG + '_' + DATE_STAMP + '.h5', 
    smpl_stats_f  = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv', 
    coldata_f     = qc_dir   + '/coldata_dt_all_samples_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    integration_f = zoom_dir + '/{zoom_name}/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  params:
    zoom_int_reduction = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['INT_REDUCTION'],
    zoom_int_n_dims    = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['INT_N_DIMS'],
    zoom_int_cl_method = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['INT_CL_METHOD'], 
    zoom_int_theta     = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['INT_THETA'],
    zoom_int_res_ls    = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['INT_RES_LS']
  threads: 8
  retries: RETRIES 
  resources:
    mem_mb   = lambda wildcards, attempt: attempt * MB_RUN_INTEGRATION
  conda: 
    '../envs/rlibs.yaml'
  shell:
    """
    # run harmony
    Rscript -e "source('scripts/integration.R'); source('scripts/ambient.R'); source('scripts/zoom.R'); 
      run_zoom_integration( 
        hvg_mat_f        = '{input.hvg_mat_f}', 
        smpl_stats_f     = '{input.smpl_stats_f}', 
        coldata_f        = '{input.coldata_f}', 
        demux_type       = '{DEMUX_TYPE}', 
        exclude_mito     = '{EXCLUDE_MITO}', 
        reduction        = '{params.zoom_int_reduction}',
        n_dims           = {params.zoom_int_n_dims}, 
        cl_method        = '{params.zoom_int_cl_method}', 
        theta            = {params.zoom_int_theta}, 
        res_ls_concat    = '{params.zoom_int_res_ls}', 
        integration_f    = '{output.integration_f}', 
        batch_var        = '{BATCH_VAR}', 
        n_cores          = {threads})"
    """


rule zoom_run_marker_genes:
  input:
    sces_yaml_f    = int_dir  + '/sce_clean_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml',
    integration_f  = zoom_dir + '/{zoom_name}/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    pb_f      = zoom_dir + '/{zoom_name}/pb_' + FULL_TAG + '_' + '{mkr_sel_res}_' + DATE_STAMP + '.rds',
    mkrs_f    = zoom_dir + '/{zoom_name}/pb_marker_genes_' + FULL_TAG + '_' + '{mkr_sel_res}_' + DATE_STAMP + '.txt.gz',
    pb_hvgs_f = zoom_dir + '/{zoom_name}/pb_hvgs_' + FULL_TAG + '_' + '{mkr_sel_res}_' + DATE_STAMP + '.txt.gz', 
    **get_zoom_conditional_outputs(SPECIES)
  params: 
    fgsea_args = lambda wildcards, output: ", ".join([
        f"fgsea_go_bp_f = '{output.get('fgsea_go_bp_f', '')}'",
        f"fgsea_go_cc_f = '{output.get('fgsea_go_cc_f', '')}'",
        f"fgsea_go_mf_f = '{output.get('fgsea_go_mf_f', '')}'",
        f"fgsea_paths_f = '{output.get('fgsea_paths_f', '')}'",
        f"fgsea_hlmk_f = '{output.get('fgsea_hlmk_f', '')}',"
    ]), 
    zoom_mkr_sel_res     = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['MKR_SEL_RES'],
    zoom_mkr_min_cl_size = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['MKR_MIN_CL_SIZE'], 
    zoom_mkr_min_cells   = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['MKR_MIN_CELLS'], 
    zoom_mkr_not_ok_re   = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['MKR_NOT_OK_RE'],
    zoom_mkr_min_cpm_go  = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['MKR_MIN_CPM_GO'],
    zoom_mkr_max_zero_p  = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['MKR_MAX_ZERO_P'], 
    zoom_mkr_gsea_cut    = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['MKR_GSEA_CUT']
  threads: 8
  retries: RETRIES
  resources:
    mem_mb = lambda wildcards, attempt: attempt * MB_RUN_MARKER_GENES
  conda:
    '../envs/rlibs.yaml'
  shell:
    """
    Rscript -e "source('scripts/utils.R'); source('scripts/marker_genes.R'); calculate_marker_genes(
      integration_f = '{input.integration_f}', 
      sces_yaml_f   = '{input.sces_yaml_f}',
      pb_f          = '{output.pb_f}',
      mkrs_f        = '{output.mkrs_f}',
      pb_hvgs_f     = '{output.pb_hvgs_f}',
      {params.fgsea_args}
      species       = '{SPECIES}',
      gtf_dt_f      = '{AF_GTF_DT_F}',
      gsea_dir      = '{MKR_GSEA_DIR}',
      sel_res       = '{params.zoom_mkr_sel_res}',
      min_cl_size   = {params.zoom_mkr_min_cl_size},
      min_cells     = {params.zoom_mkr_min_cells},
      not_ok_re     = '{params.zoom_mkr_not_ok_re}',
      min_cpm_go    = {params.zoom_mkr_min_cpm_go},
      max_zero_p    = {params.zoom_mkr_max_zero_p},
      gsea_cut      = {params.zoom_mkr_gsea_cut},
      zoom          = 'True', 
      n_cores       = {threads})"
    """


rule zoom_make_subset_sces:
  input:
    integration_f = zoom_dir + '/{zoom_name}/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    smpl_stats_f  = zoom_dir + '/{zoom_name}/zoom_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv',
    sces_yaml_f   = int_dir + '/sce_clean_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml'
  output:
    clean_sce_f = zoom_dir + '/{zoom_name}/sce_objects/sce_cells_clean_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  params:
    zoom_lbls_f     = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['LABELS_F'],
    zoom_lbls_var   = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['LABELS_VAR'],
    zoom_lbls       = lambda wildcards: ','.join(ZOOM_PARAMS_DICT[wildcards.zoom_name]['LABELS'])
  threads: 1
  retries: RETRIES
  resources:
    mem_mb   = lambda wildcards, attempt: attempt * MB_MAKE_SUBSET_SCES
  conda: 
    '../envs/rlibs.yaml'
  shell:
    """
    Rscript -e "source('scripts/zoom.R');
     make_subset_sces(
     sel_s          = '{wildcards.sample}',
     clean_sce_f    = '{output.clean_sce_f}',
     integration_f  = '{input.integration_f}',
     smpl_stats_f   = '{input.smpl_stats_f}',
     sces_yaml_f    = '{input.sces_yaml_f}',
     subset_f       = '{params.zoom_lbls_f}',
     subset_col     = '{params.zoom_lbls_var}',
     subset_str     = '{params.zoom_lbls}')"
    """


# render_html_zoom
rule render_html_zoom:
  input:
    r_utils_f           = code_dir + '/utils.R',
    r_hvgs_f            = code_dir + '/hvgs.R', 
    r_int_f             = code_dir + '/integration.R',
    r_mkr_f             = code_dir + '/marker_genes.R',
    zoom_int_f          = zoom_dir + '/{zoom_name}/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    zoom_pb_f           = zoom_dir + '/{zoom_name}/pb_' + FULL_TAG + '_' + '{mkr_sel_res}_' + DATE_STAMP + '.rds',
    zoom_mkrs_f         = zoom_dir + '/{zoom_name}/pb_marker_genes_' + FULL_TAG + '_' + '{mkr_sel_res}_' + DATE_STAMP + '.txt.gz',
    zoom_mkrs_hvgs_f    = zoom_dir + '/{zoom_name}/pb_hvgs_' + FULL_TAG + '_' + '{mkr_sel_res}_' + DATE_STAMP + '.txt.gz',
    zoom_hvgs_f         = zoom_dir + '/{zoom_name}/hvg_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    zoom_empty_gs_f     = zoom_dir + '/{zoom_name}/edger_empty_genes_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    pb_empty_f          = pb_dir + '/pb_empties_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
    **get_zoom_conditional_outputs(SPECIES)
  output:
    rmd_f       = rmd_dir + '/' + SHORT_TAG + '_zoom' + '_{zoom_name}_{mkr_sel_res}.Rmd',
    html_f      = docs_dir + '/' + SHORT_TAG + '_zoom' + '_{zoom_name}_{mkr_sel_res}.html'
  params:
    meta_vars   = ','.join(METADATA_VARS), 
    fgsea_args  = lambda wildcards, input: ", ".join([
        f"fgsea_go_bp_f = '{input.get('fgsea_go_bp_f', '')}'",
        f"fgsea_go_cc_f = '{input.get('fgsea_go_cc_f', '')}'",
        f"fgsea_go_mf_f = '{input.get('fgsea_go_mf_f', '')}'",
        f"fgsea_paths_f = '{input.get('fgsea_paths_f', '')}'",
        f"fgsea_hlmk_f = '{input.get('fgsea_hlmk_f', '')}',"
    ]), 
    zoom_mkr_sel_res     = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['MKR_SEL_RES'],
    zoom_mkr_min_cpm_mkr = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['MKR_MIN_CPM_MKR'], 
    zoom_mkr_min_cells   = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['MKR_MIN_CELLS'],  
    zoom_mkr_not_ok_re   = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['MKR_NOT_OK_RE'], 
    zoom_mkr_gsea_cut    = lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['MKR_GSEA_CUT'], 
    zoom_custom_mkr_names= lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['CUSTOM_MKR_NAMES'], 
    zoom_custom_mkr_paths= lambda wildcards: ZOOM_PARAMS_DICT[wildcards.zoom_name]['CUSTOM_MKR_PATHS']
  threads: 1
  retries: RETRIES 
  conda: 
    '../envs/rlibs.yaml'
  resources:
    mem_mb =  lambda wildcards, attempt: attempt * 8192
  shell:
    """

    template_f=$(realpath resources/rmd_templates/zoom.Rmd.template)
    rule="zoom"

    Rscript --vanilla -e "source('scripts/render_reports.R'); \
    render_reports(
    rule_name         = '$rule', 
    proj_dir          = '{PROJ_DIR}', 
    temp_f            =  '$template_f', 
    rmd_f             = '{output.rmd_f}', 
    YOUR_NAME         = '{YOUR_NAME}', 
    AFFILIATION       = '{AFFILIATION}', 
    PROJ_DIR          = '{PROJ_DIR}', 
    SHORT_TAG         = '{SHORT_TAG}', 
    DATE_STAMP        = '{DATE_STAMP}', 
    threads           = {threads},
    zoom_dir          = '{zoom_dir}', 
    zoom_name         = '{wildcards.zoom_name}', 
    meta_f            = '{METADATA_F}', 
    meta_vars_ls      = '{params.meta_vars}', 
    gtf_dt_f          = '{AF_GTF_DT_F}', 
    int_f             = '{input.zoom_int_f}', 
    pb_f              = '{input.zoom_pb_f}', 
    mkrs_f            = '{input.zoom_mkrs_f}', 
    mkrs_hvgs_f       = '{input.zoom_mkrs_hvgs_f}',
    hvgs_f            = '{input.zoom_hvgs_f}',
    empty_gs_f        = '{input.zoom_empty_gs_f}', 
    pb_empty_f        = '{input.pb_empty_f}', 
    {params.fgsea_args}
    CUSTOM_MKR_NAMES  = '{params.zoom_custom_mkr_names}',
    CUSTOM_MKR_PATHS  = '{params.zoom_custom_mkr_paths}',
    MKR_SEL_RES       = '{params.zoom_mkr_sel_res}', 
    MKR_NOT_OK_RE     = '{params.zoom_mkr_not_ok_re}', 
    MKR_MIN_CPM_MKR   = {params.zoom_mkr_min_cpm_mkr}, 
    MKR_MIN_CELLS     = {params.zoom_mkr_min_cells}, 
    MKR_GSEA_CUT      = {params.zoom_mkr_gsea_cut}, 
    SPECIES           = '{SPECIES}')"

    """

