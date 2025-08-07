
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
    zoom_lbls       = lambda wildcards: ','.join(ZOOM_PARAMS_DICT[wildcards.zoom_name]['LABELS']),
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

# rule zoom_one_zoom:
#   input:
#     sce_all_f   = sce_dir + '/sce_cells_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
#     dbl_f       = dbl_dir + '/scDblFinder_combined_outputs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
#     hmny_f      = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
#   output:
#     zoom_sce_sub_f      = zoom_dir + '/{zoom_name}/' + 'zoom_sce_clean_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.rds',
#     zoom_hmny_f         = zoom_dir + '/{zoom_name}/' + 'zoom_integrated_dt_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.txt.gz',
#     zoom_pb_f           = zoom_dir + '/{zoom_name}/' + 'zoom_pb_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.rds',
#     zoom_mkrs_f         = zoom_dir + '/{zoom_name}/' + 'zoom_pb_marker_genes_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.txt.gz',
#     zoom_hvgs_f         = zoom_dir + '/{zoom_name}/' + 'zoom_pb_hvgs_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.txt.gz',    
#     zoom_fgsea_go_bp_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_go_bp_' + DATE_STAMP +'.txt.gz',
#     zoom_fgsea_go_cc_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_go_cc_' + DATE_STAMP +'.txt.gz',
#     zoom_fgsea_go_mf_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_go_mf_' + DATE_STAMP +'.txt.gz',
#     zoom_fgsea_paths_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_paths_' + DATE_STAMP +'.txt.gz',
#     zoom_fgsea_hlmk_f   = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_hlmk_' + DATE_STAMP +'.txt.gz',
#     zoom_imputed_f      = zoom_dir + '/{zoom_name}/' + 'zoom_imputed_dt_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.txt.gz'
#   params:
#     zoom_sel_cls      = lambda wildcards: ' '.join(ZOOM_SPEC_LS[wildcards.zoom_name]['sel_cls']),
#     zoom_res          = lambda wildcards: ZOOM_SPEC_LS[wildcards.zoom_name]['zoom_res'],
#     zoom_n_hvgs       = lambda wildcards: ZOOM_SPEC_LS[wildcards.zoom_name]['n_hvgs'],
#     zoom_n_dims       = lambda wildcards: ZOOM_SPEC_LS[wildcards.zoom_name]['n_dims'],
#     zoom_min_n_sample = lambda wildcards: ZOOM_SPEC_LS[wildcards.zoom_name]['min_n_sample'],
#     zoom_min_n_cl     = lambda wildcards: ZOOM_SPEC_LS[wildcards.zoom_name]['min_n_cl'],
#     zoom_n_train      = lambda wildcards: ZOOM_SPEC_LS[wildcards.zoom_name]['n_train']
#   threads: 4
#   retries: RETRIES
#   conda:
#     '../envs/rlibs.yaml'
#   resources:
#     mem_mb      = lambda wildcards, attempt: attempt * MB_ZOOM_RUN_ZOOM
#   shell:
#     """
#     Rscript -e "\
#     source('scripts/utils.R'); source('scripts/integration.R'); \
#     source('scripts/marker_genes.R'); source('scripts/zoom.R'); \
#     zoom_integrate_within_group(\
#       '{FULL_TAG}', '{DATE_STAMP}', '{zoom_dir}', \
#       '{input.hmny_f}', '{input.sce_all_f}', '{input.dbl_f}', \
#       '{SPECIES}', '{AF_GTF_DT_F}', \
#       {MKR_SEL_RES}, '{INT_EXC_REGEX}', {INT_DBL_RES}, {INT_DBL_CL_PROP}, {INT_THETA}, \
#       '{MKR_NOT_OK_RE}', '{MKR_GSEA_DIR}', {MKR_MIN_CPM_GO}, {MKR_MAX_ZERO_P}, \
#       {MKR_GSEA_CUT}, {MKR_MIN_CELLS}, '{wildcards.zoom_name}', '{params.zoom_sel_cls}', \
#       '{params.zoom_res}', {params.zoom_n_hvgs}, {params.zoom_n_dims}, {params.zoom_min_n_sample}, 
#       {params.zoom_min_n_cl}, {params.zoom_n_train}, n_cores = {threads})"

#     """
