# do labelling
rule zoom_one_zoom:
  input:
    sce_all_f   = sce_dir + '/sce_cells_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    dbl_f       = dbl_dir + '/scDblFinder_combined_outputs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    hmny_f      = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    zoom_sce_sub_f      = zoom_dir + '/{zoom_name}/' + 'zoom_sce_clean_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.rds',
    zoom_hmny_f         = zoom_dir + '/{zoom_name}/' + 'zoom_integrated_dt_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.txt.gz',
    zoom_pb_f           = zoom_dir + '/{zoom_name}/' + 'zoom_pb_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.rds',
    zoom_mkrs_f         = zoom_dir + '/{zoom_name}/' + 'zoom_pb_marker_genes_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.txt.gz',
    zoom_hvgs_f         = zoom_dir + '/{zoom_name}/' + 'zoom_pb_hvgs_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.txt.gz',    
    zoom_fgsea_go_bp_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_go_bp_' + DATE_STAMP +'.txt.gz',
    zoom_fgsea_go_cc_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_go_cc_' + DATE_STAMP +'.txt.gz',
    zoom_fgsea_go_mf_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_go_mf_' + DATE_STAMP +'.txt.gz',
    zoom_fgsea_paths_f  = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_paths_' + DATE_STAMP +'.txt.gz',
    zoom_fgsea_hlmk_f   = zoom_dir + '/{zoom_name}/' + 'zoom_fgsea_' + FULL_TAG + '_{zoom_name}_{zoom_res}_hlmk_' + DATE_STAMP +'.txt.gz',
    zoom_imputed_f      = zoom_dir + '/{zoom_name}/' + 'zoom_imputed_dt_' + FULL_TAG + '_{zoom_name}_{zoom_res}_' + DATE_STAMP +'.txt.gz'
  params:
    zoom_sel_cls      = lambda wildcards: ' '.join(ZOOM_SPEC_LS[wildcards.zoom_name]['sel_cls']),
    zoom_res          = lambda wildcards: ZOOM_SPEC_LS[wildcards.zoom_name]['zoom_res'],
    zoom_n_dims       = lambda wildcards: ZOOM_SPEC_LS[wildcards.zoom_name]['n_dims'],
    zoom_n_hvgs       = lambda wildcards: ZOOM_SPEC_LS[wildcards.zoom_name]['n_hvgs'],
    zoom_min_n_sample = lambda wildcards: ZOOM_SPEC_LS[wildcards.zoom_name]['min_n_sample'],
    zoom_min_n_cl     = lambda wildcards: ZOOM_SPEC_LS[wildcards.zoom_name]['min_n_cl'],
    zoom_n_train      = lambda wildcards: ZOOM_SPEC_LS[wildcards.zoom_name]['n_train']
  threads: 4
  retries: RETRIES
  conda:
    '../envs/rlibs.yml'
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * MB_ZOOM_RUN_ZOOM
  shell:
    """
  
    ZOOM_SEL_CLS={params.zoom_sel_cls} 
    ZOOM_RES={params.zoom_res} 
    ZOOM_N_DIMS={params.zoom_n_dims} 
    ZOOM_N_HVGS={params.zoom_n_hvgs}
    ZOOM_MIN_N_SAMPLE={params.zoom_min_n_sample}
    ZOOM_MIN_N_CL={params.zoom_min_n_cl} 
    ZOOM_N_TRAIN={params.zoom_n_train} 

  
    Rscript -e "\
    source('scripts/utils.R'); source('scripts/integration.R'); \
    source('scripts/marker_genes.R'); source('scripts/zoom.R'); \
    zoom_integrate_within_group(\
      '{FULL_TAG}', '{DATE_STAMP}', '{zoom_dir}', \
      '{input.hmny_f}', '{input.sce_all_f}', '{input.dbl_f}', \
      '{SPECIES}', '{AF_GTF_DT_F}', \
      {INT_SEL_RES}, '{INT_EXC_REGEX}', {INT_DBL_RES}, {INT_DBL_CL_PROP}, {INT_THETA}, \
      '{MKR_NOT_OK_RE}', '{MKR_GSEA_DIR}', {MKR_MIN_CPM_GO}, {MKR_MAX_ZERO_P}, {MKR_GSEA_CUT}, {MKR_MIN_CELLS}, \
      '{wildcards.zoom_name}', '$ZOOM_SEL_CLS', $ZOOM_RES, $ZOOM_N_HVGS, $ZOOM_N_DIMS, \
      $ZOOM_MIN_N_SAMPLE, $ZOOM_MIN_N_CL, $ZOOM_N_TRAIN, n_cores = {threads})"

    """
