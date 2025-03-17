# snakemake rule for calculating marker genes

rule run_marker_genes:
  input:
    sce_clean_f   = int_dir + '/sce_clean_' + FULL_TAG + '_' + DATE_STAMP + '.rds'
  output:
    pb_f          = mkr_dir + '/pb_'              + FULL_TAG + f'_{MKR_SEL_RES}_' + DATE_STAMP + '.rds',
    mkrs_f        = mkr_dir + '/pb_marker_genes_' + FULL_TAG + f'_{MKR_SEL_RES}_' + DATE_STAMP + '.txt.gz',
    pb_hvgs_f     = mkr_dir + '/pb_hvgs_'         + FULL_TAG + f'_{MKR_SEL_RES}_' + DATE_STAMP + '.txt.gz',
    fgsea_go_bp_f = mkr_dir + '/fgsea_'           + FULL_TAG + f'_{MKR_SEL_RES}_' + 'go_bp_' + DATE_STAMP + '.txt.gz',
    fgsea_go_cc_f = mkr_dir + '/fgsea_'           + FULL_TAG + f'_{MKR_SEL_RES}_' + 'go_cc_' + DATE_STAMP + '.txt.gz',
    fgsea_go_mf_f = mkr_dir + '/fgsea_'           + FULL_TAG + f'_{MKR_SEL_RES}_' + 'go_mf_' + DATE_STAMP + '.txt.gz',
    fgsea_paths_f = mkr_dir + '/fgsea_'           + FULL_TAG + f'_{MKR_SEL_RES}_' + 'paths_' + DATE_STAMP + '.txt.gz',
    fgsea_hlmk_f  = mkr_dir + '/fgsea_'           + FULL_TAG + f'_{MKR_SEL_RES}_' + 'hlmk_' + DATE_STAMP + '.txt.gz'
  threads: 8
  retries: RETRIES
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * MB_RUN_MARKER_GENES
  conda: 
    '../envs/rlibs.yml'
  shell:
    """
    # save sce object
    Rscript -e "\
      source('scripts/utils.R'); \
      source('scripts/marker_genes.R'); \
      calculate_marker_genes(
        sce_clean_f   = '{input.sce_clean_f}', \
        pb_f          = '{output.pb_f}', \
        mkrs_f        = '{output.mkrs_f}', \
        pb_hvgs_f     = '{output.pb_hvgs_f}', \
        fgsea_go_bp_f = '{output.fgsea_go_bp_f}', \
        fgsea_go_cc_f = '{output.fgsea_go_cc_f}', \
        fgsea_go_mf_f = '{output.fgsea_go_mf_f}', \
        fgsea_paths_f = '{output.fgsea_paths_f}', \
        fgsea_hlmk_f  = '{output.fgsea_hlmk_f}', \
        species       = '{SPECIES}', \
        gtf_dt_f      = '{AF_GTF_DT_F}', \
        gsea_dir      = '{MKR_GSEA_DIR}', \
        sel_res       = {MKR_SEL_RES}, \
        exc_regex     = '{INT_EXC_REGEX}', \
        min_cl_size   = {MKR_MIN_CL_SIZE}, \
        min_cells     = {MKR_MIN_CELLS}, \
        not_ok_re     = '{MKR_NOT_OK_RE}', \
        min_cpm_go    = {MKR_MIN_CPM_GO}, \
        max_zero_p    = {MKR_MAX_ZERO_P}, \
        gsea_cut      = {MKR_GSEA_CUT}, \
        n_cores       = {threads})"
    """
