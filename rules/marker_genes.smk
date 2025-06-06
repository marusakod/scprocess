# snakemake rule for calculating marker genes

def get_conditional_outputs(species):
    if species in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020']:
        return {
            'fgsea_go_bp_f': mkr_dir + '/fgsea_' + FULL_TAG + f'_{INT_SEL_RES}_go_bp_' + DATE_STAMP + '.txt.gz',
            'fgsea_go_cc_f': mkr_dir + '/fgsea_' + FULL_TAG + f'_{INT_SEL_RES}_go_cc_' + DATE_STAMP + '.txt.gz',
            'fgsea_go_mf_f': mkr_dir + '/fgsea_' + FULL_TAG + f'_{INT_SEL_RES}_go_mf_' + DATE_STAMP + '.txt.gz',
            'fgsea_paths_f': mkr_dir + '/fgsea_' + FULL_TAG + f'_{INT_SEL_RES}_paths_' + DATE_STAMP + '.txt.gz',
            'fgsea_hlmk_f':  mkr_dir + '/fgsea_' + FULL_TAG + f'_{INT_SEL_RES}_hlmk_' + DATE_STAMP + '.txt.gz'
        }
    else:
        return {}

rule run_marker_genes:
  input:
    sces_yaml_f    = qc_dir  + '/sce_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml',
    integration_f  = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    pb_f      = mkr_dir + '/pb_' + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.rds',
    mkrs_f    = mkr_dir + '/pb_marker_genes_' + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.txt.gz',
    pb_hvgs_f = mkr_dir + '/pb_hvgs_' + FULL_TAG + f'_{INT_SEL_RES}_' + DATE_STAMP + '.txt.gz',
    **get_conditional_outputs(SPECIES)
  params:
    fgsea_args = lambda wildcards, output: ", ".join([
        f"fgsea_go_bp_f = '{output.get('fgsea_go_bp_f', '')}'",
        f"fgsea_go_cc_f = '{output.get('fgsea_go_cc_f', '')}'",
        f"fgsea_go_mf_f = '{output.get('fgsea_go_mf_f', '')}'",
        f"fgsea_paths_f = '{output.get('fgsea_paths_f', '')}'",
        f"fgsea_hlmk_f = '{output.get('fgsea_hlmk_f', '')}',"
    ])
  threads: 8
  retries: RETRIES
  resources:
    mem_mb = lambda wildcards, attempt: attempt * MB_RUN_MARKER_GENES
  conda: '../envs/rlibs.yml'
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
        sel_res       = {INT_SEL_RES},
        min_cl_size   = {MKR_MIN_CL_SIZE},
        min_cells     = {MKR_MIN_CELLS},
        not_ok_re     = '{MKR_NOT_OK_RE}',
        min_cpm_go    = {MKR_MIN_CPM_GO},
        max_zero_p    = {MKR_MAX_ZERO_P},
        gsea_cut      = {MKR_GSEA_CUT},
        n_cores       = {threads})"
    """
