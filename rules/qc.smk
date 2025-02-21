# snakemake rule for doing QC on sce object

rule run_qc:
  input:
    sce_all_f   = sce_dir + '/sce_cells_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
    dbl_f       = dbl_dir + '/scDblFinder_combined_outputs_' + FULL_TAG +'_' + DATE_STAMP + '.txt.gz'
  output:
    qc_f        = qc_dir  + '/qc_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
    keep_f      = qc_dir  + '/keep_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  threads: 8
  retries: RETRIES 
  resources:
    mem_mb      =   lambda wildcards, attempt: attempt * MB_RUN_QC
  conda:
    '../envs/rlibs.yml'
  shell:
    """
    Rscript -e "source('scripts/SampleQC.R'); \
      main_qc(sce_f = '{input.sce_all_f}', dbl_f = '{input.dbl_f}', dbl_smpl_var = '{SAMPLE_VAR}', \
        hard_min_counts = {QC_HARD_MIN_COUNTS}, hard_min_feats = {QC_HARD_MIN_FEATS}, \
        hard_max_mito = {QC_HARD_MAX_MITO}, min_counts = {QC_MIN_COUNTS}, 
        min_feats = {QC_MIN_FEATS}, min_mito = {QC_MIN_MITO}, max_mito = {QC_MAX_MITO}, \
        min_splice = {QC_MIN_SPLICE}, max_splice = {QC_MAX_SPLICE}, \
        min_cells = {QC_MIN_CELLS}, filter_bender = '{QC_FILTER_BENDER}', \
        amb_method = '{AMBIENT_METHOD}', qc_f = '{output.qc_f}', keep_f = '{output.keep_f}')"

    """

