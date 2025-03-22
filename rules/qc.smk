# snakemake rule for doing QC on sce object

# rule run_qc:
#   input:
#     sce_all_f   = sce_dir + '/sce_cells_all_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
#     dbl_f       = dbl_dir + '/scDblFinder_combined_outputs_' + FULL_TAG +'_' + DATE_STAMP + '.txt.gz'
#   output:
#     qc_f        = qc_dir  + '/qc_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz',
#     keep_f      = qc_dir  + '/keep_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
#   threads: 8
#   retries: RETRIES 
#   resources:
#     mem_mb      =   lambda wildcards, attempt: attempt * MB_RUN_QC
#   conda:
#     '../envs/rlibs.yml'
#   shell:
#     """
#     Rscript -e "source('scripts/SampleQC.R'); \
#       main_qc(sce_f = '{input.sce_all_f}', dbl_f = '{input.dbl_f}', dbl_smpl_var = '{SAMPLE_VAR}', \
#         hard_min_counts = {QC_HARD_MIN_COUNTS}, hard_min_feats = {QC_HARD_MIN_FEATS}, \
#         hard_max_mito = {QC_HARD_MAX_MITO}, min_counts = {QC_MIN_COUNTS}, 
#         min_feats = {QC_MIN_FEATS}, min_mito = {QC_MIN_MITO}, max_mito = {QC_MAX_MITO}, \
#         min_splice = {QC_MIN_SPLICE}, max_splice = {QC_MAX_SPLICE}, \
#         min_cells = {QC_MIN_CELLS}, filter_bender = '{QC_FILTER_BENDER}', \
#         amb_method = '{AMBIENT_METHOD}', qc_f = '{output.qc_f}', keep_f = '{output.keep_f}')"

#     """



# one rule to get all qc metrics, calc doublets, create filtered sc object, save dt with QC metrics for ALL cells
# should save files per sample_id not per pool
# should save files only for samples that user doesn't want to exclude
# should save empty files for samples that don't have enough cells
# should save a csv file with sample_id and n_cells, keep columns (to keep track of all samples that were removed because not enough cells) 


# rule run_qc:
#   input:
#     smpl_stats_f = amb_dir + '/ambient_sample_statistics_' + DATE_STAMP + '.txt',
#     amb_yaml_f   = expand(amb_dir + '/ambient_{sample}/ambient_{sample}_' + DATE_STAMP + '_output_paths.yaml', sample=runs),
#     demux_f      = expand(sce_dir + '/sce_cells_htos_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.rds', sample=runs) 
#       if DEMUX_TYPE == 'af' else ([DEMUX_F] if DEMUX_TYPE == 'custom' else [])
#   output:
#     sce_f = expand(sce_dir + '/qc_{sample}/sce_cells_clean_{sample}_' + DATE_STAMP + '.rds', sample=SAMPLES if DEMUX_TYPE != '' else runs),
#     qc_f  = expand(sce_dir + '/qc_{sample}/qc_dt_{sample}_' + SHORT_TAG + '_' + DATE_STAMP + '.txt.gz', sample=SAMPLES if DEMUX_TYPE != '' else runs),
#     dimred_f = expand(dbl_dir + '/dbl_{sample}/scDblFinder_{sample}_dimreds_' + SHORT_TAG + '_' + DATE_STAMP + '.txt.gz', sample=runs),
#     dbl_f    = expand(dbl_dir + '/dbl_{sample}/scDblFinder_{sample}_outputs_' + SHORT_TAG + '_' + DATE_STAMP + '.txt.gz', sample=runs)
#   threads: 4
#   retries: RETRIES
#   resources:
#     mem_mb = lambda wildcards, attempt: attempt * MB_RUN_QC
#   conda:
#     '../envs/rlibs.yml'
#   shell:
#     """
#     # save sce object
#     Rscript -e "source('scripts/SampleQC.R'); \
#         main_qc( \
#           sel_sample     = '{wildcards.sample}', \
#           meta_f         = '{METADATA_F}', \
#           amb_yaml_f     = '{input.amb_yaml_f}', \
#           sample_stats_f = '{output.smpl_stats_f}', \
#           demux_f        = '{input.demux_f}', \
#           gtf_dt_f       = '{AF_GTF_DT_F}', \
#           ambient_method = '{AMBIENT_METHOD}', \
#           sce_f          = '{output.sce_f}', \
#           rowdata_f      = '{output.rowdata_f}', \
#           dbl_dimred_f   = '{output.dimred_f}', \
#           qc_f           = '{output.qc_f}', \
#           hard_min_counts= {QC_HARD_MIN_COUNTS}, \
#           hard_min_feats = {QC_HARD_MIN_FEATS}, \
#           hard_max_mito  = {QC_HARD_MAX_MITO}, \
#           min_counts     = {QC_MIN_COUNTS}, \
#           min_feats      = {QC_MIN_FEATS}, \
#           min_mito       = {QC_MIN_MITO}, \
#           max_mito       = {QC_MAX_MITO}, \
#           min_splice     = {QC_MIN_SPLICE}, \
#           max_splice     = {QC_MAX_SPLICE}, \
#           min_cells      = {QC_MIN_CELLS}, \
#           sample_var     = '{SAMPLE_VAR}', \
#           demux_type     = '{DEMUX_TYPE}', \
#           dbl_min_feats  = {DBL_MIN_FEATS})"
#     """






def get_samples_per_run():
    metadata = pd.read_csv(METADATA_F)
    sample_mapping = metadata.groupby("pool_id")["sample_id"].apply(list).to_dict()
    return sample_mapping

SAMPLE_MAPPING = get_samples_per_run()


rule run_qc:
  input:
    smpl_stats_f = amb_dir + '/ambient_sample_statistics_' + DATE_STAMP + '.txt',
    amb_yaml_f   = lambda wildcards: amb_dir + f'/ambient_{wildcards.run}/ambient_{wildcards.run}_' + DATE_STAMP + '_output_paths.yaml',
    demux_f      = lambda wildcards: expand(
        sce_dir + '/hto_{sample}/sce_cells_htos_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.rds', 
        sample=SAMPLE_MAPPING[wildcards.run]) 
      if DEMUX_TYPE == 'af' else ([DEMUX_F] if DEMUX_TYPE == 'custom' else [])
  output:
    sce_f = lambda wildcards: expand(
        sce_dir + '/qc_{sample}/sce_cells_clean_{sample}_' + DATE_STAMP + '.rds', 
        sample=SAMPLE_MAPPING[wildcards.run] if DEMUX_TYPE != '' else [wildcards.run]),
    qc_f  = lambda wildcards: expand(
        sce_dir + '/qc_{sample}/qc_dt_{sample}_' + SHORT_TAG + '_' + DATE_STAMP + '.txt.gz', 
        sample=SAMPLE_MAPPING[wildcards.run] if DEMUX_TYPE != '' else [wildcards.run]),
  params:
    samples = lambda wildcards: SAMPLE_MAPPING[wildcards.run] if DEMUX_TYPE != '' else [wildcards.run],
    sce_f_str = lambda wildcards: ",".join(
        expand(
            sce_dir + '/qc_{sample}/sce_cells_clean_{sample}_' + DATE_STAMP + '.rds',
            sample=SAMPLE_MAPPING[wildcards.run] if DEMUX_TYPE != '' else [wildcards.run]
        )
    ),
    qc_f_str = lambda wildcards: ",".join(
        expand(
            sce_dir + '/qc_{sample}/qc_dt_{sample}_' + SHORT_TAG + '_' + DATE_STAMP + '.txt.gz',
            sample=SAMPLE_MAPPING[wildcards.run] if DEMUX_TYPE != '' else [wildcards.run]
        )
    )
  threads: 4
  retries: RETRIES
  resources:
    mem_mb = lambda wildcards, attempt: attempt * MB_RUN_QC
  conda:
    '../envs/rlibs.yml'
  shell:
    """
    Rscript -e "source('scripts/SampleQC.R'); \
        test_qc( \
          sel_sample     = '{wildcards.run}', \
          meta_f         = '{METADATA_F}', \
          amb_yaml_f     = '{input.amb_yaml_f}', \
          sample_stats_f = '{input.smpl_stats_f}', \
          demux_f        = '{input.demux_f}', \
          gtf_dt_f       = '{AF_GTF_DT_F}', \
          ambient_method = '{AMBIENT_METHOD}', \
          sce_f          = '{params.sce_f_str}', \
          qc_f           = '{params.qc_f_str}', \
          hard_min_counts= {QC_HARD_MIN_COUNTS}, \
          hard_min_feats = {QC_HARD_MIN_FEATS}, \
          hard_max_mito  = {QC_HARD_MAX_MITO}, \
          min_counts     = {QC_MIN_COUNTS}, \
          min_feats      = {QC_MIN_FEATS}, \
          min_mito       = {QC_MIN_MITO}, \
          max_mito       = {QC_MAX_MITO}, \
          min_splice     = {QC_MIN_SPLICE}, \
          max_splice     = {QC_MAX_SPLICE}, \
          min_cells      = {QC_MIN_CELLS}, \
          sample_var     = '{SAMPLE_VAR}', \
          demux_type     = '{DEMUX_TYPE}', \
          dbl_min_feats  = {DBL_MIN_FEATS})"
    """

# rule that takes all temporary qc files and combines them into a single files (temporary files should be removed)

# rule that takes the entire qc file, determines which samples should be excluded and writes a checkpoint csv file.

