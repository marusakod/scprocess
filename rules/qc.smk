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

# get output file paths as string
def get_qc_files_str(run, SAMPLE_MAPPING, qc_dir, FULL_TAG, DATE_STAMP):
  if SAMPLE_MAPPING is None:
    sce_str = f"{qc_dir}/qc_{run}/sce_cells_clean_{run}_{FULL_TAG}_{DATE_STAMP}.rds"
    qc_str  = f"{qc_dir}/qc_{run}/qc_dt_{run}_{FULL_TAG}_{DATE_STAMP}.txt.gz"
  else:
    sce_fs_ls = []
    qc_fs_ls  = []
    for s in SAMPLE_MAPPING[run]:
      sce_fs_ls.append(f"{qc_dir}/qc_{s}/sce_cells_clean_{s}_{FULL_TAG}_{DATE_STAMP}.rds")
      qc_fs_ls.append(f"{qc_dir}/qc_{s}/qc_dt_{s}_{FULL_TAG}_{DATE_STAMP}.txt.gz")

    sce_str = ','.join(sce_fs_ls)
    qc_str =  ','.join(qc_fs_ls)
  
  return sce_str, qc_str
    


rule run_qc:
  input:
    smpl_stats_f = amb_dir + '/ambient_sample_statistics_' + DATE_STAMP + '.txt',
    amb_yaml_f   = amb_dir + '/ambient_{run}/ambient_{run}_' + DATE_STAMP + '_output_paths.yaml',
    demux_f      = (demux_dir + '/demux_{run}/sce_cells_htos_{run}_' + FULL_TAG + '_' + DATE_STAMP + '.rds') if DEMUX_TYPE == 'af' else ([DEMUX_F] if DEMUX_TYPE == 'custom' else [])
  output:
    sce_f = qc_dir  + '/qc_{sample}/sce_cells_clean_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.rds', sample =lambda wildcards: SAMPLE_MAPPING.get(wildcards.run, [wildcards.run] if SAMPLE_MAPPING else [wildcards.run]),
    qc_f  = qc_dir  + '/qc_{sample}/qc_dt_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', sample =lambda wildcards: SAMPLE_MAPPING.get(wildcards.run, [wildcards.run] if SAMPLE_MAPPING else [wildcards.run]), 
    dbl_f = dbl_dir + '/dbl_{run}/scDblFinder_{run}_outputs_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  params:
    sce_fs_str = lambda wildcards: get_qc_files_str(wildcards.run, SAMPLE_MAPPING, qc_dir, FULL_TAG, DATE_STAMP)[0],
    qc_fs_str  = lambda wildcards: get_qc_files_str(wildcards.run, SAMPLE_MAPPING, qc_dir, FULL_TAG, DATE_STAMP)[1]
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
          dbl_f          = '{output.dbl_f}', \
          gtf_dt_f       = '{AF_GTF_DT_F}', \
          ambient_method = '{AMBIENT_METHOD}', \
          sce_f          = '{params.sce_fs_str}', \
          qc_f           = '{params.qc_fs_str}', \
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

