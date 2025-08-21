# load modules
import yaml
import pandas as pd
import os
import re
import glob
from snakemake.utils import validate, min_version

# do labelling
rule get_xgboost_labels:
  input:
    sces_yaml_f        = int_dir  + '/sce_clean_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml',
    integration_f      = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
    qc_sample_stats_f  = qc_dir + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
  output:
    hvg_mat_f   = lbl_dir + '/hvg_mat_for_labelling_' + LBL_GENE_VAR + '_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    guesses_f   = lbl_dir + '/cell_annotations_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  threads: 4
  retries: RETRIES 
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * MB_LABEL_CELLTYPES
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_label_celltypes/get_xgboost_labels_' + DATE_STAMP + '.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell:
    """
    # save sce object
    Rscript -e "source('scripts/label_celltypes.R'); source('scripts/integration.R'); \
    label_celltypes_with_xgboost(
      xgb_f              = '{LBL_XGB_F}', 
      allow_f            = '{LBL_XGB_CLS_F}', 
      sces_yaml_f        = '{input.sces_yaml_f}',
      integration_f      = '{input.integration_f}',
      qc_sample_stats_f  = '{input.qc_sample_stats_f}', 
      hvg_mat_f          = '{output.hvg_mat_f}',
      guesses_f          = '{output.guesses_f}',
      exclude_mito       = '{EXCLUDE_MITO}', 
      sel_res            = {MKR_SEL_RES}, 
      gene_var           = '{LBL_GENE_VAR}',
      min_pred           = {LBL_MIN_PRED},
      min_cl_prop        = {LBL_MIN_CL_PROP},
      min_cl_size        = {LBL_MIN_CL_SIZE},
      n_cores            = {threads})"
    """

# if not LBL_SCE_SUBSETS:
#     localrules: lbl_save_subset_sces  
#     rule lbl_save_subset_sces:
#         message:
#             "No subsets specified. Skipping subset generation."
#         output:
#             subsets_df = f'{lbl_dir}/sce_subset_specifications_{FULL_TAG}_{DATE_STAMP}.csv',
#             sce_ls = []
#         shell:
#           """
#           touch {output.subsets_df}
#           """
# else:
#     rule lbl_save_subset_sces:
#         input:
#             sce_clean_f = f'{int_dir}/sce_clean_{FULL_TAG}_{DATE_STAMP}.rds',
#             guesses_f   = f'{lbl_dir}/cell_annotations_{FULL_TAG}_{DATE_STAMP}.txt.gz'
#         output:
#             subsets_df  = f'{lbl_dir}/sce_subset_specifications_{FULL_TAG}_{DATE_STAMP}.csv',
#             sce_ls      = expand(
#                 f"{lbl_dir}/sce_subset_{FULL_TAG}_{{s}}_{DATE_STAMP}.rds",
#                 s=[*LBL_SCE_SUBSETS]
#             )
#         threads: 4
#         retries: RETRIES
#         params:
#             sub_names = ' '.join([*LBL_SCE_SUBSETS])
#         conda:
#             '../envs/rlibs.yaml'
#         resources:
#             mem_mb = lambda wildcards, attempt: attempt * MB_LBL_SAVE_SUBSET_SCES
#         shell:
#             """
#             # Make dataframe with subset specifications
#             python3 -c "
# import pandas as pd
# LBL_SCE_SUBSETS = {LBL_SCE_SUBSETS}
# df = pd.concat([pd.DataFrame({{'subset_name': k, 'guess': v}}) for k, v in LBL_SCE_SUBSETS.items()])
# df.to_csv('{output.subsets_df}', index=False)
# "

#             # Save SCE objects
#             Rscript -e "source('scripts/label_celltypes.R'); \
#                 save_subset_sces(sce_f = '{input.sce_clean_f}', \
#                                  guesses_f = '{input.guesses_f}', \
#                                  sel_res_cl = '{LBL_SEL_RES_CL}', \
#                                  subset_df_f = '{output.subsets_df}', \
#                                  subset_names_concat = '{params.sub_names}', \
# 				 custom_labels_f = '{CUSTOM_LABELS_F}', \
#                                  sce_ls_concat = '{output.sce_ls}', \
#                                  allowed_cls_f = '{LBL_XGB_CLS_F}', \
#                                  n_cores = {threads})"
#             """

