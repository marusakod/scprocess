# load modules
import yaml
import pandas as pd
import os
import re
import glob
from snakemake.utils import validate, min_version

# do labelling
rule lbl_label_celltypes:
  input:
    sce_clean_f = int_dir + '/sce_clean_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    harmony_f   = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    hvg_mat_f   = lbl_dir + '/hvg_mat_for_labelling_' + LBL_GENE_VAR + '_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
    guesses_f   = lbl_dir + '/xgboost_guesses_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  threads: 4
  resources:
    mem_mb      = 8192
  conda: 
    '../envs/rlibs.yml'
  shell:
    """
    # save sce object
    Rscript -e "source('scripts/label_celltypes.R'); \
    label_celltypes_with_xgboost(xgb_f = '{LBL_XGB_F}', 
      sce_f = '{input.sce_clean_f}', harmony_f = '{input.harmony_f}', \
      hvg_mat_f = '{output.hvg_mat_f}', guesses_f = '{output.guesses_f}', \
      gene_var = '{LBL_GENE_VAR}', min_pred = {LBL_MIN_PRED}, min_cl_prop = {LBL_MIN_CL_PROP}, \
      min_cl_size = {LBL_MIN_CL_SIZE}, n_cores = {threads})"
    """


# save subsets
rule lbl_save_subset_sces:
  input:
    sce_clean_f = f'{int_dir}/sce_clean_{FULL_TAG}_{DATE_STAMP}.rds',
    guesses_f   = f'{lbl_dir}/xgboost_guesses_{FULL_TAG}_{DATE_STAMP}.txt.gz'
  output:
    subsets_df  = f'{lbl_dir}/sce_subset_specifications_{FULL_TAG}_{DATE_STAMP}.csv',
    sce_ls      = expand( [lbl_dir +'/sce_subset_' + FULL_TAG + '_{s}_' + DATE_STAMP + '.rds'], 
      s = None if LBL_SCE_SUBSETS is None else [*LBL_SCE_SUBSETS] )
  threads: 4
  params:
    sub_names = ' '.join([*LBL_SCE_SUBSETS])
  conda:
    '../envs/rlibs.yml'
  resources:
    mem_mb      = 8192
  shell:
    """
    # make dataframe with subset specifications
    python3 -c "
    LBL_SCE_SUBSETS = {LBL_SCE_SUBSETS}
    df = pd.concat([ pd.DataFrame({{'subset_name': k, 'guess': v}}) for k, v in LBL_SCE_SUBSETS.items() ])
    df.to_csv({output.subsets_df}, index = False)
    "
      # save sce object
      Rscript -e "source('scripts/label_celltypes.R'); \
      save_subset_sces(sce_f = '{input.sce_clean_f}', \
        guesses_f = '{input.guesses_f}', sel_res_cl = '{LBL_SEL_RES_CL}', \
        subset_df_f = '{output.subsets_df}', subset_names_concat = '{params.sub_names}', \
        sce_ls_concat = '{output.sce_ls}', allowed_cls_f = '{LBL_XGB_CLS_F}', \
        n_cores = {threads})"
    
    """
