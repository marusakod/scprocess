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
