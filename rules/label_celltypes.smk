# load modules
import yaml
import pandas as pd
import os
import re
import glob
from snakemake.utils import validate, min_version

# # do labelling
# rule get_xgboost_labels:
#   input:
#     sces_yaml_f       = int_dir  + '/sce_clean_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml',
#     integration_f     = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz', 
#     qc_sample_stats_f = qc_dir + '/qc_sample_statistics_' + FULL_TAG + '_' + DATE_STAMP + '.csv'
#   output:
#     hvg_mat_f         = lbl_dir + '/hvg_mat_for_labelling_' + 'gene_id' + '_' + FULL_TAG + '_' + DATE_STAMP + '.rds',
#     guesses_f         = lbl_dir + '/cell_annotations_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
#   params:
#     lbl_xgb_f         = [ config['label_celltypes']['lbl_xgb_f'] if 'label_celltypes' in config else [] ],
#     lbl_xgb_cls_f     = [ config['label_celltypes']['lbl_xgb_cls_f'] if 'label_celltypes' in config else [] ],
#     exclude_mito      = [ config['qc']['exclude_mito'] if 'label_celltypes' in config else [] ],
#     mkr_sel_res       = [ config['marker_genes']['mkr_sel_res'] if 'label_celltypes' in config else [] ],
#     lbl_gene_var      = [ config['label_celltypes']['lbl_gene_var'] if 'label_celltypes' in config else [] ],
#     lbl_min_pred      = [ config['label_celltypes']['lbl_min_pred'] if 'label_celltypes' in config else [] ],
#     lbl_min_cl_prop   = [ config['label_celltypes']['lbl_min_cl_prop'] if 'label_celltypes' in config else [] ],
#     lbl_min_cl_size   = [ config['label_celltypes']['lbl_min_cl_size'] if 'label_celltypes' in config else [] ]
#   threads: 4
#   retries: config['resources']['retries']
#   resources:
#     mem_mb      = lambda wildcards, attempt: attempt * config['resources']['gb_label_celltypes'] * MB_PER_GB
#   benchmark:
#     benchmark_dir + '/' + SHORT_TAG + '_label_celltypes/get_xgboost_labels_' + DATE_STAMP + '.benchmark.txt'
#   conda: 
#     '../envs/rlibs.yaml'
#   shell: """
#     # save sce object
#     Rscript -e "source('scripts/label_celltypes.R'); source('scripts/integration.R'); \
#     label_celltypes_with_xgboost(
#       xgb_f              = '{params.lbl_xgb_f}', 
#       allow_f            = '{params.lbl_xgb_cls_f}', 
#       sces_yaml_f        = '{input.sces_yaml_f}',
#       integration_f      = '{input.integration_f}',
#       qc_sample_stats_f  = '{input.qc_sample_stats_f}', 
#       hvg_mat_f          = '{output.hvg_mat_f}',
#       guesses_f          = '{output.guesses_f}',
#       exclude_mito       = '{params.exclude_mito}', 
#       sel_res            =  {params.mkr_sel_res}, 
#       gene_var           = 'gene_id',
#       min_pred           =  {params.lbl_min_pred},
#       min_cl_prop        =  {params.lbl_min_cl_prop},
#       min_cl_size        =  {params.lbl_min_cl_size},
#       n_cores            =  {threads}
#     )"
#     """


rule make_tmp_mtx_file:
  input:
    sces_yaml_f   = int_dir + '/sce_clean_paths_' + FULL_TAG + '_' + DATE_STAMP + '.yaml'
  output:
    mtx_f         = temp(lbl_dir + '/tmp_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.mtx'),
    cells_f       = temp(lbl_dir + '/tmp_cells_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'),
    genes_f       = temp(lbl_dir + '/tmp_genes_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz')
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb        = lambda wildcards, attempt: attempt * config['resources']['gb_label_celltypes'] * MB_PER_GB
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    # save sce object
    Rscript -e "source('scripts/label_celltypes.R'); \
    save_sce_to_mtx(
      sces_yaml_f     = '{input.sces_yaml_f}',
      sel_sample      = '{wildcards.sample}',
      mtx_f           = '{output.mtx_f}',
      cells_f         = '{output.cells_f}',
      genes_f         = '{output.genes_f}'
    )"
    """


# do labelling with celltypist
rule run_celltypist:
  input:
    mtx_f         = lbl_dir + '/tmp_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.mtx',
    cells_f       = lbl_dir + '/tmp_cells_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz',
    genes_f       = lbl_dir + '/tmp_genes_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  output:
    pred_f        = temp(lbl_dir + '/tmp_celltypist_{model}_' + FULL_TAG + '_' + DATE_STAMP + '_{sample}.csv.gz')
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * config['resources']['gb_label_celltypes'] * MB_PER_GB
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_label_celltypes/label_celltypist_{model}_' + DATE_STAMP + '_{sample}.benchmark.txt'
  conda: 
    '../envs/celltypist.yaml'
  shell:"""
    python3 scripts/label_celltypes.py celltypist_one_sample \
      {wildcards.sample} {wildcards.model} \
      --mtx_f     {input.mtx_f} \
      --cells_f   {input.cells_f} \
      --genes_f   {input.genes_f} \
      --pred_f    {output.pred_f}
    """


# do labelling with celltypist
rule merge_celltypist:
  input:
    pred_fs       = expand(lbl_dir + '/tmp_celltypist_{model}_' + FULL_TAG + '_' + DATE_STAMP + '_{sample}.csv.gz', 
      sample = SAMPLES, allow_missing = True),
    integration_f = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
  output:
    pred_out_f    = lbl_dir + '/celltypist_{model}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  params:
    pred_fs_ls    = input.pred_fs,
    hi_res_cl     = config['label_celltypes']['lbl_hi_res_cl'],
    min_cl_size   = config['label_celltypes']['lbl_min_cl_size'],
    min_cl_prop   = config['label_celltypes']['lbl_min_cl_prop']
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb        = lambda wildcards, attempt: attempt * config['resources']['gb_label_celltypes'] * MB_PER_GB
  conda: 
    '../envs/celltypist.yaml'
  shell:"""
    python3 scripts/label_celltypes.py aggregate_predictions \
      {params.pred_fs_ls} \
      --int_f           {input.integration_f} \
      --hi_res_cl       {params.hi_res_cl} \
      --min_cl_size     {params.min_cl_size} \
      --min_cl_prop     {params.min_cl_prop} \
      --agg_f           {output.pred_out_f}
    """
