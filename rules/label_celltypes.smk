# load modules
import yaml
import pandas as pd
import os
import re
import glob
from snakemake.utils import validate, min_version

# prep counts matrix for each sample
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
    pred_f        = temp(lbl_dir + '/tmp_labels_celltypist_model_{model}_' + FULL_TAG + '_' + DATE_STAMP + '_{sample}.csv.gz')
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * config['resources']['gb_label_celltypes'] * MB_PER_GB
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_label_celltypes/labels_celltypist_{model}_' + DATE_STAMP + '_{sample}.benchmark.txt'
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


# do labelling
rule run_scprocess_labeller:
  input:
    mtx_f     = lbl_dir + '/tmp_counts_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.mtx',
    cells_f   = lbl_dir + '/tmp_cells_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz',
    genes_f   = lbl_dir + '/tmp_genes_{sample}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
  output:
    pred_f    = temp(lbl_dir + '/tmp_labels_scprocess_model_{model}_' + FULL_TAG + '_' + DATE_STAMP + '_{sample}.csv.gz')
  params:
    xgb_f     = lambda wildcards: [ entry['xgb_f'] for entry in LABELLER_PARAMS 
      if (entry['labeller'] == "scprocess") and (entry['model'] == wildcards.model) ],
    xgb_cls_f = lambda wildcards: [ entry['xgb_cls_f'] for entry in LABELLER_PARAMS 
      if (entry['labeller'] == "scprocess") and (entry['model'] == wildcards.model) ]
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb      = lambda wildcards, attempt: attempt * config['resources']['gb_label_celltypes'] * MB_PER_GB
  benchmark:
    benchmark_dir + '/' + SHORT_TAG + '_label_celltypes/labels_scprocess_{model}_' + DATE_STAMP + '_{sample}.benchmark.txt'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    # save sce object
    Rscript -e "source('scripts/label_celltypes.R'); source('scripts/integration.R'); \
    label_with_xgboost_one_sample(
      sel_sample  = '{wildcards.sample}', 
      model_name  = '{wildcards.model}', 
      xgb_f       = '{params.xgb_f}', 
      xgb_cls_f   = '{params.xgb_cls_f}', 
      mtx_f       = '{input.mtx_f}',
      cells_f     = '{input.cells_f}',
      genes_f     = '{input.genes_f}',
      pred_f      = '{output.pred_f}'
    )"
    """

# only do this bit if other parts are finished
qc_stats_f  = pathlib.Path(f"{qc_dir}/qc_sample_statistics_{FULL_TAG}_{DATE_STAMP}.csv")
if ('label_celltypes' in config) & qc_stats_f.is_file():

  def parse_merge_labels_parameters(LABELLER_PARAMS, labeller, model):
    # get the one that matches
    this_entry  = [ entry for entry in LABELLER_PARAMS 
      if ((entry['labeller'] == labeller) and (entry['model'] == model)) ]

    # check exactly one match
    if len(this_entry) != 1:
      raise ValueError("only one entry should match this")

    return this_entry[0]

  # load up qc outputs
  qc_stats    = pl.read_csv(qc_stats_f).filter( pl.col('bad_sample') == False )

  # do labelling with celltypist
  rule merge_labels:
    input:
      pred_fs       = expand(lbl_dir + '/tmp_labels_{labeller}_model_{model}_' + FULL_TAG + '_' + DATE_STAMP + '_{sample}.csv.gz', 
        sample = qc_stats["sample_id"].to_list(), allow_missing = True),
      integration_f = int_dir + '/integrated_dt_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
    output:
      pred_out_f    = lbl_dir + '/labels_{labeller}_model_{model}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
    params:
      pred_fs_ls    = input.pred_fs,
      hi_res_cl     = lambda wildcards: parse_merge_labels_parameters(LABELLER_PARAMS, wildcards.labeller, wildcards.model)["hi_res_cl"],
      min_cl_size   = lambda wildcards: parse_merge_labels_parameters(LABELLER_PARAMS, wildcards.labeller, wildcards.model)["min_cl_size"],
      min_cl_prop   = lambda wildcards: parse_merge_labels_parameters(LABELLER_PARAMS, wildcards.labeller, wildcards.model)["min_cl_prop"]
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
