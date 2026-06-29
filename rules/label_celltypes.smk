# load modules
import yaml
import pandas as pd
import os
import re
import glob
from snakemake.utils import validate, min_version


wildcard_constraints:
  labeller = "celltypist|scprocess|custom"

# do labelling with celltypist
rule run_celltypist:
  input:
    adata_f       = f'{int_dir}/anndata_cells_clean_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5ad'
  output:
    pred_f        = temp(f'{lbl_dir}/tmp_labels_celltypist_model_{{model}}_{FULL_TAG}_{DATE_STAMP}_{{batch}}.csv.gz')
  params:
    batch_var     = BATCH_VAR
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_celltypist', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_celltypist', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/label_celltypes/run_celltypist_{{model}}_{{batch}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/label_celltypes/run_celltypist_{{model}}_{{batch}}_{DATE_STAMP}.log'
  conda: 
    '../envs/celltypist.yaml'
  shell:"""
    exec &>> {log}
    
    python3 scripts/label_celltypes.py celltypist_one_batch \
      {wildcards.batch} {params.batch_var} {wildcards.model} \
      --adata_f   {input.adata_f} \
      --pred_f    {output.pred_f}
    """


# do labelling with xgboost
rule run_scprocess_labeller:
  input:
    adata_f   = f'{int_dir}/anndata_cells_clean_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5ad'
  output:
    pred_f    = temp(f'{lbl_dir}/tmp_labels_scprocess_model_{{model}}_{FULL_TAG}_{DATE_STAMP}_{{batch}}.csv.gz')
  params:
    xgb_f     = lambda wildcards: [ entry['xgb_f'] for entry in LABELLER_PARAMS 
      if (entry['labeller'] == "scprocess") and (entry['model'] == wildcards.model) ],
    xgb_cls_f = lambda wildcards: [ entry['xgb_cls_f'] for entry in LABELLER_PARAMS 
      if (entry['labeller'] == "scprocess") and (entry['model'] == wildcards.model) ], 
    batch_var = BATCH_VAR
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_scprocess_labeller', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'run_scprocess_labeller', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/label_celltypes/run_scprocess_labeller_{{model}}_{{batch}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/label_celltypes/run_scprocess_labeller_{{model}}_{{batch}}_{DATE_STAMP}.log'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    exec &>> {log}

    # save sce object
    Rscript -e "source('scripts/label_celltypes.R'); source('scripts/integration.R'); \
    label_with_xgboost_one_batch(
      sel_batch   = '{wildcards.batch}', 
      batch_var   = '{params.batch_var}',
      model_name  = '{wildcards.model}', 
      xgb_f       = '{params.xgb_f}', 
      xgb_cls_f   = '{params.xgb_cls_f}', 
      adata_f     = '{input.adata_f}',
      pred_f      = '{output.pred_f}'
    )"
    """

# --- Merge per-batch predictions into final labels ---

def _parse_merge_labels_parameters(LABELLER_PARAMS, labeller, model):
  this_entry  = [ entry for entry in LABELLER_PARAMS
    if ((entry['labeller'] == labeller) and (entry['model'] == model)) ]
  if len(this_entry) != 1:
    raise ValueError("only one entry should match this")
  return this_entry[0]


def _get_good_batch_labels(wildcards):
  """Defer batch list resolution until the QC checkpoint completes."""
  qc_stats_f = checkpoints.get_qc_sample_statistics.get().output.qc_stats_f
  bad_var = 'bad_' + BATCH_VAR
  good_batches = pl.read_csv(qc_stats_f).filter(pl.col(bad_var) == False)[BATCH_VAR].to_list()
  return expand(
    f'{lbl_dir}/tmp_labels_{wildcards.labeller}_model_{wildcards.model}_{FULL_TAG}_{DATE_STAMP}_{{batch}}.csv.gz',
    batch = good_batches
  )


rule merge_labels:
  input:
    pred_fs       = _get_good_batch_labels,
    integration_f = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    pred_out_f    = f'{lbl_dir}/labels_{{labeller}}_model_{{model}}_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  params:
    pred_fs_ls    = lambda wildcards, input: input.pred_fs,
    hi_res_cl     = lambda wildcards: _parse_merge_labels_parameters(LABELLER_PARAMS, wildcards.labeller, wildcards.model)["hi_res_cl"],
    min_cl_prop   = lambda wildcards: _parse_merge_labels_parameters(LABELLER_PARAMS, wildcards.labeller, wildcards.model)["min_cl_prop"],
    batch_var     = BATCH_VAR
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'merge_labels', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'merge_labels', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/label_celltypes/merge_labels_{{labeller}}_{{model}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/label_celltypes/merge_labels_{{labeller}}_{{model}}_{DATE_STAMP}.log'
  conda:
    '../envs/celltypist.yaml'
  shell:"""
    exec &>> {log}

    python3 scripts/label_celltypes.py aggregate_predictions \
      {params.pred_fs_ls} \
      --int_f           {input.integration_f} \
      --hi_res_cl       {params.hi_res_cl} \
      --min_cl_prop     {params.min_cl_prop} \
      --batch_var       {params.batch_var} \
      --agg_f           {output.pred_out_f}
    """


_names_entries = [e for e in LABELLER_PARAMS if e.get('save_cluster_names_file', False)]
if _names_entries:
  _names_mkr_sel_res = config['marker_genes']['mkr_sel_res']

  rule save_cluster_names:
    input:
      labels_f      = f'{lbl_dir}/labels_{{labeller}}_model_{{model}}_{FULL_TAG}_{DATE_STAMP}.csv.gz',
      integration_f = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
    output:
      names_f       = f'{lbl_dir}/cluster_names_for_shiny_{{labeller}}_{{model}}_{FULL_TAG}_{_names_mkr_sel_res}_{DATE_STAMP}.csv'
    params:
      mkr_sel_res   = _names_mkr_sel_res
    retries: config['resources']['retries']
    log:
      f'{logs_dir}/label_celltypes/save_cluster_names_{{labeller}}_{{model}}_{DATE_STAMP}.log'
    conda:
      '../envs/scprocess_local.yaml'
    shell: """
      exec &>> {log}
      python3 scripts/label_celltypes.py save_cluster_names \
        --labels_f      {input.labels_f} \
        --integration_f {input.integration_f} \
        --mkr_sel_res   {params.mkr_sel_res} \
        --output_f      {output.names_f}
      """

