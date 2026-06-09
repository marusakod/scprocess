# load modules
import yaml
import pandas as pd
import os
import re
import glob
from snakemake.utils import validate, min_version


def _get_labeller_entry(labeller, model):
  matches = [ entry for entry in LABELLER_PARAMS
    if ((entry['labeller'] == labeller) and (entry['model'] == model)) ]
  if len(matches) != 1:
    raise ValueError(f"Expected exactly one labeller entry for {labeller}/{model}, got {len(matches)}")
  return matches[0]


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
    '../envs/label_celltypes.yaml'
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
    model_f   = lambda wildcards: _get_labeller_entry("scprocess", wildcards.model)['model_f'],
    cls_f     = lambda wildcards: _get_labeller_entry("scprocess", wildcards.model)['cls_f'],
    genes_f   = lambda wildcards: _get_labeller_entry("scprocess", wildcards.model)['genes_f'],
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
    '../envs/label_celltypes.yaml'
  shell: """
    exec &>> {log}

    python3 scripts/label_celltypes.py xgboost_one_batch \
      {wildcards.batch} {params.batch_var} {wildcards.model} \
      --adata_f   {input.adata_f} \
      --model_f   {params.model_f} \
      --cls_f     {params.cls_f} \
      --genes_f   {params.genes_f} \
      --pred_f    {output.pred_f}
    """

# only do this bit if other parts are finished
qc_stats_f  = pathlib.Path(f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv')
if ('label_celltypes' in config) & qc_stats_f.is_file():

  # load up qc outputs
  bad_var     = 'bad_' + BATCH_VAR
  qc_stats    = pl.read_csv(qc_stats_f).filter( pl.col(bad_var) == False )

  # do labelling with celltypist
  rule merge_labels:
    input:
      pred_fs       = expand(f'{lbl_dir}/tmp_labels_{{labeller}}_model_{{model}}_{FULL_TAG}_{DATE_STAMP}_{{batch}}.csv.gz', 
        batch = qc_stats[BATCH_VAR].to_list(), allow_missing = True),
      integration_f = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz'
    output:
      pred_out_f    = f'{lbl_dir}/labels_{{labeller}}_model_{{model}}_{FULL_TAG}_{DATE_STAMP}.csv.gz'
    params:
      pred_fs_ls    = input.pred_fs,
      hi_res_cl     = lambda wildcards: _get_labeller_entry(wildcards.labeller, wildcards.model)["hi_res_cl"],
      min_cl_prop   = lambda wildcards: _get_labeller_entry(wildcards.labeller, wildcards.model)["min_cl_prop"],
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
      '../envs/label_celltypes.yaml'
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

