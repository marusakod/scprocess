# load modules
import os
import sys
import re
import glob
import pandas as pd
import polars as pl
import yaml
import warnings
from snakemake.utils import validate, min_version

# import utils
scprocess_dir = pathlib.Path(config.pop('scprocess_dir'))
sys.path.append(str(scprocess_dir / 'scripts'))
from scprocess_utils import *
from zoom import *


# define some things
proj_schema_f   = scprocess_dir / "resources/schemas/config.schema.json"
zoom_schema_f   = scprocess_dir / "resources/schemas/zoom.schema.json"
scdata_dir      = pathlib.Path(os.getenv('SCPROCESS_DATA_DIR'))
# check config
config          = check_config(config, proj_schema_f, scdata_dir, scprocess_dir)

# get lists of parameters
LIB_PARAMS, LIB_VAR = get_lib_parameters(config, scdata_dir)
LIBS                = list(LIB_PARAMS.keys())
RUN_PARAMS, RUN_VAR = get_run_parameters(config, scdata_dir, LIB_VAR, LIBS)
RUNS                = list(RUN_PARAMS.keys())
BATCH_PARAMS, BATCH_VAR, SAMPLES = get_batch_parameters(config, RUNS, scdata_dir)
BATCHES             = list(BATCH_PARAMS.keys())
RUNS_TO_BATCHES, RUNS_TO_SAMPLES, _ = get_runs_to_batches(config, RUNS, BATCHES, BATCH_VAR, LIBS)
RESOURCE_PARAMS     = prep_resource_params(config, proj_schema_f, scprocess_dir, LIB_PARAMS, BATCHES)
LABELLER_PARAMS     = get_labeller_parameters(config, proj_schema_f, scdata_dir)

# get zoom parameters
ZOOM_PARAMS         = get_zoom_parameters(config, zoom_schema_f, scdata_dir)
ZOOMS               = list(ZOOM_PARAMS.keys())

# unpack some variables that we use a lot
PROJ_DIR        = config['project']['proj_dir']
FULL_TAG        = config['project']['full_tag']
SHORT_TAG       = config['project']['short_tag']
DATE_STAMP      = config['project']['date_stamp']
GENOME_REF      = config['project'].get('probe_set', config['project'].get('ref_txome', ''))

# specify locations
benchmark_dir = f"{PROJ_DIR}/.resources"
logs_dir      = f"{PROJ_DIR}/.log"
code_dir      = f"{PROJ_DIR}/code"
amb_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_ambient"
qc_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_qc"
int_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_integration"
pb_dir        = f"{PROJ_DIR}/output/{SHORT_TAG}_pseudobulk"
hvg_dir       = f"{PROJ_DIR}/output/{SHORT_TAG}_hvg"
rmd_dir       = f"{PROJ_DIR}/analysis"
docs_dir      = f"{PROJ_DIR}/public"
zoom_dir      = f"{PROJ_DIR}/output/{SHORT_TAG}_zoom"

# define zoom marker outputs
zoom_mkr_report_outs = [
  file
  for zoom_name, mkr_sel_res, do_gsea in zip(
    ZOOMS,
    [ZOOM_PARAMS[zoom]['marker_genes']['mkr_sel_res'] for zoom in ZOOMS],
    [ZOOM_PARAMS[zoom]['marker_genes']['mkr_do_gsea'] for zoom in ZOOMS]
  )
  for file in (
    [
      '%s/%s/pb_%s_%s_%s_%s.rds' % (zoom_dir, zoom_name, FULL_TAG, zoom_name, mkr_sel_res, DATE_STAMP),
      '%s/%s/pb_marker_genes_%s_%s_%s_%s.csv.gz' % (zoom_dir, zoom_name, FULL_TAG, zoom_name, mkr_sel_res, DATE_STAMP),
      '%s/%s/pb_hvgs_%s_%s_%s_%s.csv.gz' % (zoom_dir, zoom_name, FULL_TAG, zoom_name, mkr_sel_res, DATE_STAMP),
      '%s/%s_zoom_%s_%s.Rmd' % (rmd_dir, SHORT_TAG, zoom_name, mkr_sel_res),
      '%s/%s_zoom_%s_%s.html' % (docs_dir, SHORT_TAG, zoom_name, mkr_sel_res)
    ] + (
      [
      '%s/%s/fgsea_%s_%s_%s_go_bp_%s.csv.gz' % (zoom_dir, zoom_name, FULL_TAG, zoom_name, mkr_sel_res, DATE_STAMP),
      '%s/%s/fgsea_%s_%s_%s_go_cc_%s.csv.gz' % (zoom_dir, zoom_name, FULL_TAG, zoom_name, mkr_sel_res, DATE_STAMP),
      '%s/%s/fgsea_%s_%s_%s_go_mf_%s.csv.gz' % (zoom_dir, zoom_name, FULL_TAG, zoom_name, mkr_sel_res, DATE_STAMP)
      ] if do_gsea and (GENOME_REF in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020', 'human_v1', 'mouse_v1', 'human_v2', 'mouse_v2'])
        else []
      )
  )
]

zooms_to_save_sce     = [ zoom_name for zoom_name in ZOOMS if ZOOM_PARAMS[zoom_name]['zoom']['save_subset_sces']]
zooms_to_save_anndata = [ zoom_name for zoom_name in ZOOMS if ZOOM_PARAMS[zoom_name]['zoom']['save_subset_anndata']]

ZOOM_OUT_MAP = {}
for zoom_name in ZOOMS:
  ZOOM_OUT_MAP[zoom_name] = {}
  for b in BATCHES:
    ZOOM_OUT_MAP[zoom_name][b] = {}
    if zoom_name in zooms_to_save_sce:
      ZOOM_OUT_MAP[zoom_name][b]["sce"]   = f"{zoom_dir}/{zoom_name}/sce_objects/sce_cells_clean_{b}_{FULL_TAG}_{zoom_name}_{DATE_STAMP}.rds"
    if zoom_name in zooms_to_save_anndata:
      ZOOM_OUT_MAP[zoom_name][b]["adata"] = f"{zoom_dir}/{zoom_name}/anndata_objects/anndata_cells_clean_{b}_{FULL_TAG}_{zoom_name}_{DATE_STAMP}.h5ad"

zoom_all_subset_fs = [
  path
  for name_dict in ZOOM_OUT_MAP.values()
  for b_dict in name_dict.values()
  for path in b_dict.values()
]

zoom_cluster_names_fs = [
  f'{zoom_dir}/{zoom_name}/cluster_names_for_shiny_{FULL_TAG}_{zoom_name}_{ZOOM_PARAMS[zoom_name]["marker_genes"]["mkr_sel_res"]}_{DATE_STAMP}.csv'
  for zoom_name in ZOOMS
  if ZOOM_PARAMS[zoom_name]['zoom']['save_cluster_names_file']
]

zoom_xgb_outs = []
for zoom_name in ZOOMS:
  if 'train_xgboost' in ZOOM_PARAMS[zoom_name]:
    _zxgb = ZOOM_PARAMS[zoom_name]['train_xgboost']
    _zxgb_ref_tag = _zxgb['ref_tag']
    _zxgb_mkr_sel_res = ZOOM_PARAMS[zoom_name]['marker_genes']['mkr_sel_res']
    zoom_xgb_outs.extend([
      f'{zoom_dir}/{zoom_name}/train_xgboost/{_zxgb_ref_tag}_xgboost_model.json',
      f'{zoom_dir}/{zoom_name}/train_xgboost/{_zxgb_ref_tag}_predictions.csv.gz',
      f'{docs_dir}/{SHORT_TAG}_zoom_{zoom_name}_{_zxgb_mkr_sel_res}.html',
    ])

rule zoom:
  input:
    # zoom sample qc
    expand(f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv',
      zoom_name = ZOOMS),
    # zoom pseudobulks and empties
    expand(f'{zoom_dir}/{{zoom_name}}/pb_cells_{BATCH_VAR}_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.rds',
      zoom_name = ZOOMS),
    expand(f'{zoom_dir}/{{zoom_name}}/edger_empty_genes_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz', 
      zoom_name = ZOOMS), 
    # zoom hvgs
    expand(f'{zoom_dir}/{{zoom_name}}/hvg_paths_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv', 
      zoom_name = ZOOMS),  
    expand(f'{zoom_dir}/{{zoom_name}}/standardized_variance_stats_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz', 
      zoom_name = ZOOMS), 
    expand(f'{zoom_dir}/{{zoom_name}}/hvg_dt_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz', 
      zoom_name = ZOOMS), 
    expand(f'{zoom_dir}/{{zoom_name}}/top_hvgs_counts_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.h5', 
      zoom_name = ZOOMS), 
    # zoom integration
    expand(f'{zoom_dir}/{{zoom_name}}/integrated_dt_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz', 
      zoom_name = ZOOMS),
    # zoom marker genes, fgsea and html report
    zoom_mkr_report_outs, 
    # zoom sce and anndata subsets (optional)
    zoom_all_subset_fs,
    # suggested cluster names for Shiny (optional)
    zoom_cluster_names_fs,


if zoom_xgb_outs:
  rule zoom_train_xgboost_all:
    input:
      zoom_xgb_outs



localrules: zoom_make_tmp_pb_cells_df, zoom_make_hvg_df, zoom_merge_group_mean_var, zoom_merge_stats_for_std_variance, zoom_copy_train_xgboost_r


rule zoom_check_clusters:
  input:
    labels_f    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['_original_labels_f']
  output:
    check_f     = f'{zoom_dir}/{{zoom_name}}/zoom_clusters_check_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.ok'
  params:
    labels_col  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'],
    sel_labels  = lambda wildcards: ','.join(ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels'])
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'zoom_check_clusters', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'zoom_check_clusters', 'time', attempt)
  log:
    f'{logs_dir}/zoom/zoom_check_clusters_{{zoom_name}}_{DATE_STAMP}.log'
  conda:
    '../envs/scprocess_local.yaml'
  shell: """
    exec &>> {log}

    PYTHONPATH='{scprocess_dir}/scripts' python3 -c "from zoom import zoom_check_clusters; \
      zoom_check_clusters('{input.labels_f}', '{params.labels_col}', '{params.sel_labels}', '{output.check_f}')"
    """


rule zoom_save_cluster_names:
  """Generate zoom cluster names from the retained parent cell-type predictions."""
  input:
    labels_f      = f'{zoom_dir}/{{zoom_name}}/filtered_labels_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz',
    integration_f = f'{zoom_dir}/{{zoom_name}}/integrated_dt_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  output:
    names_f       = f'{zoom_dir}/{{zoom_name}}/cluster_names_for_shiny_{FULL_TAG}_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.csv'
  retries: config['resources']['retries']
  log:
    f'{logs_dir}/zoom/zoom_save_cluster_names_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.log'
  conda:
    '../envs/scprocess_local.yaml'
  shell: """
    exec &>> {log}
    python3 {scprocess_dir}/scripts/label_celltypes.py save_cluster_names \
      --labels_f      {input.labels_f} \
      --integration_f {input.integration_f} \
      --mkr_sel_res   {wildcards.mkr_sel_res} \
      --output_f      {output.names_f}
    """


rule zoom_filter_cells_qc:
  input:
    labels_f    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['_original_labels_f'],
    check_f     = f'{zoom_dir}/{{zoom_name}}/zoom_clusters_check_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.ok',
    qc_all_f    = f'{qc_dir}/qc_all_samples_{FULL_TAG}_{DATE_STAMP}.csv.gz'
  output:
    filtered_f  = f'{zoom_dir}/{{zoom_name}}/filtered_labels_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  params:
    labels_col    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'],
    sel_labels    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels'],
    qc_min_counts = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_min_counts'),
    qc_max_counts = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_max_counts'),
    qc_min_feats  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_min_feats'),
    qc_max_feats  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_max_feats'),
    qc_min_mito   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_min_mito'),
    qc_max_mito   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_max_mito'),
    qc_min_splice = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_min_splice'),
    qc_max_splice = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_max_splice')
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'zoom_filter_cells_qc', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'zoom_filter_cells_qc', 'time', attempt)
  log:
    f'{logs_dir}/zoom/zoom_filter_cells_qc_{{zoom_name}}_{DATE_STAMP}.log'
  run:
    import sys
    with open(str(log), "a") as f:
      sys.stdout = f
      sys.stderr = f

      filter_zoom_labels_by_qc(
        labels_f      = str(input.labels_f),
        qc_all_f      = input.qc_all_f,
        output_f      = output.filtered_f,
        labels_col    = params.labels_col,
        sel_labels    = params.sel_labels,
        qc_min_counts = params.qc_min_counts,
        qc_max_counts = params.qc_max_counts,
        qc_min_feats  = params.qc_min_feats,
        qc_max_feats  = params.qc_max_feats,
        qc_min_mito   = params.qc_min_mito,
        qc_max_mito   = params.qc_max_mito,
        qc_min_splice = params.qc_min_splice,
        qc_max_splice = params.qc_max_splice
      )


rule zoom_aggregate_labels:
  """Re-aggregate retained naive predictions using the zoom's high-resolution clusters."""
  input:
    labels_f      = f'{zoom_dir}/{{zoom_name}}/filtered_labels_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz',
    integration_f = f'{zoom_dir}/{{zoom_name}}/integrated_dt_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  output:
    labels_f      = f'{zoom_dir}/{{zoom_name}}/aggregated_labels_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  params:
    hi_res_cl     = lambda wildcards: get_zoom_labeller_entry(
      ZOOM_PARAMS, LABELLER_PARAMS, wildcards.zoom_name)['hi_res_cl'],
    min_cl_prop   = lambda wildcards: get_zoom_labeller_entry(
      ZOOM_PARAMS, LABELLER_PARAMS, wildcards.zoom_name)['min_cl_prop'],
    label_map_f   = lambda wildcards: get_zoom_labeller_entry(
      ZOOM_PARAMS, LABELLER_PARAMS, wildcards.zoom_name).get('label_map_f', ''),
    batch_var     = BATCH_VAR
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input,
      'zoom_aggregate_labels', 'memory', attempt, csv_rule = 'merge_labels'),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input,
      'zoom_aggregate_labels', 'time', attempt, csv_rule = 'merge_labels')
  log:
    f'{logs_dir}/zoom/zoom_aggregate_labels_{{zoom_name}}_{DATE_STAMP}.log'
  benchmark:
    f'{benchmark_dir}/zoom/zoom_aggregate_labels_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/label_celltypes.yaml'
  shell: """
    exec &>> {log}

    python3 {scprocess_dir}/scripts/label_celltypes.py aggregate_predictions \
      {input.labels_f} \
      --int_f           {input.integration_f} \
      --hi_res_cl       {params.hi_res_cl} \
      --min_cl_prop     {params.min_cl_prop} \
      --batch_var       {params.batch_var} \
      --agg_f           {output.labels_f} \
      $( [ "{params.label_map_f}" != "" ] && echo "--label_map_f {params.label_map_f}" )
    """


rule zoom_copy_train_xgboost_r:
  input:
    src = scprocess_dir / "scripts" / "train_xgboost.R"
  output:
    dst = f"{code_dir}/train_xgboost.R"
  shell: "cp {input.src} {output.dst}"


rule zoom_get_sample_statistics:
  input:
    qc_stats_f      = f'{qc_dir}/qc_{BATCH_VAR}_statistics_{FULL_TAG}_{DATE_STAMP}.csv',
    zoom_lbls_f     = f'{zoom_dir}/{{zoom_name}}/filtered_labels_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  output:
    zoom_stats_f    = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv'
  params:
    zoom_lbls_col   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'],
    zoom_lbls       = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels'],
    batch_var       = BATCH_VAR,
    batches         = BATCHES,
    zoom_min_n_smpl = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc']['qc_min_cells'],
    ambient_method  = config['ambient']['ambient_method']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'zoom_get_sample_statistics', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'zoom_get_sample_statistics', 'time', attempt)
  log:
    f'{logs_dir}/zoom/zoom_get_sample_statistics_{{zoom_name}}_{DATE_STAMP}.log'
  run:
    import sys
    with open(str(log), "a") as f:
      rows = []
      sys.stdout = f
      sys.stderr = f

      zoom_stats_df   = extract_zoom_sample_statistics(input.qc_stats_f, input.zoom_lbls_f,
        params.zoom_lbls_col, params.zoom_lbls, params.batches, params.batch_var,
        params.zoom_min_n_smpl, params.ambient_method)
      zoom_stats_df.write_csv(output.zoom_stats_f)


# pseudobulks and empties
rule zoom_make_one_pb_cells:
  input:
    batch_lu_f    = f'{pb_dir}/runs_to_batches_{FULL_TAG}_{DATE_STAMP}.csv',
    filt_counts_f = lambda wildcards: get_filtered_counts_file(wildcards.run, config['ambient']['ambient_method'], amb_dir, DATE_STAMP),
    h5_paths_f    = f'{hvg_dir}/hvg_paths_{FULL_TAG}_{DATE_STAMP}.csv',
    coldata_f     = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    qc_stats_f    = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv',
    zoom_lbls_f   = f'{zoom_dir}/{{zoom_name}}/filtered_labels_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  output:
    pb_cells_f    = temp(f'{zoom_dir}/{{zoom_name}}/tmp_pb_cells_{{zoom_name}}_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds')
  params:
    run_var       = RUN_VAR,
    batch_var     = BATCH_VAR,
    zoom_lbls_col = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'],
    zoom_lbls     = lambda wildcards: ','.join(ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels'])
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_make_one_pb_cells', 'memory', attempt, wildcards.run, csv_rule = 'make_one_pb_cells'),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input,
      'zoom_make_one_pb_cells', 'time', attempt, wildcards.run, csv_rule = 'make_one_pb_cells')
  benchmark:
    f'{benchmark_dir}/zoom/zoom_make_one_pb_cells_{{zoom_name}}_{{run}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/zoom/zoom_make_one_pb_cells_{{zoom_name}}_{{run}}_{DATE_STAMP}.log'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    exec &>> {log}

    Rscript -e "source('{scprocess_dir}/scripts/utils.R'); source('{scprocess_dir}/scripts/pseudobulk_and_empties.R'); \
    make_pb_cells(
      sel_run     = '{wildcards.run}',
      batch_lu_f  = '{input.batch_lu_f}',
      qc_stats_f  = '{input.qc_stats_f}',
      h5_paths_f  = '{input.h5_paths_f}', 
      coldata_f   = '{input.coldata_f}',
      run_var     = '{params.run_var}',
      batch_var   = '{params.batch_var}',
      subset_f    = '{input.zoom_lbls_f}',
      subset_col  = '{params.zoom_lbls_col}',
      subset_str  = '{params.zoom_lbls}',
      pb_cells_f  = '{output.pb_cells_f}'
    )"
    """


rule zoom_make_tmp_pb_cells_df:
  input:
    pb_cells_fs = lambda wildcards: expand(
        f'{zoom_dir}/{{zoom_name}}/tmp_pb_cells_{{zoom_name}}_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds',
        run=RUNS,
        zoom_name=[wildcards.zoom_name]  
    )
  output:
    cells_paths_f = temp(f'{zoom_dir}/{{zoom_name}}/tmp_pb_cells_paths_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv')
  params:
    run_var = RUN_VAR,
    runs    = RUNS
  log:
    f'{logs_dir}/zoom/zoom_make_tmp_pb_cells_df_{{zoom_name}}_{DATE_STAMP}.log'
  run:
    import sys
    with open(str(log), "a") as f:
      rows = []
      sys.stdout = f
      sys.stderr = f
    
      import os
      import polars as pl

      paths_df = pl.DataFrame({
        params.run_var: params.runs,
        "pb_path": input.pb_cells_fs
      })
      paths_df = paths_df.filter(
        pl.col("pb_path").map_elements(os.path.getsize, return_dtype=pl.Int64) > 0
      )
      paths_df.write_csv(output.cells_paths_f)


rule zoom_merge_pb_cells:
  input:
    rowdata_f     = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    cells_paths_f = f'{zoom_dir}/{{zoom_name}}/tmp_pb_cells_paths_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv',
    pb_cells_fs   = lambda wildcards: expand(
        f'{zoom_dir}/{{zoom_name}}/tmp_pb_cells_{{zoom_name}}_{{run}}_{FULL_TAG}_{DATE_STAMP}.rds',
        run=RUNS, zoom_name=[wildcards.zoom_name])
  output:
    pb_cells_f    = f'{zoom_dir}/{{zoom_name}}/pb_cells_{BATCH_VAR}_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.rds'
  params:
    batch_var     = BATCH_VAR
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_merge_pb_cells', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_merge_pb_cells', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/zoom/zoom_merge_pb_cells_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/zoom/zoom_merge_pb_cells_{{zoom_name}}_{DATE_STAMP}.log'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    exec &>> {log}

    Rscript -e "source('{scprocess_dir}/scripts/utils.R'); source('{scprocess_dir}/scripts/pseudobulk_and_empties.R'); \
    merge_pbs_cells( \
      cells_paths_f = '{input.cells_paths_f}', 
      rowdata_f     = '{input.rowdata_f}',
      batch_var     = '{params.batch_var}',
      pb_cells_f    = '{output.pb_cells_f}'
    )"
    """


rule zoom_calculate_ambient_genes:
  input:
    pb_empty_f      = f'{pb_dir}/pb_empties_{FULL_TAG}_{DATE_STAMP}.rds', 
    zoom_pb_f       = f'{zoom_dir}/{{zoom_name}}/pb_cells_{BATCH_VAR}_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.rds'
  output:
    zoom_empty_gs_f = f'{zoom_dir}/{{zoom_name}}/edger_empty_genes_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  params:
    zoom_fdr_thr    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['pb_empties']['ambient_genes_fdr_thr'],
    zoom_logfc_thr  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['pb_empties']['ambient_genes_logfc_thr']
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_calculate_ambient_genes', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_calculate_ambient_genes', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/zoom/zoom_calculate_ambient_genes_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/zoom/zoom_calculate_ambient_genes_{{zoom_name}}_{DATE_STAMP}.log'
  conda: 
    '../envs/rlibs.yaml'
  shell: """
    exec &>> {log}

    Rscript -e "source('{scprocess_dir}/scripts/utils.R'); source('{scprocess_dir}/scripts/pseudobulk_and_empties.R'); \
    calc_empty_genes(
      pb_cells_f = '{input.zoom_pb_f}',
      pb_empty_f = '{input.pb_empty_f}',
      fdr_thr    = {params.zoom_fdr_thr}, 
      logfc_thr  = {params.zoom_logfc_thr},
      empty_gs_f = '{output.zoom_empty_gs_f}'
    )"
    """

# highly variable genes
rule zoom_make_hvg_df:
  input:
    amb_yaml_fs = expand([f'{amb_dir}/ambient_{{run}}/ambient_{{run}}_{DATE_STAMP}_output_paths.yaml'], run=RUNS)
  output:
    hvg_paths_f = f'{zoom_dir}/{{zoom_name}}/hvg_paths_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv'
  params:
    demux_type  = config['multiplexing']['demux_type'],
    run_var     = RUN_VAR,
    runs        = RUNS,
    mapping     = RUNS_TO_BATCHES,
    batch_var   = BATCH_VAR
  log: 
    f'{logs_dir}/zoom/zoom_make_hvg_df_{{zoom_name}}_{DATE_STAMP}.log'
  run:
    import sys
    with open(str(log), "a") as f:
      rows = []
      sys.stdout = f
      sys.stderr = f
    
      hvg_df = make_hvgs_input_df( 
        params.runs, input.amb_yaml_fs, params.run_var, params.batch_var, params.mapping,
        params.demux_type, FULL_TAG, DATE_STAMP, f"{zoom_dir}/{wildcards.zoom_name}"
      )
      hvg_df.write_csv(output.hvg_paths_f)


rule zoom_make_tmp_csr_matrix:
  input:
    rowdata_f       = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    hvg_paths_f     = f'{zoom_dir}/{{zoom_name}}/hvg_paths_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv',
    smpl_stats_f    = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv',
    zoom_lbls_f     = f'{zoom_dir}/{{zoom_name}}/filtered_labels_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  output:
    clean_h5_f      = temp(expand([
      f'{zoom_dir}/{{zoom_name}}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5'
      ], batch = BATCHES, allow_missing = True))
  params:
    zoom_lbls_col   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'],
    zoom_lbls       = lambda wildcards: ','.join(ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels']),
    run_var         = RUN_VAR,
    batch_var       = BATCH_VAR,
    demux_type      = config['multiplexing']['demux_type'],
    zoom_chunk_size = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_chunk_size'],
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_make_tmp_csr_matrix', 'memory', attempt, csv_rule = 'make_tmp_csr_matrix'),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input,
      'zoom_make_tmp_csr_matrix', 'time', attempt, csv_rule = 'make_tmp_csr_matrix')
  benchmark:
    f'{benchmark_dir}/zoom/zoom_make_tmp_csr_matrix_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/zoom/zoom_make_tmp_csr_matrix_{{zoom_name}}_{DATE_STAMP}.log'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}

    python3 {scprocess_dir}/scripts/hvgs.py get_csr_counts \
      {input.hvg_paths_f} \
      {input.zoom_lbls_f} \
      "{params.zoom_lbls_col}" \
      "{input.smpl_stats_f}" \
      {input.rowdata_f} \
      {params.run_var} \
      {params.batch_var} \
      {params.demux_type} \
      --keep_vals_str "{params.zoom_lbls}" \
      --chunksize {params.zoom_chunk_size} \
      --ncores {threads}
    """


rule zoom_get_stats_for_std_variance_for_sample:
  input: 
    rowdata_f     = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    clean_h5_f    = f'{zoom_dir}/{{zoom_name}}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5',
    smpl_stats_f  = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv'
  output:
    std_var_stats_f = temp(f'{zoom_dir}/{{zoom_name}}/tmp_std_var_stats_{{batch}}_sample_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz')
  params:
    batch_var = BATCH_VAR
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_stats_for_std_variance_for_sample', 'memory', attempt, csv_rule = 'get_stats_for_std_variance_for_sample'),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input,
      'zoom_get_stats_for_std_variance_for_sample', 'time', attempt, csv_rule = 'get_stats_for_std_variance_for_sample')
  benchmark:
    f'{benchmark_dir}/zoom/zoom_get_stats_for_std_variance_for_sample_{{zoom_name}}_{{batch}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/zoom/zoom_get_stats_for_std_variance_for_sample_{{zoom_name}}_{{batch}}_{DATE_STAMP}.log'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}

    python3 {scprocess_dir}/scripts/hvgs.py calculate_std_var_stats_for_sample \
      {wildcards.batch} \
      {params.batch_var} \
      {input.smpl_stats_f} \
      {input.clean_h5_f} \
      {input.rowdata_f} \
      {output.std_var_stats_f}
    """


rule zoom_get_mean_var_for_group:
  input:
    rowdata_f     = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    clean_h5_f    = expand(
      f'{zoom_dir}/{{zoom_name}}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5',
      batch = BATCHES, allow_missing = True),
    hvg_paths_f   = f'{zoom_dir}/{{zoom_name}}/hvg_paths_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv',
    smpl_stats_f  = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv'
  output: 
    mean_var_f    = temp(f'{zoom_dir}/{{zoom_name}}/tmp_mean_var_{{group}}_group_chunk_{{chunk}}_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz')
  params:
    metadata_f          = config['project']['sample_metadata'],
    batch_var           = BATCH_VAR,
    zoom_hvg_method     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_method'], 
    zoom_hvg_chunk_size = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_chunk_size'], 
    zoom_group_var      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_metadata_split_var']
  threads: 8
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_mean_var_for_group', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_mean_var_for_group', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/zoom/zoom_get_mean_var_for_group_{{zoom_name}}_{{group}}_chunk_{{chunk}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/zoom/zoom_get_mean_var_for_group_{{zoom_name}}_{{group}}_chunk_{{chunk}}_{DATE_STAMP}.log'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}

    GROUPVAR_FLAG=""
    if [ "{params.zoom_hvg_method}" = "groups" ]; then
      GROUPVAR_FLAG="--groupvar {params.zoom_group_var}"
    fi

    python3 {scprocess_dir}/scripts/hvgs.py calculate_mean_var_for_chunk \
      {input.hvg_paths_f} \
      {input.rowdata_f} \
      {params.metadata_f} \
      {input.smpl_stats_f} \
      {output.mean_var_f} \
      {wildcards.chunk} \
      {params.zoom_hvg_method} \
      {params.batch_var} \
      --chunksize {params.zoom_hvg_chunk_size} \
      --group {wildcards.group} \
      --ncores {threads} \
      $GROUPVAR_FLAG
    """


rule zoom_merge_group_mean_var:
  input:         
    mean_var_f    = lambda wildcards: get_zoom_raw_mean_var_files(wildcards.zoom_name, zoom_dir, ZOOM_PARAMS, FULL_TAG, DATE_STAMP)
  output:
    mean_var_merged_f = temp(f'{zoom_dir}/{{zoom_name}}/means_variances_dt_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz')
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_merge_group_mean_var', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_merge_group_mean_var', 'time', attempt)
  log:
    f'{logs_dir}/zoom/zoom_get_mean_var_for_group_{{zoom_name}}_{DATE_STAMP}.log'
  run:
    import sys
    with open(str(log), "a") as f:
      rows = []
      sys.stdout = f
      sys.stderr = f
      
      merge_tmp_files(input.mean_var_f, output.mean_var_merged_f)


rule zoom_get_estimated_variances:
  input:
    mean_var_merged_f = f'{zoom_dir}/{{zoom_name}}/means_variances_dt_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  output:
    estim_vars_f      = temp(f'{zoom_dir}/{{zoom_name}}/estimated_variances_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz')
  params: 
    zoom_hvg_method   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_method'],
    batch_var         = BATCH_VAR
  threads: 1
  retries: config['resources']['retries']
  conda:
    '../envs/hvgs.yaml'
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_estimated_variances', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_estimated_variances', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/zoom/zoom_get_estimated_variances_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/zoom/zoom_get_estimated_variances_{{zoom_name}}_{DATE_STAMP}.log'
  shell: """
    exec &>> {log}

    python3 {scprocess_dir}/scripts/hvgs.py calculate_estimated_vars \
      {output.estim_vars_f} \
      {params.zoom_hvg_method} \
      {params.batch_var} \
      {input.mean_var_merged_f}
    """


rule zoom_get_stats_for_std_variance_for_group:
  input: 
    rowdata_f     = f'{qc_dir}/rowdata_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz', 
    clean_h5_fs   = expand(
      f'{zoom_dir}/{{zoom_name}}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5',
      batch = BATCHES, allow_missing = True
    ),
    estim_vars_f  = f'{zoom_dir}/{{zoom_name}}/estimated_variances_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz', 
    hvg_paths_f   = f'{zoom_dir}/{{zoom_name}}/hvg_paths_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv',
    smpl_stats_f  = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv'
  output:
    std_var_stats_f = temp(f'{zoom_dir}/{{zoom_name}}/tmp_std_var_stats_{{group}}_group_chunk_{{chunk}}_{FULL_TAG}_{DATE_STAMP}.csv.gz')
  params:
    metadata_f          = config['project']['sample_metadata'],
    batch_var       = BATCH_VAR,
    zoom_hvg_method     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_method'], 
    zoom_hvg_chunk_size = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_chunk_size'], 
    zoom_hvg_group_var  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_metadata_split_var']
  threads: 8
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_stats_for_std_variance_for_group', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_stats_for_std_variance_for_group', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/zoom/zoom_get_stats_for_std_variance_for_group_{{zoom_name}}_{{group}}_chunk_{{chunk}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/zoom/zoom_get_stats_for_std_variance_for_group_{{zoom_name}}_{{group}}_chunk_{{chunk}}_{DATE_STAMP}.log'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}

    python3 {scprocess_dir}/scripts/hvgs.py calculate_std_var_stats_for_chunk \
      {input.hvg_paths_f} \
      {input.rowdata_f} \
      {params.metadata_f} \
      {input.smpl_stats_f} \
      {output.std_var_stats_f} \
      {input.estim_vars_f} \
      {wildcards.chunk} \
      {params.zoom_hvg_method} \
      {params.batch_var} \
      --chunksize {params.zoom_hvg_chunk_size} \
      --group {wildcards.group} \
      --groupvar {params.zoom_hvg_group_var} \
      --ncores {threads} 
    """


rule zoom_merge_stats_for_std_variance:
  input:
    tmp_std_var_stats_fs = lambda wildcards: get_zoom_std_var_stats_files(wildcards.zoom_name, \
      zoom_dir, ZOOM_PARAMS, FULL_TAG, DATE_STAMP, BATCHES)
  output:
    std_var_stats_merged_f= f'{zoom_dir}/{{zoom_name}}/standardized_variance_stats_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  log:
    f'{logs_dir}/zoom/zoom_merge_stats_for_std_variance_{{zoom_name}}_{DATE_STAMP}.log'
  run:
    import sys
    with open(str(log), "a") as f:
      rows = []
      sys.stdout = f
      sys.stderr = f
    
      merge_tmp_files(input.tmp_std_var_stats_fs, output.std_var_stats_merged_f)

        
rule zoom_get_highly_variable_genes:
  input:
    std_var_stats_f = f'{zoom_dir}/{{zoom_name}}/standardized_variance_stats_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz', 
    empty_gs_fs     = f'{zoom_dir}/{{zoom_name}}/edger_empty_genes_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz' 
  output:
    hvg_f           = f'{zoom_dir}/{{zoom_name}}/hvg_dt_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  threads: 1
  retries: config['resources']['retries']
  params:
    zoom_hvg_method   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_method'],
    zoom_n_hvgs       = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_n_hvgs'],
    batch_var         = BATCH_VAR,
    zoom_exc_gs_f     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_exclude_from_file'],
    zoom_exc_ambient  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['hvg']['hvg_exclude_ambient_genes']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_get_highly_variable_genes', 'memory', attempt, csv_rule = 'get_highly_variable_genes'),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input,
      'zoom_get_highly_variable_genes', 'time', attempt, csv_rule = 'get_highly_variable_genes')
  benchmark:
    f'{benchmark_dir}/zoom/zoom_get_highly_variable_genes_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/zoom/zoom_get_highly_variable_genes_{{zoom_name}}_{DATE_STAMP}.log'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}

    NOAMBIENT_FLAG=""
    if [ "{params.zoom_exc_ambient}" = "True" ]; then
      NOAMBIENT_FLAG="--noambient"
    fi
    EXC_GS_F_FLAG=""
    if [ "{params.zoom_exc_gs_f}" != "None" ]; then
      EXC_GS_F_FLAG="--exc_gs_f {params.zoom_exc_gs_f}"
    fi

    python3 {scprocess_dir}/scripts/hvgs.py calculate_hvgs \
      {input.std_var_stats_f} \
      {output.hvg_f} \
      {params.zoom_hvg_method} \
      {params.batch_var} \
      {params.zoom_n_hvgs} \
      --empty_gs_f {input.empty_gs_fs} \
      $NOAMBIENT_FLAG \
      $EXC_GS_F_FLAG
    """


rule zoom_create_hvg_matrix:
  input: 
    clean_h5_f    = expand(
      f'{zoom_dir}/{{zoom_name}}/chunked_counts_{{batch}}_{FULL_TAG}_{DATE_STAMP}.h5',
      zoom_name = ZOOMS, batch = BATCHES
    ),
    smpl_stats_f  = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv', 
    hvg_paths_f   = f'{zoom_dir}/{{zoom_name}}/hvg_paths_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv', 
    hvg_f         = f'{zoom_dir}/{{zoom_name}}/hvg_dt_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  output:
    hvg_mat_f     = f'{zoom_dir}/{{zoom_name}}/top_hvgs_counts_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.h5'
  params:
    demux_type    = config['multiplexing']['demux_type'], 
    batch_var     = BATCH_VAR
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_create_hvg_matrix', 'memory', attempt, csv_rule = 'create_hvg_matrix'),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input,
      'zoom_create_hvg_matrix', 'time', attempt, csv_rule = 'create_hvg_matrix')
  benchmark:
    f'{benchmark_dir}/zoom/zoom_create_hvg_matrix_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/zoom/zoom_create_hvg_matrix_{{zoom_name}}_{DATE_STAMP}.log'
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}

    python3 {scprocess_dir}/scripts/hvgs.py create_hvg_matrix \
      {input.smpl_stats_f} \
      {input.hvg_paths_f} \
      {input.hvg_f} \
      {output.hvg_mat_f} \
      {params.demux_type} \
      {params.batch_var} \
      --ncores {threads}
    """


rule zoom_run_integration:
  input:
    coldata_f     = f'{qc_dir}/coldata_dt_all_cells_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    hvg_mat_f     = f'{zoom_dir}/{{zoom_name}}/top_hvgs_counts_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.h5',
    sample_qc_f   = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv'
  output:
    integration_f = f'{zoom_dir}/{{zoom_name}}/integrated_dt_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  params:
    demux_type            = config['multiplexing']['demux_type'],
    exclude_mito          = config['qc']['exclude_mito'],
    zoom_int_embedding    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_embedding'],
    zoom_int_n_dims       = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_n_dims'],
    zoom_int_cl_method    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_cl_method'],
    zoom_int_theta        = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_theta'],
    zoom_int_use_paga     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_use_paga'],
    zoom_int_paga_cl_res  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_paga_cl_res'],
    zoom_int_res_ls       = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_res_ls'],
    zoom_int_use_gpu      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_use_gpu'],
    batch_var             = BATCH_VAR,
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_run_integration', 'memory', attempt, csv_rule = 'run_integration'),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input,
      'zoom_run_integration', 'time', attempt, csv_rule = 'run_integration')
  benchmark: 
    f'{benchmark_dir}/zoom/zoom_run_integration_{{zoom_name}}_{DATE_STAMP}.benchmark.txt'
  log: 
    f'{logs_dir}/zoom/zoom_run_integration_{{zoom_name}}_{DATE_STAMP}.log'
  conda: 
    '../envs/integration.yaml'
  shell: """
    exec &>> {log}

    python3 {scprocess_dir}/scripts/integration.py run_zoom_integration \
      --hvg_mat_f     {input.hvg_mat_f} \
      --sample_qc_f   {input.sample_qc_f} \
      --coldata_f     {input.coldata_f} \
      --demux_type    {params.demux_type} \
      --exclude_mito  "{params.exclude_mito}" \
      --embedding     {params.zoom_int_embedding} \
      --n_dims        {params.zoom_int_n_dims} \
      --cl_method     {params.zoom_int_cl_method} \
      --theta         {params.zoom_int_theta} \
      --res_ls_concat "{params.zoom_int_res_ls}" \
      --integration_f {output.integration_f} \
      --batch_var     {params.batch_var} \
      $( [ "{params.zoom_int_use_paga}" == "True" ] && echo "--use-paga" ) \
      $( [ "{params.zoom_int_use_paga}" == "True" ] && echo "--paga-cl-res {params.zoom_int_paga_cl_res}" ) \
      $( [ "{params.zoom_int_use_gpu}" == "True" ] && echo "--use-gpu" )
    """


rule zoom_run_marker_genes:
  input:
    h5ads_yaml_f    = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml',
    integration_f   = f'{zoom_dir}/{{zoom_name}}/integrated_dt_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  output:
    pb_f            = f'{zoom_dir}/{{zoom_name}}/pb_{FULL_TAG}_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.rds',
    mkrs_f          = f'{zoom_dir}/{{zoom_name}}/pb_marker_genes_{FULL_TAG}_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.csv.gz',
    pb_hvgs_f       = f'{zoom_dir}/{{zoom_name}}/pb_hvgs_{FULL_TAG}_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.csv.gz'
  params:
    genome_ref          = GENOME_REF,
    af_gtf_dt_f         = config['mapping_af']['gene_info_f'],
    batch_var           = BATCH_VAR,
    zoom_mkr_sel_res     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_sel_res'],
    zoom_mkr_min_cl_size = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cl_size'], 
    zoom_mkr_min_cells   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cells']
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_run_marker_genes', 'memory', attempt, csv_rule = 'run_marker_genes'),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input,
      'zoom_run_marker_genes', 'time', attempt, csv_rule = 'run_marker_genes')
  benchmark:
    f'{benchmark_dir}/zoom/zoom_run_marker_genes_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/zoom/zoom_run_marker_genes_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.log'
  conda:
    '../envs/rlibs.yaml'
  shell: """
    exec &>> {log}

    Rscript -e "source('{scprocess_dir}/scripts/utils.R'); source('{scprocess_dir}/scripts/marker_genes.R'); calculate_marker_genes(
      integration_f = '{input.integration_f}', 
      h5ads_yaml_f  = '{input.h5ads_yaml_f}',
      pb_f          = '{output.pb_f}',
      mkrs_f        = '{output.mkrs_f}',
      pb_hvgs_f     = '{output.pb_hvgs_f}',
      gtf_dt_f      = '{params.af_gtf_dt_f}',
      sel_res       = '{params.zoom_mkr_sel_res}',
      min_cl_size   =  {params.zoom_mkr_min_cl_size},
      min_cells     =  {params.zoom_mkr_min_cells},
      zoom          = 'True', 
      batch_var     = '{params.batch_var}', 
      n_cores       =  {threads})"
    """


rule zoom_run_fgsea:
  input:
    mkrs_f        = f'{zoom_dir}/{{zoom_name}}/pb_marker_genes_{FULL_TAG}_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.csv.gz'
  output:
    fgsea_go_bp_f = f'{zoom_dir}/{{zoom_name}}/fgsea_{FULL_TAG}_{{zoom_name}}_{{mkr_sel_res}}_go_bp_{DATE_STAMP}.csv.gz', 
    fgsea_go_cc_f = f'{zoom_dir}/{{zoom_name}}/fgsea_{FULL_TAG}_{{zoom_name}}_{{mkr_sel_res}}_go_cc_{DATE_STAMP}.csv.gz',
    fgsea_go_mf_f = f'{zoom_dir}/{{zoom_name}}/fgsea_{FULL_TAG}_{{zoom_name}}_{{mkr_sel_res}}_go_mf_{DATE_STAMP}.csv.gz'
  params:
    genome_ref           = GENOME_REF,
    mkr_gsea_dir         = config['marker_genes']['mkr_gsea_dir'],
    zoom_mkr_min_cpm_go  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cpm_go'],
    zoom_mkr_max_zero_p  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_max_zero_p'],
    zoom_mkr_gsea_cut    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_gsea_cut'], 
    zoom_mkr_not_ok_re   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_not_ok_re'],
    zoom_mkr_gsea_var    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_gsea_var']
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 
      'zoom_run_fgsea', 'memory', attempt, csv_rule = 'run_fgsea'),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input,
      'zoom_run_fgsea', 'time', attempt, csv_rule = 'run_fgsea')
  benchmark:
    f'{benchmark_dir}/zoom/zoom_run_fgsea_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/zoom/zoom_run_fgsea_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.log'
  conda: '../envs/rlibs.yaml'
  shell:"""
    exec &>> {log}

    Rscript -e "source('{scprocess_dir}/scripts/utils.R'); source('{scprocess_dir}/scripts/fgsea.R'); run_fgsea(
      mkrs_f        = '{input.mkrs_f}', 
      fgsea_go_bp_f = '{output.fgsea_go_bp_f}', 
      fgsea_go_cc_f = '{output.fgsea_go_cc_f}', 
      fgsea_go_mf_f = '{output.fgsea_go_mf_f}', 
      genome_ref    = '{params.genome_ref}',
      gsea_dir      = '{params.mkr_gsea_dir}',
      min_cpm_go    = {params.zoom_mkr_min_cpm_go},
      max_zero_p    = {params.zoom_mkr_max_zero_p},
      gsea_cut      = {params.zoom_mkr_gsea_cut},
      not_ok_re     = '{params.zoom_mkr_not_ok_re}',
      gsea_var      = '{params.zoom_mkr_gsea_var}',
      n_cores       =  {threads})"
    """

rule zoom_make_subsets:
  input:
    h5ads_yaml_f  = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml',
    integration_f = f'{zoom_dir}/{{zoom_name}}/integrated_dt_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz',
    smpl_stats_f  = f'{zoom_dir}/{{zoom_name}}/zoom_{BATCH_VAR}_statistics_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv',
    zoom_lbls_f   = f'{zoom_dir}/{{zoom_name}}/filtered_labels_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  output:
    f"{zoom_dir}/{{zoom_name}}/{{prefix}}_objects/{{prefix}}_cells_clean_{{batch}}_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.{{ext}}"
  params:
    batch_var     = BATCH_VAR,
    zoom_lbls_col = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'],
    zoom_lbls     = lambda wildcards: ','.join(ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels']),
    save_sce      = lambda wildcards: "TRUE" if wildcards.zoom_name in zooms_to_save_sce else "FALSE",
    save_adata    = lambda wildcards: "TRUE" if wildcards.zoom_name in zooms_to_save_anndata else "FALSE",
    sce_path      = lambda wildcards: ZOOM_OUT_MAP[wildcards.zoom_name][wildcards.batch].get("sce", ""),
    adata_path    = lambda wildcards: ZOOM_OUT_MAP[wildcards.zoom_name][wildcards.batch].get("adata", "")
  threads: 1
  benchmark:
    f'{benchmark_dir}/zoom/zoom_make_subsets_{{zoom_name}}_{{prefix}}_{{ext}}_{{batch}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/zoom/zoom_make_subsets_{{zoom_name}}_{{prefix}}_{{ext}}_{{batch}}_{DATE_STAMP}.log'
  resources:
    mem_mb  = lambda w, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'zoom_make_subsets', 'memory', attempt),
    runtime = lambda w, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'zoom_make_subsets', 'time', attempt)
  conda: '../envs/rlibs.yaml'
  shell:"""
    exec &>> {log}

    Rscript -e "source('{scprocess_dir}/scripts/zoom.R');
    make_subset_objects(
      sel_b         = '{wildcards.batch}',
      batch_var     = '{params.batch_var}',
      smpl_stats_f  = '{input.smpl_stats_f}',
      h5ads_yaml_f  = '{input.h5ads_yaml_f}',
      subset_f      = '{input.zoom_lbls_f}',
      subset_col    = '{params.zoom_lbls_col}',
      subset_str    = '{params.zoom_lbls}',
      integration_f = '{input.integration_f}',
      save_sce      = {params.save_sce},
      subset_sce_f  = '{params.sce_path}',
      save_adata    = {params.save_adata},
      subset_h5ad_f = '{params.adata_path}'
    )"
    """


rule zoom_train_xgboost:
  input:
    annots_f    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost']['annots_f'],
    cluster_csv = f'{zoom_dir}/{{zoom_name}}/integrated_dt_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz',
    h5ads_yaml  = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml',
  output:
    model_f  = f'{zoom_dir}/{{zoom_name}}/train_xgboost/{{ref_tag}}_xgboost_model.json',
    cls_f    = f'{zoom_dir}/{{zoom_name}}/train_xgboost/{{ref_tag}}_allowed_cls.csv',
    genes_f  = f'{zoom_dir}/{{zoom_name}}/train_xgboost/{{ref_tag}}_selected_genes.txt',
    imp_f    = f'{zoom_dir}/{{zoom_name}}/train_xgboost/{{ref_tag}}_gene_importance.csv',
    preds_f  = f'{zoom_dir}/{{zoom_name}}/train_xgboost/{{ref_tag}}_predictions.csv.gz',
    pb_f     = f'{zoom_dir}/{{zoom_name}}/train_xgboost/{{ref_tag}}_pseudobulk.h5ad',
  params:
    output_dir           = lambda wildcards: f'{zoom_dir}/{wildcards.zoom_name}/train_xgboost',
    batch_var            = BATCH_VAR,
    int_res_ls           = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_res_ls'],
    label_map_f          = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('label_map_f') or '',
    refine_labels        = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('refine_labels', True),
    purity_threshold     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('purity_threshold', 0.65),
    n_cells_per_type     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('n_cells_per_type', 1000),
    min_cells_per_type   = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('min_cells_per_type', 20),
    min_cells_expressed  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('min_cells_expressed', 10),
    gene_exclude_re      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('gene_exclude_re', '(lincRNA|lncRNA|pseudogene|antisense)'),
    seed                 = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('seed', 42),
    use_gpu              = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('use_gpu', False),
    pass1_subsample      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('pass1_subsample', 0.632),
    pass1_colsample_bytree = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('pass1_colsample_bytree', 0.1),
    pass1_learning_rate  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('pass1_learning_rate', 0.1),
    pass1_nrounds        = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('pass1_nrounds', 300),
    pass1_early_stopping = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('pass1_early_stopping', 10),
    pass2_colsample_bytree = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('pass2_colsample_bytree', 0.5),
    pass2_learning_rate  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('pass2_learning_rate', 0.05),
    pass2_nrounds        = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('pass2_nrounds', 500),
    pass2_early_stopping = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('pass2_early_stopping', 10),
    gain_threshold       = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('gain_threshold', 0.9),
    min_genes            = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('min_genes', 100),
    max_genes            = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['train_xgboost'].get('max_genes', 3000),
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input,
      'zoom_train_xgboost', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input,
      'zoom_train_xgboost', 'time', attempt)
  log:
    f'{logs_dir}/zoom/zoom_train_xgboost_{{zoom_name}}_{{ref_tag}}_{DATE_STAMP}.log'
  benchmark:
    f'{benchmark_dir}/zoom/zoom_train_xgboost_{{zoom_name}}_{{ref_tag}}_{DATE_STAMP}.benchmark.txt'
  conda:
    '../envs/label_celltypes.yaml'
  shell: """
    exec &>> {log}
    python3 {scprocess_dir}/scripts/train_xgboost.py train \
      --annots_f          {input.annots_f} \
      --cluster_csv       {input.cluster_csv} \
      --h5ads_yaml        {input.h5ads_yaml} \
      --ref_tag           {wildcards.ref_tag} \
      --output_dir        {params.output_dir} \
      --batch_var         {params.batch_var} \
      --int_res_ls        "{params.int_res_ls}" \
      --refine_labels     {params.refine_labels} \
      --purity_threshold  {params.purity_threshold} \
      --n_cells_per_type  {params.n_cells_per_type} \
      --min_cells_per_type {params.min_cells_per_type} \
      --min_cells_expressed {params.min_cells_expressed} \
      --gene_exclude_re  "{params.gene_exclude_re}" \
      --seed              {params.seed} \
      --n_cores           {threads} \
      --pass1_subsample      {params.pass1_subsample} \
      --pass1_colsample_bytree {params.pass1_colsample_bytree} \
      --pass1_learning_rate  {params.pass1_learning_rate} \
      --pass1_nrounds        {params.pass1_nrounds} \
      --pass1_early_stopping {params.pass1_early_stopping} \
      --pass2_colsample_bytree {params.pass2_colsample_bytree} \
      --pass2_learning_rate  {params.pass2_learning_rate} \
      --pass2_nrounds        {params.pass2_nrounds} \
      --pass2_early_stopping {params.pass2_early_stopping} \
      --gain_threshold    {params.gain_threshold} \
      --min_genes         {params.min_genes} \
      --max_genes         {params.max_genes} \
      $( [ "{params.label_map_f}" != "" ] && echo "--label_map_f {params.label_map_f}" ) \
      $( [ "{params.use_gpu}" == "True" ] && echo "--use_gpu" )
    """


# render_html_zoom
rule render_html_zoom:
  input:
    unpack(lambda wildcards: get_zoom_conditional_fgsea_files(GENOME_REF, zoom_dir,
        FULL_TAG, DATE_STAMP, ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_do_gsea'])),
    unpack(lambda wildcards: get_zoom_conditional_xgboost_files(
        ZOOM_PARAMS[wildcards.zoom_name], zoom_dir, wildcards.zoom_name, code_dir)),
    rmd_template_f        = str(scprocess_dir / "resources/rmd_templates/zoom.Rmd.template"),
    r_utils_f             = f'{code_dir}/utils.R',
    r_hvgs_f              = f'{code_dir}/hvgs.R',
    r_int_f               = f'{code_dir}/integration.R',
    r_mkr_f               = f'{code_dir}/marker_genes.R',
    r_label_celltypes_f   = str(scprocess_dir / "scripts/label_celltypes.R"),
    metadata_f            = config['project']['sample_metadata'],
    qc_f                  = f'{qc_dir}/qc_all_samples_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    zoom_lbls_f           = f'{zoom_dir}/{{zoom_name}}/filtered_labels_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz',
    zoom_report_lbls_f    = lambda wildcards: get_zoom_report_labels_f(
      ZOOM_PARAMS, zoom_dir, FULL_TAG, DATE_STAMP, wildcards.zoom_name),
    pb_empty_f            = f'{pb_dir}/pb_empties_{FULL_TAG}_{DATE_STAMP}.rds',
    zoom_cell_hvgs_f      = f'{zoom_dir}/{{zoom_name}}/hvg_dt_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz',
    zoom_int_f            = f'{zoom_dir}/{{zoom_name}}/integrated_dt_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz',
    zoom_pb_f             = f'{zoom_dir}/{{zoom_name}}/pb_{FULL_TAG}_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.rds',
    zoom_pb_hvgs_f        = f'{zoom_dir}/{{zoom_name}}/pb_hvgs_{FULL_TAG}_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.csv.gz',
    zoom_mkrs_f           = f'{zoom_dir}/{{zoom_name}}/pb_marker_genes_{FULL_TAG}_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.csv.gz',
    zoom_empty_gs_f       = f'{zoom_dir}/{{zoom_name}}/edger_empty_genes_{FULL_TAG}_{{zoom_name}}_{DATE_STAMP}.csv.gz'
  output:
    rmd_f                 = f'{rmd_dir}/{SHORT_TAG}_zoom_{{zoom_name}}_{{mkr_sel_res}}.Rmd',
    html_f                = f'{docs_dir}/{SHORT_TAG}_zoom_{{zoom_name}}_{{mkr_sel_res}}.html'
  params:
    your_name             = config['project']['your_name'],
    affiliation           = config['project']['affiliation'],
    short_tag             = config['project']['short_tag'],
    date_stamp            = config['project']['date_stamp'],
    proj_dir              = config['project']['proj_dir'],
    genome_ref            = GENOME_REF,
    zoom_dir              = zoom_dir,
    batch_var             = BATCH_VAR,
    zoom_original_lbls_f  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['_original_labels_f'],
    zoom_lbls_col         = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_col'],
    zoom_label_hi_res_cl  = lambda wildcards: (
      get_zoom_labeller_entry(
        ZOOM_PARAMS, LABELLER_PARAMS, wildcards.zoom_name)['hi_res_cl']
      if ZOOM_PARAMS[wildcards.zoom_name]['zoom']['labels_source'] in ['celltypist', 'scprocess']
      else ''),
    zoom_sel_labels        = lambda wildcards: ','.join(ZOOM_PARAMS[wildcards.zoom_name]['zoom']['sel_labels']),
    zoom_qc_min_counts    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_min_counts', ''),
    zoom_qc_max_counts    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_max_counts', ''),
    zoom_qc_min_feats     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_min_feats', ''),
    zoom_qc_max_feats     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_max_feats', ''),
    zoom_qc_min_mito      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_min_mito', ''),
    zoom_qc_max_mito      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_max_mito', ''),
    zoom_qc_min_splice    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_min_splice', ''),
    zoom_qc_max_splice    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['qc'].get('qc_max_splice', ''),
    meta_vars             = ','.join(config['project']['metadata_vars']),
    fgsea_args            = lambda wildcards, input: ", ".join([
      f"fgsea_go_bp_f = '{input.get('fgsea_go_bp_f', '')}'",
      f"fgsea_go_cc_f = '{input.get('fgsea_go_cc_f', '')}'",
      f"fgsea_go_mf_f = '{input.get('fgsea_go_mf_f', '')}'"
    ]),
    af_gtf_dt_f           = config['mapping_af']['gene_info_f'],
    zoom_int_res_ls       = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['integration']['int_res_ls'],
    zoom_mkr_sel_res      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_sel_res'],
    zoom_mkr_min_cpm_mkr  = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cpm_mkr'],
    zoom_mkr_min_cells    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_min_cells'],
    zoom_mkr_not_ok_re    = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_not_ok_re'],
    zoom_mkr_do_gsea      = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_do_gsea'],
    zoom_mkr_gsea_var     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_gsea_var'],
    zoom_mkr_gsea_cut     = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['mkr_gsea_cut'],
    zoom_custom_mkr_names = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['custom_mkr_names'],
    zoom_custom_mkr_paths = lambda wildcards: ZOOM_PARAMS[wildcards.zoom_name]['marker_genes']['custom_mkr_paths'],
    do_xgboost            = lambda wildcards: 'train_xgboost' in ZOOM_PARAMS[wildcards.zoom_name],
    xgb_args              = lambda wildcards, input: ", ".join([
      f"xgb_predictions_f = '{input.get('xgb_preds_f', '')}'",
      f"xgb_importance_f  = '{input.get('xgb_imp_f', '')}'",
      f"xgb_pseudobulk_f  = '{input.get('xgb_pb_f', '')}'",
      f"xgb_has_coarse    = '{'true' if ZOOM_PARAMS[wildcards.zoom_name].get('train_xgboost', {}).get('label_map_f') else 'false'}'",
      f"xgb_min_cells     = '{ZOOM_PARAMS[wildcards.zoom_name].get('train_xgboost', {}).get('min_cells_expressed', 10)}'"
    ])
  threads: 1
  retries: config['resources']['retries']
  conda:
    '../envs/rlibs.yaml'
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input,
      'render_html_zoom', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input,
      'render_html_zoom', 'time', attempt)
  benchmark:
    f'{benchmark_dir}/zoom/render_html_zoom_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.benchmark.txt'
  log:
    f'{logs_dir}/zoom/render_html_zoom_{{zoom_name}}_{{mkr_sel_res}}_{DATE_STAMP}.log'
  shell: """
    exec &>> {log}

    cp {scprocess_dir}/scripts/utils.R {input.r_utils_f}
    cp {scprocess_dir}/scripts/SampleQC.R $(dirname {input.r_utils_f})/qc.R

    template_f=$(realpath {input.rmd_template_f})
    rule="zoom"

    Rscript --vanilla -e "source('{scprocess_dir}/scripts/render_htmls.R'); \
    render_html(
      rule_name         = '$rule',
      proj_dir          = '{PROJ_DIR}',
      temp_f            =  '$template_f',
      rmd_f             = '{output.rmd_f}',
      your_name         = '{params.your_name}',
      affiliation       = '{params.affiliation}',
      short_tag         = '{params.short_tag}',
      date_stamp        = '{params.date_stamp}',
      threads           =  {threads},
      zoom_dir          = '{params.zoom_dir}',
      zoom_name         = '{wildcards.zoom_name}',
      metadata_f        = '{input.metadata_f}',
      meta_vars_ls      = '{params.meta_vars}',
      gtf_dt_f          = '{params.af_gtf_dt_f}',
      qc_f              = '{input.qc_f}',
      r_label_celltypes_f = '{input.r_label_celltypes_f}',
      zoom_lbls_f       = '{input.zoom_lbls_f}',
      zoom_report_lbls_f = '{input.zoom_report_lbls_f}',
      zoom_original_lbls_f = '{params.zoom_original_lbls_f}',
      zoom_lbls_col     = '{params.zoom_lbls_col}',
      zoom_label_hi_res_cl = '{params.zoom_label_hi_res_cl}',
      zoom_sel_labels   = '{params.zoom_sel_labels}',
      zoom_qc_min_counts = '{params.zoom_qc_min_counts}',
      zoom_qc_max_counts = '{params.zoom_qc_max_counts}',
      zoom_qc_min_feats  = '{params.zoom_qc_min_feats}',
      zoom_qc_max_feats  = '{params.zoom_qc_max_feats}',
      zoom_qc_min_mito   = '{params.zoom_qc_min_mito}',
      zoom_qc_max_mito   = '{params.zoom_qc_max_mito}',
      zoom_qc_min_splice = '{params.zoom_qc_min_splice}',
      zoom_qc_max_splice = '{params.zoom_qc_max_splice}',
      cell_hvgs_f       = '{input.zoom_cell_hvgs_f}',
      int_f             = '{input.zoom_int_f}',
      pb_f              = '{input.zoom_pb_f}',
      mkrs_f            = '{input.zoom_mkrs_f}',
      pb_hvgs_f         = '{input.zoom_pb_hvgs_f}',
      empty_gs_f        = '{input.zoom_empty_gs_f}',
      pb_empty_f        = '{input.pb_empty_f}',
      {params.fgsea_args},
      int_res_ls        = '{params.zoom_int_res_ls}',
      custom_mkr_names  = '{params.zoom_custom_mkr_names}',
      custom_mkr_paths  = '{params.zoom_custom_mkr_paths}',
      mkr_sel_res       = '{params.zoom_mkr_sel_res}',
      mkr_not_ok_re     = '{params.zoom_mkr_not_ok_re}',
      mkr_min_cpm_mkr   =  {params.zoom_mkr_min_cpm_mkr},
      mkr_min_cells     =  {params.zoom_mkr_min_cells},
      mkr_gsea_var      = '{params.zoom_mkr_gsea_var}',
      mkr_gsea_cut      =  {params.zoom_mkr_gsea_cut},
      ref_txome         = '{params.genome_ref}',
      batch_var         = '{params.batch_var}',
      do_gsea           = '{params.zoom_mkr_do_gsea}',
      do_xgboost        = '{params.do_xgboost}',
      {params.xgb_args}
    )"
    """
