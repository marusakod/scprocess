# join.smk — standalone Snakemake workflow for scprocess join
#
# Integrates outputs from multiple completed scprocess projects into a single
# joint analysis (HVG ranking → matrix assembly → integration → marker genes →
# optional GSEA). Called via `scprocess join join.yaml`.
#
# Usage:
#   scprocess join join.yaml
#   scprocess join join.yaml -n            # dry run
#   scprocess join join.yaml --unlock

import os
import sys
import pathlib
import yaml
import polars as pl

sys.path.append('scripts')
from scprocess_utils import (check_join_config, get_resources, prep_resource_params,
  get_join_project_parameters, get_join_source_labels_f, get_join_batch_sources)

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

scprocess_dir = pathlib.Path(config.pop('scprocess_dir'))
scdata_dir    = pathlib.Path(os.getenv('SCPROCESS_DATA_DIR'))
join_schema_f = scprocess_dir / "resources/schemas/join.schema.json"

# validate config, apply defaults, and check shiny parameters
config = check_join_config(config, join_schema_f, scdata_dir)

# resource parameters (no ML model for join — uses schema defaults + user overrides)
lm_f            = scprocess_dir / "resources/snakemake/resources_lm_params_2025-12-16.csv"
RESOURCE_PARAMS = prep_resource_params(config, join_schema_f, lm_f)

# ---------------------------------------------------------------------------
# Project and derived parameters
# ---------------------------------------------------------------------------

JOIN_PROJ = get_join_project_parameters(config)

# unpack project parameters
JOIN_PROJECT_IDS  = JOIN_PROJ['project_ids']
PROJECT_CFGS      = JOIN_PROJ['project_cfgs']
VAR_STATS_FS      = JOIN_PROJ['var_stats_fs']
H5ADS_YAML_FS     = JOIN_PROJ['h5ads_yaml_fs']
INTEGRATED_FS     = JOIN_PROJ['integrated_fs']
SAMPLE_META_FS    = JOIN_PROJ['sample_meta_fs']
JOIN_BATCH_KEYS   = JOIN_PROJ['batch_keys']

# join metadata and paths
JOIN_NAME     = config['join']['name']
JOIN_DIR      = pathlib.Path(config['join']['proj_dir'])
JOIN_TAG      = f"{JOIN_NAME}_join"
DATE_STAMP    = config['join']['date_stamp']
REF_TXOME     = config['join']['ref_txome']
join_int_dir  = str(JOIN_DIR / f"output/{JOIN_TAG}")
join_mkr_dir  = join_int_dir
logs_dir      = str(JOIN_DIR / ".log/join")
benchmark_dir = str(JOIN_DIR / ".resources")

# integration
INT_N_DIMS     = config['integration']['int_n_dims']
INT_RES_LS     = config['integration']['int_res_ls']
INT_PCA_METHOD = config['integration']['int_pca_method']

# marker genes
MKR_SEL_RES = config['marker_genes']['mkr_sel_res']
N_HVGS      = config['hvg']['hvg_n_hvgs']
GSEA_TXOMES = ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020']
DO_GSEA     = config['marker_genes']['mkr_do_gsea'] and (REF_TXOME in GSEA_TXOMES)

# label_celltypes (optional)
LABELLER_PARAMS = config.get('label_celltypes', [])
DO_LABEL = len(LABELLER_PARAMS) > 0

# ---------------------------------------------------------------------------
# Output file paths
# ---------------------------------------------------------------------------

joint_hvgs_f        = f"{join_int_dir}/joint_hvgs_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
joint_counts_f      = f"{join_int_dir}/joint_counts_hvgs_{JOIN_TAG}_{DATE_STAMP}.h5"
joint_coldata_f     = f"{join_int_dir}/joint_coldata_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
joint_sample_meta_f = f"{join_int_dir}/joint_sample_meta_{JOIN_TAG}_{DATE_STAMP}.csv"
joint_pca_f         = f"{join_int_dir}/joint_pca_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
joint_integration_f = f"{join_int_dir}/integrated_dt_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
joint_h5ads_yaml_f  = f"{join_int_dir}/h5ads_clean_paths_{JOIN_TAG}_{DATE_STAMP}.yaml"
h5ads_dir           = f"{join_int_dir}/h5ads"

pb_f        = f"{join_mkr_dir}/pb_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.rds"
mkrs_f      = f"{join_mkr_dir}/pb_marker_genes_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv.gz"
pb_hvgs_f   = f"{join_mkr_dir}/pb_hvgs_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv.gz"
fgsea_bp_f  = f"{join_mkr_dir}/fgsea_{JOIN_TAG}_{MKR_SEL_RES}_go_bp_{DATE_STAMP}.csv.gz"
fgsea_cc_f  = f"{join_mkr_dir}/fgsea_{JOIN_TAG}_{MKR_SEL_RES}_go_cc_{DATE_STAMP}.csv.gz"
fgsea_mf_f  = f"{join_mkr_dir}/fgsea_{JOIN_TAG}_{MKR_SEL_RES}_go_mf_{DATE_STAMP}.csv.gz"

join_lbl_dir = join_int_dir
label_fs = [
  f"{join_lbl_dir}/labels_{e['labeller']}_model_{e['model']}_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
  for e in LABELLER_PARAMS
]
_names_entries = [e for e in LABELLER_PARAMS if e.get('save_cluster_names_file', False)]
cluster_names_fs = [
  f"{join_lbl_dir}/cluster_names_for_shiny_{e['labeller']}_{e['model']}_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv"
  for e in _names_entries
]

docs_dir  = str(JOIN_DIR / "public")
rmd_dir   = str(JOIN_DIR / "analysis")
code_dir  = str(JOIN_DIR / "code")
html_f    = f"{docs_dir}/{JOIN_TAG}.html"
rmd_f     = f"{rmd_dir}/{JOIN_TAG}.Rmd"

# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------

rule all:
  input:
    joint_integration_f,
    joint_h5ads_yaml_f,
    mkrs_f,
    pb_hvgs_f,
    *([fgsea_bp_f, fgsea_cc_f, fgsea_mf_f] if DO_GSEA else []),
    *label_fs,
    *cluster_names_fs,
    html_f


rule join_select_hvgs:
  """Select joint HVGs using mean-rank aggregation across projects."""
  input:
    var_stats_fs = VAR_STATS_FS
  output:
    joint_hvgs_f = joint_hvgs_f
  params:
    project_ids  = " ".join(JOIN_PROJECT_IDS),
    n_hvgs       = N_HVGS
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_select_hvgs', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_select_hvgs', 'time', attempt)
  log:
    f"{logs_dir}/join_select_hvgs_{JOIN_TAG}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_select_hvgs_{JOIN_TAG}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    mkdir -p {join_int_dir}
    python3 scripts/join.py select_joint_hvgs \
      --var_stats_fs  {input.var_stats_fs} \
      --project_ids   {params.project_ids} \
      --n_hvgs        {params.n_hvgs} \
      --out_f         {output.joint_hvgs_f}
    """


rule join_build_matrix:
  """Assemble joint HVG count matrix from per-project h5ads."""
  input:
    joint_hvgs_f  = joint_hvgs_f,
    h5ads_yaml_fs = H5ADS_YAML_FS,
    integrated_fs = INTEGRATED_FS
  output:
    joint_counts_f = joint_counts_f
  params:
    project_ids   = " ".join(JOIN_PROJECT_IDS)
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_build_matrix', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_build_matrix', 'time', attempt)
  log:
    f"{logs_dir}/join_build_matrix_{JOIN_TAG}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_build_matrix_{JOIN_TAG}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 scripts/join.py build_joint_matrix \
      --joint_hvgs_f        {input.joint_hvgs_f} \
      --h5ads_yaml_fs       {input.h5ads_yaml_fs} \
      --project_ids         {params.project_ids} \
      --integrated_dt_fs    {input.integrated_fs} \
      --out_h5_f            {output.joint_counts_f}
    """


rule join_build_coldata:
  """Build joint coldata and sample metadata from matrix barcodes."""
  input:
    joint_counts_f = joint_counts_f,
    integrated_fs  = INTEGRATED_FS,
    sample_meta_fs = SAMPLE_META_FS
  output:
    joint_coldata_f     = joint_coldata_f,
    joint_sample_meta_f = joint_sample_meta_f
  params:
    project_ids = " ".join(JOIN_PROJECT_IDS)
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_build_coldata', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_build_coldata', 'time', attempt)
  log:
    f"{logs_dir}/join_build_coldata_{JOIN_TAG}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_build_coldata_{JOIN_TAG}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 scripts/join.py build_joint_coldata \
      --h5_f                {input.joint_counts_f} \
      --project_ids         {params.project_ids} \
      --integrated_dt_fs    {input.integrated_fs} \
      --sample_meta_fs      {input.sample_meta_fs} \
      --out_coldata_f       {output.joint_coldata_f} \
      --out_sample_meta_f   {output.joint_sample_meta_f}
    """


rule join_pca:
  """Compute PCA on joint matrix using BPCells disk-backed streaming SVD."""
  input:
    counts_h5_f = joint_counts_f
  output:
    pca_f = joint_pca_f
  params:
    n_dims = INT_N_DIMS
  threads: 8
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_pca', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_pca', 'time', attempt)
  log:
    f"{logs_dir}/join_pca_{JOIN_TAG}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_pca_{JOIN_TAG}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/bpcells_pca.yaml'
  shell: """
    exec &>> {log}
    export OPENBLAS_NUM_THREADS={threads}
    export MKL_NUM_THREADS={threads}
    export OMP_NUM_THREADS={threads}

    # Copy HDF5 to local disk for fast repeated reads during SVD
    LOCAL_H5=$(mktemp /tmp/join_pca_XXXXXX.h5)
    trap "rm -f $LOCAL_H5" EXIT
    cp {input.counts_h5_f} $LOCAL_H5

    Rscript -e "source('scripts/join_pca.R'); run_join_pca(
      counts_h5_f = '$LOCAL_H5',
      n_dims      =  {params.n_dims},
      out_pca_f   = '{output.pca_f}')"
    """


rule join_integration:
  """Run Harmony integration on the joint HVG matrix."""
  input:
    hvg_mat_f    = joint_counts_f,
    coldata_f    = joint_coldata_f,
    sample_qc_f  = joint_sample_meta_f,
    pca_f        = joint_pca_f if INT_PCA_METHOD == 'bpcells' else []
  output:
    integration_f = joint_integration_f
  params:
    embedding         = config['integration']['int_embedding'],
    n_dims            = INT_N_DIMS,
    cl_method         = config['integration']['int_cl_method'],
    theta_concat      = config['integration']['int_theta'],
    batch_var_concat  = config['integration']['int_batch_var'],
    res_ls_concat     = config['integration']['int_res_ls'],
    use_paga          = config['integration']['int_use_paga'],
    paga_cl_res       = config['integration']['int_paga_cl_res'],
    int_use_gpu       = config['integration']['int_use_gpu'],
    pca_method        = INT_PCA_METHOD
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_integration', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_integration', 'time', attempt)
  log:
    f"{logs_dir}/join_integration_{JOIN_TAG}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_integration_{JOIN_TAG}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/integration.yaml'
  shell: """
    exec &>> {log}

    set +u
    USE_GPU_FLAG=""
    if [ "{params.int_use_gpu}" == "True" ]; then
      if [ -n "$CUDA_VISIBLE_DEVICES" ]; then
        echo "running on GPU"
        USE_GPU_FLAG="--use-gpu"
      else
        echo "GPU requested but no GPU available, running on CPU"
      fi
    fi
    set -u

    PCA_FLAG=""
    if [ "{params.pca_method}" == "bpcells" ]; then
      PCA_FLAG="--precomputed_pca_f {input.pca_f}"
    fi

    python3 scripts/integration.py run_zoom_integration \
      --hvg_mat_f        {input.hvg_mat_f} \
      --sample_qc_f      {input.sample_qc_f} \
      --coldata_f        {input.coldata_f} \
      --demux_type       none \
      --exclude_mito     False \
      --embedding        {params.embedding} \
      --n_dims           {params.n_dims} \
      --cl_method        {params.cl_method} \
      --theta_concat     "{params.theta_concat}" \
      --batch_var_concat "{params.batch_var_concat}" \
      --res_ls_concat    "{params.res_ls_concat}" \
      --integration_f    {output.integration_f} \
      $(if [ "{params.use_paga}" == "True" ]; then echo "--use-paga"; fi) \
      $(if [ "{params.use_paga}" == "True" ]; then echo "--paga-cl-res {params.paga_cl_res}"; fi) \
      $USE_GPU_FLAG \
      $PCA_FLAG
    """


rule join_build_h5ads_yaml:
  """Create symlinks and the joint h5ads YAML manifest."""
  input:
    h5ads_yaml_fs = H5ADS_YAML_FS
  output:
    joint_h5ads_yaml_f = joint_h5ads_yaml_f,
    h5ad_symlinks = [f"{h5ads_dir}/{bk}.h5ad" for bk in JOIN_BATCH_KEYS]
  params:
    project_ids = " ".join(JOIN_PROJECT_IDS),
    h5ads_dir   = h5ads_dir
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_build_h5ads_yaml', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_build_h5ads_yaml', 'time', attempt)
  log:
    f"{logs_dir}/join_build_h5ads_yaml_{JOIN_TAG}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_build_h5ads_yaml_{JOIN_TAG}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 scripts/join.py build_join_h5ads_yaml \
      --h5ads_yaml_fs  {input.h5ads_yaml_fs} \
      --project_ids    {params.project_ids} \
      --h5ads_dir      {params.h5ads_dir} \
      --out_yaml_f     {output.joint_h5ads_yaml_f}
    """


rule join_marker_genes:
  """Pseudobulk marker gene detection on the joint integration."""
  input:
    h5ads_yaml_f  = joint_h5ads_yaml_f,
    integration_f = joint_integration_f
  output:
    pb_f      = pb_f,
    mkrs_f    = mkrs_f,
    pb_hvgs_f = pb_hvgs_f
  params:
    gene_info_f = config['marker_genes']['gene_info_f'],
    sel_res     = MKR_SEL_RES,
    min_cl_size = config['marker_genes']['mkr_min_cl_size'],
    min_cells   = config['marker_genes']['mkr_min_cells'],
    batch_var   = "sample_id"
  threads: 8
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_marker_genes', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_marker_genes', 'time', attempt)
  log:
    f"{logs_dir}/join_marker_genes_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_marker_genes_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/rlibs_bpcells.yaml'
  shell: """
    exec &>> {log}
    Rscript -e "source('scripts/utils.R'); source('scripts/marker_genes.R'); calculate_marker_genes(
      integration_f = '{input.integration_f}',
      h5ads_yaml_f  = '{input.h5ads_yaml_f}',
      pb_f          = '{output.pb_f}',
      mkrs_f        = '{output.mkrs_f}',
      pb_hvgs_f     = '{output.pb_hvgs_f}',
      gtf_dt_f      = '{params.gene_info_f}',
      sel_res       = '{params.sel_res}',
      min_cl_size   =  {params.min_cl_size},
      min_cells     =  {params.min_cells},
      zoom          = 'True',
      batch_var     = '{params.batch_var}',
      n_cores       =  {threads},
      use_bpcells   = 'True')"
    """


rule join_fgsea:
  """GSEA on join marker genes (runs only for supported transcriptomes)."""
  input:
    mkrs_f = mkrs_f
  output:
    fgsea_go_bp_f = fgsea_bp_f,
    fgsea_go_cc_f = fgsea_cc_f,
    fgsea_go_mf_f = fgsea_mf_f
  params:
    ref_txome   = REF_TXOME,
    gsea_dir    = config['marker_genes']['mkr_gsea_dir'],
    min_cpm_go  = config['marker_genes']['mkr_min_cpm_go'],
    max_zero_p  = config['marker_genes']['mkr_max_zero_p'],
    gsea_cut    = config['marker_genes']['mkr_gsea_cut'],
    not_ok_re   = config['marker_genes']['mkr_not_ok_re'],
    gsea_var    = config['marker_genes']['mkr_gsea_var']
  threads: 8
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_fgsea', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_fgsea', 'time', attempt)
  log:
    f"{logs_dir}/join_fgsea_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_fgsea_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/rlibs.yaml'
  shell: """
    exec &>> {log}
    Rscript -e "source('scripts/utils.R'); source('scripts/fgsea.R'); run_fgsea(
      mkrs_f        = '{input.mkrs_f}',
      fgsea_go_bp_f = '{output.fgsea_go_bp_f}',
      fgsea_go_cc_f = '{output.fgsea_go_cc_f}',
      fgsea_go_mf_f = '{output.fgsea_go_mf_f}',
      genome_ref    = '{params.ref_txome}',
      gsea_dir      = '{params.gsea_dir}',
      min_cpm_go    = {params.min_cpm_go},
      max_zero_p    = {params.max_zero_p},
      gsea_cut      = {params.gsea_cut},
      not_ok_re     = '{params.not_ok_re}',
      gsea_var      = '{params.gsea_var}',
      n_cores       =  {threads})"
    """


wildcard_constraints:
  labeller = "celltypist|scprocess|custom"

def _parse_join_label_params(labeller, model):
  matches = [e for e in LABELLER_PARAMS
    if e['labeller'] == labeller and e['model'] == model]
  if len(matches) != 1:
    raise ValueError(f"Expected 1 match for {labeller}/{model}, got {len(matches)}")
  return matches[0]

def _join_merge_labels_inputs(wildcards):
  fresh_batches, reuse_label_fs    get_join_batch_sources(
    PROJECT_CFGS,     JN_PROJECT_IDS, JOIN_BATCH_KEYS, wildcards.labeller, wildcards.model)
  pred_fs = expand(
    f"{join_lbl_dir}/tmp_labels_{{labeller}}_model_{{model}}_{JOIN_TAG}_{DATE_STAMP}_{{batch}}.csv.gz",
    batch=fresh_batches, allow_missing=True)
  if reuse_label_fs:
    pred_fs.append(f" {join_lbl_dir}/tmp_labels_{wildcards.labeller}_model_{wildcards.model}_{JOIN_TAG}_{DATE_STAMP}_reused.csv.gz")
  return pred_fs


rule join_celltypist:
  """Run CellTypist on each batch h5ad for the join integration."""
  input:
    adata_f = f"{join_int_dir}/h5ads/{{batch}}.h5ad"
  output:
    pred_f  = temp(f"{join_lbl_dir}/tmp_labels_celltypist_model_{{model}}_{JOIN_TAG}_{DATE_STAMP}_{{batch}}.csv.gz")
  threads: 4
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_celltypist', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_celltypist', 'time', attempt)
  log:
    f"{logs_dir}/join_celltypist_{{model}}_{{batch}}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_celltypist_{{model}}_{{batch}}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/celltypist.yaml'
  shell: """
    exec &>> {log}
    python3 scripts/label_celltypes.py celltypist_one_batch \
      {wildcards.batch} sample_id {wildcards.model} \
      --adata_f   {input.adata_f} \
      --pred_f    {output.pred_f}
    """


rule join_scprocess_labeller:
  """Run scprocess XGBoost labeller on each batch h5ad."""
  input:
    adata_f = f"{join_int_dir}/h5ads/{{batch}}.h5ad"
  output:
    pred_f  = temp(f"{join_lbl_dir}/tmp_labels_scprocess_model_{{model}}_{JOIN_TAG}_{DATE_STAMP}_{{batch}}.csv.gz")
  params:
    xgb_f     = lambda wildcards: _parse_join_label_params('scprocess', wildcards.model)['xgb_f'],
    xgb_cls_f = lambda wildcards: _parse_join_label_params('scprocess', wildcards.model)['xgb_cls_f']
  threads: 1
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_scprocess_labeller', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_scprocess_labeller', 'time', attempt)
  log:
    f"{logs_dir}/join_scprocess_labeller_{{model}}_{{batch}}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_scprocess_labeller_{{model}}_{{batch}}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/rlibs.yaml'
  shell: """
    exec &>> {log}
    Rscript -e "source('scripts/label_celltypes.R'); source('scripts/integration.R'); \\
    label_with_xgboost_one_batch(
      sel_batch   = '{wildcards.batch}',
      batch_var   = 'sample_id',
      model_name  = '{wildcards.model}',
      xgb_f       = '{params.xgb_f}',
      xgb_cls_f   = '{params.xgb_cls_f}',
      adata_f     = '{input.adata_f}',
      pred_f      = '{output.pred_f}'
    )"
    """


rule join_extract_labels:
  """Extract naive predictions from source project labels for cells in the join."""
  input:
    source_labels_fs   lambda wildcards: get_join_batch_sources(PROJECT_CFGS, JN_PROJECT_IDS, JOIN_BATCH_KEYS, wildcards.labeller, wildcards.model)[1],
        integration_f    = joint_integration_f
  output:
    pred_f = temp(f"{join_lbl_dir}/tmp_labels_{{labeller}}_model_{{model}}_{JOIN_TAG}_{DATE_STAMP}_reused.csv.gz")
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_extract_labels', 'memory', attempt),
    run time = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_extract_labels', 'time', attempt)
  log:
    f"{logs_dir}/join_extract_labels_{{labeller}}_{{model}}_{DATE_STAMP}.log"
  conda:
    '../envs/scprocess_local.yaml'
  shell: """
    exec &>> {log}
    python3 scripts/label_celltypes.py extract_naive_predictions \
      {input.source_labels_fs} \
      --int_f     {input.integration_f} \
      --model     {wildcards.model} \
      --pred_f    {output.pred_f}
    """


rule join_merge_labels:
  """Aggregate per-batch predictions by majority voting."""
  input:
    pred_fs       = _join_merge_labels_inputs,
    integration_f = joint_integration_f
  output:
    pred_out_f    = f"{join_lbl_dir}/labels_{{labeller}}_model_{{model}}_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
  params:
    hi_res_cl   = lambda wildcards: _parse_join_label_params(wildcards.labeller, wildcards.model)['hi_res_cl'],
    min_cl_prop = lambda wildcards: _parse_join_label_params(wildcards.labeller, wildcards.model)['min_cl_prop']
  threads: 4
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_merge_labels', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_merge_labels', 'time', attempt)
  log:
    f"{logs_dir}/join_merge_labels_{{labeller}}_{{model}}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_merge_labels_{{labeller}}_{{model}}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/celltypist.yaml'
  shell: """
    exec &>> {log}
    python3 scripts/label_celltypes.py aggregate_predictions \
      {input.pred_fs} \
      --int_f           {input.integration_f} \
      --hi_res_cl       {params.hi_res_cl} \
      --min_cl_prop     {params.min_cl_prop} \
      --batch_var       sample_id \
      --agg_f           {output.pred_out_f}
    """


rule join_save_cluster_names:
  """Generate cluster_name CSV from aggregated labels at the marker gene resolution."""
  input:
    labels_f      = f"{join_lbl_dir}/labels_{{labeller}}_model_{{model}}_{JOIN_TAG}_{DATE_STAMP}.csv.gz",
    integration_f = joint_integration_f
  output:
    names_f       = f"{join_lbl_dir}/cluster_names_for_shiny_{{labeller}}_{{model}}_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv"
  params:
    mkr_sel_res   = MKR_SEL_RES
  log:
    f"{logs_dir}/join_save_cluster_names_{{labeller}}_{{model}}_{DATE_STAMP}.log"
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


rule join_render_html:
  """Render Rmd report and HTML for the join integration."""
  input:
    integration_f = joint_integration_f,
    mkrs_f        = mkrs_f,
    pb_hvgs_f     = pb_hvgs_f,
    pb_f          = pb_f,
    fgsea_files   = [fgsea_bp_f, fgsea_cc_f, fgsea_mf_f] if DO_GSEA else [],
    label_files   = label_fs,
    cluster_names = cluster_names_fs
  output:
    r_utils_f     = f"{code_dir}/utils.R",
    r_int_f       = f"{code_dir}/integration.R",
    r_mkr_f       = f"{code_dir}/marker_genes.R",
    r_fgsea_f     = f"{code_dir}/fgsea.R",
    r_lbl_f       = f"{code_dir}/label_celltypes.R",
    rmd_f         = rmd_f,
    html_f        = html_f
  params:
    your_name     = config['join']['your_name'],
    affiliation   = config['join']['affiliation'],
    join_name     = JOIN_NAME,
    join_tag      = JOIN_TAG,
    join_int_dir  = join_int_dir,
    join_mkr_dir  = join_mkr_dir,
    ref_txome     = REF_TXOME,
    mkr_sel_res   = MKR_SEL_RES,
    int_res_ls    = config['integration']['int_res_ls'],
    metadata_vars    = " ".join(config['join'].get('metadata_vars', [])),
    custom_mkr_names = config['marker_genes']['custom_mkr_names'],
    custom_mkr_paths = config['marker_genes']['custom_mkr_paths'],
    label_f_ls       = ' '.join(label_fs),
    labeller_ls      = ' '.join(e['labeller']  for e in LABELLER_PARAMS) if DO_LABEL else '',
    model_ls         = ' '.join(e['model']     for e in LABELLER_PARAMS) if DO_LABEL else '',
    hi_res_cl_ls     = ' '.join(e['hi_res_cl'] for e in LABELLER_PARAMS) if DO_LABEL else '',
    min_cl_prop_ls   = ' '.join(str(e['min_cl_prop']) for e in LABELLER_PARAMS) if DO_LABEL else '',
    mkr_min_cpm_mkr  = config['marker_genes']['mkr_min_cpm_mkr'],
    mkr_min_cells    = config['marker_genes']['mkr_min_cells'],
    mkr_gsea_cut     = config['marker_genes']['mkr_gsea_cut'],
    integration_f    = joint_integration_f,
    sample_meta_f    = joint_sample_meta_f,
    mkrs_f           = mkrs_f,
    pb_hvgs_f        = pb_hvgs_f,
    pb_f             = pb_f,
    fgsea_go_bp_f    = fgsea_bp_f if DO_GSEA else '',
    fgsea_go_cc_f    = fgsea_cc_f if DO_GSEA else '',
    fgsea_go_mf_f    = fgsea_mf_f if DO_GSEA else '',
    proj_dir         = str(JOIN_DIR),
    scprocess_dir    = str(scprocess_dir),
    date_stamp       = DATE_STAMP
  threads: 1
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_render_html', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_render_html', 'time', attempt)
  log:
    f"{logs_dir}/join_render_html_{JOIN_TAG}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_render_html_{JOIN_TAG}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/rlibs.yaml'
  shell: """
    exec &>> {log}

    # copy R code over
    echo "copying relevant R files over"
    cp scripts/utils.R        {output.r_utils_f}
    cp scripts/integration.R  {output.r_int_f}
    cp scripts/marker_genes.R {output.r_mkr_f}
    cp scripts/fgsea.R        {output.r_fgsea_f}
    cp scripts/label_celltypes.R {output.r_lbl_f}

    # define rule and template
    template_f=$(realpath resources/rmd_templates/join.Rmd.template)
    rule="join"

    # rendering html
    Rscript --vanilla -e "source('scripts/render_htmls.R'); \\
    render_html(
      rule_name     = '$rule',
      temp_f        = '$template_f',
      rmd_f         = '{output.rmd_f}',
      proj_dir      = '{params.proj_dir}',
      your_name     = '{params.your_name}',
      affiliation   = '{params.affiliation}',
      join_name     = '{params.join_name}',
      join_tag      = '{params.join_tag}',
      join_int_dir  = '{params.join_int_dir}',
      join_mkr_dir  = '{params.join_mkr_dir}',
      ref_txome     = '{params.ref_txome}',
      mkr_sel_res   =  {params.mkr_sel_res},
      int_res_ls    = '{params.int_res_ls}',
      metadata_vars    = '{params.metadata_vars}',
      custom_mkr_names = '{params.custom_mkr_names}',
      custom_mkr_paths = '{params.custom_mkr_paths}',
      label_f_ls       = '{params.label_f_ls}',
      labeller_ls      = '{params.labeller_ls}',
      model_ls         = '{params.model_ls}',
      hi_res_cl_ls     = '{params.hi_res_cl_ls}',
      min_cl_prop_ls   = '{params.min_cl_prop_ls}',
      mkr_min_cpm_mkr  =  {params.mkr_min_cpm_mkr},
      mkr_min_cells    =  {params.mkr_min_cells},
      mkr_gsea_cut     =  {params.mkr_gsea_cut},
      integration_f    = '{params.integration_f}',
      sample_meta_f    = '{params.sample_meta_f}',
      mkrs_f           = '{params.mkrs_f}',
      pb_hvgs_f        = '{params.pb_hvgs_f}',
      pb_f             = '{params.pb_f}',
      fgsea_go_bp_f    = '{params.fgsea_go_bp_f}',
      fgsea_go_cc_f    = '{params.fgsea_go_cc_f}',
      fgsea_go_mf_f    = '{params.fgsea_go_mf_f}',
      scprocess_dir    = '{params.scprocess_dir}',
      date_stamp       = '{params.date_stamp}'
    )"
    """
