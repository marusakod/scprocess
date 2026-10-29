# join.smk — standalone Snakemake workflow for scprocess join
#
# Integrates outputs from multiple completed scprocess projects into a single
# joint analysis (HVG ranking → matrix assembly → integration → marker genes →
# optional GSEA). Called via `scprocess join join.yaml`.


import os
import sys
import pathlib
import yaml
import polars as pl

scprocess_dir = pathlib.Path(config.pop('scprocess_dir'))
sys.path.append(str(scprocess_dir / 'scripts'))
from scprocess_utils import (check_join_config, get_resources, prep_resource_params,
  get_join_project_parameters, get_join_source_labels_f, get_join_batch_sources,
  get_paga_resolution, get_xgboost_spec_inputs, get_join_labeller_parameters)

# Prevent underscore-containing model names from being partially captured as the
# labeller in join label, log, benchmark, and cluster-name file patterns.
wildcard_constraints:
  labeller = "celltypist|scprocess"

# setup
scdata_dir    = pathlib.Path(os.getenv('SCPROCESS_DATA_DIR'))
join_schema_f = scprocess_dir / "resources/schemas/join.schema.json"

# validate config, apply defaults, and check shiny parameters
JOIN_ZOOM_PARAMS = globals().get('JOIN_ZOOM_PARAMS', {})
config = config if JOIN_ZOOM_PARAMS else check_join_config(config, join_schema_f, scdata_dir)

# Track these options only for zooms. Even an empty-string parameter would
# invalidate completed standard joins when upgrading to the zoom-capable rules.
JOIN_REQUIRE_ALL_CELLS = '--require-all-cells' if JOIN_ZOOM_PARAMS else ''
JOIN_RETAIN_METADATA = '--retain-metadata' if JOIN_ZOOM_PARAMS else ''

# resource parameters (no ML model for join — uses schema defaults + user overrides)
RESOURCE_PARAMS = prep_resource_params(config, join_schema_f, scprocess_dir)

# ---------------------------------------------------------------------------
# Project and derived parameters
# ---------------------------------------------------------------------------

JOIN_PROJ = JOIN_ZOOM_PARAMS.get('project') or get_join_project_parameters(config)

# unpack project parameters
JOIN_PROJECT_IDS  = JOIN_PROJ['project_ids']
PROJECT_CFGS      = JOIN_PROJ['project_cfgs']
HVG_FS            = JOIN_PROJ['hvg_fs']
H5ADS_YAML_FS     = JOIN_PROJ['h5ads_yaml_fs']
INTEGRATED_FS     = JOIN_PROJ['integrated_fs']
SAMPLE_META_FS    = JOIN_PROJ['sample_meta_fs']
QC_ALL_FS         = JOIN_PROJ['qc_all_fs']
# Prefer merged real-cell diagnostics; legacy projects may have only per-run scores.
import shlex
JOIN_DOUBLET_INPUTS = []
if not JOIN_ZOOM_PARAMS:
  for _pid in JOIN_PROJECT_IDS:
    _pcfg = PROJECT_CFGS[_pid]['project']
    _ddir = pathlib.Path(_pcfg['proj_dir']) / 'output' / f"{_pcfg['short_tag']}_doublet_id"
    _merged = _ddir / f"scDblFinder_real_all_runs_{_pcfg['full_tag']}_{_pcfg['date_stamp']}.csv.gz"
    _paths = [_merged] if _merged.is_file() else sorted(_ddir.glob(
      f"dbl_*/scDblFinder_*_outputs_{_pcfg['full_tag']}_{_pcfg['date_stamp']}.csv.gz"))
    JOIN_DOUBLET_INPUTS.extend((_pid, str(path)) for path in _paths)
JOIN_DOUBLET_ARGS = ('--doublet_project_ids ' + shlex.join([pid for pid, path in JOIN_DOUBLET_INPUTS]) +
                    ' --doublet_fs ' + shlex.join([path for pid, path in JOIN_DOUBLET_INPUTS])) if JOIN_DOUBLET_INPUTS else ''

JOIN_BATCH_KEYS   = JOIN_PROJ['batch_keys']

# join metadata and paths
JOIN_NAME     = config['join']['name']
JOIN_DIR      = pathlib.Path(config['join']['proj_dir'])
JOIN_TAG      = JOIN_ZOOM_PARAMS.get('tag', f"{JOIN_NAME}_join")
DATE_STAMP    = config['join']['date_stamp']
GENOME_REF    = config['join'].get('probe_set', config['join'].get('ref_txome', ''))
join_int_dir  = JOIN_ZOOM_PARAMS.get('out_dir', str(JOIN_DIR / f"output/{JOIN_TAG}"))
join_mkr_dir  = join_int_dir
logs_dir      = JOIN_ZOOM_PARAMS.get('logs_dir', str(JOIN_DIR / ".log/join"))
benchmark_dir = JOIN_ZOOM_PARAMS.get('benchmark_dir', str(JOIN_DIR / ".resources"))
marker_edger_bp_script = str(scprocess_dir / "scripts/marker_genes_edger_bp.R")

# integration
INT_N_DIMS     = config['integration']['int_n_dims']
INT_RES_LS     = config['integration']['int_res_ls']
INT_PCA_METHOD = config['integration']['int_pca_method']
INT_BATCH_VAR  = config['integration']['int_batch_var']
INT_CYCLE_MODE = config['integration'].get('int_cell_cycle_correction')
DO_CYCLE_CORRECTION = INT_CYCLE_MODE is not None

# marker genes
_mkr_cfg    = config['marker_genes']
MKR_SEL_RES = _mkr_cfg['mkr_sel_res']
N_HVGS      = config['hvg']['hvg_n_hvgs']
GENE_INFO_F = _mkr_cfg['gene_info_f']
GSEA_DIR    = _mkr_cfg['mkr_gsea_dir']
MKR_CUSTOM_NAMES = _mkr_cfg['custom_mkr_names']
MKR_CUSTOM_PATHS = _mkr_cfg['custom_mkr_paths']

GSEA_REFS   = ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020',
               'human_v1', 'human_v2', 'mouse_v1', 'mouse_v2']
DO_GSEA     = _mkr_cfg['mkr_do_gsea'] and (GENOME_REF in GSEA_REFS)

# Load optional annotation settings and validate installed model metadata.
LABELLER_PARAMS = get_join_labeller_parameters(config, join_schema_f, scdata_dir)
DO_LABEL        = bool(LABELLER_PARAMS)

# output paths
joint_hvgs_f        = f"{join_int_dir}/joint_hvgs_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
joint_label_counts_f      = f"{join_int_dir}/joint_counts_hvgs_{JOIN_TAG}_{DATE_STAMP}.h5"
joint_bpcells_counts_dir  = f"{join_int_dir}/joint_counts_hvgs_{JOIN_TAG}_{DATE_STAMP}.bpcells"
joint_coldata_f     = f"{join_int_dir}/joint_coldata_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
joint_sample_meta_f = f"{join_int_dir}/joint_sample_meta_{JOIN_TAG}_{DATE_STAMP}.csv"
joint_qc_f          = f"{join_int_dir}/joint_qc_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
joint_pca_f         = f"{join_int_dir}/joint_pca_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
joint_cycle_regression_f = f"{join_int_dir}/cell_cycle_regression_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
joint_cycle_genes_f = f"{join_int_dir}/cell_cycle_regression_genes_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
JOIN_CYCLE_ARGS = (
  f", tricycle_f = '{joint_coldata_f}', "
  f"correction_centering = '{INT_CYCLE_MODE}', "
  f"harmonics = {config['integration']['int_cell_cycle_harmonics']}, "
  f"ridge_lambda = {config['integration']['int_cell_cycle_ridge_lambda']}, "
  f"out_regression_f = '{joint_cycle_regression_f}', "
  f"out_regression_genes_f = '{joint_cycle_genes_f}'"
) if DO_CYCLE_CORRECTION else ''
joint_integration_f = f"{join_int_dir}/integrated_dt_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
joint_h5ads_yaml_f  = f"{join_int_dir}/h5ads_clean_paths_{JOIN_TAG}_{DATE_STAMP}.yaml"
h5ads_dir           = f"{join_int_dir}/h5ads"

pb_dir      = f"{join_mkr_dir}/pb_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.bpcells"
pb_prepared_coldata_f = f"{join_mkr_dir}/pb_prepared_coldata_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv.gz"
pb_prepared_rowdata_f = f"{join_mkr_dir}/pb_prepared_rowdata_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv.gz"
mkrs_f      = f"{join_mkr_dir}/pb_marker_genes_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv.gz"
pb_hvgs_f   = f"{join_mkr_dir}/pb_hvgs_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv.gz"
pb_plot_data_f = f"{join_mkr_dir}/pb_plot_data_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv.gz"
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
code_dir  = str(pathlib.Path(join_int_dir) / 'code') if JOIN_ZOOM_PARAMS else str(JOIN_DIR / "code")
html_f    = f"{docs_dir}/{JOIN_TAG}.html"
rmd_f     = f"{rmd_dir}/{JOIN_TAG}.Rmd"

YOUR_NAME   = config['join']['your_name']
AFFILIATION = config['join']['affiliation']

INT_RES_LS_STR = ' '.join(str(r) for r in INT_RES_LS)

# train_xgboost (optional)
DO_TRAIN_XGB = 'train_xgboost' in config
if DO_TRAIN_XGB:
  _join_xgb_cfg = config['train_xgboost']
  join_xgb_dir  = str(JOIN_DIR / f"output/{JOIN_NAME}_train_xgboost")
  _join_xgb_ref_tag = _join_xgb_cfg['xgboost_tag']
  _join_xgb_model_f    = f"{join_xgb_dir}/{_join_xgb_ref_tag}_xgboost_model.json"
  _join_xgb_fulldata_f = f"{join_xgb_dir}/{_join_xgb_ref_tag}_fulldata_predictions.csv.gz"

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


rule label_celltypes:
  input:
    label_fs,
    cluster_names_fs


if DO_TRAIN_XGB:
  rule train_xgboost:
    input:
      _join_xgb_model_f,
      _join_xgb_fulldata_f,
      f"{rmd_dir}/{JOIN_TAG}_train_xgboost.Rmd",
      f"{docs_dir}/{JOIN_TAG}_train_xgboost.html"


rule join_select_hvgs:
  input:
    hvg_fs = HVG_FS
  output:
    joint_hvgs_f = joint_hvgs_f
  params:
    project_ids  = "zoom" if JOIN_ZOOM_PARAMS else " ".join(JOIN_PROJECT_IDS),
    n_hvgs       = N_HVGS
  retries: config['resources']['retries']
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
    python3 {scprocess_dir}/scripts/join.py select_joint_hvgs \
      --hvg_fs       {input.hvg_fs} \
      --project_ids  {params.project_ids} \
      --n_hvgs       {params.n_hvgs} \
      --out_f        {output.joint_hvgs_f}
    """


rule join_build_matrix:
  """Assemble joint HVG count matrix from per-project h5ads."""
  input:
    joint_hvgs_f  = joint_hvgs_f,
    h5ads_yaml_fs = H5ADS_YAML_FS,
    integrated_fs = INTEGRATED_FS
  output:
    joint_label_counts_f = joint_label_counts_f
  params:
    project_ids   = " ".join(JOIN_PROJECT_IDS),
    **({'require_all_cells': JOIN_REQUIRE_ALL_CELLS} if JOIN_ZOOM_PARAMS else {})
  retries: config['resources']['retries']
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
    python3 {scprocess_dir}/scripts/join.py build_joint_matrix \
      --joint_hvgs_f        {input.joint_hvgs_f} \
      --h5ads_yaml_fs       {input.h5ads_yaml_fs} \
      --project_ids         {params.project_ids} \
      --integrated_dt_fs    {input.integrated_fs} \
      --out_h5_f            {output.joint_label_counts_f} {JOIN_REQUIRE_ALL_CELLS}
    """


rule join_build_coldata:
  """Build cell-level joint coldata from matrix barcodes and integrations."""
  input:
    joint_label_counts_f = joint_label_counts_f,
    integrated_fs = INTEGRATED_FS
  output:
    joint_coldata_f = joint_coldata_f
  params:
    project_ids = " ".join(JOIN_PROJECT_IDS)
  retries: config['resources']['retries']
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
    python3 {scprocess_dir}/scripts/join.py build_joint_coldata \
      --h5_f                {input.joint_label_counts_f} \
      --project_ids         {params.project_ids} \
      --integrated_dt_fs    {input.integrated_fs} \
      --out_coldata_f       {output.joint_coldata_f}
    """


rule join_build_sample_metadata:
  """Combine descriptive source sample metadata for reports and apps."""
  input:
    sample_meta_fs = SAMPLE_META_FS
  output:
    sample_meta_f = joint_sample_meta_f
  params:
    project_ids = " ".join(JOIN_PROJECT_IDS)
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_build_sample_metadata', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_build_sample_metadata', 'time', attempt)
  log:
    f"{logs_dir}/join_build_sample_metadata_{JOIN_TAG}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_build_sample_metadata_{JOIN_TAG}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 {scprocess_dir}/scripts/join.py build_joint_sample_metadata \
      --project_ids    {params.project_ids} \
      --sample_meta_fs {input.sample_meta_fs} \
      --out_f          {output.sample_meta_f}
    """


rule join_build_qc:
  """Combine source QC metrics for cells included in the joint analysis."""
  input:
    qc_fs     = QC_ALL_FS,
    doublet_fs = [path for pid, path in JOIN_DOUBLET_INPUTS],
    coldata_f = joint_coldata_f
  output:
    qc_f = joint_qc_f
  params:
    project_ids = " ".join(JOIN_PROJECT_IDS),
    doublet_args = JOIN_DOUBLET_ARGS
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_build_qc', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_build_qc', 'time', attempt)
  log:
    f"{logs_dir}/join_build_qc_{JOIN_TAG}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_build_qc_{JOIN_TAG}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 {scprocess_dir}/scripts/join.py build_joint_qc \
      --qc_fs       {input.qc_fs} \
      --project_ids {params.project_ids} \
      --coldata_f   {input.coldata_f} \
      --out_f       {output.qc_f} \
      {params.doublet_args}
    """


rule join_build_bpcells_matrix:
  """Convert joint HDF5 counts to persistent native BPCells storage."""
  input:
    counts_h5_f = joint_label_counts_f
  output:
    counts_bp_dir = directory(joint_bpcells_counts_dir)
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_build_bpcells_matrix', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_build_bpcells_matrix', 'time', attempt)
  log:
    f"{logs_dir}/join_build_bpcells_matrix_{JOIN_TAG}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_build_bpcells_matrix_{JOIN_TAG}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/bpcells_pca.yaml'
  shell: """
    exec &>> {log}
    export OPENBLAS_NUM_THREADS={threads}
    export MKL_NUM_THREADS={threads}
    export OMP_NUM_THREADS={threads}

    Rscript --vanilla -e "source('{scprocess_dir}/scripts/join_pca.R'); write_join_bpcells_matrix(
      counts_h5_f = '{input.counts_h5_f}',
      out_dir     = '{output.counts_bp_dir}')"
    """


rule join_pca:
  """Compute centered and scaled PCA from native BPCells storage."""
  input:
    counts_bp_dir = joint_bpcells_counts_dir,
    tricycle_f = joint_coldata_f if DO_CYCLE_CORRECTION else []
  output:
    pca_f = joint_pca_f,
    regression_f = joint_cycle_regression_f if DO_CYCLE_CORRECTION else [],
    regression_genes_f = joint_cycle_genes_f if DO_CYCLE_CORRECTION else []
  params:
    n_dims = INT_N_DIMS,
    # Only track correction parameters when enabled; preserve ordinary join caches.
    **({'cycle_args': JOIN_CYCLE_ARGS} if DO_CYCLE_CORRECTION else {})
  threads: 8
  retries: config['resources']['retries']
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

    # Copy compact native BPCells storage to local disk for repeated SVD reads.
    LOCAL_BP=$(mktemp -d /tmp/scprocess_join_pca_XXXXXX)
    trap 'rm -rf -- "$LOCAL_BP"' EXIT
    cp -a {input.counts_bp_dir}/. "$LOCAL_BP"/

    Rscript --vanilla -e "source('{scprocess_dir}/scripts/bpcells_pca.R'); source('{scprocess_dir}/scripts/join_pca.R'); run_join_pca(
      counts_bp_dir = '$LOCAL_BP',
      n_dims        =  {params.n_dims},
      n_threads     =  {threads},
      out_pca_f     = '{output.pca_f}'{JOIN_CYCLE_ARGS})"
    """


rule join_integration:
  input:
    hvg_mat_f    = joint_label_counts_f,
    coldata_f    = joint_coldata_f,
    pca_f        = joint_pca_f if INT_PCA_METHOD == 'bpcells' else []
  output:
    integration_f = joint_integration_f
  params:
    embedding         = config['integration']['int_embedding'],
    n_dims            = INT_N_DIMS,
    umap_n_iter = config['integration']['int_umap_n_iter'],
    cl_method         = config['integration']['int_cl_method'],
    theta_concat      = config['integration']['int_theta'],
    batch_var_concat  = config['integration']['int_batch_var'],
    res_ls_concat     = config['integration']['int_res_ls'],
    use_paga          = config['integration']['int_use_paga'],
    paga_cl_res       = get_paga_resolution(config['integration']),
    int_use_gpu       = config['integration']['int_use_gpu'],
    pca_method        = INT_PCA_METHOD,
    **({'retain_metadata': JOIN_RETAIN_METADATA} if JOIN_ZOOM_PARAMS else {})
  retries: config['resources']['retries']
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

    python3 {scprocess_dir}/scripts/integration.py run_zoom_integration \
      --hvg_mat_f        {input.hvg_mat_f} \
      --coldata_f        {input.coldata_f} \
      --demux_type       none \
      --exclude_mito     False \
      --embedding        {params.embedding} \
      --n_dims           {params.n_dims} \
      --umap_n_iter {params.umap_n_iter} \
      --cl_method        {params.cl_method} \
      --theta_concat     "{params.theta_concat}" \
      --batch_var_concat "{params.batch_var_concat}" \
      --res_ls_concat    "{params.res_ls_concat}" \
      --integration_f    {output.integration_f} \
      {JOIN_RETAIN_METADATA} \
      $( [ "{params.use_paga}" == "True" ] && echo "--use-paga" ) \
      $( [ "{params.use_paga}" == "True" ] && echo "--paga-cl-res {params.paga_cl_res}" ) \
      $( [ "{params.int_use_gpu}" == "True" ] && echo "--use-gpu" ) \
      $( [ "{params.pca_method}" == "bpcells" ] && echo "--precomputed_pca_f {input.pca_f}" )
    """


rule join_build_h5ads_yaml:
  input:
    h5ads_yaml_fs = H5ADS_YAML_FS
  output:
    joint_h5ads_yaml_f = joint_h5ads_yaml_f,
    # Zoom preparation omits batches with no retained cells. Their identities
    # are not known at DAG construction, so track the generated link directory.
    h5ad_symlinks = directory(h5ads_dir) if JOIN_ZOOM_PARAMS else [f"{h5ads_dir}/{bk}.h5ad" for bk in JOIN_BATCH_KEYS]
  params:
    project_ids = " ".join(JOIN_PROJECT_IDS),
    h5ads_dir   = h5ads_dir
  retries: config['resources']['retries']
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
    python3 {scprocess_dir}/scripts/join.py build_join_h5ads_yaml \
      --h5ads_yaml_fs  {input.h5ads_yaml_fs} \
      --project_ids    {params.project_ids} \
      --h5ads_dir      {params.h5ads_dir} \
      --out_yaml_f     {output.joint_h5ads_yaml_f}
    """


rule join_make_pseudobulks:
  input:
    h5ads_yaml_f  = joint_h5ads_yaml_f,
    integration_f = joint_integration_f,
    marker_script = marker_edger_bp_script
  output:
    pb_dir = directory(pb_dir)
  params:
    sel_res     = MKR_SEL_RES,
    min_cl_size = config['marker_genes']['mkr_min_cl_size'],
    batch_var   = "sample_id"
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_make_pseudobulks', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_make_pseudobulks', 'time', attempt)
  log:
    f"{logs_dir}/join_make_pseudobulks_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_make_pseudobulks_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/rlibs_bpcells.yaml'
  shell: """
    exec &>> {log}
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    Rscript --vanilla -e "source('{input.marker_script}'); make_pseudobulks_bpcells(
      integration_f = '{input.integration_f}',
      h5ads_yaml_f  = '{input.h5ads_yaml_f}',
      pb_dir        = '{output.pb_dir}',
      sel_res       = '{params.sel_res}',
      min_cl_size   =  {params.min_cl_size},
      batch_var     = '{params.batch_var}',
      n_cores       =  {threads})"
    """


rule join_marker_genes_prep:
  input:
    pb_dir        = pb_dir,
    marker_script = marker_edger_bp_script
  output:
    prepared_coldata_f = pb_prepared_coldata_f,
    prepared_rowdata_f = pb_prepared_rowdata_f
  params:
    min_cells        = config['marker_genes']['mkr_min_cells'],
    adjust_project_id = str(config['marker_genes']['mkr_adjust_project_id']).upper()
  threads: 8
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_marker_genes_prep', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_marker_genes_prep', 'time', attempt)
  log:
    f"{logs_dir}/join_marker_genes_prep_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_marker_genes_prep_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/rlibs_bpcells.yaml'
  shell: """
    exec &>> {log}
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    Rscript --vanilla -e "source('{input.marker_script}'); prepare_pseudobulks_bpcells(
      pb_dir             = '{input.pb_dir}',
      prepared_coldata_f = '{output.prepared_coldata_f}',
      prepared_rowdata_f = '{output.prepared_rowdata_f}',
      min_cells          =  {params.min_cells},
      adjust_project_id  =  {params.adjust_project_id},
      n_cores            =  {threads})"
    """


rule join_calc_hvgs:
  input:
    pb_dir             = pb_dir,
    prepared_coldata_f = pb_prepared_coldata_f,
    marker_script      = marker_edger_bp_script
  output:
    pb_hvgs_f = pb_hvgs_f
  params:
    gene_info_f = GENE_INFO_F
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_calc_hvgs', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_calc_hvgs', 'time', attempt)
  log:
    f"{logs_dir}/join_calc_hvgs_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_calc_hvgs_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/rlibs_bpcells.yaml'
  shell: """
    exec &>> {log}
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    Rscript --vanilla -e "source('{input.marker_script}'); calculate_pseudobulk_hvgs_bpcells(
      pb_dir             = '{input.pb_dir}',
      prepared_coldata_f = '{input.prepared_coldata_f}',
      pb_hvgs_f          = '{output.pb_hvgs_f}',
      biotypes_dt        = load_gene_biotypes_bpcells('{params.gene_info_f}'),
      n_cores            =  {threads})"
    """


rule join_marker_genes_calc:
  input:
    pb_dir             = pb_dir,
    prepared_coldata_f = pb_prepared_coldata_f,
    prepared_rowdata_f = pb_prepared_rowdata_f,
    marker_script      = marker_edger_bp_script
  output:
    mkrs_f = mkrs_f
  params:
    gene_info_f       = GENE_INFO_F,
    adjust_project_id = str(config['marker_genes']['mkr_adjust_project_id']).upper(),
    dispersion_use_clusters = str(config['marker_genes']['mkr_dispersion_use_clusters']).upper()
  threads: config['resources']['n_join_marker_genes_calc']
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_marker_genes_calc', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_marker_genes_calc', 'time', attempt)
  log:
    f"{logs_dir}/join_marker_genes_calc_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_marker_genes_calc_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/rlibs_bpcells.yaml'
  shell: """
    exec &>> {log}
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    export OMP_NUM_THREADS=1

    Rscript --vanilla -e "source('{input.marker_script}'); calculate_marker_genes_bpcells(
      pb_dir             = '{input.pb_dir}',
      prepared_coldata_f = '{input.prepared_coldata_f}',
      prepared_rowdata_f = '{input.prepared_rowdata_f}',
      mkrs_f             = '{output.mkrs_f}',
      biotypes_dt        = load_gene_biotypes_bpcells('{params.gene_info_f}'),
      adjust_project_id  =  {params.adjust_project_id},
      dispersion_use_clusters = {params.dispersion_use_clusters},
      n_cores            =  {threads})"
    """


rule join_marker_plot_data:
  input:
    pb_dir             = pb_dir,
    prepared_coldata_f = pb_prepared_coldata_f,
    mkrs_f             = mkrs_f,
    pb_hvgs_f          = pb_hvgs_f
  output:
    pb_plot_data_f = pb_plot_data_f
  params:
    min_cpm_mkr = config['marker_genes']['mkr_min_cpm_mkr']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_marker_plot_data', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_marker_plot_data', 'time', attempt)
  log:
    f"{logs_dir}/join_marker_plot_data_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.log"
  conda:
    '../envs/rlibs_bpcells.yaml'
  shell: """
    exec &>> {log}
    Rscript --vanilla -e "source('{scprocess_dir}/scripts/utils.R'); source('{scprocess_dir}/scripts/marker_genes_edger_bp.R'); source('{scprocess_dir}/scripts/marker_genes.R'); make_report_logcpms_bpcells(
      pb_dir             = '{input.pb_dir}',
      prepared_coldata_f = '{input.prepared_coldata_f}',
      mkrs_f             = '{input.mkrs_f}',
      pb_hvgs_f          = '{input.pb_hvgs_f}',
      out_f              = '{output.pb_plot_data_f}',
      min_cpm_mkr        = {params.min_cpm_mkr})"
    """


rule join_fgsea:
  input:
    mkrs_f = mkrs_f
  output:
    fgsea_go_bp_f = fgsea_bp_f,
    fgsea_go_cc_f = fgsea_cc_f,
    fgsea_go_mf_f = fgsea_mf_f
  params:
    genome_ref  = GENOME_REF,
    gsea_dir    = GSEA_DIR,
    min_cpm_go  = _mkr_cfg['mkr_min_cpm_go'],
    max_zero_p  = _mkr_cfg['mkr_max_zero_p'],
    gsea_cut    = _mkr_cfg['mkr_gsea_cut'],
    not_ok_re   = _mkr_cfg['mkr_not_ok_re'],
    gsea_var    = _mkr_cfg['mkr_gsea_var']
  threads: 8
  retries: config['resources']['retries']
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
    Rscript --vanilla -e "source('{scprocess_dir}/scripts/utils.R'); source('{scprocess_dir}/scripts/fgsea.R'); run_fgsea(
      mkrs_f        = '{input.mkrs_f}',
      fgsea_go_bp_f = '{output.fgsea_go_bp_f}',
      fgsea_go_cc_f = '{output.fgsea_go_cc_f}',
      fgsea_go_mf_f = '{output.fgsea_go_mf_f}',
      genome_ref    = '{params.genome_ref}',
      gsea_dir      = '{params.gsea_dir}',
      min_cpm_go    = {params.min_cpm_go},
      max_zero_p    = {params.max_zero_p},
      gsea_cut      = {params.gsea_cut},
      not_ok_re     = '{params.not_ok_re}',
      gsea_var      = '{params.gsea_var}',
      n_cores       =  {threads})"
    """


def _get_labeller_entry(labeller, model):
  matches = [e for e in LABELLER_PARAMS
    if e['labeller'] == labeller and e['model'] == model]
  if len(matches) != 1:
    raise ValueError(f"Expected exactly one labeller entry for {labeller}/{model}, got {len(matches)}")
  return matches[0]

def _join_merge_labels_inputs(wildcards):
  fresh_batches, reuse_label_fs = get_join_batch_sources(
    PROJECT_CFGS, JOIN_PROJECT_IDS, JOIN_BATCH_KEYS, wildcards.labeller, wildcards.model)
  pred_fs = expand(
    f"{join_lbl_dir}/tmp_labels_{{labeller}}_model_{{model}}_{JOIN_TAG}_{DATE_STAMP}_{{batch}}.csv.gz",
    batch=fresh_batches, allow_missing=True)
  if reuse_label_fs:
    pred_fs.append(f"{join_lbl_dir}/tmp_labels_{wildcards.labeller}_model_{wildcards.model}_{JOIN_TAG}_{DATE_STAMP}_reused.csv.gz")
  return pred_fs


rule join_celltypist:
  """Run CellTypist on each batch h5ad for the join integration."""
  input:
    adata_f = f"{join_int_dir}/h5ads/{{batch}}.h5ad"
  output:
    pred_f  = temp(f"{join_lbl_dir}/tmp_labels_celltypist_model_{{model}}_{JOIN_TAG}_{DATE_STAMP}_{{batch}}.csv.gz")
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_celltypist', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_celltypist', 'time', attempt)
  log:
    f"{logs_dir}/join_celltypist_{{model}}_{{batch}}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_celltypist_{{model}}_{{batch}}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/label_celltypes.yaml'
  shell: """
    exec &>> {log}
    python3 {scprocess_dir}/scripts/label_celltypes.py celltypist_one_batch \
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
    model_f = lambda wildcards: _get_labeller_entry('scprocess', wildcards.model)['model_f'],
    cls_f   = lambda wildcards: _get_labeller_entry('scprocess', wildcards.model)['cls_f'],
    genes_f = lambda wildcards: _get_labeller_entry('scprocess', wildcards.model)['genes_f']
  threads: 1
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_scprocess_labeller', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_scprocess_labeller', 'time', attempt)
  log:
    f"{logs_dir}/join_scprocess_labeller_{{model}}_{{batch}}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_scprocess_labeller_{{model}}_{{batch}}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/label_celltypes.yaml'
  shell: """
    exec &>> {log}
    python3 {scprocess_dir}/scripts/label_celltypes.py xgboost_one_batch \
      {wildcards.batch} sample_id {wildcards.model} \
      --adata_f   {input.adata_f} \
      --model_f   {params.model_f} \
      --cls_f     {params.cls_f} \
      --genes_f   {params.genes_f} \
      --pred_f    {output.pred_f}
    """


rule join_extract_labels:
  """Extract naive predictions from source project labels for cells in the join."""
  input:
    source_labels_fs = lambda wildcards: get_join_batch_sources(PROJECT_CFGS, JOIN_PROJECT_IDS, JOIN_BATCH_KEYS, wildcards.labeller, wildcards.model)[1],
    integration_f    = joint_integration_f
  output:
    pred_f = temp(f"{join_lbl_dir}/tmp_labels_{{labeller}}_model_{{model}}_{JOIN_TAG}_{DATE_STAMP}_reused.csv.gz")
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_extract_labels', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_extract_labels', 'time', attempt)
  log:
    f"{logs_dir}/join_extract_labels_{{labeller}}_{{model}}_{DATE_STAMP}.log"
  conda:
    '../envs/label_celltypes.yaml'
  shell: """
    exec &>> {log}
    python3 {scprocess_dir}/scripts/label_celltypes.py extract_naive_predictions \
      {input.source_labels_fs} \
      --int_f     {input.integration_f} \
      --model     {wildcards.model} \
      --pred_f    {output.pred_f}
    """


rule join_merge_labels:
  """Aggregate per-batch predictions by majority voting."""
  input:
    pred_fs       = _join_merge_labels_inputs,
    integration_f = joint_integration_f,
    label_map_f   = lambda wildcards: (
      _get_labeller_entry(wildcards.labeller, wildcards.model).get('label_map_f') or [])
  output:
    pred_out_f    = f"{join_lbl_dir}/labels_{{labeller}}_model_{{model}}_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
  params:
    hi_res_cl   = lambda wildcards: _get_labeller_entry(wildcards.labeller, wildcards.model)['hi_res_cl'],
    min_cl_prop = lambda wildcards: _get_labeller_entry(wildcards.labeller, wildcards.model)['min_cl_prop'],

    # An absent optional input becomes an empty argument for the CLI.
    label_map_f = lambda wildcards, input: str(input.label_map_f)
  threads: 4
  retries: config['resources']['retries']
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_merge_labels', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_merge_labels', 'time', attempt)
  log:
    f"{logs_dir}/join_merge_labels_{{labeller}}_{{model}}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_merge_labels_{{labeller}}_{{model}}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/label_celltypes.yaml'
  shell: """
    exec &>> {log}
    python3 {scprocess_dir}/scripts/label_celltypes.py aggregate_predictions \
      {input.pred_fs} \
      --int_f           {input.integration_f} \
      --hi_res_cl       {params.hi_res_cl} \
      --min_cl_prop     {params.min_cl_prop} \
      --batch_var       sample_id \
      --label_map_f     {params.label_map_f:q} \
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
    '../envs/label_celltypes.yaml'
  shell: """
    exec &>> {log}
    python3 {scprocess_dir}/scripts/label_celltypes.py save_cluster_names \
      --labels_f      {input.labels_f} \
      --integration_f {input.integration_f} \
      --mkr_sel_res   {params.mkr_sel_res} \
      --output_f      {output.names_f}
    """


rule join_render_html:
  input:
    xgb_specs = ([f'{join_xgb_dir}/{_join_xgb_ref_tag}_label_specs.csv', f'{join_xgb_dir}/{_join_xgb_ref_tag}_coarse_label_specs.csv'] if DO_TRAIN_XGB else []),
    integration_f = joint_integration_f,
    qc_f          = joint_qc_f,
    sample_meta_f = joint_sample_meta_f,
    mkrs_f        = mkrs_f,
    pb_hvgs_f     = pb_hvgs_f,
    pb_plot_data_f = pb_plot_data_f,
    fgsea_files   = [fgsea_bp_f, fgsea_cc_f, fgsea_mf_f] if DO_GSEA else [],
    label_files   = label_fs,
    cluster_names = cluster_names_fs,
    xgb_files = ([f'{join_xgb_dir}/{_join_xgb_ref_tag}_model_card.json',
      _join_xgb_fulldata_f, str(scprocess_dir / 'scripts/train_xgboost.R')] if DO_TRAIN_XGB else []),
    report_helper_f = str(scprocess_dir / "scripts/integration.R"),
    rmd_template_f = str(scprocess_dir / "resources/rmd_templates/join.Rmd.template")
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
    ref_txome     = GENOME_REF,
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
    qc_f             = joint_qc_f,
    mkrs_f           = mkrs_f,
    pb_hvgs_f        = pb_hvgs_f,
    pb_plot_data_f   = pb_plot_data_f,
    fgsea_go_bp_f    = fgsea_bp_f if DO_GSEA else '',
    fgsea_go_cc_f    = fgsea_cc_f if DO_GSEA else '',
    fgsea_go_mf_f    = fgsea_mf_f if DO_GSEA else '',
    proj_dir         = str(JOIN_DIR),
    scprocess_dir    = str(scprocess_dir),
    date_stamp       = DATE_STAMP,
    do_xgboost       = DO_TRAIN_XGB,
    xgb_predictions_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_predictions.csv.gz' if DO_TRAIN_XGB else '',
    xgb_importance_f  = f'{join_xgb_dir}/{_join_xgb_ref_tag}_gene_importance.csv' if DO_TRAIN_XGB else '',
    xgb_model_card_f  = f'{join_xgb_dir}/{_join_xgb_ref_tag}_model_card.json' if DO_TRAIN_XGB else '',
    xgb_pseudobulk_f  = f'{join_xgb_dir}/{_join_xgb_ref_tag}_pseudobulk.h5ad' if DO_TRAIN_XGB else '',
    xgb_has_coarse    = ("true" if _join_xgb_cfg.get('label_map_f') else "false") if DO_TRAIN_XGB else 'false',
    xgb_min_cells     = _join_xgb_cfg.get('min_cells_expressed', 10) if DO_TRAIN_XGB else 10
  threads: 1
  retries: config['resources']['retries']
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
    cp {scprocess_dir}/scripts/utils.R        {output.r_utils_f}
    cp {scprocess_dir}/scripts/integration.R  {output.r_int_f}
    cp {scprocess_dir}/scripts/marker_genes.R {output.r_mkr_f}
    cp {scprocess_dir}/scripts/fgsea.R        {output.r_fgsea_f}
    cp {scprocess_dir}/scripts/label_celltypes.R {output.r_lbl_f}

    # define rule and template
    template_f=$(realpath {input.rmd_template_f})
    rule="join"

    # rendering html
    Rscript --vanilla -e "source('{scprocess_dir}/scripts/render_htmls.R'); \\
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
      qc_f             = '{params.qc_f}',
      sample_meta_f    = '{input.sample_meta_f}',
      mkrs_f           = '{params.mkrs_f}',
      pb_hvgs_f        = '{params.pb_hvgs_f}',
      pb_plot_data_f   = '{params.pb_plot_data_f}',
      fgsea_go_bp_f    = '{params.fgsea_go_bp_f}',
      fgsea_go_cc_f    = '{params.fgsea_go_cc_f}',
      fgsea_go_mf_f    = '{params.fgsea_go_mf_f}',
      scprocess_dir    = '{params.scprocess_dir}',
      date_stamp       = '{params.date_stamp}',
      do_xgboost       = '{params.do_xgboost}',
      xgb_predictions_f = '{params.xgb_predictions_f}',
      xgb_importance_f  = '{params.xgb_importance_f}',
      xgb_model_card_f  = '{params.xgb_model_card_f}',
      xgb_pseudobulk_f  = '{params.xgb_pseudobulk_f}',
      xgb_has_coarse    = '{params.xgb_has_coarse}',
      xgb_min_cells     = '{params.xgb_min_cells}'
    )"
    """


if DO_TRAIN_XGB:

  rule join_train_xgboost_prepare:
    input:
      annots_f = _join_xgb_cfg['annots_f'],
      integration_f = joint_integration_f,
      h5ads_yaml = joint_h5ads_yaml_f,
      label_map = [_join_xgb_cfg['label_map_f']] if _join_xgb_cfg.get('label_map_f') else [],
      planner = str(scprocess_dir / 'scripts/train_xgboost.py'),
    output:
      training_metadata_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_training_metadata.csv.gz',
      coverage_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_split_coverage.csv',
      counts_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_label_counts.csv',
      metadata_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_training_metadata.csv.gz' + '.json',
    params:
      excluded_labels                  = _join_xgb_cfg['excluded_labels'],
      label_map_f                      = _join_xgb_cfg.get('label_map_f') or '',
      cluster_majority_vote_resolution = _join_xgb_cfg.get('cluster_majority_vote_resolution', ''),
      cluster_majority_vote_min_purity = _join_xgb_cfg['cluster_majority_vote_min_purity'],
      test_group_var                   = _join_xgb_cfg.get('test_group_var') or '',
      test_groups                      = _join_xgb_cfg.get('test_groups', []),
      validation_group_var             = _join_xgb_cfg['validation_group_var'],
      validation_fraction              = _join_xgb_cfg['validation_fraction'],
      validation_variable_combinations = _join_xgb_cfg['validation_variable_combinations'],
      sampling_variable_combinations   = _join_xgb_cfg['sampling_variable_combinations'],
      train_cells_per_label             = _join_xgb_cfg['train_cells_per_label'],
      validation_cells_per_label        = _join_xgb_cfg['validation_cells_per_label'],
      min_cells_per_label               = _join_xgb_cfg['min_cells_per_label'],
      seed                             = _join_xgb_cfg['seed'],
      batch_var = 'sample_id',
      is_join   = 1,
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb = config['resources']['gb_join_train_xgboost_prepare'] * 1024,
      runtime = config['resources']['mins_join_train_xgboost_prepare']
    log: f'{join_xgb_dir}/{_join_xgb_ref_tag}_training_metadata.csv.gz' + '.log'
    conda: '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/train_xgboost.py prepare \
        --annots_f {input.annots_f:q} --integration_f {input.integration_f:q} \
        --h5ads_yaml {input.h5ads_yaml:q} --batch_var {params.batch_var:q} \
        --is_join {params.is_join} \
        --excluded_labels {params.excluded_labels:q} \
        --label_map_f={params.label_map_f:q} \
        --cluster_majority_vote_resolution={params.cluster_majority_vote_resolution:q} \
        --cluster_majority_vote_min_purity {params.cluster_majority_vote_min_purity:q} \
        --test_group_var={params.test_group_var:q} \
        --test_groups {params.test_groups:q} \
        --validation_group_var {params.validation_group_var:q} \
        --validation_fraction {params.validation_fraction:q} \
        --validation_variable_combinations {params.validation_variable_combinations:q} \
        --sampling_variable_combinations {params.sampling_variable_combinations:q} \
        --train_cells_per_label {params.train_cells_per_label:q} \
        --validation_cells_per_label {params.validation_cells_per_label:q} \
        --min_cells_per_label {params.min_cells_per_label:q} \
        --seed {params.seed:q} --training_metadata_f {output.training_metadata_f:q} \
        --coverage_f {output.coverage_f:q} --counts_f {output.counts_f:q}
      """


  rule join_train_xgboost_train:
    input:
      training_metadata_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_training_metadata.csv.gz',
      h5ads_yaml = joint_h5ads_yaml_f,
      trainer = str(scprocess_dir / "scripts/train_xgboost.py"),
    output:
      model_f   = f'{join_xgb_dir}/{_join_xgb_ref_tag}_xgboost_model.json',
      cls_f     = f'{join_xgb_dir}/{_join_xgb_ref_tag}_allowed_cls.csv',
      genes_f   = f'{join_xgb_dir}/{_join_xgb_ref_tag}_selected_genes.txt',
      imp_f     = f'{join_xgb_dir}/{_join_xgb_ref_tag}_gene_importance.csv',
      preds_f   = f'{join_xgb_dir}/{_join_xgb_ref_tag}_predictions.csv.gz',
      pb_f      = f'{join_xgb_dir}/{_join_xgb_ref_tag}_pseudobulk.h5ad',
      model_card_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_model_card.json',
    params:
      output_dir = join_xgb_dir,
      xgboost_tag            = _join_xgb_cfg['xgboost_tag'],
      use_gpu                = int(_join_xgb_cfg['use_gpu']),
      gene_exclude_re        = _join_xgb_cfg.get('gene_exclude_re') or '',
      min_cells_expressed    = _join_xgb_cfg['min_cells_expressed'],
      seed                   = _join_xgb_cfg['seed'],
      pass1_subsample        = _join_xgb_cfg['pass1_subsample'],
      pass1_colsample_bytree = _join_xgb_cfg['pass1_colsample_bytree'],
      pass1_learning_rate    = _join_xgb_cfg['pass1_learning_rate'],
      pass1_nrounds          = _join_xgb_cfg['pass1_nrounds'],
      pass1_early_stopping   = _join_xgb_cfg['pass1_early_stopping'],
      gain_threshold         = _join_xgb_cfg['gain_threshold'],
      min_genes              = _join_xgb_cfg['min_genes'],
      max_genes              = _join_xgb_cfg['max_genes'],
      pass2_colsample_bytree = _join_xgb_cfg['pass2_colsample_bytree'],
      pass2_learning_rate    = _join_xgb_cfg['pass2_learning_rate'],
      pass2_nrounds          = _join_xgb_cfg['pass2_nrounds'],
      pass2_early_stopping   = _join_xgb_cfg['pass2_early_stopping'],
      label_map_f            = _join_xgb_cfg.get('label_map_f') or '',
    threads: 8
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_train_xgboost_train', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_train_xgboost_train', 'time', attempt)
    log:
      f"{logs_dir}/join_train_xgboost_{_join_xgb_ref_tag}_{DATE_STAMP}.log"
    benchmark:
      f"{benchmark_dir}/join_train_xgboost_{_join_xgb_ref_tag}_{DATE_STAMP}.benchmark.txt"
    conda:
      '../envs/train_xgboost_gpu.yaml' if _join_xgb_cfg['use_gpu'] else '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/train_xgboost.py train \
        --training_metadata_f {input.training_metadata_f:q} --h5ads_yaml {input.h5ads_yaml:q} \
        --xgboost_tag {params.xgboost_tag:q} \
        --use_gpu {params.use_gpu:q} \
        --gene_exclude_re={params.gene_exclude_re:q} \
        --min_cells_expressed {params.min_cells_expressed:q} \
        --seed {params.seed:q} \
        --pass1_subsample {params.pass1_subsample:q} \
        --pass1_colsample_bytree {params.pass1_colsample_bytree:q} \
        --pass1_learning_rate {params.pass1_learning_rate:q} \
        --pass1_nrounds {params.pass1_nrounds:q} \
        --pass1_early_stopping {params.pass1_early_stopping:q} \
        --gain_threshold {params.gain_threshold:q} \
        --min_genes {params.min_genes:q} \
        --max_genes {params.max_genes:q} \
        --pass2_colsample_bytree {params.pass2_colsample_bytree:q} \
        --pass2_learning_rate {params.pass2_learning_rate:q} \
        --pass2_nrounds {params.pass2_nrounds:q} \
        --pass2_early_stopping {params.pass2_early_stopping:q} \
        --label_map_f={params.label_map_f:q} --output_dir {params.output_dir:q} \
        --n_cores {threads}
      """

  rule join_train_xgboost_predict:
    input:
      adata_f = f"{join_int_dir}/h5ads/{{batch}}.h5ad",
      model_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_xgboost_model.json',
      cls_f   = f'{join_xgb_dir}/{_join_xgb_ref_tag}_allowed_cls.csv',
      genes_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_selected_genes.txt',
    output:
      pred_f  = temp(f"{join_xgb_dir}/tmp_fulldata_predict_{_join_xgb_ref_tag}_{{batch}}.csv.gz")
    params:
      ref_tag = _join_xgb_ref_tag
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_train_xgboost_predict', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_train_xgboost_predict', 'time', attempt)
    log:
      f"{logs_dir}/join_train_xgboost_predict_{_join_xgb_ref_tag}_{{batch}}_{DATE_STAMP}.log"
    benchmark:
      f"{benchmark_dir}/join_train_xgboost_predict_{_join_xgb_ref_tag}_{{batch}}_{DATE_STAMP}.benchmark.txt"
    conda:
      '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/label_celltypes.py xgboost_one_batch \
        {wildcards.batch} batch {params.ref_tag} \
        --adata_f   {input.adata_f} \
        --model_f   {input.model_f} \
        --cls_f     {input.cls_f} \
        --genes_f   {input.genes_f} \
        --pred_f    {output.pred_f}
      """


  rule join_train_xgboost_combine:
    input:
      pred_fs          = expand(
        f"{join_xgb_dir}/tmp_fulldata_predict_{_join_xgb_ref_tag}_{{batch}}.csv.gz",
        batch=JOIN_BATCH_KEYS),


      training_metadata_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_training_metadata.csv.gz',
      integration_f = joint_integration_f,
    output:
      fulldata_f = _join_xgb_fulldata_f
    params:
      ref_tag      = _join_xgb_ref_tag,
      label_map_f  = _join_xgb_cfg.get('label_map_f') or '',
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_train_xgboost_combine', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_train_xgboost_combine', 'time', attempt)
    log:
      f"{logs_dir}/join_train_xgboost_combine_{_join_xgb_ref_tag}_{DATE_STAMP}.log"
    benchmark:
      f"{benchmark_dir}/join_train_xgboost_combine_{_join_xgb_ref_tag}_{DATE_STAMP}.benchmark.txt"
    conda:
      '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/train_xgboost.py combine \
        {input.pred_fs:q} --training_metadata_f {input.training_metadata_f:q} --integration_f {input.integration_f:q} \
        --output_f {output.fulldata_f:q}
      """


  rule join_train_xgboost_specs:
    input:
      cls_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_allowed_cls.csv',
      counts_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_label_counts.csv',
      predictions_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_fulldata_predictions.csv.gz',
      specifications = get_xgboost_spec_inputs(_join_xgb_cfg),
      script_f = str(scprocess_dir / "scripts/train_xgboost.py"),
    output:
      label_specs_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_label_specs.csv',
      coarse_label_specs_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_coarse_label_specs.csv',
      coarse_evaluation_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_coarse_evaluation.csv',
      coarse_per_label_evaluation_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_coarse_per_label_evaluation.csv',
    params:
      xgboost_tag          = _join_xgb_cfg['xgboost_tag'],
      label_specs_f        = _join_xgb_cfg.get('label_specs_f') or '',
      coarse_label_specs_f = _join_xgb_cfg.get('coarse_label_specs_f') or '',
      label_map_f          = _join_xgb_cfg.get('label_map_f') or '',
      output_dir = str(pathlib.Path(f'{join_xgb_dir}/{_join_xgb_ref_tag}').parent),
    threads: 1
    resources:
      mem_mb = 2048,
      runtime = 5
    log: f'{join_xgb_dir}/{_join_xgb_ref_tag}_specs.log'
    conda: '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 {scprocess_dir}/scripts/train_xgboost.py specs \
        --cls_f {input.cls_f:q} --counts_f {input.counts_f:q} \
        --predictions_f {input.predictions_f:q} --xgboost_tag {params.xgboost_tag:q} \
        --label_specs_f={params.label_specs_f:q} \
        --coarse_label_specs_f={params.coarse_label_specs_f:q} \
        --label_map_f={params.label_map_f:q} \
        --output_dir {params.output_dir:q}
      """


  rule join_render_html_train_xgboost:
    input:
      model_card_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_model_card.json',
      train_script_f = str(scprocess_dir / 'scripts/train_xgboost.R'),
      template_f = str(scprocess_dir / 'resources/rmd_templates/train_xgboost.Rmd.template'),
      label_specs_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_label_specs.csv',
      coarse_label_specs_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_coarse_label_specs.csv',
      preds_f              = f'{join_xgb_dir}/{_join_xgb_ref_tag}_predictions.csv.gz',
      imp_f                = f'{join_xgb_dir}/{_join_xgb_ref_tag}_gene_importance.csv',
      pb_f                 = f'{join_xgb_dir}/{_join_xgb_ref_tag}_pseudobulk.h5ad',
      label_counts_f       = f'{join_xgb_dir}/{_join_xgb_ref_tag}_label_counts.csv',
      fulldata_preds_f     = _join_xgb_fulldata_f
    output:
      rmd_f         = f"{rmd_dir}/{JOIN_TAG}_train_xgboost.Rmd",
      html_f        = f"{docs_dir}/{JOIN_TAG}_train_xgboost.html"
    params:
      your_name            = config['join'].get('your_name', ''),
      affiliation          = config['join'].get('affiliation', ''),
      short_tag            = JOIN_TAG,
      ref_tag              = _join_xgb_ref_tag,
      proj_dir             = str(JOIN_DIR),
      output_dir           = join_xgb_dir,
      predictions_f        = f'{join_xgb_dir}/{_join_xgb_ref_tag}_predictions.csv.gz',
      importance_f         = f'{join_xgb_dir}/{_join_xgb_ref_tag}_gene_importance.csv',
      pseudobulk_f         = f'{join_xgb_dir}/{_join_xgb_ref_tag}_pseudobulk.h5ad',
      label_counts_f       = f'{join_xgb_dir}/{_join_xgb_ref_tag}_label_counts.csv',
      fulldata_predictions_f = _join_xgb_fulldata_f,
      integration_f        = joint_integration_f,
      has_coarse           = "true" if _join_xgb_cfg.get('label_map_f') else "false",
      min_cells            = _join_xgb_cfg.get('min_cells_expressed', 10),
    threads: 1
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_render_html_train_xgboost', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_render_html_train_xgboost', 'time', attempt)
    log:
      f"{logs_dir}/join_render_html_train_xgboost_{_join_xgb_ref_tag}_{DATE_STAMP}.log"
    benchmark:
      f"{benchmark_dir}/join_render_html_train_xgboost_{_join_xgb_ref_tag}_{DATE_STAMP}.benchmark.txt"
    conda:
      '../envs/rlibs.yaml'
    shell: """
      exec &>> {log}
      cp {scprocess_dir}/scripts/utils.R          {params.proj_dir}/code/utils.R
      cp {scprocess_dir}/scripts/marker_genes.R   {params.proj_dir}/code/marker_genes.R
      cp {scprocess_dir}/scripts/train_xgboost.R  {params.proj_dir}/code/train_xgboost.R

      template_f={scprocess_dir}/resources/rmd_templates/train_xgboost.Rmd.template
      rule="train_xgboost"

      Rscript --vanilla -e "source('{scprocess_dir}/scripts/render_htmls.R'); \\
        render_html(
          rule_name       = '$rule',
          proj_dir        = '{params.proj_dir}',
          temp_f          = '$template_f',
          rmd_f           = '{output.rmd_f}',
          your_name       = '{params.your_name}',
          affiliation     = '{params.affiliation}',
          short_tag       = '{params.short_tag}',
          ref_tag         = '{params.ref_tag}',
          predictions_f   = '{params.predictions_f}',
          importance_f    = '{params.importance_f}',
          model_card_f    = '{input.model_card_f}',
          pseudobulk_f    = '{params.pseudobulk_f}',
          label_counts_f         = '{params.label_counts_f}',
          fulldata_predictions_f = '{params.fulldata_predictions_f}',
          integration_f          = '{params.integration_f}',
          has_coarse             = '{params.has_coarse}',
          min_cells              = '{params.min_cells}',
          n_cores                = '{threads}'
        )"
      """
