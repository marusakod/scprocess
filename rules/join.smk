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
import json

sys.path.append('scripts')
from scprocess_utils import (check_join_config, get_resources, prep_resource_params,
  get_join_project_parameters, get_join_source_labels_f, get_join_batch_sources)

# setup
scprocess_dir = pathlib.Path(config.pop('scprocess_dir'))
scdata_dir    = pathlib.Path(os.getenv('SCPROCESS_DATA_DIR'))
join_schema_f = scprocess_dir / "resources/schemas/join.schema.json"

# validate config, apply defaults, and check shiny parameters
config = check_join_config(config, join_schema_f, scdata_dir)

# resource parameters (no ML model for join — uses schema defaults + user overrides)
RESOURCE_PARAMS = prep_resource_params(config, join_schema_f, scprocess_dir)

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
GENOME_REF    = config['join'].get('probe_set', config['join'].get('ref_txome', ''))
join_int_dir  = str(JOIN_DIR / f"output/{JOIN_TAG}")
join_mkr_dir  = join_int_dir
logs_dir      = str(JOIN_DIR / ".log/join")
benchmark_dir = str(JOIN_DIR / ".resources")

# integration
INT_N_DIMS     = config['integration']['int_n_dims']
INT_RES_LS     = config['integration']['int_res_ls']
INT_PCA_METHOD = config['integration']['int_pca_method']
INT_BATCH_VAR  = config['integration']['int_batch_var']

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

# label_celltypes (optional; validated and defaults applied by check_join_config)
_lbl_cfg = config.get('label_celltypes', [])
DO_LABEL = len(_lbl_cfg) > 0
if DO_LABEL:
  with open(join_schema_f) as _f:
    _join_schema = json.load(_f)
  _lbl_schema_props = _join_schema['properties']['label_celltypes']['items']['properties']
  for entry in _lbl_cfg:
    for key, prop in _lbl_schema_props.items():
      if key not in entry and 'default' in prop:
        entry[key] = prop['default']

  _typist_ls_f  = scdata_dir / 'celltypist/celltypist_models.csv'
  _mdls_typist  = pl.read_csv(_typist_ls_f)['model'].to_list() if _typist_ls_f.is_file() else []
  _xgb_csv_f    = scdata_dir / 'xgboost' / 'xgboost_models.csv'
  _xgb_df       = pl.read_csv(_xgb_csv_f) if _xgb_csv_f.is_file() else pl.DataFrame()
  _mdls_scproc  = _xgb_df['model'].to_list() if len(_xgb_df) > 0 else []

  for entry in _lbl_cfg:
    if entry['labeller'] == 'celltypist':
      if entry['model'] not in _mdls_typist:
        raise ValueError(f"CellTypist model '{entry['model']}' not found. Valid: {', '.join(_mdls_typist)}")
    elif entry['labeller'] == 'scprocess':
      if entry['model'] not in _mdls_scproc:
        raise ValueError(f"scprocess model '{entry['model']}' not found. Valid: {', '.join(_mdls_scproc)}")
      row = _xgb_df.filter(pl.col('model') == entry['model']).row(0, named=True)
      entry['model_f'] = row['model_f']
      entry['cls_f']   = row['cls_f']
      entry['genes_f'] = row['genes_f']
      for key in ['model_f', 'cls_f', 'genes_f']:
        if not pathlib.Path(entry[key]).is_file():
          raise FileNotFoundError(f"file {entry[key]} doesn't exist; consider (re)running scprocess setup")

  LABELLER_PARAMS = _lbl_cfg
else:
  LABELLER_PARAMS = []

# output paths
joint_hvgs_f        = f"{join_int_dir}/joint_hvgs_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
joint_label_counts_f      = f"{join_int_dir}/joint_counts_hvgs_{JOIN_TAG}_{DATE_STAMP}.h5"
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

YOUR_NAME   = config['join']['your_name']
AFFILIATION = config['join']['affiliation']

INT_RES_LS_STR = ' '.join(str(r) for r in INT_RES_LS)

# train_xgboost (optional)
DO_TRAIN_XGB = 'train_xgboost' in config
if DO_TRAIN_XGB:
  _join_xgb_cfg = config['train_xgboost']
  join_xgb_dir  = str(JOIN_DIR / f"output/{JOIN_NAME}_train_xgboost")
  _join_xgb_ref_tag = _join_xgb_cfg['ref_tag']
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


if DO_TRAIN_XGB:
  rule train_xgboost:
    input:
      _join_xgb_model_f,
      _join_xgb_fulldata_f,
      f"{rmd_dir}/{JOIN_TAG}_train_xgboost.Rmd",
      f"{docs_dir}/{JOIN_TAG}_train_xgboost.html"


rule join_select_hvgs:
  input:
    var_stats_fs = VAR_STATS_FS
  output:
    joint_hvgs_f = joint_hvgs_f
  params:
    project_ids  = " ".join(JOIN_PROJECT_IDS),
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
    joint_label_counts_f = joint_label_counts_f
  params:
    project_ids   = " ".join(JOIN_PROJECT_IDS)
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
    python3 scripts/join.py build_joint_matrix \
      --joint_hvgs_f        {input.joint_hvgs_f} \
      --h5ads_yaml_fs       {input.h5ads_yaml_fs} \
      --project_ids         {params.project_ids} \
      --integrated_dt_fs    {input.integrated_fs} \
      --out_h5_f            {output.joint_label_counts_f}
    """


rule join_build_coldata:
  """Build joint coldata and sample metadata from matrix barcodes."""
  input:
    joint_label_counts_f = joint_label_counts_f,
    integrated_fs  = INTEGRATED_FS,
    sample_meta_fs = SAMPLE_META_FS
  output:
    joint_coldata_f     = joint_coldata_f,
    joint_sample_meta_f = joint_sample_meta_f
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
    python3 scripts/join.py build_joint_coldata \
      --h5_f                {input.joint_label_counts_f} \
      --project_ids         {params.project_ids} \
      --integrated_dt_fs    {input.integrated_fs} \
      --sample_meta_fs      {input.sample_meta_fs} \
      --out_coldata_f       {output.joint_coldata_f} \
      --out_sample_meta_f   {output.joint_sample_meta_f}
    """


rule join_pca:
  """Compute PCA on joint matrix using BPCells disk-backed streaming SVD."""
  input:
    counts_h5_f = joint_label_counts_f
  output:
    pca_f = joint_pca_f
  params:
    n_dims = INT_N_DIMS
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
  input:
    hvg_mat_f    = joint_label_counts_f,
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
  input:
    h5ads_yaml_fs = H5ADS_YAML_FS
  output:
    joint_h5ads_yaml_f = joint_h5ads_yaml_f,
    h5ad_symlinks = [f"{h5ads_dir}/{bk}.h5ad" for bk in JOIN_BATCH_KEYS]
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
    python3 scripts/join.py build_join_h5ads_yaml \
      --h5ads_yaml_fs  {input.h5ads_yaml_fs} \
      --project_ids    {params.project_ids} \
      --h5ads_dir      {params.h5ads_dir} \
      --out_yaml_f     {output.joint_h5ads_yaml_f}
    """


rule join_marker_genes:
  input:
    h5ads_yaml_f  = joint_h5ads_yaml_f,
    integration_f = joint_integration_f
  output:
    pb_f      = pb_f,
    mkrs_f    = mkrs_f,
    pb_hvgs_f = pb_hvgs_f
  params:
    gene_info_f = GENE_INFO_F,
    sel_res     = MKR_SEL_RES,
    min_cl_size = config['marker_genes']['mkr_min_cl_size'],
    min_cells   = config['marker_genes']['mkr_min_cells'],
    batch_var   = "sample_id"
  threads: 8
  retries: config['resources']['retries']
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
    Rscript -e "source('scripts/utils.R'); source('scripts/fgsea.R'); run_fgsea(
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
    python3 scripts/label_celltypes.py xgboost_one_batch \
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
    '../envs/hvgs.yaml'
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
    hi_res_cl   = lambda wildcards: _get_labeller_entry(wildcards.labeller, wildcards.model)['hi_res_cl'],
    min_cl_prop = lambda wildcards: _get_labeller_entry(wildcards.labeller, wildcards.model)['min_cl_prop']
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
    '../envs/hvgs.yaml'
  shell: """
    exec &>> {log}
    python3 scripts/label_celltypes.py save_cluster_names \
      --labels_f      {input.labels_f} \
      --integration_f {input.integration_f} \
      --mkr_sel_res   {params.mkr_sel_res} \
      --output_f      {output.names_f}
    """


rule join_render_html:
  input:
    integration_f = joint_integration_f,
    mkrs_f        = mkrs_f,
    pb_hvgs_f     = pb_hvgs_f,
    pb_f          = pb_f,
    fgsea_files   = [fgsea_bp_f, fgsea_cc_f, fgsea_mf_f] if DO_GSEA else [],
    label_files   = label_fs,
    cluster_names = cluster_names_fs,
    xgb_files     = [],
    rmd_template_f = "resources/rmd_templates/join.Rmd.template"
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
    sample_meta_f    = joint_sample_meta_f,
    mkrs_f           = mkrs_f,
    pb_hvgs_f        = pb_hvgs_f,
    pb_f             = pb_f,
    fgsea_go_bp_f    = fgsea_bp_f if DO_GSEA else '',
    fgsea_go_cc_f    = fgsea_cc_f if DO_GSEA else '',
    fgsea_go_mf_f    = fgsea_mf_f if DO_GSEA else '',
    proj_dir         = str(JOIN_DIR),
    scprocess_dir    = str(scprocess_dir),
    date_stamp       = DATE_STAMP,
    do_xgboost       = DO_TRAIN_XGB,
    xgb_predictions_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_predictions.csv.gz' if DO_TRAIN_XGB else '',
    xgb_importance_f  = f'{join_xgb_dir}/{_join_xgb_ref_tag}_gene_importance.csv' if DO_TRAIN_XGB else '',
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
    cp scripts/utils.R        {output.r_utils_f}
    cp scripts/integration.R  {output.r_int_f}
    cp scripts/marker_genes.R {output.r_mkr_f}
    cp scripts/fgsea.R        {output.r_fgsea_f}
    cp scripts/label_celltypes.R {output.r_lbl_f}

    # define rule and template
    template_f=$(realpath {input.rmd_template_f})
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
      date_stamp       = '{params.date_stamp}',
      do_xgboost       = '{params.do_xgboost}',
      xgb_predictions_f = '{params.xgb_predictions_f}',
      xgb_importance_f  = '{params.xgb_importance_f}',
      xgb_pseudobulk_f  = '{params.xgb_pseudobulk_f}',
      xgb_has_coarse    = '{params.xgb_has_coarse}',
      xgb_min_cells     = '{params.xgb_min_cells}'
    )"
    """


if DO_TRAIN_XGB:

  rule join_train_xgboost_train:
    input:
      annots_f    = _join_xgb_cfg['annots_f'],
      cluster_csv = joint_integration_f,
      h5ads_yaml  = joint_h5ads_yaml_f,
    output:
      model_f   = f'{join_xgb_dir}/{_join_xgb_ref_tag}_xgboost_model.json',
      cls_f     = f'{join_xgb_dir}/{_join_xgb_ref_tag}_allowed_cls.csv',
      genes_f   = f'{join_xgb_dir}/{_join_xgb_ref_tag}_selected_genes.txt',
      imp_f     = f'{join_xgb_dir}/{_join_xgb_ref_tag}_gene_importance.csv',
      preds_f   = f'{join_xgb_dir}/{_join_xgb_ref_tag}_predictions.csv.gz',
      pb_f      = f'{join_xgb_dir}/{_join_xgb_ref_tag}_pseudobulk.h5ad',
      label_counts_f  = f'{join_xgb_dir}/{_join_xgb_ref_tag}_label_counts.csv',
    params:
      ref_tag              = _join_xgb_ref_tag,
      output_dir           = join_xgb_dir,
      batch_var            = INT_BATCH_VAR if not isinstance(INT_BATCH_VAR, list) else INT_BATCH_VAR[0],
      int_res_ls           = INT_RES_LS,
      label_map_f          = _join_xgb_cfg.get('label_map_f') or '',
      refine_labels        = _join_xgb_cfg.get('refine_labels', True),
      purity_threshold     = _join_xgb_cfg.get('purity_threshold', 0.65),
      n_cells_per_type     = _join_xgb_cfg.get('n_cells_per_type', 1000),
      min_cells_per_type   = _join_xgb_cfg.get('min_cells_per_type', 20),
      min_cells_expressed  = _join_xgb_cfg.get('min_cells_expressed', 10),
      gene_exclude_re      = _join_xgb_cfg.get('gene_exclude_re', '(lincRNA|lncRNA|pseudogene|antisense)'),
      seed                 = _join_xgb_cfg.get('seed', 42),
      use_gpu              = _join_xgb_cfg.get('use_gpu', False),
      pass1_subsample      = _join_xgb_cfg.get('pass1_subsample', 0.632),
      pass1_colsample_bytree = _join_xgb_cfg.get('pass1_colsample_bytree', 0.1),
      pass1_learning_rate  = _join_xgb_cfg.get('pass1_learning_rate', 0.1),
      pass1_nrounds        = _join_xgb_cfg.get('pass1_nrounds', 300),
      pass1_early_stopping = _join_xgb_cfg.get('pass1_early_stopping', 10),
      pass2_colsample_bytree = _join_xgb_cfg.get('pass2_colsample_bytree', 0.5),
      pass2_learning_rate  = _join_xgb_cfg.get('pass2_learning_rate', 0.05),
      pass2_nrounds        = _join_xgb_cfg.get('pass2_nrounds', 500),
      pass2_early_stopping = _join_xgb_cfg.get('pass2_early_stopping', 10),
      gain_threshold       = _join_xgb_cfg.get('gain_threshold', 0.9),
      min_genes            = _join_xgb_cfg.get('min_genes', 100),
      max_genes            = _join_xgb_cfg.get('max_genes', 3000),
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
      '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 scripts/train_xgboost.py train \
        --annots_f          {input.annots_f} \
        --cluster_csv       {input.cluster_csv} \
        --h5ads_yaml        {input.h5ads_yaml} \
        --ref_tag           {params.ref_tag} \
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
      python3 scripts/label_celltypes.py xgboost_one_batch \
        {wildcards.batch} batch {params.ref_tag} \
        --adata_f   {input.adata_f} \
        --model_f   {input.model_f} \
        --cls_f     {input.cls_f} \
        --genes_f   {input.genes_f} \
        --pred_f    {output.pred_f}
      """


  rule join_train_xgboost_aggregate:
    input:
      pred_fs          = expand(
        f"{join_xgb_dir}/tmp_fulldata_predict_{_join_xgb_ref_tag}_{{batch}}.csv.gz",
        batch=JOIN_BATCH_KEYS),
      subsample_preds_f = f'{join_xgb_dir}/{_join_xgb_ref_tag}_predictions.csv.gz',
      annots_f          = _join_xgb_cfg['annots_f'],
    output:
      fulldata_f = _join_xgb_fulldata_f
    params:
      ref_tag      = _join_xgb_ref_tag,
      label_map_f  = _join_xgb_cfg.get('label_map_f') or '',
    retries: config['resources']['retries']
    resources:
      mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_train_xgboost_aggregate', 'memory', attempt),
      runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'join_train_xgboost_aggregate', 'time', attempt)
    log:
      f"{logs_dir}/join_train_xgboost_aggregate_{_join_xgb_ref_tag}_{DATE_STAMP}.log"
    benchmark:
      f"{benchmark_dir}/join_train_xgboost_aggregate_{_join_xgb_ref_tag}_{DATE_STAMP}.benchmark.txt"
    conda:
      '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 scripts/train_xgboost.py aggregate_fulldata \
        {input.pred_fs} \
        --annots_f          {input.annots_f} \
        --subsample_preds_f {input.subsample_preds_f} \
        --output_f          {output.fulldata_f} \
        $( [ "{params.label_map_f}" != "" ] && echo "--label_map_f {params.label_map_f}" )
      """


  rule join_render_html_train_xgboost:
    input:
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
      cp scripts/utils.R          {params.proj_dir}/code/utils.R
      cp scripts/marker_genes.R   {params.proj_dir}/code/marker_genes.R
      cp scripts/train_xgboost.R  {params.proj_dir}/code/train_xgboost.R

      template_f=$(realpath resources/rmd_templates/train_xgboost.Rmd.template)
      rule="train_xgboost"

      Rscript --vanilla -e "source('scripts/render_htmls.R'); \\
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
          pseudobulk_f    = '{params.pseudobulk_f}',
          label_counts_f         = '{params.label_counts_f}',
          fulldata_predictions_f = '{params.fulldata_predictions_f}',
          integration_f          = '{params.integration_f}',
          has_coarse             = '{params.has_coarse}',
          min_cells              = '{params.min_cells}',
          n_cores                = '{threads}'
        )"
      """

