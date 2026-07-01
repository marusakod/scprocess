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
import json

sys.path.append('scripts')
from scprocess_utils import check_join_config

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

scprocess_dir = pathlib.Path(config.pop('scprocess_dir'))
scdata_dir    = pathlib.Path(os.getenv('SCPROCESS_DATA_DIR'))
join_schema_f = scprocess_dir / "resources/schemas/join.schema.json"

# validate config, apply defaults, and check shiny parameters
config = check_join_config(config, join_schema_f)

# ---------------------------------------------------------------------------
# Project config loading
# ---------------------------------------------------------------------------

JOIN_PROJECT_IDS = list(config['projects'].keys())
_project_cfgs    = {}

for _pid in JOIN_PROJECT_IDS:
  _cfg_f = config['projects'][_pid]['config']
  with open(_cfg_f) as _f:
    _project_cfgs[_pid] = yaml.safe_load(_f)

def _proj_dir(pid):
  return pathlib.Path(_project_cfgs[pid]['project']['proj_dir'])

def _proj_short_tag(pid):
  return _project_cfgs[pid]['project']['short_tag']

def _proj_full_tag(pid):
  return _project_cfgs[pid]['project']['full_tag']

def _proj_date(pid):
  return _project_cfgs[pid]['project']['date_stamp']

def _proj_hvg_dir(pid):
  return _proj_dir(pid) / f"output/{_proj_short_tag(pid)}_hvg"

def _proj_int_dir(pid):
  return _proj_dir(pid) / f"output/{_proj_short_tag(pid)}_integration"

def _proj_var_stats_f(pid):
  zoom_name = config['projects'][pid].get('zoom_name')
  if zoom_name:
    zoom_dir = _proj_dir(pid) / f"output/{_proj_short_tag(pid)}_zoom"
    return zoom_dir / zoom_name / f"standardized_variance_stats_{_proj_full_tag(pid)}_{zoom_name}_{_proj_date(pid)}.csv.gz"
  return _proj_hvg_dir(pid) / f"standardized_variance_stats_{_proj_full_tag(pid)}_{_proj_date(pid)}.csv.gz"

def _proj_h5ads_yaml_f(pid):
  return _proj_int_dir(pid) / f"h5ads_clean_paths_{_proj_full_tag(pid)}_{_proj_date(pid)}.yaml"

def _proj_integrated_dt_f(pid):
  zoom_name = config['projects'][pid].get('zoom_name')
  if zoom_name:
    zoom_dir = _proj_dir(pid) / f"output/{_proj_short_tag(pid)}_zoom"
    return zoom_dir / zoom_name / f"integrated_dt_{_proj_full_tag(pid)}_{zoom_name}_{_proj_date(pid)}.csv.gz"
  return _proj_int_dir(pid) / f"integrated_dt_{_proj_full_tag(pid)}_{_proj_date(pid)}.csv.gz"

def _proj_sample_meta_f(pid):
  meta_f = pathlib.Path(_project_cfgs[pid]['project']['sample_metadata'])
  if not meta_f.is_absolute():
    meta_f = _proj_dir(pid) / meta_f
  return meta_f

VAR_STATS_FS   = [str(_proj_var_stats_f(pid)) for pid in JOIN_PROJECT_IDS]
H5ADS_YAML_FS  = [str(_proj_h5ads_yaml_f(pid)) for pid in JOIN_PROJECT_IDS]
INTEGRATED_FS  = [str(_proj_integrated_dt_f(pid)) for pid in JOIN_PROJECT_IDS]
SAMPLE_META_FS = [str(_proj_sample_meta_f(pid)) for pid in JOIN_PROJECT_IDS]

# Pre-compute joint batch keys from per-project h5ads YAMLs (existence
# already validated by check_join_config in scprocess_utils.py).
_JOIN_BATCH_KEYS = []
for _pid, _h5yaml in zip(JOIN_PROJECT_IDS, H5ADS_YAML_FS):
  with open(_h5yaml) as _fh:
    _h5paths = yaml.safe_load(_fh)
  for _bk in _h5paths:
    _JOIN_BATCH_KEYS.append(f"{_pid}_{_bk}")

# ---------------------------------------------------------------------------
# Derived constants
# ---------------------------------------------------------------------------

JOIN_NAME  = config['join']['name']
JOIN_DIR   = pathlib.Path(config['join']['proj_dir'])
DATE_STAMP = config['join']['date_stamp']
GENOME_REF = config['join'].get('probe_set', config['join'].get('ref_txome', ''))

JOIN_TAG = f"{JOIN_NAME}_join"
MKRS_TAG = f"{JOIN_NAME}_marker_genes"

join_int_dir  = str(JOIN_DIR / f"output/{JOIN_TAG}")
join_mkr_dir  = str(JOIN_DIR / f"output/{MKRS_TAG}")
logs_dir      = str(JOIN_DIR / ".log")
benchmark_dir = str(JOIN_DIR / ".resources")

# integration
_int_cfg        = config.get('integration', {})
INT_EMBEDDING   = _int_cfg.get('int_embedding', 'harmony')
INT_BATCH_VAR   = _int_cfg.get('int_batch_var', 'sample_id')
INT_THETA_RAW   = _int_cfg.get('int_theta', 0.1)
INT_N_DIMS      = _int_cfg.get('int_n_dims', 50)
INT_CL_METHOD   = _int_cfg.get('int_cl_method', 'leiden')
INT_USE_PAGA    = _int_cfg.get('int_use_paga', True)
INT_PAGA_CL_RES = _int_cfg.get('int_paga_cl_res', 0.2)
INT_RES_LS      = _int_cfg.get('int_res_ls', [0.1, 0.2, 0.5, 1, 2])
INT_USE_GPU     = _int_cfg.get('int_use_gpu', True)

INT_RES_LS_CONCAT = " ".join(str(r) for r in INT_RES_LS)

# build concat strings for potentially-list batch_var and theta
if isinstance(INT_BATCH_VAR, list):
  INT_BATCH_VAR_CONCAT = " ".join(INT_BATCH_VAR)
  INT_BATCH_IS_LIST    = True
else:
  INT_BATCH_VAR_CONCAT = INT_BATCH_VAR
  INT_BATCH_IS_LIST    = False

if isinstance(INT_THETA_RAW, list):
  INT_THETA_CONCAT  = " ".join(str(t) for t in INT_THETA_RAW)
  INT_THETA_IS_LIST = True
else:
  INT_THETA_CONCAT  = str(INT_THETA_RAW)
  INT_THETA_IS_LIST = False

# marker genes
_mkr_cfg        = config.get('marker_genes', {})
MKR_SEL_RES     = _mkr_cfg.get('mkr_sel_res',    0.2)
MKR_MIN_CL_SIZE = _mkr_cfg.get('mkr_min_cl_size', 100)
MKR_MIN_CELLS   = _mkr_cfg.get('mkr_min_cells',   10)
MKR_NOT_OK_RE   = _mkr_cfg.get('mkr_not_ok_re',   '(lincRNA|lncRNA|pseudogene|antisense)')
MKR_MIN_CPM_MKR = _mkr_cfg.get('mkr_min_cpm_mkr', 50)
MKR_MIN_CPM_GO  = _mkr_cfg.get('mkr_min_cpm_go',  1)
MKR_MAX_ZERO_P  = _mkr_cfg.get('mkr_max_zero_p',  0.5)
MKR_DO_GSEA     = _mkr_cfg.get('mkr_do_gsea',     True)
MKR_GSEA_CUT    = _mkr_cfg.get('mkr_gsea_cut',    0.1)
MKR_GSEA_VAR    = _mkr_cfg.get('mkr_gsea_var',    'z_score')

def _get_custom_mkr_strings(genesets, proj_dir, data_dir):
  """Convert mkr_custom_genesets list to comma-separated name/path strings."""
  if not genesets:
    return '', ''
  names, paths = [], []
  for g in genesets:
    names.append(g['name'])
    if 'file' in g:
      p = pathlib.Path(g['file'])
      if not p.is_absolute():
        p = proj_dir / p
    else:
      p = data_dir / 'marker_genes' / f"{g['name']}.csv"
    paths.append(str(p))
  return ','.join(names), ','.join(paths)

MKR_CUSTOM_NAMES, MKR_CUSTOM_PATHS = _get_custom_mkr_strings(
  _mkr_cfg.get('mkr_custom_genesets', []), JOIN_DIR, scdata_dir)

GSEA_REFS   = ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020',
               'human_v1', 'human_v2', 'mouse_v1', 'mouse_v2']
DO_GSEA     = MKR_DO_GSEA and (GENOME_REF in GSEA_REFS)

# HVG
N_HVGS = config.get('hvg', {}).get('hvg_n_hvgs', 2000)

# metadata vars (space-separated string for join.py)
METADATA_VARS_STR = " ".join(config['join'].get('metadata_vars', []))

# GTF file from index_parameters.csv (needed for marker genes)
_idx_params_f = scdata_dir / 'index_parameters.csv'
_idx_params   = pl.read_csv(_idx_params_f)
GENE_INFO_F = _idx_params.filter(pl.col('reference') == GENOME_REF)['gene_info_f'][0]
GSEA_DIR = str(scdata_dir / 'gmt_pathways')

# label_celltypes (optional; validated and defaults applied by check_join_config)
_lbl_cfg = config.get('label_celltypes', [])
DO_LABEL = len(_lbl_cfg) > 0
if DO_LABEL:
  # apply schema defaults
  with open(join_schema_f) as _f:
    _join_schema = json.load(_f)
  _lbl_schema_props = _join_schema['properties']['label_celltypes']['items']['properties']
  for entry in _lbl_cfg:
    for key, prop in _lbl_schema_props.items():
      if key not in entry and 'default' in prop:
        entry[key] = prop['default']

  # validate models and resolve paths
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

# ---------------------------------------------------------------------------
# Output file paths
# ---------------------------------------------------------------------------

joint_hvgs_f        = f"{join_int_dir}/joint_hvgs_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
joint_counts_f      = f"{join_int_dir}/joint_counts_{JOIN_TAG}_{DATE_STAMP}.h5"
joint_coldata_f     = f"{join_int_dir}/joint_coldata_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
joint_sample_meta_f = f"{join_int_dir}/joint_sample_meta_{JOIN_TAG}_{DATE_STAMP}.csv"
joint_integration_f = f"{join_int_dir}/integrated_dt_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
joint_h5ads_yaml_f  = f"{join_int_dir}/h5ads_clean_paths_{JOIN_TAG}_{DATE_STAMP}.yaml"
h5ads_dir           = f"{join_int_dir}/h5ads"

pb_f        = f"{join_mkr_dir}/pb_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.rds"
mkrs_f      = f"{join_mkr_dir}/pb_marker_genes_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv.gz"
pb_hvgs_f   = f"{join_mkr_dir}/pb_hvgs_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv.gz"
fgsea_bp_f  = f"{join_mkr_dir}/fgsea_{JOIN_TAG}_{MKR_SEL_RES}_go_bp_{DATE_STAMP}.csv.gz"
fgsea_cc_f  = f"{join_mkr_dir}/fgsea_{JOIN_TAG}_{MKR_SEL_RES}_go_cc_{DATE_STAMP}.csv.gz"
fgsea_mf_f  = f"{join_mkr_dir}/fgsea_{JOIN_TAG}_{MKR_SEL_RES}_go_mf_{DATE_STAMP}.csv.gz"

join_lbl_dir = str(JOIN_DIR / f"output/{JOIN_NAME}_label_celltypes")
if DO_LABEL:
  label_fs = [
    f"{join_lbl_dir}/labels_{e['labeller']}_model_{e['model']}_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
    for e in LABELLER_PARAMS
  ]
else:
  label_fs = []

docs_dir  = str(JOIN_DIR / "public")
rmd_dir   = str(JOIN_DIR / "analysis")
code_dir  = str(JOIN_DIR / "code")
html_f    = f"{docs_dir}/{JOIN_TAG}.html"
rmd_f     = f"{rmd_dir}/{JOIN_TAG}.Rmd"

YOUR_NAME   = config['join']['your_name']
AFFILIATION = config['join']['affiliation']

INT_RES_LS_STR = ' '.join(str(r) for r in INT_RES_LS)

# ---------------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------------

# train_xgboost (optional)
DO_TRAIN_XGB = 'train_xgboost' in config
if DO_TRAIN_XGB:
  _join_xgb_cfg = config['train_xgboost']
  join_xgb_dir  = str(JOIN_DIR / f"output/{JOIN_NAME}_train_xgboost")
  _join_xgb_ref_tag = _join_xgb_cfg['ref_tag']
  _join_xgb_model_f = f"{join_xgb_dir}/{_join_xgb_ref_tag}_xgboost_model.json"

rule all:
  input:
    joint_integration_f,
    joint_h5ads_yaml_f,
    mkrs_f,
    pb_hvgs_f,
    *([fgsea_bp_f, fgsea_cc_f, fgsea_mf_f] if DO_GSEA else []),
    *label_fs,
    *([_join_xgb_model_f] if DO_TRAIN_XGB else []),
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
  log:
    f"{logs_dir}/join_select_hvgs_{JOIN_TAG}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_select_hvgs_{JOIN_TAG}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/scprocess_local.yaml'
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
  """Assemble joint count matrix and coldata from per-project h5ads."""
  input:
    joint_hvgs_f  = joint_hvgs_f,
    h5ads_yaml_fs = H5ADS_YAML_FS,
    integrated_fs = INTEGRATED_FS,
    sample_meta_fs = SAMPLE_META_FS
  output:
    joint_counts_f      = joint_counts_f,
    joint_coldata_f     = joint_coldata_f,
    joint_sample_meta_f = joint_sample_meta_f
  params:
    project_ids   = " ".join(JOIN_PROJECT_IDS),
    metadata_vars = METADATA_VARS_STR
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
      --sample_meta_fs      {input.sample_meta_fs} \
      --metadata_vars       "{params.metadata_vars}" \
      --out_h5_f            {output.joint_counts_f} \
      --out_coldata_f       {output.joint_coldata_f} \
      --out_sample_meta_f   {output.joint_sample_meta_f}
    """


rule join_integration:
  """Run Harmony integration on the joint HVG matrix."""
  input:
    hvg_mat_f    = joint_counts_f,
    coldata_f    = joint_coldata_f,
    sample_qc_f  = joint_sample_meta_f
  output:
    integration_f = joint_integration_f
  params:
    embedding         = INT_EMBEDDING,
    n_dims            = INT_N_DIMS,
    cl_method         = INT_CL_METHOD,
    theta_concat      = INT_THETA_CONCAT,
    batch_var_concat  = INT_BATCH_VAR_CONCAT,
    res_ls_concat     = INT_RES_LS_CONCAT,
    use_paga          = INT_USE_PAGA,
    paga_cl_res       = INT_PAGA_CL_RES,
    int_use_gpu       = INT_USE_GPU
  resources:
    mem_mb      = 64 * 1024
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
    else
      echo "running on CPU"
    fi
    set -u

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
      $USE_GPU_FLAG
    """


rule join_build_h5ads_yaml:
  """Create symlinks and the joint h5ads YAML manifest."""
  input:
    h5ads_yaml_fs = H5ADS_YAML_FS
  output:
    joint_h5ads_yaml_f = joint_h5ads_yaml_f,
    h5ad_symlinks = [f"{h5ads_dir}/{bk}.h5ad" for bk in _JOIN_BATCH_KEYS]
  params:
    project_ids = " ".join(JOIN_PROJECT_IDS),
    h5ads_dir   = h5ads_dir
  log:
    f"{logs_dir}/join_build_h5ads_yaml_{JOIN_TAG}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_build_h5ads_yaml_{JOIN_TAG}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/scprocess_local.yaml'
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
    gtf_dt_f    = GENE_INFO_F,
    sel_res     = MKR_SEL_RES,
    min_cl_size = MKR_MIN_CL_SIZE,
    min_cells   = MKR_MIN_CELLS,
    batch_var   = "sample_id"
  threads: 8
  resources:
    mem_mb      = 64 * 1024
  log:
    f"{logs_dir}/join_marker_genes_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.log"
  benchmark:
    f"{benchmark_dir}/join_marker_genes_{JOIN_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.benchmark.txt"
  conda:
    '../envs/rlibs.yaml'
  shell: """
    exec &>> {log}
    Rscript -e "source('scripts/utils.R'); source('scripts/marker_genes.R'); calculate_marker_genes(
      integration_f = '{input.integration_f}',
      h5ads_yaml_f  = '{input.h5ads_yaml_f}',
      pb_f          = '{output.pb_f}',
      mkrs_f        = '{output.mkrs_f}',
      pb_hvgs_f     = '{output.pb_hvgs_f}',
      gtf_dt_f      = '{params.gtf_dt_f}',
      sel_res       = '{params.sel_res}',
      min_cl_size   =  {params.min_cl_size},
      min_cells     =  {params.min_cells},
      zoom          = 'True',
      batch_var     = '{params.batch_var}',
      n_cores       =  {threads})"
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
    genome_ref  = GENOME_REF,
    gsea_dir    = GSEA_DIR,
    min_cpm_go  = MKR_MIN_CPM_GO,
    max_zero_p  = MKR_MAX_ZERO_P,
    gsea_cut    = MKR_GSEA_CUT,
    not_ok_re   = MKR_NOT_OK_RE,
    gsea_var    = MKR_GSEA_VAR
  threads: 8
  resources:
    mem_mb      = 16 * 1024
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


if DO_LABEL:
  def _get_labeller_entry(labeller, model):
    matches = [e for e in LABELLER_PARAMS
      if e['labeller'] == labeller and e['model'] == model]
    if len(matches) != 1:
      raise ValueError(f"Expected exactly one labeller entry for {labeller}/{model}, got {len(matches)}")
    return matches[0]

  rule join_celltypist:
    """Run CellTypist on each batch h5ad for the join integration."""
    input:
      adata_f = f"{join_int_dir}/h5ads/{{batch}}.h5ad"
    output:
      pred_f  = temp(f"{join_lbl_dir}/tmp_labels_celltypist_model_{{model}}_{JOIN_TAG}_{DATE_STAMP}_{{batch}}.csv.gz")
    threads: 4
    resources:
      mem_mb = 16 * 1024
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
    resources:
      mem_mb = 16 * 1024
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


  rule join_merge_labels:
    """Aggregate per-batch CellTypist predictions by majority voting."""
    input:
      pred_fs       = lambda wildcards: expand(
        f"{join_lbl_dir}/tmp_labels_{{labeller}}_model_{{model}}_{JOIN_TAG}_{DATE_STAMP}_{{batch}}.csv.gz",
        batch=_JOIN_BATCH_KEYS, allow_missing=True),
      integration_f = joint_integration_f
    output:
      pred_out_f    = f"{join_lbl_dir}/labels_{{labeller}}_model_{{model}}_{JOIN_TAG}_{DATE_STAMP}.csv.gz"
    params:
      hi_res_cl   = lambda wildcards: _get_labeller_entry(wildcards.labeller, wildcards.model)['hi_res_cl'],
      min_cl_prop = lambda wildcards: _get_labeller_entry(wildcards.labeller, wildcards.model)['min_cl_prop']
    threads: 4
    resources:
      mem_mb = 16 * 1024
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


rule join_render_html:
  """Render Rmd report and HTML for the join integration."""
  input:
    integration_f = joint_integration_f,
    mkrs_f        = mkrs_f,
    pb_hvgs_f     = pb_hvgs_f,
    pb_f          = pb_f,
    fgsea_files   = [fgsea_bp_f, fgsea_cc_f, fgsea_mf_f] if DO_GSEA else [],
    label_files   = label_fs,
    xgb_files     = [f'{join_xgb_dir}/{_join_xgb_ref_tag}_predictions.csv.gz',
                     f'{join_xgb_dir}/{_join_xgb_ref_tag}_gene_importance.csv',
                     f'{join_xgb_dir}/{_join_xgb_ref_tag}_pseudobulk.h5ad'] if DO_TRAIN_XGB else []
  output:
    r_utils_f     = f"{code_dir}/utils.R",
    r_int_f       = f"{code_dir}/integration.R",
    r_mkr_f       = f"{code_dir}/marker_genes.R",
    r_fgsea_f     = f"{code_dir}/fgsea.R",
    r_lbl_f       = f"{code_dir}/label_celltypes.R",
    rmd_f         = rmd_f,
    html_f        = html_f
  params:
    your_name     = YOUR_NAME,
    affiliation   = AFFILIATION,
    join_name     = JOIN_NAME,
    join_tag      = JOIN_TAG,
    join_int_dir  = join_int_dir,
    join_mkr_dir  = join_mkr_dir,
    ref_txome     = GENOME_REF,
    mkr_sel_res   = MKR_SEL_RES,
    int_res_ls    = INT_RES_LS_STR,
    metadata_vars    = METADATA_VARS_STR,
    custom_mkr_names = MKR_CUSTOM_NAMES,
    custom_mkr_paths = MKR_CUSTOM_PATHS,
    label_f_ls       = ' '.join(label_fs),
    labeller_ls      = ' '.join(e['labeller']  for e in _lbl_cfg) if DO_LABEL else '',
    model_ls         = ' '.join(e['model']     for e in _lbl_cfg) if DO_LABEL else '',
    hi_res_cl_ls     = ' '.join(e['hi_res_cl'] for e in _lbl_cfg) if DO_LABEL else '',
    min_cl_prop_ls   = ' '.join(str(e['min_cl_prop']) for e in _lbl_cfg) if DO_LABEL else '',
    mkr_min_cpm_mkr  = MKR_MIN_CPM_MKR,
    mkr_min_cells    = MKR_MIN_CELLS,
    mkr_gsea_cut     = MKR_GSEA_CUT,
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
  resources:
    mem_mb = 16 * 1024
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

  rule join_train_xgboost:
    input:
      annots_f    = _join_xgb_cfg['annots_f'],
      cluster_csv = joint_integration_f,
      h5ads_yaml  = joint_h5ads_yaml_f,
    output:
      model_f  = f'{join_xgb_dir}/{_join_xgb_ref_tag}_xgboost_model.json',
      cls_f    = f'{join_xgb_dir}/{_join_xgb_ref_tag}_allowed_cls.csv',
      genes_f  = f'{join_xgb_dir}/{_join_xgb_ref_tag}_selected_genes.txt',
      imp_f    = f'{join_xgb_dir}/{_join_xgb_ref_tag}_gene_importance.csv',
      preds_f  = f'{join_xgb_dir}/{_join_xgb_ref_tag}_predictions.csv.gz',
      pb_f     = f'{join_xgb_dir}/{_join_xgb_ref_tag}_pseudobulk.h5ad',
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
    resources:
      mem_mb  = 32 * 1024,
      runtime = 180
    log:
      f"{logs_dir}/join_train_xgboost_{_join_xgb_ref_tag}_{DATE_STAMP}.log"
    benchmark:
      f"{benchmark_dir}/join_train_xgboost_{_join_xgb_ref_tag}_{DATE_STAMP}.benchmark.txt"
    conda:
      '../envs/label_celltypes.yaml'
    shell: """
      exec &>> {log}
      python3 scripts/train_xgboost.py \
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

  rule train_xgboost:
    input:
      _join_xgb_model_f
