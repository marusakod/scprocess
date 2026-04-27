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
import json
import yaml
import jsonschema
import polars as pl

sys.path.append('scripts')

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

scprocess_dir = pathlib.Path(config.pop('scprocess_dir'))
scdata_dir    = pathlib.Path(os.getenv('SCPROCESS_DATA_DIR'))
schema_f      = scprocess_dir / "resources/schemas/join.schema.json"

# validate join config
with open(schema_f) as _f:
  _join_schema = json.load(_f)
_validator = jsonschema.Draft202012Validator(_join_schema)
_errors = sorted(_validator.iter_errors(config), key=lambda e: e.path)
if _errors:
  raise ValueError("join.yaml validation errors:\n" +
    "\n".join(f"  {list(e.path)}: {e.message}" for e in _errors))

# apply schema defaults
def _apply_defaults(cfg, schema_props):
  for key, prop in schema_props.items():
    if key not in cfg and 'default' in prop:
      cfg[key] = prop['default']
    elif key in cfg and prop.get('type') == 'object' and 'properties' in prop:
      _apply_defaults(cfg[key], prop['properties'])

for section in ['hvg', 'integration', 'marker_genes']:
  if section not in config:
    config[section] = {}
  _apply_defaults(config[section], _join_schema['properties'].get(section, {}).get('properties', {}))

# ---------------------------------------------------------------------------
# Project config loading
# ---------------------------------------------------------------------------

JOIN_PROJECT_IDS = list(config['projects'].keys())
_project_cfgs    = {}

for _pid in JOIN_PROJECT_IDS:
  _cfg_f = config['projects'][_pid]['config']
  with open(_cfg_f) as _f:
    _project_cfgs[_pid] = yaml.safe_load(_f)
  _pref = _project_cfgs[_pid]['project']['ref_txome']
  if _pref != config['join']['ref_txome']:
    raise ValueError(
      f"Project '{_pid}' ref_txome={_pref!r} does not match join ref_txome={config['join']['ref_txome']!r}"
    )

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
  return _proj_hvg_dir(pid) / f"standardized_variance_stats_{_proj_full_tag(pid)}_{_proj_date(pid)}.csv.gz"

def _proj_h5ads_yaml_f(pid):
  return _proj_int_dir(pid) / f"h5ads_clean_paths_{_proj_full_tag(pid)}_{_proj_date(pid)}.yaml"

def _proj_integrated_dt_f(pid):
  return _proj_int_dir(pid) / f"integrated_dt_{_proj_full_tag(pid)}_{_proj_date(pid)}.csv.gz"

def _proj_sample_meta_f(pid):
  return pathlib.Path(_project_cfgs[pid]['project']['sample_metadata'])

VAR_STATS_FS   = [str(_proj_var_stats_f(pid)) for pid in JOIN_PROJECT_IDS]
H5ADS_YAML_FS  = [str(_proj_h5ads_yaml_f(pid)) for pid in JOIN_PROJECT_IDS]
INTEGRATED_FS  = [str(_proj_integrated_dt_f(pid)) for pid in JOIN_PROJECT_IDS]
SAMPLE_META_FS = [str(_proj_sample_meta_f(pid)) for pid in JOIN_PROJECT_IDS]

# ---------------------------------------------------------------------------
# Derived constants
# ---------------------------------------------------------------------------

JOIN_NAME  = config['join']['name']
JOIN_DIR   = pathlib.Path(config['join']['proj_dir'])
DATE_STAMP = config['join']['date_stamp']
REF_TXOME  = config['join']['ref_txome']

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

GSEA_TXOMES = ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020']
DO_GSEA     = MKR_DO_GSEA and (REF_TXOME in GSEA_TXOMES)

# HVG
N_HVGS = config.get('hvg', {}).get('hvg_n_hvgs', 2000)

# metadata vars (space-separated string for join.py)
METADATA_VARS_STR = " ".join(config['join'].get('metadata_vars', []))

# GTF file from index_parameters.csv (needed for marker genes)
_idx_params_f = scdata_dir / 'index_parameters.csv'
_idx_params   = pl.read_csv(_idx_params_f)
GTF_DT_F = _idx_params.filter(pl.col('ref_txome') == REF_TXOME)['gtf_txt_f'][0]
GSEA_DIR = str(scdata_dir / 'gmt_pathways')

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

rule all:
  input:
    joint_integration_f,
    joint_h5ads_yaml_f,
    mkrs_f,
    pb_hvgs_f,
    *([fgsea_bp_f, fgsea_cc_f, fgsea_mf_f] if DO_GSEA else []),
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
    embedding       = INT_EMBEDDING,
    n_dims          = INT_N_DIMS,
    cl_method       = INT_CL_METHOD,
    theta_concat    = INT_THETA_CONCAT,
    batch_var_concat = INT_BATCH_VAR_CONCAT,
    res_ls_concat   = INT_RES_LS_CONCAT,
    use_paga        = INT_USE_PAGA,
    paga_cl_res     = INT_PAGA_CL_RES,
    int_use_gpu     = INT_USE_GPU
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
    joint_h5ads_yaml_f = joint_h5ads_yaml_f
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
    gtf_dt_f    = GTF_DT_F,
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


if DO_GSEA:
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
        ref_txome     = '{params.ref_txome}',
        gsea_dir      = '{params.gsea_dir}',
        min_cpm_go    = {params.min_cpm_go},
        max_zero_p    = {params.max_zero_p},
        gsea_cut      = {params.gsea_cut},
        not_ok_re     = '{params.not_ok_re}',
        gsea_var      = '{params.gsea_var}',
        n_cores       =  {threads})"
      """


rule join_render_html:
  """Render Rmd report and HTML for the join integration."""
  input:
    integration_f = joint_integration_f,
    mkrs_f        = mkrs_f,
    pb_hvgs_f     = pb_hvgs_f,
    pb_f          = pb_f,
    fgsea_files   = [fgsea_bp_f, fgsea_cc_f, fgsea_mf_f] if DO_GSEA else []
  output:
    r_utils_f     = f"{code_dir}/utils.R",
    r_int_f       = f"{code_dir}/integration.R",
    r_mkr_f       = f"{code_dir}/marker_genes.R",
    r_fgsea_f     = f"{code_dir}/fgsea.R",
    rmd_f         = rmd_f,
    html_f        = html_f
  params:
    your_name     = YOUR_NAME,
    affiliation   = AFFILIATION,
    join_name     = JOIN_NAME,
    join_tag      = JOIN_TAG,
    join_int_dir  = join_int_dir,
    join_mkr_dir  = join_mkr_dir,
    ref_txome     = REF_TXOME,
    mkr_sel_res   = MKR_SEL_RES,
    int_res_ls    = INT_RES_LS_STR,
    metadata_vars = METADATA_VARS_STR,
    proj_dir      = str(JOIN_DIR),
    scprocess_dir = str(scprocess_dir),
    date_stamp    = DATE_STAMP
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
      metadata_vars = '{params.metadata_vars}',
      scprocess_dir = '{params.scprocess_dir}',
      date_stamp    = '{params.date_stamp}'
    )"
    """
