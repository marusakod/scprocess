# load modules
import os
import sys
import json
import pathlib
import yaml

sys.path.append('scripts')
from scprocess_utils import *

# define some things
scprocess_dir = pathlib.Path(config.pop('scprocess_dir'))
scdata_dir    = pathlib.Path(os.getenv('SCPROCESS_DATA_DIR'))
schema_f      = scprocess_dir / "resources/schemas/config.schema.json"
zoom_schema_f = scprocess_dir / "resources/schemas/zoom.schema.json"
join_schema_f = scprocess_dir / "resources/schemas/join.schema.json"

# detect config type: project (normal) or join
_is_join = 'join' in config

# references that support GSEA (ref_txome for poly-A, probe_set for Flex)
_GSEA_REFS = {'human_2024', 'human_2020', 'mouse_2024', 'mouse_2020',
  'human_v1', 'human_v2', 'mouse_v1', 'mouse_v2'}

if _is_join:
  # --- validate and unpack join config ---
  config = check_join_config(config, join_schema_f, scdata_dir)
  RESOURCE_PARAMS = prep_resource_params(config, join_schema_f, scprocess_dir)

  JOIN_NAME  = config['join']['name']
  PROJ_DIR   = pathlib.Path(config['join']['proj_dir'])
  DATE_STAMP = config['join']['date_stamp']
  GENOME_REF = config['join'].get('probe_set', config['join'].get('ref_txome', ''))
  JOIN_TAG   = f"{JOIN_NAME}_join"
  FULL_TAG   = JOIN_TAG
  SHORT_TAG  = JOIN_NAME
  MKR_SEL_RES = config.get('marker_genes', {}).get('mkr_sel_res', 0.2)

  int_dir   = str(PROJ_DIR / f"output/{JOIN_TAG}")
  mkr_dir   = str(PROJ_DIR / f"output/{JOIN_TAG}")
  zoom_dir  = ''   # not used for join
  logs_dir  = str(PROJ_DIR / ".log")
  docs_dir  = str(PROJ_DIR / "public")

  _shiny_cfg      = config.get('shiny', {})
  _sample_meta_f  = str(PROJ_DIR / f"output/{JOIN_TAG}/joint_sample_meta_{JOIN_TAG}_{DATE_STAMP}.csv")
  _metadata_vars_ls = _shiny_cfg.get('metadata_vars', config['join'].get('metadata_vars', []))
  _metadata_vars  = ','.join(_metadata_vars_ls)

  ZOOM_PARAMS = {}
  ZOOMS       = []

  def _gsea_supported():
    return GENOME_REF in _GSEA_REFS

else:
  # --- validate and unpack normal project config ---
  config = check_config(config, schema_f, scdata_dir, scprocess_dir)
  RESOURCE_PARAMS = prep_resource_params(config, schema_f, scprocess_dir)

  PROJ_DIR   = config['project']['proj_dir']
  FULL_TAG   = config['project']['full_tag']
  SHORT_TAG  = config['project']['short_tag']
  DATE_STAMP = config['project']['date_stamp']
  GENOME_REF = config['project'].get('probe_set', config['project'].get('ref_txome', ''))
  MKR_SEL_RES = config['marker_genes']['mkr_sel_res']

  int_dir   = f"{PROJ_DIR}/output/{SHORT_TAG}_integration"
  mkr_dir   = f"{PROJ_DIR}/output/{SHORT_TAG}_marker_genes"
  zoom_dir  = f"{PROJ_DIR}/output/{SHORT_TAG}_zoom"
  logs_dir  = f"{PROJ_DIR}/.log"
  docs_dir  = f"{PROJ_DIR}/public"

  _shiny_cfg     = config.get('shiny', {})
  _sample_meta_f = config['project']['sample_metadata']
  _metadata_vars_ls = _shiny_cfg.get('metadata_vars', config['project'].get('metadata_vars', []))
  _metadata_vars = ','.join(_metadata_vars_ls)

  ZOOM_PARAMS = get_zoom_parameters(config, zoom_schema_f, scdata_dir)
  ZOOMS       = list(ZOOM_PARAMS.keys())

  def _gsea_supported():
    return GENOME_REF in _GSEA_REFS

def _main_has_gsea():
  return config.get('marker_genes', {}).get('mkr_do_gsea', True) and _gsea_supported()

def _zoom_has_gsea(zoom_name):
  return (ZOOM_PARAMS[zoom_name]['marker_genes'].get('mkr_do_gsea', True)
          and _gsea_supported())

def _resolve_optional_path(val, proj_dir):
  """Resolve val relative to proj_dir if not absolute; return '' if unset."""
  if not val:
    return ''
  p = pathlib.Path(val)
  if not p.is_absolute():
    p = pathlib.Path(proj_dir) / p
  return str(p)

def _resolve_var_names(metadata_vars_ls, shiny_cfg):
  """Build display names list from metadata_labels (dict) or var_names (legacy array)."""
  labels = shiny_cfg.get('metadata_labels', {})
  if labels:
    return [labels.get(v, v) for v in metadata_vars_ls]
  if 'var_names' in shiny_cfg:
    return shiny_cfg['var_names']
  return list(metadata_vars_ls)

_home_md_f        = _resolve_optional_path(_shiny_cfg.get('home_md'),        PROJ_DIR)
_annotation_csv_f = _resolve_optional_path(_shiny_cfg.get('annotation_csv'), PROJ_DIR)
_main_app_tag     = get_shiny_app_tag(config)
_main_shiny_dir   = f'{docs_dir}/shiny_{_main_app_tag}'


# ---- helper: fgsea inputs for main rule (conditional) --------------------
def _main_fgsea_inputs(wildcards):
  if _main_has_gsea():
    return {
      'fgsea_go_bp_f': f'{mkr_dir}/fgsea_{FULL_TAG}_{MKR_SEL_RES}_go_bp_{DATE_STAMP}.csv.gz',
      'fgsea_go_cc_f': f'{mkr_dir}/fgsea_{FULL_TAG}_{MKR_SEL_RES}_go_cc_{DATE_STAMP}.csv.gz',
      'fgsea_go_mf_f': f'{mkr_dir}/fgsea_{FULL_TAG}_{MKR_SEL_RES}_go_mf_{DATE_STAMP}.csv.gz',
    }
  return {}


# ---- helper: fgsea inputs for zoom rule (conditional) --------------------
def _zoom_fgsea_inputs(wildcards):
  zoom_name   = wildcards.zoom_name
  mkr_sel_res = ZOOM_PARAMS[zoom_name]['marker_genes']['mkr_sel_res']
  if _zoom_has_gsea(zoom_name):
    return {
      'fgsea_go_bp_f': f'{zoom_dir}/{zoom_name}/fgsea_{FULL_TAG}_{zoom_name}_{mkr_sel_res}_go_bp_{DATE_STAMP}.csv.gz',
      'fgsea_go_cc_f': f'{zoom_dir}/{zoom_name}/fgsea_{FULL_TAG}_{zoom_name}_{mkr_sel_res}_go_cc_{DATE_STAMP}.csv.gz',
      'fgsea_go_mf_f': f'{zoom_dir}/{zoom_name}/fgsea_{FULL_TAG}_{zoom_name}_{mkr_sel_res}_go_mf_{DATE_STAMP}.csv.gz',
    }
  return {}


# ---- helper: resolve per-zoom shiny optional path -------------------------
def _zoom_optional_path(zoom_name, key):
  val = ZOOM_PARAMS[zoom_name].get('shiny', {}).get(key, '')
  return _resolve_optional_path(val, PROJ_DIR)


def _zoom_metadata_vars_ls(zoom_name):
  return ZOOM_PARAMS[zoom_name].get('shiny', {}).get('metadata_vars', _metadata_vars_ls)


# ---- aggregate rule: build all zoom shiny apps ----------------------------
rule build_all_zoom_shiny_apps:
  input:
    expand(f'{docs_dir}/shiny_zoom_{{zoom_name}}/.shiny_built_{DATE_STAMP}',
           zoom_name = ZOOMS)


rule build_shiny_app:
  input:
    unpack(_main_fgsea_inputs),
    **({'home_md_f':        _home_md_f}        if _home_md_f        else {}),
    **({'annotation_csv_f': _annotation_csv_f} if _annotation_csv_f else {}),
    h5ads_yaml_f  = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml',
    sample_meta_f = _sample_meta_f,
    integration_f = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    mkrs_f        = f'{mkr_dir}/pb_marker_genes_{FULL_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv.gz',
    pb_hvgs_f     = f'{mkr_dir}/pb_hvgs_{FULL_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv.gz',
  output:
    sentinel_f    = f'{_main_shiny_dir}/.shiny_built_{DATE_STAMP}'
  params:
    scprocess_dir = str(scprocess_dir),
    deploy_dir    = _main_shiny_dir,
    date_stamp    = DATE_STAMP,
    app_tag       = _main_app_tag,
    mkr_sel_res   = MKR_SEL_RES,
    ref_txome     = GENOME_REF,
    metadata_vars = _metadata_vars,
    app_title     = _shiny_cfg.get('app_title', SHORT_TAG),
    email         = _shiny_cfg.get('email', ''),
    keyword       = _shiny_cfg.get('keyword', 'cells'),
    default_gene  = _shiny_cfg.get('default_gene', ''),
    n_keep        = int(_shiny_cfg.get('n_keep', 30000)),
    var_names     = ','.join(_resolve_var_names(_metadata_vars_ls, _shiny_cfg)),
    metadata_combns        = json.dumps(_shiny_cfg.get('metadata_combns', [])).replace('"', '\\"'),
    home_md_f         = _home_md_f,
    annotation_csv_f  = _annotation_csv_f,
    cluster_palette   = _shiny_cfg.get('cluster_palette', ''),
    metadata_palettes = json.dumps(_shiny_cfg.get('metadata_palettes', {})).replace('"', '\\"'),
    fgsea_go_bp_f     = lambda wildcards, input: getattr(input, 'fgsea_go_bp_f', ''),
    fgsea_go_cc_f     = lambda wildcards, input: getattr(input, 'fgsea_go_cc_f', ''),
    fgsea_go_mf_f     = lambda wildcards, input: getattr(input, 'fgsea_go_mf_f', ''),
  threads: 4
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'build_shiny_app', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'build_shiny_app', 'time', attempt)
  conda: '../envs/shiny.yaml'
  log:   f'{logs_dir}/shiny/build_shiny_app_{_main_app_tag}_{DATE_STAMP}.log'
  shell: """
    exec &>> {log}
    mkdir -p {params.deploy_dir}
    mkdir -p $(dirname {log})

    Rscript --vanilla -e "
      source('scripts/shiny.R')
      make_shiny_app_scprocess(
        integration_f = '{input.integration_f}',
        h5ads_yaml_f  = '{input.h5ads_yaml_f}',
        sample_meta_f = '{input.sample_meta_f}',
        mkrs_f        = '{input.mkrs_f}',
        pb_hvgs_f     = '{input.pb_hvgs_f}',
        fgsea_bp_f    = '{params.fgsea_go_bp_f}',
        fgsea_cc_f    = '{params.fgsea_go_cc_f}',
        fgsea_mf_f    = '{params.fgsea_go_mf_f}',
        deploy_dir    = '{params.deploy_dir}',
        scprocess_dir = '{params.scprocess_dir}',
        app_tag       = '{params.app_tag}',
        date_stamp    = '{params.date_stamp}',
        mkr_sel_res   = '{params.mkr_sel_res}',
        ref_txome     = '{params.ref_txome}',
        metadata_vars = '{params.metadata_vars}',
        app_title     = '{params.app_title}',
        email         = '{params.email}',
        keyword       = '{params.keyword}',
        default_gene  = '{params.default_gene}',
        n_keep        = {params.n_keep},
        var_names     = '{params.var_names}',
        metadata_combns        = '{params.metadata_combns}',
        home_md_f         = '{params.home_md_f}',
        annotation_csv_f  = '{params.annotation_csv_f}',
        cluster_palette   = '{params.cluster_palette}',
        metadata_palettes = '{params.metadata_palettes}',
        n_cores           = {threads}
      )
    "
    touch {output.sentinel_f}
  """


rule build_zoom_shiny_app:
  input:
    unpack(_zoom_fgsea_inputs),
    h5ads_yaml_f  = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml',
    sample_meta_f = _sample_meta_f,
    integration_f = lambda wc: f'{zoom_dir}/{wc.zoom_name}/integrated_dt_{FULL_TAG}_{wc.zoom_name}_{DATE_STAMP}.csv.gz',
    mkrs_f        = lambda wc: (f'{zoom_dir}/{wc.zoom_name}/pb_marker_genes_{FULL_TAG}_{wc.zoom_name}'
                                f'_{ZOOM_PARAMS[wc.zoom_name]["marker_genes"]["mkr_sel_res"]}_{DATE_STAMP}.csv.gz'),
    pb_hvgs_f     = lambda wc: (f'{zoom_dir}/{wc.zoom_name}/pb_hvgs_{FULL_TAG}_{wc.zoom_name}'
                                f'_{ZOOM_PARAMS[wc.zoom_name]["marker_genes"]["mkr_sel_res"]}_{DATE_STAMP}.csv.gz'),
  output:
    sentinel_f    = f'{docs_dir}/shiny_zoom_{{zoom_name}}/.shiny_built_{DATE_STAMP}'
  params:
    scprocess_dir    = str(scprocess_dir),
    deploy_dir       = lambda wc: f'{docs_dir}/shiny_zoom_{wc.zoom_name}',
    date_stamp       = DATE_STAMP,
    app_tag          = lambda wc: f'{SHORT_TAG}_{wc.zoom_name}',
    mkr_sel_res      = lambda wc: ZOOM_PARAMS[wc.zoom_name]['marker_genes']['mkr_sel_res'],
    ref_txome        = GENOME_REF,
    app_title        = lambda wc: ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('app_title',
                         f'{SHORT_TAG} — {wc.zoom_name}'),
    email            = lambda wc: ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('email', ''),
    keyword          = lambda wc: ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('keyword', 'cells'),
    default_gene     = lambda wc: ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('default_gene', ''),
    n_keep           = lambda wc: int(ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('n_keep', 30000)),
    metadata_vars    = lambda wc: ','.join(_zoom_metadata_vars_ls(wc.zoom_name)),
    var_names        = lambda wc: ','.join(
                         _resolve_var_names(_zoom_metadata_vars_ls(wc.zoom_name),
                           ZOOM_PARAMS[wc.zoom_name].get('shiny', {}))),
    metadata_combns       = lambda wc: json.dumps(
                         ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('metadata_combns', [])).replace('"', '\\"'),
    home_md_f        = lambda wc: _zoom_optional_path(wc.zoom_name, 'home_md'),
    annotation_csv_f = lambda wc: _zoom_optional_path(wc.zoom_name, 'annotation_csv'),
    cluster_palette  = lambda wc: ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('cluster_palette', ''),
    metadata_palettes = lambda wc: json.dumps(
                          ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('metadata_palettes', {})).replace('"', '\\"'),
    fgsea_go_bp_f    = lambda wildcards, input: getattr(input, 'fgsea_go_bp_f', ''),
    fgsea_go_cc_f    = lambda wildcards, input: getattr(input, 'fgsea_go_cc_f', ''),
    fgsea_go_mf_f    = lambda wildcards, input: getattr(input, 'fgsea_go_mf_f', ''),
  threads: 4
  resources:
    mem_mb  = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'build_zoom_shiny_app', 'memory', attempt),
    runtime = lambda wildcards, attempt, input: get_resources(RESOURCE_PARAMS, rules, input, 'build_zoom_shiny_app', 'time', attempt)
  conda: '../envs/shiny.yaml'
  log:   f'{logs_dir}/shiny/build_zoom_shiny_app_{{zoom_name}}_{DATE_STAMP}.log'
  shell: """
    exec &>> {log}
    mkdir -p {params.deploy_dir}
    mkdir -p $(dirname {log})

    Rscript --vanilla -e "
      source('scripts/shiny.R')
      make_shiny_app_scprocess(
        integration_f = '{input.integration_f}',
        h5ads_yaml_f  = '{input.h5ads_yaml_f}',
        sample_meta_f = '{input.sample_meta_f}',
        mkrs_f        = '{input.mkrs_f}',
        pb_hvgs_f     = '{input.pb_hvgs_f}',
        fgsea_bp_f    = '{params.fgsea_go_bp_f}',
        fgsea_cc_f    = '{params.fgsea_go_cc_f}',
        fgsea_mf_f    = '{params.fgsea_go_mf_f}',
        deploy_dir    = '{params.deploy_dir}',
        scprocess_dir = '{params.scprocess_dir}',
        app_tag       = '{params.app_tag}',
        date_stamp    = '{params.date_stamp}',
        mkr_sel_res   = '{params.mkr_sel_res}',
        ref_txome     = '{params.ref_txome}',
        metadata_vars = '{params.metadata_vars}',
        app_title     = '{params.app_title}',
        email         = '{params.email}',
        keyword       = '{params.keyword}',
        default_gene  = '{params.default_gene}',
        n_keep        = {params.n_keep},
        var_names     = '{params.var_names}',
        metadata_combns        = '{params.metadata_combns}',
        home_md_f         = '{params.home_md_f}',
        annotation_csv_f  = '{params.annotation_csv_f}',
        cluster_palette   = '{params.cluster_palette}',
        metadata_palettes = '{params.metadata_palettes}',
        n_cores           = {threads}
      )
    "
    touch {output.sentinel_f}
  """
