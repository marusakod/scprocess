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

# validate and populate config with defaults
config = check_config(config, schema_f, scdata_dir, scprocess_dir)

# unpack frequently used variables
PROJ_DIR   = config['project']['proj_dir']
FULL_TAG   = config['project']['full_tag']
SHORT_TAG  = config['project']['short_tag']
DATE_STAMP = config['project']['date_stamp']

# pipeline output directories (read-only inputs to this rule)
int_dir    = f"{PROJ_DIR}/output/{SHORT_TAG}_integration"
mkr_dir    = f"{PROJ_DIR}/output/{SHORT_TAG}_marker_genes"
zoom_dir   = f"{PROJ_DIR}/output/{SHORT_TAG}_zoom"
logs_dir   = f"{PROJ_DIR}/.log"
docs_dir   = f"{PROJ_DIR}/public"

# marker gene resolution (needed to construct file paths)
MKR_SEL_RES = config['marker_genes']['mkr_sel_res']

# optional shiny config section
_shiny_cfg  = config.get('shiny', {})

# load zoom parameters (empty dict if no zoom config)
ZOOM_PARAMS = get_zoom_parameters(config, zoom_schema_f, scdata_dir)
ZOOMS       = list(ZOOM_PARAMS.keys())

# reference transcriptomes that support GSEA
_GSEA_TXOMES = {'human_2024', 'human_2020', 'mouse_2024', 'mouse_2020'}

def _gsea_supported():
  """True if the project ref_txome supports GSEA."""
  return config['project']['ref_txome'] in _GSEA_TXOMES

def _main_has_gsea():
  return config['marker_genes'].get('mkr_do_gsea', True) and _gsea_supported()

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

_home_md_f        = _resolve_optional_path(_shiny_cfg.get('home_md'),        PROJ_DIR)
_annotation_csv_f = _resolve_optional_path(_shiny_cfg.get('annotation_csv'), PROJ_DIR)


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
    integration_f = f'{int_dir}/integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz',
    mkrs_f        = f'{mkr_dir}/pb_marker_genes_{FULL_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv.gz',
    pb_hvgs_f     = f'{mkr_dir}/pb_hvgs_{FULL_TAG}_{MKR_SEL_RES}_{DATE_STAMP}.csv.gz',
  output:
    sentinel_f    = f'{docs_dir}/shiny/.shiny_built_{DATE_STAMP}'
  params:
    scprocess_dir = str(scprocess_dir),
    deploy_dir    = f'{docs_dir}/shiny',
    sample_meta_f = config['project']['sample_metadata'],
    date_stamp    = DATE_STAMP,
    app_tag       = SHORT_TAG,
    mkr_sel_res   = MKR_SEL_RES,
    ref_txome     = config['project']['ref_txome'],
    metadata_vars = ','.join(config['project'].get('metadata_vars', [])),
    app_title     = _shiny_cfg.get('app_title', SHORT_TAG),
    email         = _shiny_cfg.get('email', ''),
    keyword       = _shiny_cfg.get('keyword', 'cells'),
    default_gene  = _shiny_cfg.get('default_gene', ''),
    n_keep        = int(_shiny_cfg.get('n_keep', 30000)),
    var_names     = ','.join(_shiny_cfg.get('var_names',
                      config['project'].get('metadata_vars', []))),
    var_combns        = json.dumps(_shiny_cfg.get('var_combns', [])),
    home_md_f         = _home_md_f,
    annotation_csv_f  = _annotation_csv_f,
    fgsea_go_bp_f     = lambda wildcards, input: getattr(input, 'fgsea_go_bp_f', ''),
    fgsea_go_cc_f     = lambda wildcards, input: getattr(input, 'fgsea_go_cc_f', ''),
    fgsea_go_mf_f     = lambda wildcards, input: getattr(input, 'fgsea_go_mf_f', ''),
  threads: 4
  resources:
    mem_mb  = 64 * MB_PER_GB,
    runtime = 30
  conda: '../envs/shiny.yaml'
  log:   f'{logs_dir}/shiny/build_shiny_app_{DATE_STAMP}.log'
  shell: """
    exec &>> {log}
    mkdir -p {params.deploy_dir}
    mkdir -p $(dirname {log})

    Rscript --vanilla -e "
      source('scripts/shiny.R')
      make_shiny_app_scprocess(
        integration_f = '{input.integration_f}',
        h5ads_yaml_f  = '{input.h5ads_yaml_f}',
        sample_meta_f = '{params.sample_meta_f}',
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
        var_combns       = '{params.var_combns}',
        home_md_f        = '{params.home_md_f}',
        annotation_csv_f = '{params.annotation_csv_f}',
        n_cores          = {threads}
      )
    "
    touch {output.sentinel_f}
  """


rule build_zoom_shiny_app:
  input:
    unpack(_zoom_fgsea_inputs),
    h5ads_yaml_f  = f'{int_dir}/h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml',
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
    sample_meta_f    = config['project']['sample_metadata'],
    date_stamp       = DATE_STAMP,
    app_tag          = lambda wc: f'{SHORT_TAG}_{wc.zoom_name}',
    mkr_sel_res      = lambda wc: ZOOM_PARAMS[wc.zoom_name]['marker_genes']['mkr_sel_res'],
    ref_txome        = config['project']['ref_txome'],
    metadata_vars    = ','.join(config['project'].get('metadata_vars', [])),
    app_title        = lambda wc: ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('app_title',
                         f'{SHORT_TAG} — {wc.zoom_name}'),
    email            = lambda wc: ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('email', ''),
    keyword          = lambda wc: ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('keyword', 'cells'),
    default_gene     = lambda wc: ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('default_gene', ''),
    n_keep           = lambda wc: int(ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('n_keep', 30000)),
    var_names        = lambda wc: ','.join(
                         ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('var_names',
                         config['project'].get('metadata_vars', []))),
    var_combns       = lambda wc: json.dumps(
                         ZOOM_PARAMS[wc.zoom_name].get('shiny', {}).get('var_combns', [])),
    home_md_f        = lambda wc: _zoom_optional_path(wc.zoom_name, 'home_md'),
    annotation_csv_f = lambda wc: _zoom_optional_path(wc.zoom_name, 'annotation_csv'),
    fgsea_go_bp_f    = lambda wildcards, input: getattr(input, 'fgsea_go_bp_f', ''),
    fgsea_go_cc_f    = lambda wildcards, input: getattr(input, 'fgsea_go_cc_f', ''),
    fgsea_go_mf_f    = lambda wildcards, input: getattr(input, 'fgsea_go_mf_f', ''),
  threads: 4
  resources:
    mem_mb  = 64 * MB_PER_GB,
    runtime = 30
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
        sample_meta_f = '{params.sample_meta_f}',
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
        var_combns       = '{params.var_combns}',
        home_md_f        = '{params.home_md_f}',
        annotation_csv_f = '{params.annotation_csv_f}',
        n_cores          = {threads}
      )
    "
    touch {output.sentinel_f}
  """
