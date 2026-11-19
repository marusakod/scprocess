"""Exercise the HTO / ADT rules of hto.smk in a Snakemake dry run."""
import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]

POOL_REFS = {
  'pool1': {'hto_ids': ['Hashtag-1', 'Hashtag-2'], 'adt_ids': ['CD3', 'CD19']},
  'pool2': {'hto_ids': ['Hashtag-2', 'Hashtag-3'], 'adt_ids': ['CD3', 'CD19']}}
NO_ADT_REFS = {pool: {'hto_ids': refs['hto_ids'], 'adt_ids': []} for pool, refs in POOL_REFS.items()}


def _dry_run(tmp_path, pool_refs, targets):
  # files the rules read or that earlier stages would have made
  for name in ('features.csv', 'samples.csv', 'af/af_pool1/chemistry_statistics.yaml',
               'af/af_pool2/chemistry_statistics.yaml', 'ambient/stats.csv', 'ambient/pool1.yaml',
               'ambient/pool2.yaml'):
    path = tmp_path / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()

  config = {
    'multiplexing': {'demux_type': 'hto', 'feature_ref': str(tmp_path / 'features.csv'), 'seurat_quantile': 0.99},
    'mapping_af': {'alevin_fry_home': '/af_home', 'wl_lu_f': '/wl_lu.csv'},
    'project': {'sample_metadata': str(tmp_path / 'samples.csv')},
    'ambient': {'ambient_method': 'none'},
    'resources': {'retries': 0, 'n_run_mapping': 2}}
  lib_params = {pool: {'multiplexing': {'where': 'local', 'R1_fs': ['r1.fq'], 'R2_fs': ['r2.fq'],
    'af_chemistry': '10xv3'}} for pool in pool_refs}

  snakefile = tmp_path / 'Snakefile'
  snakefile.write_text(f'''
config = {config!r}
scprocess_dir = {str(ROOT)!r}
FULL_TAG, DATE_STAMP, RUN_VAR = "test", "date", "pool_id"
af_dir, logs_dir, benchmark_dir, amb_dir, demux_dir = "af", "logs", "benchmarks", "ambient", "demux"
af_rna_dir = ""
HTO_POOL_REFS = {pool_refs!r}
LIB_PARAMS = {lib_params!r}
RESOURCE_PARAMS = {{}}
def get_resources(*args):
    return 1000
include: {str(ROOT / 'rules/hto.smk')!r}
''')
  env = {k: v for k, v in os.environ.items() if not k.startswith('SNAKEMAKE_')}
  env['XDG_CACHE_HOME'] = str(tmp_path / 'cache')
  result = subprocess.run([
    sys.executable, '-m', 'snakemake', '--snakefile', str(snakefile), '--directory', str(tmp_path),
    '--dry-run', '--printshellcmds', '--cores', '2', '--executor', 'local', '--profile', 'none',
    '--workflow-profile', 'none', '--nolock'] + targets, text=True, capture_output=True, timeout=60, env=env)
  output = result.stdout + result.stderr
  assert result.returncode == 0, output
  return output


def test_each_pool_gets_its_own_index_with_its_htos_and_all_adts(tmp_path):
  output = _dry_run(tmp_path, POOL_REFS, ['af/hto_ref/pool1/hto_index', 'af/hto_ref/pool2/hto_index'])
  assert output.count('mapping.py build_hto_index') == 2
  assert '--hto_ids       "Hashtag-1,Hashtag-2"' in output
  assert '--hto_ids       "Hashtag-2,Hashtag-3"' in output
  assert output.count('--adt_ids       "CD3,CD19"') == 2
  assert '--out_dir       "af/hto_ref/pool1"' in output
  assert '--out_dir       "af/hto_ref/pool2"' in output


def test_index_without_adts_has_an_empty_adt_list(tmp_path):
  output = _dry_run(tmp_path, NO_ADT_REFS, ['af/hto_ref/pool1/hto_index'])
  assert '--hto_ids       "Hashtag-1,Hashtag-2"' in output
  assert '--adt_ids       ""' in output


def test_hto_mapping_uses_the_pool_index_and_no_rna_orientation(tmp_path):
  output = _dry_run(tmp_path, POOL_REFS, ['af/af_pool1/hto/af_quant/'])
  assert 'hto_ref/pool1/hto_index' in output
  assert 'map_fastqs_to_counts pool1' in output
  assert '--what            "hto"' in output
  assert '--exp_ori' not in output
  assert 'selected_ori' not in output


def test_hto_counts_keep_only_the_hto_rows_of_the_pool_index(tmp_path):
  output = _dry_run(tmp_path, POOL_REFS, ['af/af_pool1/hto/af_hto_counts_mat.h5'])
  assert "keep_features = get_pool_feature_ids('af/hto_ref/pool1/features.csv', 'hto')" in output
  assert 'rule save_alevin_adt_to_h5:' not in output


def test_adt_counts_are_a_separate_target_from_the_adt_rows_of_the_index(tmp_path):
  output = _dry_run(tmp_path, POOL_REFS, ['af/af_pool1/hto/af_adt_counts_mat.h5'])
  assert 'save_alevin_adt_h5(' in output
  assert "features_f = 'af/hto_ref/pool1/features.csv'" in output
  assert 'rule save_alevin_adt_to_h5:' in output
  assert 'rule save_alevin_hto_to_h5:' not in output
