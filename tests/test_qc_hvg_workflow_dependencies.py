import json
import os
import re
from pathlib import Path
import subprocess
import sys

import polars as pl
import pytest

from scripts.scprocess_utils import (
  _get_default_config_from_schema, prep_resource_params, get_resources,
  write_qc_thresholds, extract_qc_sample_statistics,
  _check_qc_parameters, _check_hvg_parameters,
)
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]


def run_snakemake(tmp_path, target, *options):
  env = {key: value for key, value in os.environ.items() if not key.startswith('SNAKEMAKE_')}
  env['XDG_CACHE_HOME'] = str(tmp_path / 'cache')
  result = subprocess.run([
    sys.executable, '-m', 'snakemake', '--snakefile', str(tmp_path / 'Snakefile'),
    '--directory', str(tmp_path), '--cores', '1', '--executor', 'local',
    '--profile', 'none', '--workflow-profile', 'none', '--nolock', target,
    *options], env=env, capture_output=True, text=True, timeout=45)
  output = result.stdout + result.stderr
  assert result.returncode == 0, output
  return output


def setup_rules(tmp_path, rule_file, method='sample', exclude=False, custom=False):
  schema = json.loads((ROOT / 'resources/schemas/config.schema.json').read_text())
  config = _get_default_config_from_schema(schema, {})
  config['project'] = {'sample_metadata': str(tmp_path / 'metadata.csv')}
  config['mapping_af'] = {'gene_info_f': str(tmp_path / 'genes.tsv'), 'af_mito_str': 'MT-'}
  config['multiplexing']['demux_type'] = 'none'
  config['ambient']['ambient_method'] = 'none'
  config['resources']['retries'] = 0
  config = _check_qc_parameters(config)
  config['hvg'].update(hvg_method=method, hvg_group_names=['g'],
                       hvg_num_chunks=1, hvg_metadata_split_var='group' if method == 'groups' else None)
  (tmp_path / 'metadata.csv').write_text('sample_id,group\ns1,g\ns2,g\n')
  (tmp_path / 'genes.tsv').write_text('gene_id\tsymbol\ng1\tG1\n')
  if exclude:
    config['hvg']['hvg_exclude_from_file'] = str(tmp_path / 'exclude.csv')
    (tmp_path / 'exclude.csv').write_text('gene_id\ng1\n')
  if custom:
    config['project']['custom_sample_params'] = str(tmp_path / 'custom.yaml')
    (tmp_path / 'custom.yaml').write_text('s1:\n  qc_min_counts: 20\n')
  fields = ['hvg/hvg_paths_test_date.csv', 'qc/coldata_dt_all_cells_test_date.csv.gz',
    'qc/rowdata_dt_test_date.csv.gz', 'qc/qc_sample_id_statistics_test_date.csv',
    'hvg/standardized_variance_stats_test_date.csv.gz']
  for run in ('s1','s2'):
    fields += [f'ambient/ambient_{run}/ambient_{run}_date_output_paths.yaml',
      f'ambient/ambient_{run}/uncorrected_{run}_date_filtered.h5',
      f'hvg/chunked_counts_{run}_test_date.h5', f'hvg/tmp_std_var_stats_{run}_test_date.csv.gz']
  fields += ['hvg/tmp_std_var_stats_g_chunk_0_test_date.csv.gz',
    'hvg/tmp_mean_var_g_chunk_0_test_date.csv.gz', 'hvg/estimated_variances_test_date.csv.gz',
    'hvg/means_variances_dt_test_date.csv.gz']
  for name in fields:
    path = tmp_path / name
    path.parent.mkdir(parents=True,exist_ok=True)
    path.touch()
  (tmp_path / 'Snakefile').write_text(f'''
import sys
sys.path.insert(0, {str(ROOT / 'scripts')!r})
from scprocess_utils import *
config = {config!r}
scprocess_dir = {str(ROOT)!r}
FULL_TAG, DATE_STAMP = 'test', 'date'
hvg_dir, qc_dir, amb_dir, empty_dir = 'hvg', 'qc', 'ambient', 'empty'
af_dir, af_rna_dir, dbl_dir, demux_dir = 'af', '', 'doublets', 'demux'
logs_dir, benchmark_dir = 'logs', 'benchmarks'
RUNS = BATCHES = ['s1', 's2']
RUN_VAR = BATCH_VAR = 'sample_id'
RUNS_TO_BATCHES = {{run: [run] for run in RUNS}}
BATCH_PARAMS = {{run: {{'qc': {{'qc_min_counts': 20, 'qc_max_counts': None}}}} for run in RUNS}}
CAN_CALC_AMBIENT_GENES = False
RESOURCE_PARAMS = {{}}
def get_resources(*args, **kwargs):
  return 1000
include: {str(ROOT / 'rules' / rule_file)!r}
''')
  return config


@pytest.mark.parametrize('method', ['sample','all','groups'])
@pytest.mark.parametrize('exclude', [False,True])
def test_hvg_optional_exclusion_and_metadata_inputs(tmp_path, method, exclude):
  setup_rules(tmp_path, 'hvgs.smk', method, exclude)
  target = 'hvg/hvg_dt_test_date.csv.gz'
  output = run_snakemake(tmp_path, target, '--dry-run', '--printshellcmds')
  flag_block = re.search(r'if \[ "([^"]*)" != "" \]; then\s*EXC_GS_F_FLAG="[^"]*"\s*fi', output)
  assert flag_block, output
  effective = subprocess.check_output([
    'bash', '-c', 'EXC_GS_F_FLAG=""; ' + flag_block[0] + '\nprintf "%s" "$EXC_GS_F_FLAG"'], text=True)
  assert effective == ('--exc_gs_f ' + str(tmp_path / 'exclude.csv') if exclude else '')
  if method != 'sample':
    # Force a statistics rule to expose its actual metadata file input.
    stats = run_snakemake(tmp_path, 'hvg/tmp_mean_var_g_chunk_0_test_date.csv.gz',
                          '--dry-run', '--force', '--printshellcmds')
    assert str(tmp_path / 'metadata.csv') in stats
  if exclude:
    (tmp_path / target).touch()
    run_snakemake(tmp_path, target, '--touch', '--force')
    assert 'Nothing to be done' in run_snakemake(tmp_path, target, '--dry-run')
    (tmp_path / 'exclude.csv').write_text('gene_id\ng2\n')
    rerun = run_snakemake(tmp_path, target, '--dry-run')
    assert 'rule get_highly_variable_genes:' in rerun


@pytest.mark.parametrize('custom',[False,True])
def test_qc_threshold_local_rule_and_custom_file_recomputation(tmp_path,custom):
  setup_rules(tmp_path, 'qc.smk',custom=custom)
  target = 'qc/qc_thresholds_by_sample_id_test_date.csv'
  run_snakemake(tmp_path,target)
  table = pl.read_csv(tmp_path / target)
  assert table['sample_id'].to_list() == ['s1','s2']
  assert table['min_counts'].to_list() == [20,20]
  assert table['max_counts'].null_count() == 2
  assert 'Nothing to be done' in run_snakemake(tmp_path,target,'--dry-run')
  if custom:
    (tmp_path / 'custom.yaml').write_text('s1:\n  qc_min_counts: 30\n')
    assert 'rule make_qc_thresholds_csv:' in run_snakemake(tmp_path,target,'--dry-run')


def test_moved_qc_statistics_preserves_counts_and_minimum_checks(tmp_path):
  qc = tmp_path / 'qc.csv'
  cuts = tmp_path / 'cuts.csv'
  pl.DataFrame({'sample_id':['s1','s1','s2'],'keep':[True,False,True]}).write_csv(qc)
  pl.DataFrame({'sample_id':['s1','s2','s3'],'min_cells':[1,2,1]}).write_csv(cuts)
  result = extract_qc_sample_statistics(None,qc,cuts,{'ambient':{'ambient_method':'none'}},
                                       ['s1','s2','s3'],{},'sample_id','sample_id')
  assert result.sort('sample_id')['n_cells'].to_list() == [1,1,0]
  assert result.sort('sample_id')['bad_sample_id'].to_list() == [False,True,True]


def test_hvg_validation_preserves_absent_exclusion_key(tmp_path):
  config = setup_rules(tmp_path, 'hvgs.smk')
  assert 'hvg_exclude_from_file' not in _check_hvg_parameters(config)['hvg']


def test_doublet_merge_resources_follow_schema_and_overrides():
  schema = ROOT / 'resources/schemas/config.schema.json'
  name = 'merge_doublet_diagnostics'
  rules = SimpleNamespace(**{name: object()})
  params = prep_resource_params({},schema,ROOT)
  assert get_resources(params,rules,SimpleNamespace(),name,'memory',1) == 4096
  assert get_resources(params,rules,SimpleNamespace(),name,'memory',2) == 6144
  assert get_resources(params,rules,SimpleNamespace(),name,'time',1) == 15
  params = prep_resource_params({'resources':{f'gb_{name}':1,f'mins_{name}':2}},schema,ROOT)
  assert get_resources(params,rules,SimpleNamespace(),name,'memory',1) == 1024
  assert get_resources(params,rules,SimpleNamespace(),name,'time',1) == 2
