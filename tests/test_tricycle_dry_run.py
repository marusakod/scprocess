"""Exercise origin arguments in the real tricycle Snakemake rules."""
import json
import os
from pathlib import Path
import subprocess
import sys
import pytest

from scripts.scprocess_utils import _get_default_config_from_schema


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('setting', [{'cyc_density_percentile': 0}, {'cyc_density_percentile': 1.3}, {'cyc_density_percentile': 100},
  {'cyc_density_percentile': 1.3, 'cyc_origin_manual': [-3, 1]}])
@pytest.mark.parametrize('target', ['origin', 'standalone_default', 'parent_default'])
def test_tricycle_origin_arguments(tmp_path, setting, target):
  schema = json.loads((ROOT / 'resources/schemas/config.schema.json').read_text())
  config = _get_default_config_from_schema(schema, {'cell_cycle': setting})
  config['cell_cycle'].update(setting, _species='mouse')
  config['resources']['retries'] = 0
  for name in ('counts.h5', 'qc/coldata_dt_all_cells_test_date.csv.gz',
               'qc/rowdata_dt_test_date.csv.gz', 'hvg/top_hvgs_counts_test_date.h5',
               'hvg/top_hvgs_doublet_counts_test_date.h5'):
    path = tmp_path / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()
  snakefile = tmp_path / 'Snakefile'
  parent_rule = '''rule all:
  input: "cycle/tricycle_origin_test_date.csv"
''' if target == 'parent_default' else ''
  rule_path = str(ROOT / 'rules/tricycle.smk')
  include_rule = (
    f'workflow.include({rule_path!r}, overwrite_default_target=True)'
    if target == 'standalone_default' else f'include: {rule_path!r}')
  snakefile.write_text(f'''
config = {config!r}
scprocess_dir = {str(ROOT)!r}
FULL_TAG, DATE_STAMP = "test", "date"
cell_cycle_dir, qc_dir, hvg_dir = "cycle", "qc", "hvg"
logs_dir, benchmark_dir, amb_dir = "logs", "benchmarks", "ambient"
RUN_VAR, RUNS = "pool_id", ["pool1", "pool2"]
RESOURCE_PARAMS = {{}}
def get_filtered_counts_file(*args):
    return "counts.h5"
def get_resources(*args):
    return 1000
{parent_rule}
{include_rule}
''')
  env = {k: v for k, v in os.environ.items() if not k.startswith('SNAKEMAKE_')}
  env['XDG_CACHE_HOME'] = str(tmp_path / 'cache')
  result = subprocess.run([
    sys.executable, '-m', 'snakemake', '--snakefile', str(snakefile),
    '--directory', str(tmp_path), '--dry-run', '--printshellcmds',
    '--cores', '2', '--executor', 'local', '--profile', 'none',
    '--workflow-profile', 'none', '--nolock',
  ] + (['cycle/tricycle_origin_test_date.csv'] if target == 'origin' else []),
    text=True, capture_output=True, timeout=45, env=env)
  output = result.stdout + result.stderr
  assert result.returncode == 0, output
  assert '--density_percentile ' + str(setting['cyc_density_percentile']) in output
  assert ('estimate_cell_cycle_blob(' in output) == (target == 'standalone_default')
  assert 'prepare_tricycle' not in output
  assert output.count('run_column = "pool_id"') == 2
  if 'cyc_origin_manual' in setting:
    assert '--origin_manual -3 1' in output
  assert '--origin_line_fraction' not in output
  assert '--density_percentile_grid 100' in output
  assert '--target_prop' not in output
