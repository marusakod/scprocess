"""Exercise origin arguments in the real tricycle Snakemake rules."""
import json
import os
from pathlib import Path
import subprocess
import sys
import pytest

from scripts.scprocess_utils import _get_default_config_from_schema


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('setting', [{'cyc_target_prop': 0.01}, {'cyc_target_prop': 1},
  {'cyc_target_prop': 0.01, 'cyc_origin_manual': [-3, 1]}])
def test_tricycle_origin_arguments(tmp_path, setting):
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
include: {str(ROOT / 'rules/tricycle.smk')!r}
''')
  env = {k: v for k, v in os.environ.items() if not k.startswith('SNAKEMAKE_')}
  env['XDG_CACHE_HOME'] = str(tmp_path / 'cache')
  result = subprocess.run([
    sys.executable, '-m', 'snakemake', '--snakefile', str(snakefile),
    '--directory', str(tmp_path), '--dry-run', '--printshellcmds',
    '--cores', '2', '--executor', 'local', '--profile', 'none',
    '--workflow-profile', 'none', '--nolock', 'cycle/tricycle_origin_test_date.csv',
  ], text=True, capture_output=True, timeout=45, env=env)
  output = result.stdout + result.stderr
  assert result.returncode == 0, output
  assert '--target_prop ' + str(setting['cyc_target_prop']) in output
  assert 'estimate_cell_cycle_blob(' not in output
  assert 'prepare_tricycle' not in output
  assert output.count('run_column = "pool_id"') == 2
  if 'cyc_origin_manual' in setting:
    assert '--origin_manual -3 1' in output
  assert '--origin_line_fraction' not in output
