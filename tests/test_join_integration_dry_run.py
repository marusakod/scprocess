"""Exercise real join PCA/integration rules with isolated placeholder inputs."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
from scripts.scprocess_utils import _check_join_integration_parameters

ROOT = Path(__file__).resolve().parents[1]

@pytest.mark.parametrize('correction', [None, 'within_sample', 'global'])
def test_join_correction_dag_and_shell(tmp_path, correction):
    integration = {'int_pca_method': 'bpcells', 'int_embedding': 'harmony',
      'int_n_dims': 6, 'int_umap_n_iter': 2000, 'int_cl_method': 'leiden',
      'int_theta': [0.1, 0.1], 'int_batch_var': ['sample_id', 'project_id'],
      'int_res_ls': [0.5, 2], 'int_use_paga': False, 'int_use_gpu': True}
    if correction:
        integration['int_cell_cycle_correction'] = correction
    config = {'integration': integration, 'resources': {'retries': 0}}
    _check_join_integration_parameters(config)
    rules = (ROOT / 'rules/join.smk').read_text()
    # Use the actual production path/argument definitions and rule bodies.
    paths = rules[rules.index('# output paths'):rules.index('pb_dir      =')]
    bodies = rules[rules.index('rule join_pca:'):rules.index('rule join_build_h5ads_yaml:')]
    (tmp_path / 'counts.bpcells').mkdir()
    (tmp_path / 'coldata.csv.gz').touch()
    (tmp_path / 'counts.h5').touch()
    snakefile = tmp_path / 'Snakefile'
    snakefile.write_text(f'''config = {config!r}
scprocess_dir = {str(ROOT)!r}
join_int_dir, JOIN_TAG, DATE_STAMP = 'output', 'test', 'date'
INT_N_DIMS, INT_PCA_METHOD = 6, 'bpcells'
INT_CYCLE_MODE = {correction!r}
DO_CYCLE_CORRECTION = {correction is not None!r}
JOIN_ZOOM_CONTEXT, JOIN_RETAIN_METADATA = {{}}, ''
logs_dir, benchmark_dir = 'logs', 'benchmarks'
RESOURCE_PARAMS = {{}}
def get_resources(*args): return 1000
def get_paga_resolution(*args): return None
{paths}
joint_bpcells_counts_dir = 'counts.bpcells'
joint_coldata_f = 'coldata.csv.gz'
joint_label_counts_f = 'counts.h5'
# Rebind the score path after the placeholder metadata path is chosen.
JOIN_CYCLE_ARGS = JOIN_CYCLE_ARGS.replace('output/joint_coldata_test_date.csv.gz', joint_coldata_f)
{bodies}
''')
    env = {k: v for k, v in os.environ.items() if not k.startswith('SNAKEMAKE_')}
    env['XDG_CACHE_HOME'] = str(tmp_path / 'cache')
    run = subprocess.run([sys.executable, '-m', 'snakemake', '--snakefile', str(snakefile),
      '--directory', str(tmp_path), '--dry-run', '--printshellcmds', '--cores', '2',
      '--executor', 'local', '--profile', 'none', '--workflow-profile', 'none',
      '--nolock', 'output/integrated_dt_test_date.csv.gz'],
      text=True, capture_output=True, timeout=45, env=env)
    output = run.stdout + run.stderr
    assert run.returncode == 0, output
    assert 'join_pca' in output and 'join_integration' in output
    assert '--precomputed_pca_f' in output
    if correction:
        assert f"correction_centering = '{correction}'" in output
        assert "tricycle_f = 'coldata.csv.gz'" in output
        assert 'cell_cycle_regression_genes_test_date.csv.gz' in output
    else:
        assert 'correction_centering =' not in output
        assert 'cell_cycle_regression_genes_test_date.csv.gz' not in output
