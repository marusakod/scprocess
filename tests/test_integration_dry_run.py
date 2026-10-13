"""Plan the real integration rules using isolated placeholder inputs only."""
import json
import os
from pathlib import Path
import re
import subprocess
import sys

import jsonschema
import pytest

from scripts.scprocess_utils import (
    _check_cell_cycle_parameters, _get_default_config_from_schema,
)


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('paga_resolution', [None, 2, 0.5])
@pytest.mark.parametrize('backend,correction,estimate_phase', [
    ('bpcells', 'within_sample', True), ('bpcells', 'global', True),
    ('bpcells', None, True), ('scanpy', None, True),
    ('bpcells', None, False), ('scanpy', None, False),
])
def test_integration_dag_and_effective_arguments(tmp_path, backend, correction, estimate_phase, paga_resolution):
    paga = paga_resolution is not None
    schema = json.loads((ROOT / 'resources/schemas/config.schema.json').read_text())
    config = _get_default_config_from_schema(schema, {})
    config['project'] = {'ref_txome': 'mouse_2024'}
    config['multiplexing']['demux_type'] = 'none'
    config['ambient']['ambient_method'] = 'none'
    config['resources']['retries'] = 0
    config['integration'].update(int_pca_method=backend, int_umap_n_iter=1234,
                                 int_use_gpu=False, int_use_paga=paga)
    if not paga:
        config['integration'].pop('int_paga_cl_res')
    else:
        config['integration']['int_paga_cl_res'] = paga_resolution
    if estimate_phase:
        config['cell_cycle'] = {'cyc_target_prop': 0.01}
    if correction:
        config['integration']['int_cell_cycle_correction'] = correction
    config = _check_cell_cycle_parameters(config)
    jsonschema.validate(config['integration'], schema['properties']['integration'])

    for name in (
        'hvg/top_hvgs_counts_test_date.h5',
        'hvg/top_hvgs_doublet_counts_test_date.h5',
        'qc/coldata_dt_all_cells_test_date.csv.gz',
        'qc/qc_sample_id_statistics_test_date.csv',
        'qc/qc_passed_test_date.flag', 'cycle/scores.csv.gz',
    ):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()

    snakefile = tmp_path / 'Snakefile'
    snakefile.write_text(f'''
config = {config!r}
import sys
sys.path.insert(0, {str(ROOT / 'scripts')!r})
from scprocess_utils import get_paga_resolution
scprocess_dir = {str(ROOT)!r}
FULL_TAG, DATE_STAMP = "test", "date"
int_dir, hvg_dir, qc_dir = "integration", "hvg", "qc"
cell_cycle_dir, logs_dir, benchmark_dir = "cycle", "logs", "benchmarks"
amb_dir, af_dir, af_rna_dir = "ambient", "af", ""
BATCH_VAR, RUN_VAR = "sample_id", "sample_id"
RUNS, BATCHES, RUNS_TO_BATCHES = ["s1"], ["s1"], {{"s1": ["s1"]}}
CELL_CYCLE_ENABLED = {estimate_phase!r}
TRICYCLE_SCORES_F, TRICYCLE_ORIGIN_F = "cycle/scores.csv.gz", "cycle/origin.csv"
RESOURCE_PARAMS = {{}}
def get_resources(*args, **kwargs):
    return 1000
include: {str(ROOT / 'rules/integration.smk')!r}
''')
    env = {k: v for k, v in os.environ.items() if not k.startswith('SNAKEMAKE_')}
    env['XDG_CACHE_HOME'] = str(tmp_path / 'cache')
    result = subprocess.run([
        sys.executable, '-m', 'snakemake', '--snakefile', str(snakefile),
        '--directory', str(tmp_path), '--dry-run', '--printshellcmds',
        '--cores', '2', '--executor', 'local', '--profile', 'none',
        '--workflow-profile', 'none', '--nolock',
        'integration/integrated_dt_test_date.csv.gz',
    ], text=True, capture_output=True, timeout=45, env=env)
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    commands = re.findall(r'--umap_n_iter\s+(\d+)', output)
    assert commands == ['1234'] * (2 if backend == 'bpcells' else 1), output
    # Evaluate only the optional shell fragment, never the analysis command.
    # This checks that disabled PAGA does not forward a None/empty resolution.
    fragment = re.findall(r'\$\( \[ "(?:True|False)" == "True" \] && echo "--paga-cl-res [^"]+" \)', output)
    assert len(fragment) == 1, output
    effective = subprocess.check_output(['bash', '-c', 'printf "%s" ' + fragment[0]], text=True)
    assert effective == (f'--paga-cl-res{paga_resolution}' if paga else '')
    if correction:
        assert len(re.findall(r"correction_centering\s*=\s*'" + correction + "'", output)) == 2, output
        assert 'pca_with_doublets' in output and 'pca_singlets_only' in output
    else:
        if backend == 'bpcells':
            assert len(re.findall(r"correction_centering\s*=\s*''", output)) == 2, output
            assert len(re.findall(r"harmonics\s*=\s*NULL", output)) == 2, output
            assert len(re.findall(r"ridge_lambda\s*=\s*NULL", output)) == 2, output
    assert ('--tricycle_f' in output) == estimate_phase
