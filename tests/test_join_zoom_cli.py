"""CLI dispatch checks: no scheduler calls, no workflow state changes."""
import importlib.machinery
import importlib.util
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def cli():
    loader = importlib.machinery.SourceFileLoader('scprocess_test_cli', str(ROOT / 'scprocess'))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def test_multiple_zooms_run_sequentially_without_mutating_arguments(cli, tmp_path, monkeypatch):
    config = tmp_path / 'join.yaml'
    config.write_text('join: {}\n')
    monkeypatch.setattr(cli.scprocess_utils, 'get_join_zoom_specs', lambda *a: {'neurons': {}, 'glia': {}})
    calls = []
    def run(*args, **kwargs):
        calls.append((list(args[2]), kwargs))
        args[2].append('mutated-by-workflow-setup')
    monkeypatch.setattr(cli, 'run_join', run)
    extra = ['--printshellcmds']
    cli.run_join_zooms(config, ROOT, extra, True)
    assert [kw['zoom_name'] for _, kw in calls] == ['neurons', 'glia']
    assert [args for args, _ in calls] == [extra, extra]
    assert extra == ['--printshellcmds']


def test_shiny_dispatch_and_unknown_zoom(cli, tmp_path, monkeypatch):
    config = tmp_path / 'join.yaml'
    config.write_text('join: {}\n')
    monkeypatch.setattr(cli.scprocess_utils, 'get_join_zoom_specs', lambda *a: {'neurons': {}})
    run = Mock()
    monkeypatch.setattr(cli, 'run_join', run)
    cli.run_join_zooms(config, ROOT, [], True, 'neurons', shiny=True)
    assert run.call_args.kwargs == {'rule': 'configure_shiny_app', 'zoom_name': 'neurons'}
    with pytest.raises(ValueError, match='Unknown or missing'):
        cli.run_join_zooms(config, ROOT, [], True, 'absent')


@pytest.mark.parametrize('zoom,rule,snakefile', [
    (None, 'all', 'join.smk'), ('neurons', 'all', 'join_zoom.smk'),
    ('neurons', 'configure_shiny_app', 'shiny.smk')])
def test_run_join_preserves_profile_and_extra_arguments(cli, tmp_path, monkeypatch, zoom, rule, snakefile):
    config = tmp_path / 'join.yaml'
    config.write_text(yaml.safe_dump({'join': {'proj_dir': str(tmp_path)}}))
    profile = tmp_path / 'profile'
    profile.mkdir()
    (profile / 'config.yaml').write_text('executor: slurm\n')
    def setup(sc_dir, extra):
        return tmp_path, [*extra, '--workflow-profile', str(profile)], {'user': {'profile_dir': profile}}
    monkeypatch.setattr(cli.scprocess_utils, 'check_setup_before_running_scprocess', setup)
    monkeypatch.setattr(cli.scprocess_utils, 'get_conda_prefix', lambda *a: str(tmp_path / 'envs'))
    monkeypatch.setattr(cli, '_get_main_log', lambda *a: (tmp_path / 'driver.log', 'test\n'))
    monkeypatch.setattr(cli, '_clean_conda_env', lambda *a: {})
    launch = Mock(return_value=SimpleNamespace(returncode=0))
    monkeypatch.setattr(cli.subprocess, 'run', launch)
    cli.run_join(config, ROOT, ['--set-threads', 'join_marker_genes_calc=4'], True, rule, zoom)
    command = launch.call_args.args[0]
    assert f'rules/{snakefile}' in command
    assert '--workflow-profile' in command and str(profile) in command
    assert '--set-threads join_marker_genes_calc=4' in command
    assert '--slurm-logdir' in command and '-np' in command
    assert ('join_zoom_name=neurons' in command) == (zoom is not None)
    assert launch.call_args.kwargs['cwd'] == tmp_path
