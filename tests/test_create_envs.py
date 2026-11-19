"""scprocess setup --create-envs: environment coverage and Snakemake command."""
import ast
import importlib.machinery
import importlib.util
import re
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def cli():
    loader = importlib.machinery.SourceFileLoader('scprocess_test_cli', str(ROOT / 'scprocess'))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def _create_envs_list():
    text = (ROOT / 'rules/create_envs.smk').read_text()
    match = re.search(r'RULE_ENV_NAMES = (\[.*?\])', text, re.S)
    return ast.literal_eval(match.group(1))


def test_create_envs_covers_every_rule_environment():
    # environments named by conda: directives in workflow rules, e.g. '../envs/rlibs.yaml' -> 'rlibs'
    rule_envs = set()
    for smk_f in (ROOT / 'rules').glob('*.smk'):
        if smk_f.name == 'create_envs.smk':
            continue
        rule_envs.update(re.findall(r"envs/([A-Za-z0-9_]+)\.yaml", smk_f.read_text()))

    env_names = _create_envs_list()
    assert len(env_names) == len(set(env_names))
    assert set(env_names) == rule_envs
    assert all((ROOT / 'envs' / f'{name}.yaml').is_file() for name in env_names)


def test_create_envs_uses_run_conda_prefix(cli, tmp_path, monkeypatch):
    prefix = tmp_path / 'envs'
    monkeypatch.setattr(cli.scprocess_utils, 'get_conda_prefix', lambda *a: prefix)
    monkeypatch.setattr(cli, '_clean_conda_env', lambda *a: {})
    launch = Mock(return_value=SimpleNamespace(returncode=0))
    monkeypatch.setattr(cli.subprocess, 'run', launch)

    cli._create_conda_envs(ROOT, {'user': {}}, [], tmp_path, tmp_path / 'setup.log', False)
    command = launch.call_args.args[0]
    assert 'rules/create_envs.smk' in command
    assert '--conda-create-envs-only' in command
    assert f'--conda-prefix {prefix}' in command
    assert launch.call_args.kwargs['cwd'] == tmp_path

    # Snakemake would create envs despite -n, so dry runs omit --conda-create-envs-only
    cli._create_conda_envs(ROOT, {'user': {}}, ['-np'], tmp_path, tmp_path / 'setup.log', True)
    command = launch.call_args.args[0]
    assert '-np' in command and '--conda-create-envs-only' not in command
