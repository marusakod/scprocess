"""Driver log headers identify the scprocess checkout, not the caller's repository."""
import importlib.machinery
import importlib.util
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load_cli():
    loader = importlib.machinery.SourceFileLoader('scprocess_test_cli', str(ROOT / 'scprocess'))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def test_log_header_uses_scprocess_commit_from_another_git_repo(tmp_path, monkeypatch):
    # a project directory that is itself a git (e.g. workflowr) repository
    project_dir = tmp_path / 'project'
    project_dir.mkdir()
    subprocess.run(['git', 'init', '-q', str(project_dir)], check=True)
    monkeypatch.chdir(project_dir)

    cli = _load_cli()
    _, header = cli._get_main_log(project_dir / 'config.yaml', tmp_path / 'logs', True)

    scprocess_commit = subprocess.check_output(
        ['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    assert f'GIT COMMIT:   {scprocess_commit}' in header


def test_log_header_outside_git_checkout(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cli = _load_cli()

    # simulate a packaged install where git cannot describe the scprocess directory
    def no_git(*args, **kwargs):
        raise subprocess.CalledProcessError(128, args[0])
    monkeypatch.setattr(cli.subprocess, 'check_output', no_git)

    _, header = cli._get_main_log(tmp_path / 'config.yaml', tmp_path / 'logs', True)
    version = (ROOT / 'VERSION').read_text().strip()
    assert f'GIT VERSION:  {version} (not a git checkout)' in header
    assert 'GIT COMMIT:   unavailable' in header
