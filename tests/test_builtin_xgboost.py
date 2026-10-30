"""Built-in model integrity, user registration and inference provenance."""
from pathlib import Path
import hashlib
import importlib.util
import io
import json
import sys
import tarfile
import polars as pl
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import train_xgboost as models
from label_celltypes import aggregate_predictions
from scprocess_utils import get_xgboost_model_paths


@pytest.mark.parametrize('payload', [b'', b'scprocess\x00\xff', bytes(range(256)) * 4096 + b'tail'],
  ids=['empty', 'binary', 'multiple_blocks'])
def test_sha256_matches_reference(tmp_path, payload):
  path = tmp_path / 'payload.bin'
  path.write_bytes(payload)
  expected = hashlib.sha256(payload).hexdigest()
  assert models.sha256(path) == expected


@pytest.fixture
def setup_tools():
    spec = importlib.util.spec_from_file_location('scprocess_builtin_setup', ROOT / 'scripts/setup.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def installer_fixture(tmp_path, monkeypatch, setup_tools, bad_member=False):
    content = {'human_cns_xgboost_model.json': b'weights',
        'human_cns_allowed_cls.csv': b'class\na\nb\n', 'human_cns_selected_genes.txt': b'g1\n'}
    archive = tmp_path / 'model.tar.gz'
    with tarfile.open(archive, 'w:gz') as tar:
        for name, value in content.items():
            info = tarfile.TarInfo(('..' if bad_member else 'human_cns-v1') + '/' + name)
            info.size = len(value)
            tar.addfile(info, io.BytesIO(value))
    manifest = dict(name='human_cns', model_id='human_cns-v1', url=archive.as_uri(),
        archive_sha256=models.sha256(archive),
        files={name: hashlib.sha256(value).hexdigest() for name, value in content.items()})
    monkeypatch.setattr(setup_tools, 'builtin_models', lambda: {'human_cns': manifest})
    return archive, manifest


def test_install_and_resolve_identity(tmp_path, monkeypatch, setup_tools):
    _, manifest = installer_fixture(tmp_path, monkeypatch, setup_tools)
    dest = tmp_path / 'scdata/xgboost/human_cns'
    setup_tools._install_builtin_xgboost_model('human_cns', dest)
    model = dest / 'human_cns_xgboost_model.json'
    assert models.get_xgboost_model_info(model, 'human_cns') == dict(
        model_id='human_cns-v1', model_sha256=manifest['files'][model.name])
    pl.DataFrame(dict(model=['human_cns'], model_f=[str(model)],
        cls_f=[str(dest / 'human_cns_allowed_cls.csv')],
        genes_f=[str(dest / 'human_cns_selected_genes.txt')], label_map_f=[''])).write_csv(dest.parent / 'xgboost_models.csv')
    entry = get_xgboost_model_paths({'model': 'human_cns'}, tmp_path / 'scdata')
    assert entry['model_id'] == 'human_cns-v1'
    model.write_text('changed')
    with pytest.raises(ValueError, match='checksum mismatch'):
        get_xgboost_model_paths({'model': 'human_cns'}, tmp_path / 'scdata')


@pytest.mark.parametrize('failure', ['archive', 'member', 'traversal'])
def test_bad_archive_preserves_installation(tmp_path, monkeypatch, failure, setup_tools):
    _, manifest = installer_fixture(tmp_path, monkeypatch, setup_tools, bad_member=failure == 'traversal')
    if failure == 'archive':
        manifest['archive_sha256'] = '0' * 64
    if failure == 'member':
        manifest['files']['human_cns_xgboost_model.json'] = '0' * 64
    dest = tmp_path / 'installed'; dest.mkdir()
    marker = dest / 'existing'; marker.write_text('preserve')
    with pytest.raises(ValueError):
        setup_tools._install_builtin_xgboost_model('human_cns', dest)
    assert list(dest.iterdir()) == [marker] and marker.read_text() == 'preserve'


def test_user_registration_is_separate(tmp_path):
    spec = importlib.util.spec_from_file_location('scprocess_setup', ROOT / 'scripts/setup.py')
    setup = importlib.util.module_from_spec(spec); spec.loader.exec_module(setup)
    assert 'human_cns' in setup.get_xgboost_parameters({})
    source = tmp_path / 'source'; source.mkdir()
    for suffix in ['xgboost_model.json', 'allowed_cls.csv', 'selected_genes.txt']:
        (source / f'trained_{suffix}').write_text('custom')
    setup.add_xgboost_classifier(tmp_path / 'scdata', 'my_classifier', 'trained', source, False)
    model = tmp_path / 'scdata/xgboost/my_classifier/my_classifier_xgboost_model.json'
    assert models.get_xgboost_model_info(model, 'my_classifier')['model_id'] == 'my_classifier'


def test_aggregate_preserves_identity_and_empty_batches(tmp_path):
    valid = pl.DataFrame(dict(cell_id=['a'], predicted_label_naive=['x'], probability_naive=[.9],
        model=['human_cns'], model_id=['human_cns-v1'], model_sha256=['abc'], sample_id=['s']))
    first = tmp_path / 'valid.csv'; valid.write_csv(first)
    empty = tmp_path / 'empty.csv'; valid.clear().write_csv(empty)
    integration = tmp_path / 'integration.csv'
    pl.DataFrame({'cell_id': ['a'], 'sample_id': ['s'], 'RNA_snn_res.4': ['cl1']}).write_csv(integration)
    result = aggregate_predictions([first, empty], integration, 'RNA_snn_res.4', .5, 'sample_id')
    assert result['model_id'].to_list() == ['human_cns-v1']
    assert result['model_sha256'].to_list() == ['abc']



def test_registry_import_does_not_load_training_libraries():
  import subprocess
  code = """
import importlib.abc
import sys
class BlockTraining(importlib.abc.MetaPathFinder):
  def find_spec(self, fullname, path=None, target=None):
    if fullname.split('.')[0] in {'anndata', 'xgboost', 'scipy'}:
      raise AssertionError('Unexpected training import: ' + fullname)
sys.meta_path.insert(0, BlockTraining())
sys.path.insert(0, sys.argv[1])
import train_xgboost
assert 'human_cns' in train_xgboost.builtin_models()
assert train_xgboost.read_label_specs(labels=['a'])[0][0]['label'] == 'a'
import setup
assert 'human_cns' in setup.get_xgboost_parameters({})
"""
  result = subprocess.run([sys.executable, '-c', code, str(ROOT / 'scripts')],
    text = True, capture_output = True, timeout = 30)
  assert result.returncode == 0, result.stdout + result.stderr
