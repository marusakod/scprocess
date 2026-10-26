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
import xgboost_models as models
from label_celltypes import aggregate_predictions


def installer_fixture(tmp_path, monkeypatch, bad_member=False):
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
    monkeypatch.setattr(models, 'builtin_models', lambda: {'human_cns': manifest})
    return archive, manifest


def test_install_and_resolve_identity(tmp_path, monkeypatch):
    _, manifest = installer_fixture(tmp_path, monkeypatch)
    dest = tmp_path / 'scdata/xgboost/human_cns'
    models.install_builtin('human_cns', dest)
    model = dest / 'human_cns_xgboost_model.json'
    assert models.model_identity(model, 'human_cns') == dict(
        model_id='human_cns-v1', model_sha256=manifest['files'][model.name])
    pl.DataFrame(dict(model=['human_cns'], model_f=[str(model)],
        cls_f=[str(dest / 'human_cns_allowed_cls.csv')],
        genes_f=[str(dest / 'human_cns_selected_genes.txt')], label_map_f=[''])).write_csv(dest.parent / 'xgboost_models.csv')
    entry = models.resolve_xgboost_entry({'model': 'human_cns'}, tmp_path / 'scdata')
    assert entry['model_id'] == 'human_cns-v1'
    model.write_text('changed')
    with pytest.raises(ValueError, match='checksum mismatch'):
        models.resolve_xgboost_entry({'model': 'human_cns'}, tmp_path / 'scdata')


@pytest.mark.parametrize('failure', ['archive', 'member', 'traversal'])
def test_bad_archive_preserves_installation(tmp_path, monkeypatch, failure):
    _, manifest = installer_fixture(tmp_path, monkeypatch, bad_member=failure == 'traversal')
    if failure == 'archive':
        manifest['archive_sha256'] = '0' * 64
    if failure == 'member':
        manifest['files']['human_cns_xgboost_model.json'] = '0' * 64
    dest = tmp_path / 'installed'; dest.mkdir()
    marker = dest / 'existing'; marker.write_text('preserve')
    with pytest.raises(ValueError):
        models.install_builtin('human_cns', dest)
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
    assert models.model_identity(model, 'my_classifier')['model_id'] == 'my_classifier'


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
