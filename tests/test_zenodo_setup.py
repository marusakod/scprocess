"""Verified Zenodo installation, isolated updates and canonical model metadata."""
from pathlib import Path
import hashlib
import gzip
import importlib.util
import io
import json
import sys
import tarfile

import polars as pl
import pytest


ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))


@pytest.fixture
def setup_tools():
  spec = importlib.util.spec_from_file_location('zenodo_setup_tests', ROOT / 'scripts/setup.py')
  module = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(module)
  return module


@pytest.fixture
def support_archive(tmp_path, monkeypatch, setup_tools):
  contents = {'marker_genes/human_brain.csv': b'gene,cluster\nAPOE,Astrocyte\n',
    'gmt_pathways/example.gmt': b'pathway\tna\tAPOE\n'}
  archive_f = tmp_path / 'support.tar.gz'
  manifest_f = tmp_path / 'support.json'

  def write_archive(failure=None):
    with tarfile.open(archive_f, 'w:gz') as archive:
      for name, content in contents.items():
        if failure == 'missing' and name.startswith('gmt_pathways/'):
          continue
        if failure == 'traversal' and name.startswith('gmt_pathways/'):
          name = '../outside.gmt'
        member = tarfile.TarInfo(name)
        member.size = len(content)
        archive.addfile(member, io.BytesIO(content))
        if failure == 'duplicate':
          archive.addfile(member, io.BytesIO(content))
    manifest = dict(name='scprocess_data', version='test', url=archive_f.as_uri(),
      archive_sha256=setup_tools.sha256(archive_f),
      files={name: hashlib.sha256(content).hexdigest() for name, content in contents.items()})
    if failure == 'archive':
      manifest['archive_sha256'] = '0' * 64
    if failure == 'member':
      manifest['files']['marker_genes/human_brain.csv'] = '0' * 64
    manifest_f.write_text(json.dumps(manifest))
    return manifest

  write_archive()
  monkeypatch.setattr(setup_tools, 'SCPROCESS_SUPPORT_MANIFEST_F', manifest_f)
  return contents, write_archive


def test_support_install_and_repeat_preserve_unrelated_files(tmp_path, setup_tools, support_archive):
  contents, _ = support_archive
  data_dir = tmp_path / 'scdata'
  (data_dir / 'marker_genes').mkdir(parents=True)
  custom_f = data_dir / 'marker_genes/custom.csv'
  custom_f.write_text('user markers')
  before_cwd = Path.cwd()
  setup_tools.get_scprocess_data(data_dir)
  setup_tools.get_scprocess_data(data_dir)
  assert Path.cwd() == before_cwd
  assert custom_f.read_text() == 'user markers'
  for name, expected in contents.items():
    assert (data_dir / name).read_bytes() == expected
  installed = json.loads((data_dir / 'scprocess_support_manifest.json').read_text())
  assert installed['archive_sha256'] == json.loads(setup_tools.SCPROCESS_SUPPORT_MANIFEST_F.read_text())['archive_sha256']


@pytest.mark.parametrize('failure', ['archive', 'member', 'traversal', 'duplicate', 'missing'])
def test_bad_support_archive_leaves_current_data_unchanged(tmp_path, setup_tools, support_archive, failure):
  _, write_archive = support_archive
  data_dir = tmp_path / 'scdata'
  setup_tools.get_scprocess_data(data_dir)
  before = {str(path.relative_to(data_dir)): path.read_bytes()
    for path in data_dir.rglob('*') if path.is_file()}
  write_archive(failure)
  with pytest.raises(ValueError):
    setup_tools.get_scprocess_data(data_dir)
  after = {str(path.relative_to(data_dir)): path.read_bytes()
    for path in data_dir.rglob('*') if path.is_file()}
  assert before == after
  assert not (tmp_path / 'outside.gmt').exists()


def test_cellranger_download_is_independent(tmp_path, monkeypatch, setup_tools):
  calls = []
  def extract(output_dir, whitelist_f, url, version):
    Path(output_dir).mkdir(parents=True)
    calls.append((output_dir, whitelist_f, url, version))
  monkeypatch.setattr(setup_tools, 'extract_cellranger_resources', extract)
  version_f = tmp_path / 'cellranger_ref/.cellranger_version'
  setup_tools.get_cellranger_data(tmp_path, 'https://example.invalid/cellranger-9.0.1.tar.gz',
    tmp_path / 'cellranger_ref/whitelists.csv', version_f)
  assert version_f.read_text() == '9.0.1'
  assert len(calls) == 1 and calls[0][-1] == '9.0.1'
  assert not (tmp_path / 'marker_genes').exists()


def test_pinned_manifests_and_32_label_model(setup_tools):
  support = json.loads((ROOT / 'resources/scprocess_data.json').read_text())
  model = json.loads((ROOT / 'resources/xgboost/human_cns.json').read_text())
  prefix = 'https://zenodo.org/records/23803028/files/'
  assert support['url'] == prefix + 'scprocess_data_archive.tar.gz'
  assert len(support['files']) == 22
  assert all(name.startswith(('marker_genes/', 'gmt_pathways/')) for name in support['files'])
  assert model['url'] == prefix + 'human_cns-v1.tar.gz'
  assert model['n_classes'] == 32 and model['n_genes'] == 2386
  assert model['files']['human_cns_xgboost_model.json'] == '52dbe38828feb55011158ded7f6068628cc7125431ed8cb969235ade0f561cc0'
  for urls in (setup_tools.URLS_GTF_TXTS, setup_tools.URLS_ZEN_IDXS, setup_tools.URLS_PROBE_GTF_TXTS):
    assert all(url.startswith(prefix) for url in urls.values())


def test_registry_refresh_keeps_other_installed_custom_models(tmp_path, setup_tools):
  source_dir = tmp_path / 'source'
  source_dir.mkdir()
  for suffix, content in [('xgboost_model.json', 'weights'), ('allowed_cls.csv', 'class\na\n'),
      ('selected_genes.txt', 'g1\n')]:
    (source_dir / ('trained_' + suffix)).write_text(content)
  scdata_dir = tmp_path / 'scdata'
  for name in ('custom', 'configured'):
    setup_tools.add_xgboost_classifier(scdata_dir, name, 'trained', source_dir, False)
  csv_f = scdata_dir / 'xgboost/xgboost_models.csv'
  setup_tools.save_available_xgboost_csv(csv_f, scdata_dir, ['custom'])
  custom_before = pl.read_csv(csv_f).row(0, named=True)
  setup_tools.save_available_xgboost_csv(csv_f, scdata_dir, ['configured'])
  models_df = pl.read_csv(csv_f)
  assert set(models_df['model']) == {'custom', 'configured'}
  assert models_df.filter(pl.col('model') == 'custom').row(0, named=True) == custom_before


def test_model_update_retires_only_previously_declared_companions(tmp_path, monkeypatch, setup_tools):
  contents = {'human_cns_xgboost_model.json': b'new weights',
    'human_cns_allowed_cls.csv': b'class\na\n', 'human_cns_selected_genes.txt': b'g1\n'}
  archive_f = tmp_path / 'model.tar.gz'
  with tarfile.open(archive_f, 'w:gz') as archive:
    for name, content in contents.items():
      member = tarfile.TarInfo('human_cns-v1/' + name)
      member.size = len(content)
      archive.addfile(member, io.BytesIO(content))
  manifest = dict(name='human_cns', model_id='human_cns-v1', url=archive_f.as_uri(),
    archive_sha256=setup_tools.sha256(archive_f),
    files={name: hashlib.sha256(content).hexdigest() for name, content in contents.items()})
  monkeypatch.setattr(setup_tools, 'builtin_models', lambda: {'human_cns': manifest})
  dest = tmp_path / 'installed'
  dest.mkdir()
  (dest / 'old-evaluation.csv').write_text('old metrics')
  (dest / 'user-note.txt').write_text('keep')
  (dest / 'human_cns_manifest.json').write_text(json.dumps({'files': {'old-evaluation.csv': 'unused'}}))
  setup_tools._install_builtin_xgboost_model('human_cns', dest)
  assert not (dest / 'old-evaluation.csv').exists()
  assert (dest / 'user-note.txt').read_text() == 'keep'
  assert (dest / 'human_cns_xgboost_model.json').read_bytes() == b'new weights'


@pytest.mark.parametrize('stored,expected_reuse', [('current', True), ('old', False), (None, False)])
def test_join_reuse_requires_matching_weights(tmp_path, stored, expected_reuse):
  from scprocess_utils import get_join_batch_sources
  project = dict(proj_dir=str(tmp_path), short_tag='source', full_tag='source', date_stamp='2026-01-01')
  config = {'project': project, 'label_celltypes': [{'labeller': 'scprocess', 'model': 'human_cns'}]}
  labels_f = tmp_path / 'output/source_label_celltypes/labels_scprocess_model_human_cns_source_2026-01-01.csv.gz'
  labels_f.parent.mkdir(parents=True)
  labels_df = pl.DataFrame({'cell_id': ['c1'], 'predicted_label_naive': ['Neuron']})
  if stored is not None:
    labels_df = labels_df.with_columns(pl.lit(stored).alias('model_sha256'))
  with gzip.open(labels_f, 'wt') as stream:
    stream.write(labels_df.write_csv())
  fresh, reused = get_join_batch_sources({'p': config}, ['p'], ['p_s1', 'p_s2'],
    'scprocess', 'human_cns', model_sha256='current')
  assert bool(reused) == expected_reuse
  assert fresh == ([] if expected_reuse else ['p_s1', 'p_s2'])


def test_reused_extraction_keeps_identity_and_checks_every_cell(tmp_path):
  from label_celltypes import extract_naive_predictions
  labels_f = tmp_path / 'source.csv'
  labels_df = pl.DataFrame(dict(cell_id=['a', 'b'], predicted_label_naive=['Neuron', 'Neuron'],
    probability_naive=[.9, .8], model_id=['human_cns-v1'] * 2, model_sha256=['current'] * 2))
  labels_df.write_csv(labels_f)
  integration_f = tmp_path / 'integration.csv'
  pl.DataFrame({'cell_id': ['a', 'b']}).write_csv(integration_f)
  output_f = tmp_path / 'reused.csv.gz'
  extract_naive_predictions([labels_f], integration_f, 'human_cns', output_f, model_sha256='current')
  assert output_f.read_bytes().startswith(b'\x1f\x8b')
  reused_df = pl.read_csv(output_f)
  assert reused_df['model_sha256'].to_list() == ['current', 'current']
  assert reused_df['model_id'].to_list() == ['human_cns-v1', 'human_cns-v1']
  labels_df = labels_df.with_columns(pl.Series('model_sha256', ['current', 'old']))
  labels_df.write_csv(labels_f)
  with pytest.raises(ValueError, match='different model weights'):
    extract_naive_predictions([labels_f], integration_f, 'human_cns', output_f, model_sha256='current')


def test_aggregation_rejects_mixed_weight_predictions(tmp_path):
  from label_celltypes import aggregate_predictions
  labels_f = tmp_path / 'predictions.csv'
  pl.DataFrame(dict(cell_id=['a', 'b'], predicted_label_naive=['Neuron'] * 2,
    probability_naive=[.9, .8], model=['human_cns'] * 2,
    model_sha256=['current', 'old'])).write_csv(labels_f)
  integration_f = tmp_path / 'integration.csv'
  pl.DataFrame({'cell_id': ['a', 'b'], 'RNA_snn_res.4': ['cl1', 'cl1']}).write_csv(integration_f)
  with pytest.raises(ValueError, match='mixed model checksums'):
    aggregate_predictions([labels_f], integration_f, 'RNA_snn_res.4', .5, 'sample_id')
