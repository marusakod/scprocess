import json
from pathlib import Path

import jsonschema
import polars as pl
import pytest

from scripts.scprocess_utils import (
  _check_multiplexing_parameters, get_hto_pool_references, split_hto_adt_ids)


ROOT = Path(__file__).resolve().parents[1]
FEATURE_IDS = ['Hashtag-1', 'Hashtag-2', 'Hashtag-3', 'CD3', 'CD19']


def _make_config(tmp_path, regex=None, sample_hto_ids=None, pools=None):
  # two pools; pool1 uses Hashtag-1/2 and pool2 uses Hashtag-2/3, unless sample_hto_ids is given
  sample_hto_ids = sample_hto_ids or ['Hashtag-1', 'Hashtag-2', 'Hashtag-2', 'Hashtag-3']
  pools = pools or ['pool1', 'pool1', 'pool2', 'pool2']
  samples_f = tmp_path / 'samples.csv'
  pl.DataFrame({
    'sample_id': [f's{i}' for i in range(len(pools))], 'pool_id': pools, 'hto_id': sample_hto_ids,
  }).write_csv(samples_f)
  feature_ref_f = tmp_path / 'features.csv'
  pl.DataFrame({
    'hto_id': FEATURE_IDS, 'sequence': ['ACGTACGTACGTACG', 'TTGCAGTTGCAGTTG', 'CAGTCAGTCAGTCAG',
      'GGATCCGGATCCGGA', 'AATTCCAATTCCAAT']}).write_csv(feature_ref_f)
  fastq_dir = tmp_path / 'fastqs'
  fastq_dir.mkdir()

  multiplexing = {'demux_type': 'hto', 'feature_ref': str(feature_ref_f), 'fastq_dir': str(fastq_dir),
    'seurat_quantile': 0.99}
  if regex is not None:
    multiplexing['hto_id_regex'] = regex
  return {'project': {'proj_dir': tmp_path, 'sample_metadata': str(samples_f), 'short_tag': 'test'},
    'multiplexing': multiplexing}


def test_split_without_a_regex_makes_every_feature_an_hto():
  assert split_hto_adt_ids(FEATURE_IDS, '') == (FEATURE_IDS, [])
  assert split_hto_adt_ids(FEATURE_IDS, None) == (FEATURE_IDS, [])


def test_split_with_a_regex_puts_matching_features_in_hto_and_the_rest_in_adt():
  assert split_hto_adt_ids(FEATURE_IDS, '^Hashtag') == (['Hashtag-1', 'Hashtag-2', 'Hashtag-3'], ['CD3', 'CD19'])
  # the regex is searched, not anchored
  assert split_hto_adt_ids(FEATURE_IDS, 'tag-[12]$')[0] == ['Hashtag-1', 'Hashtag-2']


def test_pool_references_use_the_pool_htos_and_all_adts(tmp_path):
  config = _make_config(tmp_path, regex='^Hashtag')
  refs = get_hto_pool_references(config, ['pool1', 'pool2'])
  assert refs['pool1'] == {'hto_ids': ['Hashtag-1', 'Hashtag-2'], 'adt_ids': ['CD3', 'CD19']}
  assert refs['pool2'] == {'hto_ids': ['Hashtag-2', 'Hashtag-3'], 'adt_ids': ['CD3', 'CD19']}


def test_pool_references_without_a_regex_have_no_adts_and_only_the_pool_htos(tmp_path):
  refs = get_hto_pool_references(_make_config(tmp_path), ['pool1'])
  assert refs['pool1'] == {'hto_ids': ['Hashtag-1', 'Hashtag-2'], 'adt_ids': []}


def test_pool_references_are_empty_without_hto_demultiplexing(tmp_path):
  config = _make_config(tmp_path)
  config['multiplexing'] = {'demux_type': 'none'}
  assert get_hto_pool_references(config, ['pool1']) == {}


def test_valid_regex_passes_the_check(tmp_path):
  config = _make_config(tmp_path, regex='^Hashtag')
  assert _check_multiplexing_parameters(config)['multiplexing']['hto_id_regex'] == '^Hashtag'


def test_config_without_a_regex_passes_the_check(tmp_path):
  _check_multiplexing_parameters(_make_config(tmp_path))


def test_invalid_regex_is_an_error(tmp_path):
  with pytest.raises(ValueError, match='not a valid regular expression'):
    _check_multiplexing_parameters(_make_config(tmp_path, regex='^Hash(tag'))


def test_regex_matching_nothing_is_an_error(tmp_path):
  with pytest.raises(ValueError, match='does not match any hto_id'):
    _check_multiplexing_parameters(_make_config(tmp_path, regex='^TotalSeq'))


def test_sample_hto_id_not_matching_the_regex_is_an_error(tmp_path):
  config = _make_config(tmp_path, regex='^Hashtag', sample_hto_ids=['Hashtag-1', 'Hashtag-2', 'CD3', 'Hashtag-3'])
  with pytest.raises(ValueError, match=r"do not match hto_id_regex: \['CD3'\]"):
    _check_multiplexing_parameters(config)


def test_sample_hto_id_missing_from_the_reference_is_still_an_error(tmp_path):
  config = _make_config(tmp_path, regex='^Hashtag', sample_hto_ids=['Hashtag-1', 'Hashtag-2', 'Hashtag-7', 'Hashtag-3'])
  with pytest.raises(ValueError, match="don't match hto_id values"):
    _check_multiplexing_parameters(config)


def _multiplexing_validator():
  schema = json.loads((ROOT / 'resources/schemas/config.schema.json').read_text())
  return jsonschema.Draft7Validator(schema['properties']['multiplexing'])


def test_schema_accepts_the_regex_with_hto_demultiplexing():
  _multiplexing_validator().validate({
    'demux_type': 'hto', 'feature_ref': 'f.csv', 'fastq_dir': 'fq', 'hto_id_regex': '^Hashtag'})


@pytest.mark.parametrize('demux_type', ['none', 'ocm', 'flex'])
def test_schema_rejects_the_regex_without_hto_demultiplexing(demux_type):
  with pytest.raises(jsonschema.ValidationError):
    _multiplexing_validator().validate({'demux_type': demux_type, 'hto_id_regex': '^Hashtag'})


def test_schema_rejects_an_empty_regex():
  with pytest.raises(jsonschema.ValidationError):
    _multiplexing_validator().validate({
      'demux_type': 'hto', 'feature_ref': 'f.csv', 'fastq_dir': 'fq', 'hto_id_regex': ''})


@pytest.mark.parametrize('bad_id', ['Hashtag 1', 'Hashtag,1', 'Hashtag"1', "Hashtag'1", 'Hashtag$1', 'Hashtag;1', 'Hashtag(1)'])
def test_feature_ids_with_whitespace_commas_quotes_or_shell_characters_are_an_error(tmp_path, bad_id):
  config = _make_config(tmp_path)
  feature_ref_f = Path(config['multiplexing']['feature_ref'])
  feature_ref_df = pl.read_csv(feature_ref_f)
  feature_ref_df = feature_ref_df.with_columns(
    hto_id = pl.when(pl.col('hto_id') == 'CD19').then(pl.lit(bad_id)).otherwise(pl.col('hto_id')))
  feature_ref_df.write_csv(feature_ref_f)
  with pytest.raises(ValueError, match='whitespace, commas, quotes or shell characters'):
    _check_multiplexing_parameters(config)


@pytest.mark.parametrize('good_id', ['TotalSeqC_Hashtag_1', 'Hashtag-9', '1_CD184_CXCR4', 'CD307e.FCRL5', 'A:B'])
def test_ordinary_feature_ids_are_accepted(tmp_path, good_id):
  config = _make_config(tmp_path)
  feature_ref_f = Path(config['multiplexing']['feature_ref'])
  feature_ref_df = pl.read_csv(feature_ref_f)
  feature_ref_df = feature_ref_df.with_columns(
    hto_id = pl.when(pl.col('hto_id') == 'CD19').then(pl.lit(good_id)).otherwise(pl.col('hto_id')))
  feature_ref_df.write_csv(feature_ref_f)
  _check_multiplexing_parameters(config)
