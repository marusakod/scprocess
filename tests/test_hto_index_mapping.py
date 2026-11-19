import warnings
from pathlib import Path
from unittest.mock import patch

import polars as pl
import pytest

from scripts import mapping


ROOT = Path(__file__).resolve().parents[1]

# 15 bp features that share no 11-mer
FEATURES = {
  'Hashtag-1': 'ACGTACGTACGTACG',
  'Hashtag-2': 'TTGCAGTTGCAGTTG',
  'Hashtag-3': 'CAGTCAGTCAGTCAG',
  'CD3':       'GGATCCGGATCCGGA',
  'CD19':      'AATTCCAATTCCAAT',
}


def _feature_ref_f(tmp_path, features=FEATURES):
  feature_ref_f = tmp_path / 'features.csv'
  pl.DataFrame({'hto_id': list(features), 'sequence': list(features.values())}).write_csv(feature_ref_f)
  return feature_ref_f


def test_pool_features_are_the_pool_htos_then_all_adts_in_reference_order(tmp_path):
  feat_ref_df = pl.read_csv(_feature_ref_f(tmp_path))
  features_df = mapping._select_pool_features(feat_ref_df, ['Hashtag-3', 'Hashtag-1'], ['CD19', 'CD3'])
  assert features_df['hto_id'].to_list() == ['Hashtag-1', 'Hashtag-3', 'CD3', 'CD19']
  assert features_df['feature_type'].to_list() == ['hto', 'hto', 'adt', 'adt']


def test_pool_features_without_adts_are_only_the_pool_htos(tmp_path):
  feat_ref_df = pl.read_csv(_feature_ref_f(tmp_path))
  features_df = mapping._select_pool_features(feat_ref_df, ['Hashtag-2'], [])
  assert features_df['hto_id'].to_list() == ['Hashtag-2']


def test_pool_features_reject_an_id_missing_from_the_reference(tmp_path):
  feat_ref_df = pl.read_csv(_feature_ref_f(tmp_path))
  with pytest.raises(ValueError, match=r"not in the feature reference file: \['Hashtag-9'\]"):
    mapping._select_pool_features(feat_ref_df, ['Hashtag-9'], [])


def test_build_hto_index_writes_the_reference_files_and_runs_simpleaf_index(tmp_path):
  out_dir = tmp_path / 'hto_ref' / 'pool1'
  commands = []
  working_dirs = {}
  with patch.object(mapping.subprocess, 'run', side_effect=lambda cmd, **kw: (
      commands.append(cmd), working_dirs.update({cmd[1]: kw.get('cwd')}))):
    mapping.build_hto_index(_feature_ref_f(tmp_path), ['Hashtag-1', 'Hashtag-2'], ['CD3'],
      str(out_dir), '/af_home', 2)

  assert (out_dir / 'hto.fa').read_text() == (
    '>Hashtag-1\nACGTACGTACGTACG\n>Hashtag-2\nTTGCAGTTGCAGTTG\n>CD3\nGGATCCGGATCCGGA\n')
  assert (out_dir / 't2g_hto.tsv').read_text() == 'Hashtag-1\tHashtag-1\nHashtag-2\tHashtag-2\nCD3\tCD3\n'
  features_df = pl.read_csv(out_dir / 'features.csv')
  assert features_df['feature_type'].to_list() == ['hto', 'hto', 'adt']

  index_cmd = next(cmd for cmd in commands if cmd[:2] == ['simpleaf', 'index'])
  assert index_cmd[index_cmd.index('--kmer-length') + 1] == '11'
  assert index_cmd[index_cmd.index('--minimizer-length') + 1] == '9'
  assert index_cmd[index_cmd.index('--ref-seq') + 1] == str(out_dir / 'hto.fa')
  assert index_cmd[index_cmd.index('--output') + 1] == str(out_dir / 'hto_index')

  # piscem writes ./workdir.noindex: builds of different pools must not share a working directory
  assert working_dirs['index'] == str(out_dir)


def test_kmer_and_minimizer_lengths_are_consistent():
  assert mapping.HTO_KMER_LENGTH == 11
  assert mapping.HTO_MINIMIZER_LENGTH < mapping.HTO_KMER_LENGTH


def test_shared_kmers_are_reported_with_the_pair(tmp_path):
  # A and B share ACGTACGTACG and TACGTACGTAC
  features_df = pl.DataFrame({'hto_id': ['A', 'B'], 'sequence': ['ACGTACGTACGTACG', 'TTTTACGTACGTACG']})
  with pytest.warns(UserWarning, match=r'Features A and B share 2 11-mer\(s\)'):
    mapping._warn_shared_kmers(features_df, 11)


def test_distinct_features_give_no_shared_kmer_warning():
  features_df = pl.DataFrame({'hto_id': list(FEATURES), 'sequence': list(FEATURES.values())})
  with warnings.catch_warnings():
    warnings.simplefilter('error')
    mapping._warn_shared_kmers(features_df, 11)


def test_a_feature_shorter_than_the_kmer_is_an_error():
  features_df = pl.DataFrame({'hto_id': ['A'], 'sequence': ['ACGTACGTAC']})
  with pytest.raises(ValueError, match=r"shorter than the 11-mer length.*\['A'\]"):
    mapping._warn_shared_kmers(features_df, 11)


def _map_with_recorded_orientation(tmp_path, what, exp_ori):
  seen = {}

  def fake_resolve_chemistry(R1_fs, R2_fs, af_chemistry, exp_ori, *args):
    seen['resolve_ori'] = exp_ori
    return {'af_chemistry': '10xv3', 'exp_ori': exp_ori, 'whitelist_f': 'wl.txt'}

  def fake_quant(out_dir, R1_fs, R2_fs, threads, index_dir, chemistry, ori, t2g_f, wl_f, extra_args=None):
    seen['quant_ori'] = ori

  with patch.object(mapping, '_resolve_fastq_paths', return_value=(['r1'], ['r2'], False, None)), \
       patch.object(mapping, '_resolve_chemistry', side_effect=fake_resolve_chemistry), \
       patch.object(mapping, '_run_simpleaf_quant', side_effect=fake_quant), \
       patch.object(mapping, '_write_chemistry_stats'), \
       patch.object(mapping.pl, 'read_csv'), \
       patch.object(mapping.subprocess, 'run'):
    mapping.map_fastqs_to_counts('pool1', str(tmp_path), 'hto', what, '/af_home', 'local',
      ['r1'], ['r2'], 2, 't2g.tsv', 'index', 'wl_lu.csv', af_chemistry='10xv3', exp_ori=exp_ori,
      whitelist_f='wl.txt')
  return seen


@pytest.mark.parametrize('rna_ori', ['rc', 'fw', 'none'])
def test_hto_mapping_is_always_forward(tmp_path, rna_ori):
  seen = _map_with_recorded_orientation(tmp_path, 'hto', rna_ori)
  assert seen == {'resolve_ori': 'fw', 'quant_ori': 'fw'}


def test_rna_mapping_keeps_the_orientation_it_is_given(tmp_path):
  seen = _map_with_recorded_orientation(tmp_path, 'rna', 'rc')
  assert seen['quant_ori'] == 'rc'


def test_hto_rule_does_not_read_the_rna_orientation():
  rule_text = (ROOT / 'rules/hto.smk').read_text()
  assert 'selected_ori' not in rule_text
  assert '--exp_ori' not in rule_text
