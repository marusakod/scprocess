"""The QC rules are split so that changing QC thresholds does not rerun doublet identification."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
QC_SMK = (ROOT / 'rules' / 'qc.smk').read_text()
SAMPLE_QC_R = (ROOT / 'scripts' / 'SampleQC.R').read_text()


def _rule_block(name):
  match = re.search(rf'^rule {name}:\n(.*?)(?=^rule |^checkpoint |\Z)', QC_SMK, re.S | re.M)
  assert match, f'rule {name} not found'
  return match.group(1)


def test_the_combined_rule_is_gone():
  assert 'rule run_qc_one_run:' not in QC_SMK
  assert 'main_qc(' not in QC_SMK
  assert 'main_qc <-' not in SAMPLE_QC_R


def test_doublet_rule_does_not_depend_on_qc_thresholds():
  block = _rule_block('run_doublet_id_one_run')
  assert 'cuts_f' not in block
  assert 'qc_thresholds_by_' not in block
  assert 'qc_hard_' not in block
  assert 'main_doublet_id(' in block


def test_doublet_rule_keeps_the_doublet_settings_and_outputs():
  block = _rule_block('run_doublet_id_one_run')
  for setting in ('dbl_min_feats', 'dbl_min_counts', 'dbl_dbr_sd', 'dbl_clusters'):
    assert f'{setting}' in block
  for output in ('coldata_dbl_f', 'rowdata_f', 'dbl_f', 'dbl_artificial_f', 'dbl_summary_f'):
    assert re.search(rf'^\s+{output}\s*=', block, re.M)
  # the table of cells is a kept result: if it were temp, changing a threshold would rerun this rule
  assert not re.search(r'coldata_dbl_f\s*=\s*temp\(', block)


def test_filter_rule_reads_the_cell_table_and_thresholds_not_the_counts():
  block = _rule_block('filter_qc_one_run')
  assert 'cuts_f' in block
  assert 'coldata_dbl_f' in block
  assert 'af_h5_f' not in block
  assert 'amb_yaml_f' not in block
  assert 'main_filter_qc(' in block
  assert re.search(r'qc_f\s*=\s*temp\(', block) and re.search(r'coldata_f\s*=\s*temp\(', block)


def test_filter_rule_consumes_the_doublet_rule_output():
  doublet = _rule_block('run_doublet_id_one_run')
  filt = _rule_block('filter_qc_one_run')
  pattern = "coldata_dt_with_doublets_{{run}}_{FULL_TAG}_{DATE_STAMP}.csv.gz"
  assert pattern in doublet and pattern in filt


def test_merge_rules_still_use_the_tmp_filter_outputs():
  merge_qc = _rule_block('merge_qc')
  assert 'tmp_qc_dt_' in merge_qc and 'tmp_coldata_dt_' in merge_qc
  assert 'tmp_rowdata_dt_' in _rule_block('merge_rowdata')


def test_filter_function_takes_the_cell_table_not_an_sce():
  assert re.search(r'^\.filter_qc <- function\(coldata_in,', SAMPLE_QC_R, re.M)
  assert 'colData(sce)' not in SAMPLE_QC_R.split('.filter_qc <- function')[1].split('\n}\n')[0]
