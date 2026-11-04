"""Snakemake-side lookup of prediction chunks (base scprocess environment)."""
from pathlib import Path
import sys

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from scprocess_utils import get_good_batches, get_prediction_chunks


def test_chunks_are_restricted_to_requested_batches(tmp_path):
  chunks_f = tmp_path / 'chunks.csv'
  pl.DataFrame({'batch': ['a', 'b', 'c', '10'], 'n_cells': [5, 5, 5, 5],
    'chunk': ['chunk001', 'chunk001', 'chunk002', 'chunk003']}).write_csv(chunks_f)
  assert get_prediction_chunks(chunks_f) == {'chunk001': ['a', 'b'], 'chunk002': ['c'], 'chunk003': ['10']}

  # Only requested batches are kept and empty chunks dropped; numeric-looking
  # batch names stay strings; absent batches (no H5AD) are ignored.
  assert get_prediction_chunks(chunks_f, ['a', 'c', 'missing']) == {'chunk001': ['a'], 'chunk002': ['c']}
  assert get_prediction_chunks(chunks_f, [10]) == {'chunk003': ['10']}


def test_good_batches_exclude_failed_qc(tmp_path):
  qc_f = tmp_path / 'qc.csv'
  pl.DataFrame({'sample_id': ['s1', 's2', 's3'], 'bad_sample_id': [False, True, False]}).write_csv(qc_f)
  assert get_good_batches(qc_f, 'sample_id') == ['s1', 's3']
