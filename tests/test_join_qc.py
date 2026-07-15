import gzip
import tempfile
import unittest
from pathlib import Path

import polars as pl

from scripts.join import build_joint_qc


class TestBuildJointQc(unittest.TestCase):
  def write_csv(self, path, data):
    df = pl.DataFrame(data)
    if str(path).endswith('.gz'):
      with gzip.open(path, 'wb') as fh:
        df.write_csv(fh)
    else:
      df.write_csv(path)

  def test_filters_and_orders_qc_to_joint_coldata(self):
    with tempfile.TemporaryDirectory() as tmpdir:
      tmpdir = Path(tmpdir)
      qc_a = tmpdir / 'qc_a.csv.gz'
      qc_b = tmpdir / 'qc_b.csv.gz'
      coldata = tmpdir / 'coldata.csv.gz'
      out_f = tmpdir / 'joint_qc.csv.gz'

      self.write_csv(qc_a, {
        'cell_id': ['a1', 'a2', 'a_extra'],
        'log_counts': [2.0, 2.1, 1.0],
        'log_feats': [1.5, 1.6, 0.8],
        'logit_mito': [-3.0, -2.5, -1.0],
        'logit_spliced': [1.0, 1.2, 0.0],
      })
      self.write_csv(qc_b, {
        'cell_id': ['b1'],
        'log_counts': [2.2],
        'log_feats': [1.7],
        'logit_mito': [-2.0],
      })
      self.write_csv(coldata, {
        'cell_id': ['b1', 'a2', 'a1'],
        'sample_id': ['project_b_s2', 'project_a_s1', 'project_a_s1'],
        'project_id': ['project_b', 'project_a', 'project_a'],
      })

      build_joint_qc(
        qc_fs=[str(qc_a), str(qc_b)],
        project_ids=['project_a', 'project_b'],
        coldata_f=str(coldata),
        out_f=str(out_f),
      )

      result = pl.read_csv(out_f)
      self.assertEqual(result['cell_id'].to_list(), ['b1', 'a2', 'a1'])
      self.assertEqual(result['sample_id'].to_list(),
        ['project_b_s2', 'project_a_s1', 'project_a_s1'])
      self.assertEqual(result['project_id'].to_list(),
        ['project_b', 'project_a', 'project_a'])
      self.assertIn('logit_spliced', result.columns)
      self.assertIsNone(result.filter(pl.col('cell_id') == 'b1')['logit_spliced'][0])

  def test_errors_when_joined_cell_has_no_qc_metrics(self):
    with tempfile.TemporaryDirectory() as tmpdir:
      tmpdir = Path(tmpdir)
      qc_f = tmpdir / 'qc.csv.gz'
      coldata = tmpdir / 'coldata.csv.gz'
      out_f = tmpdir / 'joint_qc.csv.gz'
      self.write_csv(qc_f, {'cell_id': ['a1'], 'log_counts': [2.0]})
      self.write_csv(coldata, {
        'cell_id': ['missing'],
        'sample_id': ['project_a_s1'],
        'project_id': ['project_a'],
      })

      with self.assertRaisesRegex(ValueError, 'QC metrics were not found'):
        build_joint_qc(
          qc_fs=[str(qc_f)],
          project_ids=['project_a'],
          coldata_f=str(coldata),
          out_f=str(out_f),
        )


if __name__ == '__main__':
  unittest.main()
