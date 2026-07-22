import gzip
import tempfile
import unittest
from pathlib import Path

import polars as pl

from scripts.join import select_joint_hvgs


class TestJoinHvgConsensus(unittest.TestCase):
  def write_hvgs(self, path, rows):
    df = pl.DataFrame(rows)
    with gzip.open(path, 'wb') as fh:
      df.write_csv(fh)

  def test_prioritizes_project_support_then_median_final_rank(self):
    with tempfile.TemporaryDirectory() as tmpdir:
      tmpdir = Path(tmpdir)
      hvg_a = tmpdir / 'hvg_a.csv.gz'
      hvg_b = tmpdir / 'hvg_b.csv.gz'
      hvg_c = tmpdir / 'hvg_c.csv.gz'
      out_f = tmpdir / 'joint_hvgs.csv.gz'

      self.write_hvgs(hvg_a, {
        'gene_id': ['a', 'b', 'c', 'dirty'],
        'highly_variable': [True, True, True, False],
        'highly_variable_nbatches': [3, 2, 1, 0],
        'highly_variable_rank': [1.0, 2.0, 3.0, None],
      })
      self.write_hvgs(hvg_b, {
        'gene_id': ['b', 'c', 'd', 'dirty'],
        'highly_variable': [True, True, True, False],
        'highly_variable_nbatches': [3, 2, 1, 0],
        'highly_variable_rank': [1.0, 2.0, 3.0, None],
      })
      self.write_hvgs(hvg_c, {
        'gene_id': ['a', 'c', 'd', 'dirty'],
        'highly_variable': [True, True, True, False],
        'highly_variable_nbatches': [3, 2, 1, 0],
        'highly_variable_rank': [1.0, 2.0, 3.0, None],
      })

      select_joint_hvgs(
        hvg_fs=[str(hvg_a), str(hvg_b), str(hvg_c)],
        project_ids=['project_a', 'project_b', 'project_c'],
        n_hvgs=3,
        out_f=str(out_f),
      )

      result = pl.read_csv(out_f)
      self.assertEqual(result['gene_id'].to_list(), ['c', 'a', 'b'])
      self.assertEqual(result['n_projects_hvg'].to_list(), [3, 2, 2])
      self.assertEqual(result['median_project_rank'].to_list(), [2.0, 1.0, 1.5])
      self.assertNotIn('dirty', result['gene_id'].to_list())

  def test_gene_id_breaks_consensus_ties_deterministically(self):
    with tempfile.TemporaryDirectory() as tmpdir:
      tmpdir = Path(tmpdir)
      hvg_f = tmpdir / 'hvg.csv.gz'
      out_f = tmpdir / 'joint_hvgs.csv.gz'
      self.write_hvgs(hvg_f, {
        'gene_id': ['z_gene', 'a_gene'],
        'highly_variable': [True, True],
        'highly_variable_nbatches': [1, 1],
        'highly_variable_rank': [1.0, 1.0],
      })

      select_joint_hvgs(
        hvg_fs=[str(hvg_f)],
        project_ids=['project'],
        n_hvgs=2,
        out_f=str(out_f),
      )

      self.assertEqual(pl.read_csv(out_f)['gene_id'].to_list(), ['a_gene', 'z_gene'])

  def test_requires_one_hvg_file_per_project(self):
    with tempfile.TemporaryDirectory() as tmpdir:
      out_f = Path(tmpdir) / 'joint_hvgs.csv.gz'
      with self.assertRaisesRegex(ValueError, 'one HVG file per project'):
        select_joint_hvgs(
          hvg_fs=[],
          project_ids=['project'],
          n_hvgs=1,
          out_f=str(out_f),
        )


if __name__ == '__main__':
  unittest.main()
