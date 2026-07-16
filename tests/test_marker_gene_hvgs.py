import re
import unittest
from pathlib import Path


MARKER_GENES_R = Path(__file__).parents[1] / 'scripts' / 'marker_genes.R'


class TestMarkerGeneHvgs(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    source = MARKER_GENES_R.read_text()
    match = re.search(
      r'calc_hvgs_pseudobulk\s*<-\s*function\(.*?\n\}',
      source,
      flags=re.DOTALL,
    )
    if match is None:
      raise AssertionError('calc_hvgs_pseudobulk function not found')
    cls.source = match.group(0)

  def test_one_cell_pseudobulks_are_filtered_only_in_hvg_function(self):
    self.assertIn('min_hvg_cells = 2L', self.source)
    self.assertIn('vst_dt        = cpms_dt[n_cells >= min_hvg_cells]', self.source)
    self.assertIn('bulk_mat  = vst_dt', self.source)

  def test_sparse_size_factor_fallback_is_specific(self):
    self.assertIn("method = 'ratio'", self.source)
    self.assertIn("'every gene contains at least one zero'", self.source)
    self.assertIn("estimateSizeFactors(dds, type = 'poscounts')", self.source)
    self.assertIn("method = 'poscounts'", self.source)


if __name__ == '__main__':
  unittest.main()
