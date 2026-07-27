import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
MARKER_GENES_R = ROOT / 'scripts' / 'marker_genes.R'
ADAPTER_R = ROOT / 'scripts' / 'marker_genes_edger_bp.R'
JOIN_RULES = ROOT / 'rules' / 'join.smk'
EQUIVALENCE_TEST = ROOT / 'tests' / 'test_marker_genes_edger_bp.R'


class TestEdgerBpBackend(unittest.TestCase):
  def test_join_rule_loads_adapter(self):
    source = JOIN_RULES.read_text()
    self.assertIn("source('{scprocess_dir}/scripts/marker_genes_edger_bp.R')", source)
    self.assertIn("use_bpcells   = 'True'", source)

  def test_pseudobulk_counts_remain_disk_backed_for_de(self):
    source = MARKER_GENES_R.read_text()
    self.assertNotIn('as.matrix(pb_one)', source)
    self.assertIn('BPCells::write_matrix_dir(', source)
    self.assertIn('edger.bp::bp_ensure_storage_order(', source)
    self.assertIn('edger_bp_counts = edger_bp_counts', source)

  def test_adapter_uses_public_streaming_api(self):
    source = ADAPTER_R.read_text()
    self.assertIn('edger.bp::bp_check_versions(error = TRUE)', source)
    self.assertIn('edger.bp::bp_prepare_dge(', source)
    self.assertIn('edger.bp::bp_estimate_disp_stream(', source)
    self.assertIn('edger.bp::bp_glm_lrt(', source)
    self.assertIn('edger.bp::bp_glm_treat(', source)
    self.assertIn('stats::p.adjust(out$PValue, method = "BH")', source)

  def test_dense_equivalence_fixture_is_present(self):
    source = EQUIVALENCE_TEST.read_text()
    self.assertIn('estimateDisp(dge, design = dispersion_design)', source)
    self.assertIn('glmQLFit(', source)
    self.assertIn('glmTreat(fit, coef = "selectedTRUE")$table', source)
    self.assertIn('all.equal(observed$PValue, expected$PValue', source)


if __name__ == '__main__':
  unittest.main()
