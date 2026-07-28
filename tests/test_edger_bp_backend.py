import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
MARKER_GENES_R = ROOT / 'scripts' / 'marker_genes.R'
ADAPTER_R = ROOT / 'scripts' / 'marker_genes_edger_bp.R'
JOIN_RULES = ROOT / 'rules' / 'join.smk'
EQUIVALENCE_TEST = ROOT / 'tests' / 'test_marker_genes_edger_bp.R'
EDGER_BP_ENV = ROOT / 'envs' / 'rlibs_bpcells.yaml'


class TestEdgerBpBackend(unittest.TestCase):
  def test_environment_uses_public_release_channel(self):
    source = EDGER_BP_ENV.read_text()
    self.assertIn(
      'https://raw.githubusercontent.com/wmacnair/edger.bp/conda-channel',
      source,
    )
    self.assertIn('r-edger-bp=0.1.0=r45_0', source)
    self.assertNotIn('edger.bp-work', source)
    self.assertNotIn('r-edger-bp=0.0.0.9000', source)

  def test_join_rule_loads_adapter(self):
    source = JOIN_RULES.read_text()
    self.assertIn("source('{scprocess_dir}/scripts/marker_genes_edger_bp.R')", source)
    self.assertIn('rule join_make_pseudobulks:', source)
    self.assertIn('rule join_prepare_pseudobulks:', source)
    self.assertIn('rule join_calc_hvgs:', source)
    self.assertIn('rule join_marker_genes:', source)

  def test_pseudobulk_counts_remain_disk_backed_for_de(self):
    marker_source = MARKER_GENES_R.read_text()
    adapter_source = ADAPTER_R.read_text()
    self.assertNotIn('as.matrix(pb_one)', adapter_source)
    self.assertIn('calculate_marker_genes_bpcells(', marker_source)
    self.assertIn('BPCells::write_matrix_dir(', adapter_source)
    self.assertIn('file.path(pb_dir, "counts_col")', adapter_source)
    self.assertIn('file.path(pb_dir, "counts_row")', adapter_source)
    self.assertIn('edger.bp::bp_ensure_storage_order(', adapter_source)
    self.assertIn('h5ad_batches = names(h5ad_paths)', adapter_source)
    self.assertIn('ok_cells = intersect(colnames(mat), int_dt$cell_id)', adapter_source)

  def test_join_pseudobulk_is_a_directory_output(self):
    source = JOIN_RULES.read_text()
    self.assertIn('pb_dir = directory(pb_dir)', source)
    self.assertIn('.bpcells"', source)
    self.assertNotIn('pb_f          = pb_f,', source)

  def test_preparation_is_persisted_without_mutating_pseudobulk_store(self):
    source = JOIN_RULES.read_text()
    adapter_source = ADAPTER_R.read_text()
    self.assertIn('prepared_coldata_f = pb_prepared_coldata_f', source)
    self.assertIn('prepared_rowdata_f = pb_prepared_rowdata_f', source)
    self.assertIn('norm_factor = prepared$norm.factors', adapter_source)
    self.assertIn('gene_id = rownames(prepared$counts)', adapter_source)
    self.assertNotIn('marker_keep = keep', adapter_source)

  def test_variability_ranking_is_blockwise_logcpm(self):
    source = ADAPTER_R.read_text()
    self.assertIn('calc_pseudobulk_variability_bpcells <- function(', source)
    self.assertIn('blocks = .bp_row_blocks(nrow(y), block_size)', source)
    self.assertIn('block = as.matrix(y[rows, , drop = FALSE])', source)
    self.assertIn('logcpm_var = variances', source)
    self.assertIn('"tmm_logcpm_variance"', source)
    self.assertNotIn('DESeq2::vst(', source)

  def test_report_uses_compact_plot_data_not_pseudobulk(self):
    source = JOIN_RULES.read_text()
    self.assertIn('rule join_marker_plot_data:', source)
    self.assertIn('pb_plot_data_f = pb_plot_data_f', source)
    self.assertNotIn("pb_f             = '{params.pb_f}'", source)

  def test_adapter_uses_public_streaming_api(self):
    source = ADAPTER_R.read_text()
    self.assertIn('edger.bp::bp_check_versions(error = TRUE)', source)
    self.assertIn('edger.bp::bp_prepare_dge(', source)
    self.assertIn('workers = n_cores', source)
    self.assertNotIn('prepare_workers', source)
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
