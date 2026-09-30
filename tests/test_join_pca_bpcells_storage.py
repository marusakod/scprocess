import json
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]


class TestJoinPcaBpcellsStorage(unittest.TestCase):

  @classmethod
  def setUpClass(cls):
    cls.rules = (ROOT / "rules/join.smk").read_text()
    cls.script = (ROOT / "scripts/join_pca.R").read_text()
    cls.schema = json.loads(
      (ROOT / "resources/schemas/join.schema.json").read_text())

  def test_native_matrix_is_a_tracked_directory_output(self):
    self.assertIn("rule join_build_bpcells_matrix:", self.rules)
    self.assertIn("counts_bp_dir = directory(joint_bpcells_counts_dir)", self.rules)
    self.assertIn("counts_bp_dir = joint_bpcells_counts_dir", self.rules)

  def test_counts_are_bitpacked_as_uint32(self):
    self.assertIn('convert_matrix_type(mat, "uint32_t")', self.script)
    self.assertIn("write_matrix_dir(", self.script)
    self.assertIn("open_matrix_dir(counts_bp_dir)", self.script)

  def test_conversion_keeps_hdf5_input_readonly(self):
    self.assertIn('H5File$new(h5_path, mode = "r")', self.script)
    self.assertNotIn('mode = "r+"', self.script)
    self.assertIn('on.exit(h5_readonly$close_all()', self.script)
    self.assertIn('regenerate this matrix with join_build_matrix', self.script)

  def test_pca_queues_centering_and_scaling_before_threaded_svds(self):
    transform = "mat_scaled <- (mat_norm - gene_means) / gene_sds"
    solver = "svd_result <- BPCells::svds("
    self.assertIn(transform, self.script)
    self.assertIn(solver, self.script)
    self.assertLess(self.script.index(transform), self.script.index(solver))
    self.assertIn("opts = list(tol = 1e-5)", self.script)
    self.assertIn("threads = n_threads", self.script)
    self.assertIn("pc_scores <- sweep(svd_result$v", self.script)
    self.assertIn("n_threads     =  {threads}", self.rules)

  def test_schema_has_conversion_resources_and_longer_pca_default(self):
    resources = self.schema["properties"]["resources"]["properties"]
    self.assertIn("gb_join_build_bpcells_matrix", resources)
    self.assertIn("mins_join_build_bpcells_matrix", resources)
    self.assertEqual(resources["mins_join_pca"]["default"], 180)


if __name__ == "__main__":
  unittest.main()
