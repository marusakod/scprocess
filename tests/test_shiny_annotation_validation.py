import re
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
SHINY_R = ROOT / "scripts" / "shiny.R"


class TestShinyAnnotationValidation(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    source = SHINY_R.read_text()
    match = re.search(
      r"\.validate_shiny_annotation\s*<-\s*function\(.*?\n\}",
      source,
      flags=re.DOTALL,
    )
    if match is None:
      raise AssertionError(".validate_shiny_annotation function not found")
    cls.source = source
    cls.validator = match.group(0)

  def test_rejects_cluster_set_mismatches(self):
    self.assertIn("setdiff(annot_clusters, data_clusters)", self.validator)
    self.assertIn("setdiff(data_clusters, annot_clusters)", self.validator)
    self.assertIn("Present only in annotation_csv:", self.validator)
    self.assertIn("Present only in analysis data:", self.validator)
    self.assertIn("different clustering resolution or analysis run", self.validator)

  def test_rejects_invalid_or_ambiguous_rows(self):
    self.assertIn("missing or blank cluster values", self.validator)
    self.assertIn("missing or blank cluster_name values", self.validator)
    self.assertIn("duplicated cluster values", self.validator)
    self.assertIn("duplicated cluster_name values", self.validator)

  def test_validation_precedes_expensive_matrix_loading_and_cleanup(self):
    validation_call = self.source.index("annot <- .validate_shiny_annotation(")
    umap_subsampling = self.source.index('message("Downsampling UMAP to ",')
    stale_cleanup = self.source.index("stale <- list.files(")
    h5ad_loading = self.source.index('message("Reading h5ad files")')
    self.assertLess(validation_call, umap_subsampling)
    self.assertLess(validation_call, stale_cleanup)
    self.assertLess(validation_call, h5ad_loading)

  def test_full_analysis_clusters_are_not_filtered_by_marker_results(self):
    self.assertIn(
      "analysis_clusters <- sort(unique(as.character(all_meta$cluster)))",
      self.source,
    )
    self.assertNotIn(
      ".[cluster %in% unique(fread(mkrs_f)$cluster), ]",
      self.source,
    )

  def test_umap_subsampling_is_stratified_by_cluster(self):
    self.assertIn('stratify_by = "cluster"', self.source)
    self.assertIn("mandatory_idx <- vapply(group_idx", self.source)

  def test_previous_annotation_is_removed_when_no_annotation_is_configured(self):
    self.assertIn(
      "if (is.null(annot) && file.exists(annotation_dest))",
      self.source,
    )
    self.assertIn("unlink(annotation_dest)", self.source)


if __name__ == "__main__":
  unittest.main()
