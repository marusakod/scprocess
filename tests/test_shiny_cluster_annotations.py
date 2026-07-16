import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]


class TestShinyClusterAnnotations(unittest.TestCase):
  def test_fdr_positions_are_joined_by_cluster_name(self):
    source = (ROOT / "resources/shiny/utils/plots_genes.R").read_text()
    self.assertIn("anno_ypos_dt = dotplot_dt[", source)
    self.assertIn(
      'merge(anno_dt, anno_ypos_dt, by = "cluster", all = FALSE, sort = FALSE)',
      source,
    )
    self.assertNotIn("names(anno_ypos) == anno_dt$cluster", source)


if __name__ == "__main__":
  unittest.main()
