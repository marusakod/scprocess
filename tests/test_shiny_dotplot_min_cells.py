import re
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]


class TestShinyDotplotMinCells(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    cls.source = (ROOT / "resources/shiny/modules/explore_genes.R").read_text()

  def test_min_cells_is_an_advanced_option_defaulting_to_10(self):
    advanced = self.source.split('ns("gene_show_advanced"))', 1)[1]
    advanced = advanced.split("uiOutput(ns(\"gene_subset_picker\"))", 1)[0]
    match = re.search(
      r'numericInput\(\s*inputId = ns\("pb_min_cells"\).*?value\s*=\s*(\d+)',
      advanced, flags=re.S,
    )
    self.assertIsNotNone(match)
    self.assertEqual(match.group(1), "10")

  def test_dotplot_data_uses_the_input_not_a_fixed_cutoff(self):
    dotplot = self.source.split("gene_dotplot_dt <- reactive", 1)[1]
    dotplot = dotplot.split("gene_dotplot_anno_dt", 1)[0]
    self.assertIn("input$pb_min_cells", dotplot)
    self.assertNotRegex(dotplot, r"n_cells\s*>=\s*10\b")
    self.assertIn("n_cells >= min_cells", dotplot)


if __name__ == "__main__":
  unittest.main()
