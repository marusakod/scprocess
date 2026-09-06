import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]


def rule_block(source, rule_name, next_rule_name=None):
  block = source.split(f"rule {rule_name}:", 1)[1]
  if next_rule_name is not None:
    block = block.split(f"rule {next_rule_name}:", 1)[0]
  return block


class TestShinyRuleSplit(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    cls.source = (ROOT / "rules/shiny.smk").read_text()

  def test_main_cosmetic_inputs_only_feed_configure_rule(self):
    data = rule_block(self.source, "build_shiny_data", "configure_shiny_app")
    configure = rule_block(
      self.source, "configure_shiny_app", "build_zoom_shiny_data"
    )
    self.assertNotIn("annotation_csv_f = _annotation_csv_f", data)
    self.assertNotIn("home_md_f = _home_md_f", data)
    self.assertIn("annotation_csv_f  = _annotation_csv_f", configure)
    self.assertIn("home_md_f         = _home_md_f", configure)
    self.assertIn("deploy_inputs   = _shiny_deploy_inputs", configure)

  def test_zoom_cosmetic_inputs_only_feed_configure_rule(self):
    data = rule_block(
      self.source, "build_zoom_shiny_data", "configure_zoom_shiny_app"
    )
    configure = rule_block(self.source, "configure_zoom_shiny_app")
    self.assertNotIn("get_zoom_shiny_optional_inputs", data)
    self.assertIn("get_zoom_shiny_optional_inputs", configure)
    self.assertIn("deploy_inputs   = _shiny_deploy_inputs", configure)

  def test_final_rules_depend_on_data_sentinels(self):
    self.assertIn(
      "data_sentinel_f = f'{_main_shiny_dir}/.shiny_data_built_{DATE_STAMP}'",
      self.source,
    )

  def test_join_dispatch_uses_shiny_snakefile_for_configure_target(self):
    cli_source = (ROOT / "scprocess").read_text()
    self.assertIn('if rule == "configure_shiny_app":', cli_source)
    self.assertNotIn('if rule == "build_shiny_app":', cli_source)

  def test_cluster_overview_reserves_space_for_long_labels(self):
    module = (ROOT / "resources/shiny/modules/explore_genes.R").read_text()
    plots = (ROOT / "resources/shiny/utils/plots_clusters.R").read_text()
    self.assertIn('column(width = 7, uiOutput(ns("get_cluster_overview_box")))', module)
    self.assertIn("legend_width  = max(1.6, min(2.6, max_name_chars / 14))", plots)
    self.assertIn('coord_cartesian(clip = "off")', plots)
    self.assertIn("widths = c(5, legend_width)", plots)
    self.assertIn(
      "data_sentinel_f = f'{docs_dir}/shiny_zoom_{{zoom_name}}/.shiny_data_built_{DATE_STAMP}'",
      self.source,
    )


if __name__ == "__main__":
  unittest.main()
