import re
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

  def test_resource_lookups_name_defined_rules(self):
    # get_resources() raises at job submission for an undefined rule name, which a
    # dry run does not catch (e.g. the zoom data rule still naming build_zoom_shiny_app)
    defined = set(re.findall(r"^rule (\w+):", self.source, flags=re.M))
    looked_up = set(re.findall(r"get_resources\([^,]+, rules, input, '(\w+)'", self.source))
    looked_up |= set(re.findall(r"_shiny_data_resource\('(\w+)'", self.source))
    self.assertTrue(looked_up)
    self.assertEqual(looked_up - defined, set())

  def test_zoom_data_rule_keeps_old_resource_keys(self):
    data = rule_block(
      self.source, "build_zoom_shiny_data", "configure_zoom_shiny_app"
    )
    self.assertIn(
      "_shiny_data_resource('build_zoom_shiny_data', 'build_zoom_shiny_app'", data
    )

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
    self.assertIn('column(width = 7, uiOutput(ns("get_cluster_overview_box")))', module)
    # Legend wrapping and device geometry are checked by test_shiny_umap_layout.R.
    self.assertIn(
      "data_sentinel_f = f'{docs_dir}/shiny_zoom_{{zoom_name}}/.shiny_data_built_{DATE_STAMP}'",
      self.source,
    )


if __name__ == "__main__":
  unittest.main()
