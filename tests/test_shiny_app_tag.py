import json
import re
import unittest
from pathlib import Path

from scripts.scprocess_utils import get_shiny_app_tag


ROOT = Path(__file__).parents[1]


class TestShinyAppTag(unittest.TestCase):
  def test_context_specific_defaults(self):
    self.assertEqual(
      get_shiny_app_tag({"join": {"name": "joined_cells"}}),
      "joined_cells",
    )
    self.assertEqual(
      get_shiny_app_tag({"project": {"short_tag": "source_project"}}),
      "source_project",
    )

  def test_configured_tag_overrides_default(self):
    self.assertEqual(
      get_shiny_app_tag({
        "join": {"name": "joined_cells"},
        "shiny": {"app_tag": "T_NK.app-v2"},
      }),
      "T_NK.app-v2",
    )

  def test_project_and_join_schemas_use_machine_safe_pattern(self):
    for schema_name in ("config.schema.json", "join.schema.json"):
      schema = json.loads((ROOT / "resources/schemas" / schema_name).read_text())
      app_tag = schema["properties"]["shiny"]["properties"]["app_tag"]
      pattern = re.compile(app_tag["pattern"])
      self.assertIsNotNone(pattern.fullmatch("csf_T_NK.join-v2"))
      self.assertIsNone(pattern.fullmatch("csf T/NK"))

  def test_main_rule_uses_tag_for_directory_sentinel_data_and_log(self):
    source = (ROOT / "rules/shiny.smk").read_text()
    self.assertIn("_main_app_tag     = get_shiny_app_tag(config)", source)
    self.assertIn("_main_shiny_dir   = f'{docs_dir}/shiny_{_main_app_tag}'", source)
    self.assertIn("sentinel_f    = f'{_main_shiny_dir}/.shiny_built_{DATE_STAMP}'", source)
    self.assertIn("deploy_dir    = _main_shiny_dir", source)
    self.assertIn("app_tag       = _main_app_tag", source)
    self.assertIn("build_shiny_app_{_main_app_tag}_{DATE_STAMP}.log", source)


if __name__ == "__main__":
  unittest.main()
