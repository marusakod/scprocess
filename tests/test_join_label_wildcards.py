import re
import unittest
from pathlib import Path

from snakemake.io import regex_from_filepattern, update_wildcard_constraints


JOIN_SNAKEFILE = Path(__file__).parents[1] / "rules" / "join.smk"


class TestJoinLabelWildcards(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    cls.source = JOIN_SNAKEFILE.read_text()
    match = re.search(
      r'wildcard_constraints:\s+labeller\s*=\s*["\']([^"\']+)["\']',
      cls.source,
    )
    if match is None:
      raise AssertionError("join.smk must constrain the labeller wildcard")
    cls.labeller_constraint = match.group(1)

  def assert_wildcards(self, pattern, filename, labeller, model):
    constrained_pattern = update_wildcard_constraints(
      pattern,
      wildcard_constraints={},
      global_wildcard_constraints={"labeller": self.labeller_constraint},
    )
    match = re.fullmatch(regex_from_filepattern(constrained_pattern), filename)
    self.assertIsNotNone(match, filename)
    self.assertEqual(match.group("labeller"), labeller)
    self.assertEqual(match.group("model"), model)

  def test_underscore_models_remain_intact_in_all_join_label_patterns(self):
    patterns = [
      "tmp_labels_{labeller}_model_{model}_csf_T_NK_join_2026-07-15_{batch}.csv.gz",
      "tmp_labels_{labeller}_model_{model}_csf_T_NK_join_2026-07-15_reused.csv.gz",
      "labels_{labeller}_model_{model}_csf_T_NK_join_2026-07-15.csv.gz",
      "cluster_names_for_shiny_{labeller}_{model}_csf_T_NK_join_0.5_2026-07-15.csv",
      "join_extract_labels_{labeller}_{model}_2026-07-15.log",
      "join_merge_labels_{labeller}_{model}_2026-07-15.log",
      "join_merge_labels_{labeller}_{model}_2026-07-15.benchmark.txt",
      "join_save_cluster_names_{labeller}_{model}_2026-07-15.log",
    ]
    cases = [
      ("celltypist", "Immune_All_Low"),
      ("celltypist", "AIFI_L1"),
      ("scprocess", "human_brain_xgboost_v1"),
    ]

    for pattern in patterns:
      for labeller, model in cases:
        with self.subTest(pattern=pattern, labeller=labeller, model=model):
          filename = pattern.format(
            labeller=labeller,
            model=model,
            batch="project_sample_1",
          )
          self.assert_wildcards(pattern, filename, labeller, model)

  def test_constraint_rejects_unsupported_labeller_capture(self):
    pattern = "labels_{labeller}_model_{model}_csf_T_NK_join_2026-07-15.csv.gz"
    constrained_pattern = update_wildcard_constraints(
      pattern,
      wildcard_constraints={},
      global_wildcard_constraints={"labeller": self.labeller_constraint},
    )
    malformed = "labels_celltypist_Immune_All_model_Low_csf_T_NK_join_2026-07-15.csv.gz"
    self.assertIsNone(re.fullmatch(regex_from_filepattern(constrained_pattern), malformed))

  def assert_rule_uses_label_celltypes_environment(self, rule_name, next_rule_name):
    rule_block = self.source.split(f"rule {rule_name}:", 1)[1].split(
      f"rule {next_rule_name}:", 1
    )[0]
    self.assertIn("'../envs/label_celltypes.yaml'", rule_block)
    self.assertNotIn("'../envs/hvgs.yaml'", rule_block)

  def test_join_label_script_rules_use_label_celltypes_environment(self):
    self.assert_rule_uses_label_celltypes_environment(
      "join_extract_labels", "join_merge_labels"
    )
    self.assert_rule_uses_label_celltypes_environment(
      "join_save_cluster_names", "join_render_html"
    )


if __name__ == "__main__":
  unittest.main()
