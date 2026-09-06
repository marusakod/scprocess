import tempfile
import unittest
from pathlib import Path

import polars as pl
import yaml

from scripts.join import build_joint_sample_metadata
from scripts.scprocess_utils import _load_sample_metadata


ROOT = Path(__file__).parents[1]


def rule_block(source, rule_name, next_rule_name):
  return source.split(f"rule {rule_name}:", 1)[1].split(
    f"rule {next_rule_name}:", 1
  )[0]


class TestJointSampleMetadata(unittest.TestCase):
  def test_combines_mixed_schemas_and_prefixes_sample_ids(self):
    with tempfile.TemporaryDirectory() as tmpdir:
      tmpdir = Path(tmpdir)
      meta_a = tmpdir / "meta_a.csv"
      meta_b = tmpdir / "meta_b.csv"
      out_f = tmpdir / "joint_sample_meta.csv"

      pl.DataFrame({
        "sample_id": ["s1"],
        "patient_id": [101],
        "diagnosis": ["control"],
      }).write_csv(meta_a)
      pl.DataFrame({
        "sample_id": ["s2"],
        "patient_id": ["P102"],
        "therapy": ["treated"],
        "bad_sample_id": [True],
      }).write_csv(meta_b)

      build_joint_sample_metadata(
        project_ids=["project_a", "project_b"],
        sample_meta_fs=[str(meta_a), str(meta_b)],
        out_sample_meta_f=str(out_f),
      )

      result = pl.read_csv(out_f)
      self.assertEqual(
        result["sample_id"].to_list(), ["project_a_s1", "project_b_s2"]
      )
      self.assertEqual(result["project_id"].to_list(), ["project_a", "project_b"])
      self.assertEqual(result["bad_sample_id"].to_list(), [False, True])
      self.assertEqual(result.schema["patient_id"], pl.String)
      self.assertEqual(result["patient_id"].to_list(), ["101", "P102"])
      self.assertIn("diagnosis", result.columns)
      self.assertIn("therapy", result.columns)

  def test_config_metadata_loader_promotes_mixed_patient_ids_to_string(self):
    with tempfile.TemporaryDirectory() as tmpdir:
      tmpdir = Path(tmpdir)
      project_entries = {}
      for project_id, patient_id in (("consonance", 101), ("piehl", "P102")):
        project_dir = tmpdir / project_id
        project_dir.mkdir()
        metadata_f = project_dir / "metadata.csv"
        config_f = project_dir / "config.yaml"
        pl.DataFrame({
          "sample_id": [f"{project_id}_sample"],
          "patient_id": [patient_id],
        }).write_csv(metadata_f)
        with config_f.open("w") as fh:
          yaml.safe_dump({
            "project": {
              "proj_dir": str(project_dir),
              "sample_metadata": str(metadata_f),
            }
          }, fh)
        project_entries[project_id] = {"config": str(config_f)}

      result = _load_sample_metadata({
        "join": {"name": "mixed_ids"},
        "projects": project_entries,
      })

      self.assertEqual(result.schema["patient_id"], pl.String)
      self.assertEqual(result["patient_id"].to_list(), ["101", "P102"])


class TestMetadataDependencyGraph(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    cls.join = (ROOT / "rules/join.smk").read_text()
    cls.shiny = (ROOT / "rules/shiny.smk").read_text()
    cls.reports = (ROOT / "rules/render_htmls.smk").read_text()
    cls.zoom = (ROOT / "rules/zoom.smk").read_text()

  def test_join_analysis_does_not_depend_on_descriptive_metadata(self):
    coldata = rule_block(self.join, "join_build_coldata", "join_build_sample_metadata")
    integration = rule_block(self.join, "join_integration", "join_build_h5ads_yaml")
    self.assertNotIn("sample_meta", coldata)
    self.assertNotIn("sample_meta", integration)
    self.assertNotIn("sample_qc_f", integration)

  def test_join_metadata_has_its_own_rule_and_report_input(self):
    metadata = rule_block(self.join, "join_build_sample_metadata", "join_build_qc")
    report = rule_block(self.join, "join_render_html", "join_run_train_xgboost")
    self.assertIn("sample_meta_fs = SAMPLE_META_FS", metadata)
    self.assertIn("sample_meta_f = joint_sample_meta_f", metadata)
    self.assertIn("sample_meta_f = joint_sample_meta_f", report)
    self.assertIn("'{input.sample_meta_f}'", report)

  def test_shiny_apps_declare_sample_metadata_as_input(self):
    main = rule_block(self.shiny, "build_shiny_data", "configure_shiny_app")
    zoom = rule_block(
      self.shiny, "build_zoom_shiny_data", "configure_zoom_shiny_app"
    )
    for block in (main, zoom):
      self.assertIn("sample_meta_f = _sample_meta_f", block)
      self.assertIn("'{input.sample_meta_f}'", block)
      self.assertNotIn("'{params.sample_meta_f}'", block)

  def test_descriptive_reports_declare_metadata_as_input(self):
    qc = rule_block(self.reports, "render_html_qc", "render_html_integration")
    markers = rule_block(
      self.reports, "render_html_marker_genes", "render_html_label_celltypes"
    )
    zoom = self.zoom.split("rule render_html_zoom:", 1)[1]
    for block in (qc, markers, zoom):
      self.assertIn("metadata_f", block.split("output:", 1)[0])
      self.assertIn("'{input.metadata_f}'", block)


if __name__ == "__main__":
  unittest.main()
