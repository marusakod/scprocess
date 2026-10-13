import csv
import gzip
import importlib.util
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _load_merge_module():
  path = ROOT / "scripts" / "merge_doublet_diagnostics.py"
  spec = importlib.util.spec_from_file_location("merge_doublet_diagnostics", path)
  module = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(module)
  return module


def _write_csv_gz(path, rows):
  with gzip.open(path, "wt", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)


def test_doublet_schema_defaults_use_package_defaults():
  schema = json.loads(
    (ROOT / "resources" / "schemas" / "config.schema.json").read_text())
  qc = schema["properties"]["qc"]["properties"]

  assert qc["dbl_dbr_sd"]["default"] is None
  assert qc["dbl_clusters"]["default"] is None
  assert qc["dbl_dbr_sd"]["minimum"] == 0


def test_merge_doublet_diagnostics_keeps_one_header_and_all_rows(tmp_path):
  module = _load_merge_module()
  first = tmp_path / "first.csv.gz"
  second = tmp_path / "second.csv.gz"
  output = tmp_path / "combined.csv.gz"
  _write_csv_gz(
    first, [{"sample_id": "a", "ratio.k50": "0.2", "score": "0.1"}])
  _write_csv_gz(
    second, [{"sample_id": "b", "ratio.k58": "0.8", "score": "0.9"}])

  module.merge_csv_gz([first, second], output)

  with gzip.open(output, "rt", newline="") as handle:
    rows = list(csv.DictReader(handle))
  assert rows == [
    {
      "sample_id": "a", "ratio.k50": "0.2", "score": "0.1",
      "ratio.k58": "",
    },
    {
      "sample_id": "b", "ratio.k50": "", "score": "0.9",
      "ratio.k58": "0.8",
    },
  ]


def test_merge_doublet_diagnostics_cli_avoids_snakemake_serialization():
  rule = (ROOT / "rules" / "qc.smk").read_text()
  script = (ROOT / "scripts" / "merge_doublet_diagnostics.py").read_text()

  assert "script:\n    '../scripts/merge_doublet_diagnostics.py'" not in rule
  assert "python {scprocess_dir}/scripts/merge_doublet_diagnostics.py" in rule
  assert 'main(snakemake)' not in script
  assert 'parser.add_argument("--real-inputs"' in script


def test_qc_rule_passes_optional_doublet_parameters_and_outputs():
  qc_rule = (ROOT / "rules" / "qc.smk").read_text()
  sample_qc = (ROOT / "scripts" / "SampleQC.R").read_text()

  assert re.search(r"dbl_dbr_sd\s*=\s*\{params.dbl_dbr_sd\}", qc_rule)
  assert re.search(r"dbl_clusters\s*=\s*\{params.dbl_clusters\}", qc_rule)
  assert "artificial_doublets" in qc_rule
  assert "dbl_all_dt[type == 'doublet']" in sample_qc
  assert "scDblFinder.threshold" in sample_qc


def test_qc_report_combines_doublet_diagnostics_from_readable_csv_inputs():
  rule = (ROOT / "rules" / "render_htmls.smk").read_text()
  template = (
    ROOT / "resources" / "rmd_templates" / "SampleQC.Rmd.template"
  ).read_text()

  assert "scDblFinder_real_all_runs" in rule
  assert "scDblFinder_artificial_doublets_all_runs" in rule
  assert "scDblFinder_summary_all_runs" in rule
  assert "rule render_html_doublet_id:" not in rule
  assert "readRDS" not in template
  assert "source('code/doublet_id.R')" in template
  assert template.index("## Doublet identification") < template.index(
    "## Sample QC metrics"
  )
  assert "#### Linear counts{.tabset}" in template
  assert "#### Log counts{.tabset}" in template
  assert "heading_level = 5L, log_y = FALSE" in template
  assert "heading_level = 5L, log_y = TRUE" in template
  assert "Doublet-identification provenance" not in template


def test_core_workflow_has_no_standalone_doublet_report_targets():
  core = (ROOT / "rules" / "scprocess.smk").read_text()

  assert "doublet_id_data_outs" in core
  assert "_doublet_id.Rmd" not in core
  assert "_doublet_id.html" not in core
