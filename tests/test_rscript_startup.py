import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_computational_rules_do_not_load_project_r_profiles():
  bare_rscript = re.compile(r"(?<![A-Za-z0-9_-])Rscript\s+-e\b")
  offenders = []
  for rule_file in sorted((ROOT / "rules").glob("*.smk")):
    for line_number, line in enumerate(rule_file.read_text().splitlines(), 1):
      if bare_rscript.search(line):
        offenders.append(f"{rule_file.name}:{line_number}")

  assert offenders == [], (
    "Use 'Rscript --vanilla -e' so rule jobs cannot load project or user "
    f"startup profiles: {', '.join(offenders)}")
