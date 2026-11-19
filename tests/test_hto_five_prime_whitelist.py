import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import polars as pl
import yaml

from scripts import mapping


# whitelist lookup rows as in cellranger_whitelists.csv: 3v3 has HTO files, 5v3 has blanks
LOOKUP_DT = pl.DataFrame({
  "chemistry":      ["3v3", "5v3"],
  "gex_barcodes_f": ["gex_3v3.txt", "gex_5v3.txt"],
  "hto_barcodes_f": ["hto_3v3.txt", ""],
  "translation_f":  ["trans_3v3.txt", ""],
})


def _write_stats(whitelist_f):
  chem = {
    "whitelist_f": whitelist_f, "r1_length": 28, "max_overlap": 0.9,
    "cell_counts_fw": 1, "cell_counts_rc": 2, "exp_ori": "rc",
    "sample_tenx_chemistry": "5v3", "af_chemistry": "10xv3",
  }
  with tempfile.TemporaryDirectory() as tmp:
    mapping._write_chemistry_stats(tmp, "pool_1", chem, LOOKUP_DT)
    return yaml.safe_load((Path(tmp) / "chemistry_statistics.yaml").read_text())


class TestHtoFivePrimeWhitelist(unittest.TestCase):
  def test_five_prime_uses_gex_whitelist_and_no_translation(self):
    stats = _write_stats("/wl/gex_5v3.txt")
    self.assertEqual(stats["selected_hto_whitelist"], "/wl/gex_5v3.txt")
    self.assertEqual(stats["selected_translation_f"], "none")

  def test_three_prime_keeps_hto_whitelist_and_translation(self):
    stats = _write_stats("/wl/gex_3v3.txt")
    self.assertEqual(stats["selected_hto_whitelist"], "/wl/hto_3v3.txt")
    self.assertEqual(stats["selected_translation_f"], "/wl/trans_3v3.txt")

  def test_hto_validation_accepts_gex_whitelist_for_five_prime(self):
    with patch.object(mapping, "_get_r1_read_length", return_value=28):
      chem = mapping._validate_chemistry(
        ["r1.fastq.gz"], "10xv3", "rc", "/wl/gex_5v3.txt", LOOKUP_DT, "hto")
    self.assertEqual(chem["sample_tenx_chemistry"], "5v3")
    self.assertEqual(chem["whitelist_f"], "/wl/gex_5v3.txt")


if __name__ == "__main__":
  unittest.main()
