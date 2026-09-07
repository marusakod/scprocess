import ast
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).parents[1]


def _load_ordering_helper():
  source = (ROOT / "scripts/hvgs.py").read_text()
  tree = ast.parse(source)
  function = next(
    node for node in tree.body
    if isinstance(node, ast.FunctionDef)
    and node.name == "_indices_in_requested_feature_order"
  )
  namespace = {"np": np}
  exec(compile(ast.Module(body=[function], type_ignores=[]), "hvgs.py", "exec"), namespace)
  return namespace["_indices_in_requested_feature_order"]


def test_hvg_rows_are_reordered_to_match_written_feature_labels():
  ordered_indices = _load_ordering_helper()
  observed = np.array(["gene_c", "gene_a", "gene_b"])
  requested = np.array(["gene_b", "gene_c", "gene_a"])
  np.testing.assert_array_equal(ordered_indices(observed, requested), [2, 0, 1])


def test_hvg_row_ordering_rejects_missing_features():
  ordered_indices = _load_ordering_helper()
  with pytest.raises(ValueError, match="requested HVGs are absent"):
    ordered_indices(["gene_a"], ["gene_a", "gene_b"])


def test_hvg_matrix_rules_track_the_builder_script():
  rules = (ROOT / "rules/hvgs.smk").read_text()
  create_hvg = rules.split("rule create_hvg_matrix:", 1)[1].split(
    "rule create_doublets_hvg_matrix:", 1
  )[0]
  create_doublets = rules.split("rule create_doublets_hvg_matrix:", 1)[1].split(
    "if __name__", 1
  )[0]
  assert "script_f" in create_hvg
  assert "python3 {input.script_f} create_hvg_matrix" in create_hvg
  assert "script_f" in create_doublets
  assert "python3 {input.script_f} create_doublets_matrix" in create_doublets
