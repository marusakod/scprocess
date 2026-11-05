"""join.smk imports scprocess_utils names explicitly; every helper it calls must be imported."""
import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_join_imports_every_scprocess_utils_helper_it_calls():
  # Public functions defined in scprocess_utils.py
  utils_tree = ast.parse((ROOT / 'scripts' / 'scprocess_utils.py').read_text())
  utils_functions = {node.name for node in utils_tree.body
    if isinstance(node, ast.FunctionDef) and not node.name.startswith('_')}

  # Names in the explicit `from scprocess_utils import (...)` list
  join_text = (ROOT / 'rules' / 'join.smk').read_text()
  import_block = re.search(r'from scprocess_utils import \((.*?)\)', join_text, flags = re.S).group(1)
  imported = {name.strip() for name in import_block.split(',') if name.strip()}

  # Helpers called anywhere in join.smk, e.g. `get_prediction_chunks(`
  called = {name for name in utils_functions if re.search(r'\b%s\(' % name, join_text)}
  assert called <= imported, f'join.smk calls but does not import: {sorted(called - imported)}'
