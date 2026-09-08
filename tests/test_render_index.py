from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_index_zoom_parser_accepts_integer_and_decimal_resolutions():
  renderer = (ROOT / 'scripts' / 'render_htmls.R').read_text()

  assert '(.+_\\\\d+(?:\\\\.\\\\d+)?)\\\\.html$' in renderer
  assert 'if (any(is.na(params_mat[, 1])))' in renderer
  assert 'if (any(is.na(params_mat)))' not in renderer
