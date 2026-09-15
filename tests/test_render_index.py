from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_index_zoom_parser_accepts_integer_and_decimal_resolutions():
  renderer = (ROOT / 'scripts' / 'render_htmls.R').read_text()

  assert '(.+_\\\\d+(?:\\\\.\\\\d+)?)\\\\.html$' in renderer
  assert 'if (any(is.na(params_mat[, 1])))' in renderer
  assert 'if (any(is.na(params_mat)))' not in renderer


def test_index_groups_celltype_labelling_with_whole_dataset_processing():
  template = (ROOT / 'resources' / 'rmd_templates' / 'index.Rmd.template').read_text()
  renderer = (ROOT / 'scripts' / 'render_htmls.R').read_text()

  assert '## Whole dataset processing' in template
  assert '## Core steps' not in template
  assert '${lbls_title}' not in template
  assert 'Cell type annotation' not in renderer
  assert 'lbls_title' not in renderer
  assert template.index('${marker_genes_link}') < template.index('${label_celltypes_link}')
  assert template.index('${label_celltypes_link}') < template.index('${zoom_title}')
