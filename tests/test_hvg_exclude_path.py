from pathlib import Path

import polars as pl

from scripts.scprocess_utils import _check_hvg_parameters


def _make_config(proj_dir, exclude_f):
  # no ambient-gene detection and the 'all' method, so only the exclusion file is checked
  gene_info_f = proj_dir / 'gene_info.tsv'
  pl.DataFrame({'gene_id': ['G1', 'G2'], 'symbol': ['A', 'B']}).write_csv(gene_info_f, separator = '\t')

  return {
    'project'       : {'proj_dir': proj_dir},
    'multiplexing'  : {'demux_type': 'none'},
    'mapping_af'    : {'gene_info_f': gene_info_f},
    'hvg'           : {
      'hvg_method': 'all', 'hvg_exclude_ambient_genes': False, 'hvg_chunk_size': 10,
      'hvg_exclude_from_file': exclude_f},
  }


def test_relative_exclude_file_is_resolved_from_project_dir(tmp_path):
  pl.DataFrame({'symbol': ['A']}).write_csv(tmp_path / 'exclude.csv')

  config = _check_hvg_parameters(_make_config(tmp_path, 'exclude.csv'))

  assert Path(config['hvg']['hvg_exclude_from_file']) == tmp_path / 'exclude.csv'
