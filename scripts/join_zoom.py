"""Select joined-zoom cells and calculate sample-aware HVGs from clean H5ADs.

Expression is read in bounded chunks. Counts from the same biological sample
are combined across files before fitting the existing Seurat-style VST trend.
"""
import argparse
import contextlib
import gzip
import json
import math
from pathlib import Path

import polars as pl
import yaml

# CSC reads use at most 256 genes at once. The additional bound targets at most
# 1,024 × chunk_size matrix positions, with a minimum of one gene per read.
MAX_CSC_GENE_BLOCK = 256
CSC_VALUES_PER_CHUNK_CELL = 1024
MIN_VARIANCE_CELLS = 2


def prepare_sources(integration_f, coldata_f, metadata_f, labels_f, qc_f,
  selection, qc, batches, project_ids, selected_f, summary_f,
  integrated_fs, sample_meta_fs, qc_all_fs, h5ads_yaml_fs,
  int_n_dims, int_batch_var):
  """Select cells, verify their source files and write per-project input tables."""

  # Integration, coldata, labels and QC have one row per cell. Metadata has one
  # row per sample and supplies descriptive variables absent from cell tables.
  integration_df = pl.read_csv(integration_f)
  coldata_df     = pl.read_csv(coldata_f)
  metadata_df    = pl.read_csv(metadata_f)
  labels_df      = pl.read_csv(labels_f)
  qc_df          = pl.read_csv(qc_f)
  cells_df, summary_df = select_cells(
    parent_df = integration_df, coldata_df = coldata_df, metadata_df = metadata_df,
    labels_df = labels_df, qc_df = qc_df, selection = selection,
    thresholds = qc, n_dims = int_n_dims)

  # These variables must be available for every cell passed to integration.
  batch_columns = [int_batch_var] if isinstance(int_batch_var, str) else int_batch_var
  for column in batch_columns:
    if (column not in cells_df.columns or cells_df[column].null_count()
        or (cells_df[column].cast(pl.String) == '').any()):
      raise ValueError(f'Integration batch variable {column!r} is missing for selected cells')

  h5ad_lookup = _write_batch_cell_tables(cells_df, batches, project_ids)
  write_csv(cells_df, selected_f)
  write_csv(summary_df, summary_f)

  # Preserve all original QC columns, including those not used for selection
  # and only the samples/cells represented in the selected subset.
  qc_df = pl.read_csv(qc_f).join(cells_df.select('cell_id'), on = 'cell_id', how = 'semi')
  metadata_df = pl.read_csv(metadata_f).join(
    cells_df.select('sample_id').unique(), on = 'sample_id', how = 'semi')
  _write_project_tables(
    cells_df, metadata_df, qc_df, h5ad_lookup, project_ids,
    integrated_fs, sample_meta_fs, qc_all_fs, h5ads_yaml_fs)
  print(f"Retained {cells_df.height} cells across {cells_df['sample_id'].n_unique()} samples", flush = True)


def prepare(context, config):
  """Run source preparation through the existing JSON CLI interface."""
  paths   = context['paths']
  project = context['project']
  return prepare_sources(
    integration_f = paths['integration'],
    coldata_f     = paths['coldata'],
    metadata_f    = paths['metadata'],
    labels_f      = paths['labels'],
    qc_f          = paths['qc'],
    selection     = context['selection'],
    qc            = context['qc'],
    batches       = context['batches'],
    project_ids   = project['project_ids'],
    selected_f    = context['selected_f'],
    summary_f     = context['summary_f'],
    integrated_fs = project['integrated_fs'],
    sample_meta_fs = project['sample_meta_fs'],
    qc_all_fs     = project['qc_all_fs'],
    h5ads_yaml_fs = project['h5ads_yaml_fs'],
    int_n_dims    = config['integration']['int_n_dims'],
    int_batch_var = config['integration']['int_batch_var'])


def select_cells(parent_df, coldata_df, metadata_df, labels_df, qc_df, selection, thresholds, n_dims):
  """Select within the parent join and retain its IDs and descriptive metadata."""
  for cells_df, description in (
      (parent_df, 'Parent integration'), (coldata_df, 'Parent coldata'),
      (labels_df, 'Selection labels'), (qc_df, 'Parent QC')):
    unique_cells(cells_df, description)
  for column in ('sample_id', 'project_id'):
    if column not in parent_df.columns or parent_df[column].null_count():
      raise ValueError(f'Parent integration needs non-null {column}')

  # Parent coordinates and clusters take precedence over source-project values.
  extra_columns = [column for column in coldata_df.columns if column not in parent_df.columns]
  cells_df = parent_df.join(
    coldata_df.select(['cell_id', *extra_columns]),
    on = 'cell_id', how = 'left', validate = '1:1')
  if metadata_df['sample_id'].n_unique() != metadata_df.height:
    raise ValueError('Parent sample metadata has duplicate sample_id values')
  extra_columns = [column for column in metadata_df.columns if column not in cells_df.columns]
  cells_df = cells_df.join(
    metadata_df.select(['sample_id', *extra_columns]),
    on = 'sample_id', how = 'left', validate = 'm:1')

  # Check labels only among cells represented in this parent integration.
  label_column = selection['labels_col']
  if label_column not in labels_df.columns:
    raise ValueError(f'Labels column not found: {label_column}')
  labels_df = labels_df.select('cell_id', pl.col(label_column).cast(pl.String).alias('parent_label'))
  labels_df = labels_df.join(parent_df.select('cell_id'), on = 'cell_id', how = 'semi')
  labels_missing = set(selection['sel_labels']) - set(labels_df['parent_label'].drop_nulls())
  if labels_missing:
    raise ValueError(f'Requested labels absent from the parent join: {sorted(labels_missing)}')
  cells_df = cells_df.join(labels_df, on = 'cell_id', how = 'inner', validate = '1:1')
  cells_df = cells_df.filter(pl.col('parent_label').is_in(selection['sel_labels']))
  selected_counts_df = cells_df.group_by('project_id', 'sample_id').len(name = 'n_selected')

  # Record label-selected counts before applying explicit sample exclusions.
  excluded_samples = selection.get('exclude', {}).get('sample_id', [])
  samples_missing = set(excluded_samples) - set(parent_df['sample_id'])
  if samples_missing:
    raise ValueError(f'Excluded sample IDs absent from parent join: {sorted(samples_missing)}')
  cells_df = cells_df.filter(~pl.col('sample_id').is_in(excluded_samples))
  cells_df = _filter_qc_cells(cells_df, qc_df, thresholds)
  qc_counts_df = cells_df.group_by('project_id', 'sample_id').len(name = 'n_after_qc')

  # Variance requires at least two cells per sample, even if qc_min_cells is 1.
  minimum_cells = max(MIN_VARIANCE_CELLS, thresholds['qc_min_cells'])
  cells_df = cells_df.filter(pl.len().over('sample_id') >= minimum_cells)
  retained_counts_df = cells_df.group_by('project_id', 'sample_id').len(name = 'n_retained')
  summary_df = selected_counts_df.join(
    qc_counts_df, on = ['project_id', 'sample_id'], how = 'left')
  summary_df = summary_df.join(
    retained_counts_df, on = ['project_id', 'sample_id'], how = 'left')
  summary_df = summary_df.fill_null(0).sort('project_id', 'sample_id')
  if cells_df.height <= n_dims:
    raise ValueError(f'Only {cells_df.height} cells remain; need more than {n_dims} PCs')

  # Prefix parent embeddings/clusters so new zoom outputs cannot overwrite them.
  parent_columns = {
    column: 'parent_' + column for column in cells_df.columns
    if column.startswith(('RNA_snn_res.', 'pca_', 'hmny_pca_', 'harmony_', 'umap_', 'UMAP_'))
    or column in ('UMAP1', 'UMAP2')}
  if set(parent_columns.values()) & set(cells_df.columns):
    raise ValueError('Parent embedding column names conflict; recursive joined zooms are not supported')
  return cells_df.rename(parent_columns), summary_df


def _filter_qc_cells(cells_df, qc_df, thresholds):
  metric_columns = {
    'counts': 'log_counts', 'feats': 'log_feats',
    'mito': 'logit_mito', 'splice': 'logit_spliced'}
  retained_qc_df = qc_df
  for parameter, value in thresholds.items():
    if value is None or parameter == 'qc_min_cells':
      continue
    _, direction, metric = parameter.split('_', 2)
    column = metric_columns[metric]
    if column not in qc_df.columns:
      raise ValueError(f'{parameter} requires {column} in parent QC')

    # QC tables store counts/features on log10 and fractions on logit scales.
    if metric in ('counts', 'feats'):
      boundary = math.log10(value)
    elif value == 0:
      boundary = -math.inf
    elif value == 1:
      boundary = math.inf
    else:
      boundary = math.log(value / (1 - value))
    keep = pl.col(column) >= boundary if direction == 'min' else pl.col(column) <= boundary
    retained_qc_df = retained_qc_df.filter(keep)
  if set(cells_df['cell_id']) - set(qc_df['cell_id']):
    raise ValueError('Selected cells missing from parent QC')
  return cells_df.join(retained_qc_df.select('cell_id'), on = 'cell_id', how = 'semi')


def _write_batch_cell_tables(cells_df, batches, project_ids):
  import h5py

  selected_cells = set(cells_df['cell_id'])
  cell_rows = {cell: i for i, cell in enumerate(cells_df['cell_id'])}
  seen_cells = set()
  h5ad_lookup = {project_id: {} for project_id in project_ids}
  for batch in batches.values():
    h5ad_f = batch['path']

    # Every declared batch output exists, including empty/excluded batches.
    if batch.get('cells_f'):
      write_csv(cells_df.head(0).select('cell_id', 'sample_id'), batch['cells_f'])
    if Path(h5ad_f).stat().st_size == 0:
      continue
    with h5py.File(h5ad_f, 'r') as handle:
      if 'obs' not in handle or 'X' not in handle:
        continue
      obs = handle['obs']
      batch_cells = set(obs[obs.attrs.get('_index', '_index')].asstr()[:]) & selected_cells
    if batch_cells & seen_cells:
      duplicates = sorted(batch_cells & seen_cells)[:5]
      raise ValueError(f'Selected cells occur in multiple H5ADs: {duplicates}')
    project_id = batch['project_id']
    expected_cells = set(cells_df.filter(pl.col('project_id') == project_id)['cell_id'])
    if batch_cells - expected_cells:
      raise ValueError(f'Selected cell project_id disagrees with H5AD {h5ad_f}')
    seen_cells.update(batch_cells)

    # Keep selected cells in their parent-table order within each batch file.
    if batch_cells:
      batch_key = batch['key'].removeprefix(project_id + '_')
      h5ad_lookup[project_id][batch_key] = {'path': h5ad_f, 'project_id': project_id}
      if batch.get('cells_f'):
        rows = sorted(cell_rows[cell] for cell in batch_cells)
        batch_cells_df = cells_df[rows].select('cell_id', 'sample_id')
        write_csv(batch_cells_df, batch['cells_f'])
  if selected_cells != seen_cells:
    raise ValueError(f'{len(selected_cells - seen_cells)} selected cells are absent from the clean H5ADs')
  return h5ad_lookup


def _write_project_tables(cells_df, metadata_df, qc_df, h5ad_lookup, project_ids,
  integrated_fs, sample_meta_fs, qc_all_fs, h5ads_yaml_fs):
  for i, project_id in enumerate(project_ids):
    project_cells_df = cells_df.filter(pl.col('project_id') == project_id)
    prefix = project_id + '_'
    if not project_cells_df['sample_id'].str.starts_with(prefix).all():
      raise ValueError(f'Unexpected parent sample IDs for {project_id}; expected prefix {prefix}')

    # Source-project tables use local sample IDs; the selected-cell table above
    # retains project-qualified sample IDs for the joined analysis.
    project_cells_df = project_cells_df.with_columns(pl.col('sample_id').str.strip_prefix(prefix))
    project_metadata_df = metadata_df.filter(pl.col('project_id') == project_id)
    project_metadata_df = project_metadata_df.with_columns(pl.col('sample_id').str.strip_prefix(prefix))
    project_qc_df = qc_df.filter(pl.col('project_id') == project_id)
    write_csv(project_cells_df, integrated_fs[i])
    write_csv(project_metadata_df, sample_meta_fs[i])
    write_csv(project_qc_df, qc_all_fs[i])
    with open(h5ads_yaml_fs[i], 'w') as handle:
      yaml.safe_dump(h5ad_lookup[project_id], handle)


def batch_statistics(path, selected_f, chunk_size, output, trend_f = None):
  """Accumulate raw or clipped gene-count moments separately for each sample."""
  import h5py
  import numpy as np

  cells_df = pl.read_csv(selected_f, columns = ['cell_id', 'sample_id'])
  schema = {'sample_id': pl.String, 'gene_id': pl.String, 'n_cells': pl.Int64,
            'counts_sum': pl.Float64, 'squared_counts_sum': pl.Float64}
  if not cells_df.height:
    write_csv(pl.DataFrame(schema = schema), output)
    return
  trend_df = pl.read_csv(trend_f) if trend_f else None
  all_genes = np.array([], dtype = str)
  if Path(path).stat().st_size:
    with h5py.File(path, 'r') as handle:
      if 'var' in handle:
        var = handle['var']
        all_genes = var[var.attrs.get('_index', '_index')].asstr()[:]
  gene_lookup = {gene: i for i, gene in enumerate(all_genes)}
  sample_statistics = {}

  # A chunk can contain several biological samples; accumulate them separately.
  for genes, chunk_df, counts in _iter_selected_chunks(path, cells_df, chunk_size):
    gene_indices = np.array([gene_lookup[gene] for gene in genes])
    sample_cells = chunk_df.with_row_index('_chunk_row').partition_by('sample_id', as_dict = True)
    for (sample_id,), sample_cells_df in sample_cells.items():
      sample_counts = counts[sample_cells_df['_chunk_row'].to_numpy(), :].copy()
      if trend_df is not None:
        sample_trend_df = trend_df.filter(pl.col('sample_id') == sample_id)
        gene_order_df = pl.DataFrame({'gene_id': genes}).with_row_index('_gene_order')
        sample_trend_df = gene_order_df.join(sample_trend_df, on = 'gene_id', how = 'left')
        sample_trend_df = sample_trend_df.sort('_gene_order')
        if sample_trend_df['clip_val'].null_count():
          raise ValueError(f'HVG trend does not cover all genes for {sample_id}')
        clip_values = sample_trend_df['clip_val'].to_numpy()
        sample_counts.data = np.minimum(sample_counts.data, clip_values[sample_counts.indices])

      # Cell counts are per gene: CSC reads contribute one gene block at a time.
      if sample_id not in sample_statistics:
        sample_statistics[sample_id] = {
          'sample_id': sample_id, 'gene_id': all_genes,
          'n_cells': np.zeros(len(all_genes), dtype = np.int64),
          'counts_sum': np.zeros(len(all_genes)),
          'squared_counts_sum': np.zeros(len(all_genes))}
      statistics = sample_statistics[sample_id]
      statistics['n_cells'][gene_indices] += sample_counts.shape[0]
      statistics['counts_sum'][gene_indices] += np.asarray(sample_counts.sum(axis = 0)).ravel()
      statistics['squared_counts_sum'][gene_indices] += np.asarray(sample_counts.multiply(sample_counts).sum(axis = 0)).ravel()

  sample_statistics_dfs = [pl.DataFrame(values, schema = schema) for values in sample_statistics.values()]
  statistics_df = pl.concat(sample_statistics_dfs) if sample_statistics_dfs else pl.DataFrame(schema = schema)
  write_csv(statistics_df, output)


def _iter_selected_chunks(path, cells_df, chunk_size):
  """Read the stored sparse axis without AnnData's full-memory CSC fallback."""
  import h5py
  import numpy as np
  from scipy import sparse

  import anndata as ad  # only expression-reading stages need AnnData

  if Path(path).stat().st_size == 0:
    return
  with h5py.File(path, 'r') as handle:
    if 'X' not in handle or 'obs' not in handle:
      return
  adata = ad.read_h5ad(path, backed = 'r')
  try:
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
      raise ValueError(f'Duplicate cell/gene identifiers in {path}')
    row_lookup_df = pl.DataFrame({
      'cell_id': adata.obs_names.astype(str), '_row': np.arange(adata.n_obs)})
    selected_df = cells_df.join(row_lookup_df, on = 'cell_id')
    genes = np.asarray(adata.var_names, dtype = str)
    if not selected_df.height:
      return

    # CSC row indexing loads the complete matrix. Read bounded gene blocks
    # along its efficient column axis, then subset the selected cell rows.
    if getattr(adata.X, 'format', None) == 'csc':
      width = max(1, min(MAX_CSC_GENE_BLOCK,
        chunk_size * CSC_VALUES_PER_CHUNK_CELL // max(1, adata.n_obs)))
      selected_df = selected_df.sort('_row')
      for left in range(0, adata.n_vars, width):
        right = min(left + width, adata.n_vars)
        block = adata.X[:, left:right]
        counts = sparse.csr_matrix(block[selected_df['_row'].to_numpy(), :], dtype = np.float64)
        _check_counts(counts, path)
        yield genes[left:right], selected_df, counts
      return

    # CSR/dense storage can be read directly in bounded selected-cell chunks.
    for offset in range(0, selected_df.height, chunk_size):
      chunk_df = selected_df.slice(offset, chunk_size).sort('_row')
      counts = sparse.csr_matrix(adata.X[chunk_df['_row'].to_numpy(), :], dtype = np.float64)
      _check_counts(counts, path)
      yield genes, chunk_df, counts
  finally:
    adata.file.close()


def _check_counts(counts, path):
  import numpy as np

  if not np.isfinite(counts.data).all() or (counts.data < 0).any():
    raise ValueError(f'Invalid counts in {path}')


def fit_trend(files, output, batch_outputs = None):
  from hvgs import calculate_estimated_vars  # load the VST dependencies only here

  # One row per sample × gene, combining the same sample across input files.
  statistics_df = merge_moments(files).with_columns(
    (pl.col('counts_sum') / pl.col('n_cells')).alias('mean'),
    ((pl.col('squared_counts_sum') - pl.col('counts_sum')**2 / pl.col('n_cells'))
      / (pl.col('n_cells') - 1)).clip(lower_bound = 0).alias('variance'))
  if (statistics_df['n_cells'] < MIN_VARIANCE_CELLS).any():
    raise ValueError('HVG calculation needs at least two cells per sample')
  statistics_df = statistics_df.with_columns(pl.col('gene_id').alias('ensembl_id'))
  trend_df = calculate_estimated_vars(None, 'sample', 'sample_id', mean_var_df = statistics_df)
  write_csv(trend_df, output)

  # Each batch receives trends only for the samples it contains.
  if batch_outputs is not None:
    if len(files) != len(batch_outputs):
      raise ValueError('Expected one batch trend output per statistics input')
    sample_trends = trend_df.partition_by('sample_id', as_dict = True)
    for source_f, target_f in zip(files, batch_outputs):
      sample_ids = pl.read_csv(source_f, columns = ['sample_id'])['sample_id'].unique().to_list()
      batch_trend_df = (
        pl.concat([sample_trends[(sample_id,)] for sample_id in sample_ids])
        if sample_ids else trend_df.head(0))
      write_csv(batch_trend_df, target_f)


def merge_moments(files):
  # Combine overlapping samples only. Regrouping the entire atlas after each
  # input file would be quadratic in the number of batches.
  sample_statistics = {}
  for path in files:
    statistics_df = pl.read_csv(path)
    if not statistics_df.height:
      continue
    for sample_key, sample_df in statistics_df.partition_by('sample_id', as_dict = True).items():
      if sample_key in sample_statistics:
        sample_df = pl.concat([sample_statistics[sample_key], sample_df])
        sample_df = sample_df.group_by('sample_id', 'gene_id').agg(
          pl.col('n_cells').sum(), pl.col('counts_sum').sum(), pl.col('squared_counts_sum').sum())
      sample_statistics[sample_key] = sample_df
  if not sample_statistics:
    raise ValueError('No selected expression counts for HVG calculation')
  return pl.concat(list(sample_statistics.values()))


def rank_hvgs(files, trend_f, output, n_hvgs, exclude_f = None):
  from hvgs import _calculate_standardized_variance, _process_multiple_groups

  # Replace raw moments in the fitted trend with the clipped-count moments.
  moments_df = merge_moments(files)
  trend_df = pl.read_csv(trend_f).drop('counts_sum', 'squared_counts_sum')
  statistics_df = trend_df.join(
    moments_df.drop('n_cells'), on = ['sample_id', 'gene_id'], validate = '1:1')
  statistics_df = _calculate_standardized_variance(statistics_df)

  # Shared PCA uses only genes measured in every retained biological sample.
  n_samples = statistics_df['sample_id'].n_unique()
  common_genes_df = statistics_df.group_by('gene_id').agg(pl.len().alias('_n'))
  common_genes_df = common_genes_df.filter(pl.col('_n') == n_samples)
  statistics_df = statistics_df.join(common_genes_df.select('gene_id'), on = 'gene_id', how = 'semi')
  excluded_genes = pl.read_csv(exclude_f).to_series(0).to_list() if exclude_f else []
  ranked_df = _process_multiple_groups(statistics_df, 'sample_id', excluded_genes, n_hvgs)
  ranked_df = ranked_df.filter(
    pl.col('highly_variable_rank').is_not_null() & pl.col('variances_norm').is_finite())
  ranked_df = ranked_df.sort(
    ['highly_variable_nbatches', 'highly_variable_rank', 'gene_id'],
    descending = [True, False, False])
  if ranked_df.height < n_hvgs:
    raise ValueError(f'Only {ranked_df.height} eligible shared genes; requested {n_hvgs} HVGs')

  ranked_df = ranked_df.with_row_index('_rank').with_columns(
    (pl.col('_rank') < n_hvgs).alias('highly_variable')).drop('_rank')
  write_csv(ranked_df, output)


def unique_cells(df, description):
  if 'cell_id' not in df.columns or df['cell_id'].null_count():
    raise ValueError(f'{description} needs non-null cell_id values')
  if df['cell_id'].n_unique() != df.height:
    raise ValueError(f'{description} has duplicate cell_id values')


def write_csv(df, path):
  Path(path).parent.mkdir(parents = True, exist_ok = True)
  if str(path).endswith('.gz'):
    with gzip.open(path, 'wb') as handle:
      df.write_csv(handle)
  else:
    df.write_csv(path)


def main():
  parser = argparse.ArgumentParser(description = __doc__)
  subparsers = parser.add_subparsers(dest = 'command', required = True)

  # Keep the existing preparation CLI for callers outside the Snakemake rule.
  prepare_parser = subparsers.add_parser('prepare')
  for name in ('context', 'config'):
    prepare_parser.add_argument('--' + name, required = True)

  for command in ('stats', 'clip'):
    statistics_parser = subparsers.add_parser(command)
    for name in ('h5ad', 'cells', 'output'):
      statistics_parser.add_argument('--' + name, required = True)
    statistics_parser.add_argument('--trend')
    statistics_parser.add_argument('--chunk-size', type = int, required = True)

  trend_parser = subparsers.add_parser('trend')
  trend_parser.add_argument('--output', required = True)
  for name, required in (('inputs', True), ('batch-outputs', False)):
    trend_parser.add_argument('--' + name, nargs = '+', required = required)

  rank_parser = subparsers.add_parser('rank')
  for name in ('trend', 'output'):
    rank_parser.add_argument('--' + name, required = True)
  rank_parser.add_argument('--exclude')
  rank_parser.add_argument('--n-hvgs', type = int, required = True)
  rank_parser.add_argument('--inputs', nargs = '+', required = True)
  args = parser.parse_args()

  # Dispatch stages with their explicit file paths and validated settings.
  if args.command == 'prepare':
    prepare(json.loads(args.context), json.loads(args.config))
  elif args.command in ('stats', 'clip'):
    batch_statistics(path = args.h5ad, selected_f = args.cells,
      chunk_size = args.chunk_size, output = args.output, trend_f = args.trend)
  elif args.command == 'trend':
    fit_trend(files = args.inputs, output = args.output, batch_outputs = args.batch_outputs)
  else:
    rank_hvgs(files = args.inputs, trend_f = args.trend, output = args.output,
      n_hvgs = args.n_hvgs, exclude_f = args.exclude)


if 'snakemake' in globals():
  # script: carries large batch lookups through the job object, avoiding ARG_MAX.
  with open(snakemake.log[0], 'a') as log:
    with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
      prepare_sources(
        integration_f = snakemake.input.integration_f,
        coldata_f     = snakemake.input.coldata_f,
        metadata_f    = snakemake.input.metadata_f,
        labels_f      = snakemake.input.labels_f,
        qc_f          = snakemake.input.qc_f,
        selection     = snakemake.params.selection,
        qc            = snakemake.params.qc,
        batches       = snakemake.params.batches,
        project_ids   = snakemake.params.project_ids,
        selected_f    = snakemake.output.cells_f,
        summary_f     = snakemake.output.summary_f,
        integrated_fs = list(snakemake.output.integration_fs),
        sample_meta_fs = list(snakemake.output.sample_meta_fs),
        qc_all_fs     = list(snakemake.output.qc_fs),
        h5ads_yaml_fs = list(snakemake.output.h5ads_yaml_fs),
        int_n_dims    = snakemake.params.int_n_dims,
        int_batch_var = snakemake.params.int_batch_var)
elif __name__ == '__main__':
  main()
