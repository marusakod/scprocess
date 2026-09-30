"""Joined-zoom selection and sample-aware HVGs from completed clean H5ADs.

No upstream workflow is invoked. Expression is read in bounded cell chunks;
sample moments are combined across files before fitting the existing VST trend.
"""
import argparse
import gzip
import json
import math
from pathlib import Path

import polars as pl
import yaml


def write_csv(df, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    if str(path).endswith('.gz'):
        with gzip.open(path, 'wb') as fh:
            df.write_csv(fh)
    else:
        df.write_csv(path)


def unique_cells(df, description):
    if 'cell_id' not in df.columns or df['cell_id'].null_count():
        raise ValueError(f'{description} needs non-null cell_id values')
    if df['cell_id'].n_unique() != df.height:
        raise ValueError(f'{description} has duplicate cell_id values')


def select_cells(parent, coldata, metadata, labels, qc, selection, thresholds, n_dims):
    """Select strictly within the parent, preserving its IDs and descriptive metadata."""
    for frame, desc in ((parent, 'Parent integration'), (coldata, 'Parent coldata'),
                        (labels, 'Selection labels'), (qc, 'Parent QC')):
        unique_cells(frame, desc)
    for key in ('sample_id', 'project_id'):
        if key not in parent.columns or parent[key].null_count():
            raise ValueError(f'Parent integration needs non-null {key}')
    # Coldata carries metadata absent from integrated_dt. Parent integration wins
    # for coordinates and clustering, not stale source-project embeddings.
    extras = [c for c in coldata.columns if c not in parent.columns]
    cells = parent.join(coldata.select(['cell_id', *extras]), on='cell_id', how='left', validate='1:1')
    if metadata['sample_id'].n_unique() != metadata.height:
        raise ValueError('Parent sample metadata has duplicate sample_id values')
    extras = [c for c in metadata.columns if c not in cells.columns]
    cells = cells.join(metadata.select(['sample_id', *extras]), on='sample_id', how='left', validate='m:1')
    col = selection['labels_col']
    if col not in labels.columns:
        raise ValueError(f'Labels column not found: {col}')
    labels = labels.select('cell_id', pl.col(col).cast(pl.String).alias('parent_label'))
    labels = labels.join(parent.select('cell_id'), on='cell_id', how='semi')
    missing = set(selection['sel_labels']) - set(labels['parent_label'].drop_nulls())
    if missing:
        raise ValueError(f'Requested labels absent from the parent join: {sorted(missing)}')
    cells = cells.join(labels, on='cell_id', how='inner', validate='1:1')
    cells = cells.filter(pl.col('parent_label').is_in(selection['sel_labels']))
    selected_counts = cells.group_by('project_id', 'sample_id').len(name='n_selected')
    excluded = selection.get('exclude', {}).get('sample_id', [])
    missing = set(excluded) - set(parent['sample_id'])
    if missing:
        raise ValueError(f'Excluded sample IDs absent from parent join: {sorted(missing)}')
    cells = cells.filter(~pl.col('sample_id').is_in(excluded))
    metrics = {'counts': 'log_counts', 'feats': 'log_feats',
               'mito': 'logit_mito', 'splice': 'logit_spliced'}
    keep_qc = qc
    for key, value in thresholds.items():
        if value is None or key == 'qc_min_cells':
            continue
        _, direction, metric = key.split('_', 2)
        column = metrics[metric]
        if column not in qc.columns:
            raise ValueError(f'{key} requires {column} in parent QC')
        if metric in ('counts', 'feats'):
            bound = math.log10(value)
        else:
            bound = -math.inf if value == 0 else math.inf if value == 1 else math.log(value / (1-value))
        expr = pl.col(column) >= bound if direction == 'min' else pl.col(column) <= bound
        keep_qc = keep_qc.filter(expr)
    if set(cells['cell_id']) - set(qc['cell_id']):
        raise ValueError('Selected cells missing from parent QC')
    cells = cells.join(keep_qc.select('cell_id'), on='cell_id', how='semi')
    qc_counts = cells.group_by('project_id', 'sample_id').len(name='n_after_qc')
    minimum = max(2, thresholds['qc_min_cells'])  # sample variance requires >=2
    cells = cells.filter(pl.len().over('sample_id') >= minimum)
    final_counts = cells.group_by('project_id', 'sample_id').len(name='n_retained')
    summary = (selected_counts.join(qc_counts, on=['project_id', 'sample_id'], how='left')
               .join(final_counts, on=['project_id', 'sample_id'], how='left').fill_null(0)
               .sort('project_id', 'sample_id'))
    if cells.height <= n_dims:
        raise ValueError(f'Only {cells.height} cells remain; need more than {n_dims} PCs')
    # Keep parent labels/coordinates for comparison without confusing them with
    # the newly generated clusters/embeddings.
    rename = {c: 'parent_' + c for c in cells.columns
              if c.startswith(('RNA_snn_res.', 'pca_', 'hmny_pca_', 'harmony_', 'umap_', 'UMAP_'))
              or c in ('UMAP1', 'UMAP2')}
    if set(rename.values()) & set(cells.columns):
        raise ValueError('Parent embedding column names conflict; recursive joined zooms are not supported')
    return cells.rename(rename), summary


def prepare(context, config):
    import h5py  # only needed for H5AD header inspection
    paths = context['paths']
    cells, summary = select_cells(
        *(pl.read_csv(paths[key]) for key in ('integration', 'coldata', 'metadata', 'labels', 'qc')),
        context['selection'], context['qc'], config['integration']['int_n_dims'])
    batch_vars = config['integration']['int_batch_var']
    for var in ([batch_vars] if isinstance(batch_vars, str) else batch_vars):
        if var not in cells.columns or cells[var].null_count() or (cells[var].cast(pl.String) == '').any():
            raise ValueError(f'Integration batch variable {var!r} is missing for selected cells')
    selected = set(cells['cell_id'])
    cell_rows = {cell: i for i, cell in enumerate(cells['cell_id'])}
    seen = set()
    manifests = {pid: {} for pid in context['project']['project_ids']}
    for batch in context['batches'].values():
        path = batch['path']
        if batch.get('cells_f'):
            write_csv(cells.head(0).select('cell_id', 'sample_id'), batch['cells_f'])
        if Path(path).stat().st_size == 0:
            continue
        with h5py.File(path) as fh:
            if 'obs' not in fh or 'X' not in fh:
                continue
            obs = fh['obs']
            ids = set(obs[obs.attrs.get('_index', '_index')].asstr()[:]) & selected
        if ids & seen:
            raise ValueError(f'Selected cells occur in multiple H5ADs: {sorted(ids & seen)[:5]}')
        pid = batch['project_id']
        expected = set(cells.filter(pl.col('project_id') == pid)['cell_id'])
        if ids - expected:
            raise ValueError(f'Selected cell project_id disagrees with H5AD {path}')
        seen.update(ids)
        if ids:
            key = batch['key'].removeprefix(pid + '_')
            manifests[pid][key] = {'path': path, 'project_id': pid}
            if batch.get('cells_f'):
                write_csv(cells[sorted(cell_rows[cell] for cell in ids)].select('cell_id', 'sample_id'), batch['cells_f'])
    if selected != seen:
        raise ValueError(f'{len(selected - seen)} selected cells are absent from the clean H5ADs')
    write_csv(cells, context['selected_f'])
    write_csv(summary, context['summary_f'])
    qc = pl.read_csv(paths['qc']).join(cells.select('cell_id'), on='cell_id', how='semi')
    metadata = pl.read_csv(paths['metadata']).join(cells.select('sample_id').unique(), on='sample_id', how='semi')
    proj = context['project']
    for i, pid in enumerate(proj['project_ids']):
        frame = cells.filter(pl.col('project_id') == pid)
        prefix = pid + '_'
        if not frame['sample_id'].str.starts_with(prefix).all():
            raise ValueError(f'Unexpected parent sample IDs for {pid}; expected prefix {prefix}')
        frame = frame.with_columns(pl.col('sample_id').str.strip_prefix(prefix))
        write_csv(frame, proj['integrated_fs'][i])
        meta = metadata.filter(pl.col('project_id') == pid).with_columns(pl.col('sample_id').str.strip_prefix(prefix))
        write_csv(meta, proj['sample_meta_fs'][i])
        write_csv(qc.filter(pl.col('project_id') == pid), proj['qc_all_fs'][i])
        with open(proj['h5ads_yaml_fs'][i], 'w') as fh:
            yaml.safe_dump(manifests[pid], fh)
    print(f"Retained {cells.height} cells across {cells['sample_id'].n_unique()} samples", flush=True)


def iter_selected_chunks(path, cells, chunk_size):
    """Read along the stored sparse axis, avoiding AnnData's full-memory fallback."""
    import anndata as ad
    import h5py
    import numpy as np
    from scipy import sparse
    if Path(path).stat().st_size == 0:
        return
    with h5py.File(path) as fh:
        if 'X' not in fh or 'obs' not in fh:
            return  # prepare() already checks complete selected-cell coverage
    adata = ad.read_h5ad(path, backed='r')
    try:
        if not adata.obs_names.is_unique or not adata.var_names.is_unique:
            raise ValueError(f'Duplicate cell/gene identifiers in {path}')
        keep = cells.join(pl.DataFrame({'cell_id': adata.obs_names.astype(str),
                                       '_row': np.arange(adata.n_obs)}), on='cell_id')
        genes = np.asarray(adata.var_names, dtype=str)
        if not keep.height:
            return
        if getattr(adata.X, 'format', None) == 'csc':
            # CSC row indexing calls AnnData.to_memory(). Instead read bounded
            # gene blocks along its efficient column axis, then select cells.
            width = max(1, min(256, chunk_size * 1024 // max(1, adata.n_obs)))
            keep = keep.sort('_row')
            for left in range(0, adata.n_vars, width):
                right = min(left + width, adata.n_vars)
                block = adata.X[:, left:right]
                matrix = sparse.csr_matrix(block[keep['_row'].to_numpy(), :], dtype=np.float64)
                if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
                    raise ValueError(f'Invalid counts in {path}')
                yield genes[left:right], keep, matrix
            return
        for offset in range(0, keep.height, chunk_size):
            chunk = keep.slice(offset, chunk_size).sort('_row')
            matrix = sparse.csr_matrix(adata.X[chunk['_row'].to_numpy(), :], dtype=np.float64)
            if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
                raise ValueError(f'Invalid counts in {path}')
            yield genes, chunk, matrix
    finally:
        adata.file.close()


def batch_statistics(path, selected_f, chunk_size, output, trend_f=None):
    import numpy as np
    import h5py
    cells = pl.read_csv(selected_f, columns=['cell_id', 'sample_id'])
    if not cells.height:
        write_csv(pl.DataFrame(schema={'sample_id': pl.String, 'gene_id': pl.String,
            'n_cells': pl.Int64, 'counts_sum': pl.Float64, 'squared_counts_sum': pl.Float64}), output)
        return
    trend = pl.read_csv(trend_f) if trend_f else None
    all_genes = np.array([], dtype=str)
    if Path(path).stat().st_size:
        with h5py.File(path) as fh:
            if 'var' in fh:
                var = fh['var']
                all_genes = var[var.attrs.get('_index', '_index')].asstr()[:]
    gene_lookup = {gene: i for i, gene in enumerate(all_genes)}
    accum = {}
    # Only a batch's samples and selected cell chunks are held in memory.
    for genes, chunk, matrix in iter_selected_chunks(path, cells, chunk_size):
        gene_indices = np.array([gene_lookup[gene] for gene in genes])
        for (sample,), group in chunk.with_row_index('_chunk_row').partition_by('sample_id', as_dict=True).items():
            x = matrix[group['_chunk_row'].to_numpy(), :].copy()
            if trend is not None:
                t = trend.filter(pl.col('sample_id') == sample)
                order = pl.DataFrame({'gene_id': genes}).with_row_index('_gene_order')
                t = order.join(t, on='gene_id', how='left').sort('_gene_order')
                if t['clip_val'].null_count():
                    raise ValueError(f'HVG trend does not cover all genes for {sample}')
                x.data = np.minimum(x.data, t['clip_val'].to_numpy()[x.indices])
            sums = np.asarray(x.sum(axis=0)).ravel()
            squares = np.asarray(x.multiply(x).sum(axis=0)).ravel()
            if sample not in accum:
                accum[sample] = [all_genes, np.zeros(len(all_genes), dtype=np.int64),
                                 np.zeros(len(all_genes)), np.zeros(len(all_genes))]
            record = accum[sample]
            record[1][gene_indices] += x.shape[0]
            record[2][gene_indices] += sums
            record[3][gene_indices] += squares
    schema = {'sample_id': pl.String, 'gene_id': pl.String, 'n_cells': pl.Int64,
              'counts_sum': pl.Float64, 'squared_counts_sum': pl.Float64}
    frames = [pl.DataFrame({'sample_id': sample, 'gene_id': row[0], 'n_cells': row[1],
                            'counts_sum': row[2], 'squared_counts_sum': row[3]}, schema=schema)
              for sample, row in accum.items()]
    write_csv(pl.concat(frames) if frames else pl.DataFrame(schema=schema), output)


def merge_moments(files):
    # Combine only overlapping samples. Regrouping the entire atlas after each
    # file would be quadratic in the number of batches.
    samples = {}
    for path in files:
        frame = pl.read_csv(path)
        if not frame.height:
            continue
        for key, part in frame.partition_by('sample_id', as_dict=True).items():
            if key in samples:
                part = pl.concat([samples[key], part]).group_by('sample_id', 'gene_id').agg(
                    pl.col('n_cells').sum(), pl.col('counts_sum').sum(), pl.col('squared_counts_sum').sum())
            samples[key] = part
    if not samples:
        raise ValueError('No selected expression counts for HVG calculation')
    return pl.concat(list(samples.values()))


def fit_trend(files, output, batch_outputs=None):
    from hvgs import calculate_estimated_vars  # existing Seurat-style VST fit
    stats = merge_moments(files).with_columns(
        (pl.col('counts_sum') / pl.col('n_cells')).alias('mean'),
        ((pl.col('squared_counts_sum') - pl.col('counts_sum')**2 / pl.col('n_cells'))
         / (pl.col('n_cells') - 1)).clip(lower_bound=0).alias('variance'))
    if (stats['n_cells'] < 2).any():
        raise ValueError('HVG calculation needs at least two cells per sample')
    stats = stats.with_columns(pl.col('gene_id').alias('ensembl_id'))
    trend = calculate_estimated_vars(None, 'sample', 'sample_id', mean_var_df=stats)
    write_csv(trend, output)
    if batch_outputs is not None:
        if len(files) != len(batch_outputs):
            raise ValueError('Expected one batch trend output per statistics input')
        parts = trend.partition_by('sample_id', as_dict=True)
        for source, target in zip(files, batch_outputs):
            names = pl.read_csv(source, columns=['sample_id'])['sample_id'].unique().to_list()
            write_csv(pl.concat([parts[(name,)] for name in names]) if names else trend.head(0), target)


def rank_hvgs(files, trend_f, output, n_hvgs, exclude_f=None):
    from hvgs import _calculate_standardized_variance, _process_multiple_groups
    moments = merge_moments(files)
    trend = pl.read_csv(trend_f).drop('counts_sum', 'squared_counts_sum')
    stats = trend.join(moments.drop('n_cells'), on=['sample_id', 'gene_id'], validate='1:1')
    stats = _calculate_standardized_variance(stats)
    n_samples = stats['sample_id'].n_unique()
    # Only genes measured in every retained sample can enter a shared PCA.
    common = stats.group_by('gene_id').agg(pl.len().alias('_n')).filter(pl.col('_n') == n_samples)
    stats = stats.join(common.select('gene_id'), on='gene_id', how='semi')
    excluded = pl.read_csv(exclude_f).to_series(0).to_list() if exclude_f else []
    ranked = _process_multiple_groups(stats, 'sample_id', excluded, n_hvgs)
    ranked = ranked.filter(pl.col('highly_variable_rank').is_not_null() & pl.col('variances_norm').is_finite())
    ranked = ranked.sort(['highly_variable_nbatches', 'highly_variable_rank', 'gene_id'], descending=[True, False, False])
    if ranked.height < n_hvgs:
        raise ValueError(f'Only {ranked.height} eligible shared genes; requested {n_hvgs} HVGs')
    write_csv(ranked.with_row_index('_rank').with_columns(
        (pl.col('_rank') < n_hvgs).alias('highly_variable')).drop('_rank'), output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--context', required=True)
    p.add_argument('--config', required=True)
    for name in ('stats', 'clip'):
        p = sub.add_parser(name)
        p.add_argument('--h5ad', required=True)
        p.add_argument('--cells', required=True)
        p.add_argument('--chunk-size', type=int, required=True)
        p.add_argument('--output', required=True)
        p.add_argument('--trend')
    p = sub.add_parser('trend')
    p.add_argument('--inputs', nargs='+', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--batch-outputs', nargs='+')
    p = sub.add_parser('rank')
    p.add_argument('--inputs', nargs='+', required=True)
    p.add_argument('--trend', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--n-hvgs', type=int, required=True)
    p.add_argument('--exclude')
    args = parser.parse_args()
    if args.command == 'prepare':
        prepare(json.loads(args.context), json.loads(args.config))
    elif args.command in ('stats', 'clip'):
        batch_statistics(args.h5ad, args.cells, args.chunk_size, args.output, args.trend)
    elif args.command == 'trend':
        fit_trend(args.inputs, args.output, args.batch_outputs)
    else:
        rank_hvgs(args.inputs, args.trend, args.output, args.n_hvgs, args.exclude)


if 'snakemake' in globals():
    import contextlib
    with open(snakemake.log[0], 'a') as log:
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            prepare(json.loads(snakemake.params.context), json.loads(snakemake.params.effective))
elif __name__ == '__main__':
    main()
