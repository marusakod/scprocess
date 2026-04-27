import os
import sys
import gzip
import argparse
import pathlib
import yaml
import polars as pl


# ---------------------------------------------------------------------------
# Step 1: Select joint HVGs by mean rank across projects
# ---------------------------------------------------------------------------

def select_joint_hvgs(var_stats_fs, project_ids, n_hvgs, out_f):
  """
  Combine per-project standardized variance stats into a joint HVG list.

  Within each project genes are ranked by variances_norm (descending).
  Genes absent from a project receive rank = n_genes_in_project + 1.
  The mean rank across projects is computed; the top n_hvgs genes are kept.

  Parameters
  ----------
  var_stats_fs : list of str
      Paths to standardized_variance_stats_*.csv.gz, one per project.
  project_ids  : list of str
      Project IDs (parallel to var_stats_fs).
  n_hvgs       : int
      Number of HVGs to select.
  out_f        : str
      Output CSV.gz path (columns: gene_id, mean_rank, n_projects).
  """
  print(f"selecting joint HVGs (n={n_hvgs}) across {len(project_ids)} projects")

  proj_hvg_dfs = {}
  for pid, f in zip(project_ids, var_stats_fs):
    df = pl.read_csv(f).select(['gene_id', 'variances_norm']).unique('gene_id')
    proj_hvg_dfs[pid] = df

  # union of all gene ids
  all_genes = pl.Series(
    list({g for df in proj_hvg_dfs.values() for g in df['gene_id'].to_list()})
  ).alias('gene_id')
  base_df = pl.DataFrame({'gene_id': all_genes})

  # compute per-project ranks; missing genes get rank = n_genes + 1
  for pid, df in proj_hvg_dfs.items():
    n_genes = df.shape[0]
    ranked  = df.sort('variances_norm', descending=True).with_row_index('rank').select(
      ['gene_id', pl.col('rank').cast(pl.Float64) + 1]
    )
    base_df = base_df.join(ranked.rename({'rank': pid}), on='gene_id', how='left')
    base_df = base_df.with_columns(
      pl.col(pid).fill_null(n_genes + 1)
    )

  # mean rank and number of projects the gene was observed in
  rank_cols = project_ids
  mean_rank_expr    = pl.mean_horizontal([pl.col(c) for c in rank_cols])
  # a gene is "present" in project c if its rank was not the penalty rank (n_genes + 1)
  n_projects_expr   = sum(
    (pl.col(c) < len(proj_hvg_dfs[c]) + 1).cast(pl.Int32)
    for c in rank_cols
  )
  base_df = base_df.with_columns(
    mean_rank_expr.alias('mean_rank'),
    n_projects_expr.alias('n_projects')
  )

  # select top n_hvgs
  top_hvgs = (
    base_df
    .sort('mean_rank')
    .head(n_hvgs)
    .select(['gene_id', 'mean_rank', 'n_projects'])
  )

  if top_hvgs.shape[0] < n_hvgs:
    raise ValueError(
      f"Only {top_hvgs.shape[0]} genes available across projects, "
      f"but {n_hvgs} HVGs requested."
    )

  print(f"  selected {top_hvgs.shape[0]} joint HVGs")
  with gzip.open(out_f, 'wb') as fh:
    top_hvgs.write_csv(fh)
  print(f"  saved to {out_f}")


# ---------------------------------------------------------------------------
# Step 2: Build joint count matrix, coldata, and sample metadata
# ---------------------------------------------------------------------------

def _ok_cells_filter(int_dt):
  """Return a boolean filter expression for non-doublet cells.
  If is_dbl/in_dbl_cl columns are absent (e.g. zoom integrated_dt where
  doublets were already removed upstream), all cells are considered clean.
  """
  if 'is_dbl' in int_dt.columns and 'in_dbl_cl' in int_dt.columns:
    return (pl.col('is_dbl') == False) & (pl.col('in_dbl_cl') == False)
  return pl.lit(True)


def _check_sample_id_uniqueness(project_ids, integrated_dt_fs):
  """Raise ValueError if any sample_id (among non-doublet cells) appears in >1 project."""
  seen = {}
  for pid, int_f in zip(project_ids, integrated_dt_fs):
    int_dt  = pl.read_csv(int_f)
    samples = int_dt.filter(_ok_cells_filter(int_dt))['sample_id'].unique().to_list()
    for s in samples:
      if s in seen and seen[s] != pid:
        raise ValueError(
          f"sample_id '{s}' appears in both project '{seen[s]}' "
          f"and '{pid}'. Use unique sample IDs across projects, or the joint "
          f"sample_id ('{pid}_{s}') will still be unique."
        )
      seen[s] = pid


def _load_batch_hvg_matrix(h5ad_path, ok_cells, hvg_list, pid, batch_key):
  """
  Load one h5ad batch, filter to non-doublet cells, and subset to HVG genes.

  Returns (csc_sub, kept_bcs) where csc_sub is (n_hvgs x n_kept_cells) CSC,
  or None if no clean cells are present in this batch.
  """
  import numpy as np
  import h5py
  from scipy.sparse import csc_matrix, csr_matrix

  print(f"    loading {batch_key}")

  with h5py.File(h5ad_path, 'r') as f:
    obs_idx_col = f['obs'].attrs.get('_index', 'cell_id')
    barcodes    = f['obs'][obs_idx_col][:].astype(str)
    var_idx_col = f['var'].attrs.get('_index', 'gene_ids')
    features    = f['var'][var_idx_col][:].astype(str)
    data     = f['X/data'][:]
    indices  = f['X/indices'][:]
    indptr   = f['X/indptr'][:]
    n_cells  = len(barcodes)
    n_genes  = len(features)
    enc_type = f['X'].attrs.get('encoding-type', 'csr_matrix')

  # AnnData X is (n_cells, n_genes); convert to (n_genes, n_cells) CSC
  if enc_type == 'csc_matrix':
    csc_mat = csc_matrix((data, indices, indptr), shape=(n_genes, n_cells))
  else:
    csr_mat = csr_matrix((data, indices, indptr), shape=(n_cells, n_genes))
    csc_mat = csr_mat.T.tocsc()

  # filter columns to non-doublet cells
  keep_mask = np.array([bc in ok_cells for bc in barcodes])
  if keep_mask.sum() == 0:
    print(f"    no clean cells in {batch_key}, skipping")
    return None

  csc_mat  = csc_mat[:, keep_mask]
  kept_bcs = barcodes[keep_mask]

  # subset rows to HVG genes
  feat_idx    = {g: i for i, g in enumerate(features)}
  hvg_row_idx = [feat_idx[g] for g in hvg_list if g in feat_idx]
  if len(hvg_row_idx) == 0:
    raise ValueError(f"No joint HVG genes found in h5ad for {batch_key} in {pid}")
  csc_sub = csc_mat[hvg_row_idx, :]

  return csc_sub, kept_bcs.tolist()


def _build_project_coldata(int_dt, pid, smeta_dt, metadata_vars):
  """
  Build per-project coldata: filter to non-doublets, prefix sample IDs with the
  project ID, and optionally join extra metadata variables from sample_meta.
  Cell IDs are kept as-is (matching h5ad colnames) since sample IDs are
  guaranteed unique across projects by _check_sample_id_uniqueness.
  """
  proj_cells = int_dt.filter(_ok_cells_filter(int_dt)).select(
    ['cell_id', 'sample_id'] +
    [c for c in int_dt.columns if c not in ['cell_id', 'sample_id', 'project_id']])
  proj_cells = proj_cells.with_columns([
    pl.col('sample_id').map_elements(lambda x: f"{pid}_{x}", return_dtype=pl.Utf8),
    pl.lit(pid).alias('project_id')
  ])
  if metadata_vars:
    meta_cols  = ['sample_id'] + [v for v in metadata_vars if v in smeta_dt.columns]
    orig_smeta = smeta_dt.select(meta_cols).unique()
    orig_smeta = orig_smeta.with_columns(
      pl.col('sample_id').map_elements(lambda x: f"{pid}_{x}", return_dtype=pl.Utf8)
    )
    proj_cells = proj_cells.join(orig_smeta, on='sample_id', how='left')
  return proj_cells


def _load_project_data(pid, h5ads_yaml_f, int_f, smeta_f, hvg_list, metadata_vars):
  """
  Load one project: iterate batches, filter to non-doublets, build the per-project
  count matrix (n_hvgs x n_cells), coldata, and sample metadata.

  Returns (proj_mat, proj_barcodes, coldata_df, smeta_df).
  """
  from scipy.sparse import hstack

  print(f"  processing project: {pid}")

  with open(h5ads_yaml_f) as fh:
    h5ad_paths = yaml.safe_load(fh)

  int_dt   = pl.read_csv(int_f)
  ok_cells = set(
    int_dt.filter(_ok_cells_filter(int_dt))['cell_id'].to_list()
  )
  print(f"    non-doublet cells: {len(ok_cells)}")

  smeta_dt      = pl.read_csv(smeta_f)
  proj_mats     = []
  proj_barcodes = []
  for batch_key, h5ad_entry in h5ad_paths.items():
    h5ad_path = h5ad_entry if isinstance(h5ad_entry, str) else h5ad_entry['path']
    result = _load_batch_hvg_matrix(h5ad_path, ok_cells, hvg_list, pid, batch_key)
    if result is None:
      continue
    csc_sub, kept_bcs = result
    proj_mats.append(csc_sub)
    proj_barcodes.extend(kept_bcs)

  if len(proj_mats) == 0:
    raise ValueError(f"No cells loaded for project {pid}")

  proj_mat = hstack(proj_mats, format='csc') if len(proj_mats) > 1 else proj_mats[0]

  coldata_df = _build_project_coldata(int_dt, pid, smeta_dt, metadata_vars)

  smeta_df = smeta_dt.with_columns([
    pl.col('sample_id').map_elements(lambda x: f"{pid}_{x}", return_dtype=pl.Utf8),
    pl.lit(pid).alias('project_id')
  ])
  if 'bad_sample_id' not in smeta_df.columns:
    smeta_df = smeta_df.with_columns(pl.lit(False).alias('bad_sample_id'))

  return proj_mat, proj_barcodes, coldata_df, smeta_df


def _save_joint_outputs(joint_csc, hvg_list, all_barcodes, all_coldata_dfs,
                        cell_totals, all_smeta_dfs, out_h5_f, out_coldata_f,
                        out_sample_meta_f):
  """Write joint count matrix (HDF5), coldata (CSV.gz), and sample metadata (CSV)."""
  import numpy as np
  import h5py

  print(f"  saving matrix to {out_h5_f}")
  pathlib.Path(out_h5_f).parent.mkdir(parents=True, exist_ok=True)
  with h5py.File(out_h5_f, 'w') as f:
    f.create_dataset('matrix/data',          data=joint_csc.data)
    f.create_dataset('matrix/indices',       data=joint_csc.indices)
    f.create_dataset('matrix/indptr',        data=joint_csc.indptr)
    f.create_dataset('matrix/shape',         data=joint_csc.shape)
    f.create_dataset('matrix/features/name', data=np.array(hvg_list, dtype='S'))
    f.create_dataset('matrix/barcodes',      data=np.array(all_barcodes, dtype='S'))

  print("  saving coldata")
  coldata_df = pl.concat(all_coldata_dfs, how='diagonal')
  fixed_cols = ['cell_id', 'sample_id', 'project_id']
  other_cols = [c for c in coldata_df.columns if c not in fixed_cols]
  coldata_df = coldata_df.select(fixed_cols + other_cols)
  # reorder rows to match matrix barcodes and attach per-cell library size
  bc_order_df = pl.DataFrame({
    'cell_id': all_barcodes,
    '_order':  range(len(all_barcodes)),
    'total':   cell_totals,
  })
  coldata_df = coldata_df.join(bc_order_df, on='cell_id').sort('_order').drop('_order')
  with gzip.open(out_coldata_f, 'wb') as fh:
    coldata_df.write_csv(fh)

  print("  saving sample metadata")
  smeta_df = _smart_concat(all_smeta_dfs)
  smeta_df.write_csv(out_sample_meta_f)


def build_joint_matrix(joint_hvgs_f, h5ads_yaml_fs, project_ids, integrated_dt_fs,
                       sample_meta_fs, metadata_vars_str, out_h5_f, out_coldata_f,
                       out_sample_meta_f):
  """
  Assemble a joint HVG count matrix from per-project h5ads.

  Parameters
  ----------
  joint_hvgs_f        : str   Path to joint HVGs CSV.gz (ensembl_id column).
  h5ads_yaml_fs       : list  Per-project h5ads_clean_paths_*.yaml files.
  project_ids         : list  Project IDs (parallel to h5ads_yaml_fs).
  integrated_dt_fs    : list  Per-project integrated_dt_*.csv.gz (for non-doublet filter).
  sample_meta_fs      : list  Per-project sample_metadata CSV files.
  metadata_vars_str   : str   Space-separated list of metadata variable names to carry through.
  out_h5_f            : str   Output joint HVG count matrix (H5, CSC format).
  out_coldata_f       : str   Output joint coldata CSV.gz.
  out_sample_meta_f   : str   Output joint sample metadata CSV.
  """
  import numpy as np
  from scipy.sparse import hstack

  print("building joint count matrix")

  metadata_vars = metadata_vars_str.split() if metadata_vars_str.strip() else []

  hvg_df   = pl.read_csv(joint_hvgs_f)
  hvg_list = hvg_df['gene_id'].to_list()
  print(f"  joint HVGs: {len(hvg_list)}")

  _check_sample_id_uniqueness(project_ids, integrated_dt_fs)

  all_mats        = []
  all_barcodes    = []
  all_coldata_dfs = []
  all_smeta_dfs   = []

  for pid, h5ads_yaml_f, int_f, smeta_f in zip(
      project_ids, h5ads_yaml_fs, integrated_dt_fs, sample_meta_fs):
    proj_mat, proj_barcodes, coldata_df, smeta_df = _load_project_data(
      pid, h5ads_yaml_f, int_f, smeta_f, hvg_list, metadata_vars
    )
    all_mats.append(proj_mat)
    all_barcodes.extend(proj_barcodes)
    all_coldata_dfs.append(coldata_df)
    all_smeta_dfs.append(smeta_df)

  print("  concatenating matrices")
  joint_mat = hstack(all_mats, format='csc') if len(all_mats) > 1 else all_mats[0]
  joint_csc = joint_mat.tocsc()
  print(f"  joint matrix shape: {joint_csc.shape} (genes x cells)")

  if joint_csc.shape[1] != len(all_barcodes):
    raise ValueError(
      f"Matrix has {joint_csc.shape[1]} columns but {len(all_barcodes)} barcodes"
    )

  if len(all_barcodes) != len(set(all_barcodes)):
    dups = [bc for bc in all_barcodes if all_barcodes.count(bc) > 1]
    raise ValueError(f"Duplicate cell IDs after prefixing: {dups[:10]}")

  # per-cell library size; used by _normalize_hvg_mat in integration.py
  cell_totals = np.asarray(joint_csc.sum(axis=0)).flatten()

  _save_joint_outputs(
    joint_csc, hvg_list, all_barcodes, all_coldata_dfs,
    cell_totals, all_smeta_dfs, out_h5_f, out_coldata_f, out_sample_meta_f
  )

  print("done!")


def _smart_concat(dfs):
  # 1. Generate "Truth" schemas (only columns with at least one non-null value)
  truth_schemas = []
  for i, df in enumerate(dfs):
    # We find columns where null_count is less than the total number of rows
    truth_schema = {
      col: dtype 
      for col, dtype in df.schema.items()
      if df.height > 0 and df[col].null_count() < df.height
    }
    truth_schemas.append(truth_schema)

  # 2. Check consistency and aggregate to master_schema
  master_schema = {}
  for i, schema in enumerate(truth_schemas):
    for col, dtype in schema.items():
      if col in master_schema:
        # If the column exists in master, verify the types match
        if master_schema[col] != dtype:
          raise ValueError(
            f"Schema conflict for column '{col}': "
            f"Found {master_schema[col]} in previous DFs, "
            f"but DF index {i} has {dtype} (with non-null data)."
          )
      else:
        master_schema[col] = dtype

  # 3. Coerce DataFrames to follow the master_schema
  final_dfs = []
  for df in dfs:
    # Identify columns in this DF that need a cast to match the Master Truth
    casts = [
      pl.col(col).cast(master_schema[col])
      for col, dtype in df.schema.items()
      if col in master_schema and dtype != master_schema[col]
    ]
    
    if casts:
      final_dfs.append(df.with_columns(casts))
    else:
      final_dfs.append(df)

  return pl.concat(final_dfs, how='diagonal')

# ---------------------------------------------------------------------------
# Step 3: Build h5ads YAML with symlinks
# ---------------------------------------------------------------------------

def build_join_h5ads_yaml(h5ads_yaml_fs, project_ids, h5ads_dir, out_yaml_f):
  """
  Create symlinks and a joint h5ads YAML manifest.

  The output YAML uses the extended format:
    {project_id}_{batch_key}:
      path: {h5ads_dir}/{project_id}_{batch_key}.h5ad   # symlink
      project_id: {project_id}

  Parameters
  ----------
  h5ads_yaml_fs : list  Per-project h5ads_clean_paths_*.yaml files.
  project_ids   : list  Project IDs (parallel).
  h5ads_dir     : str   Directory for symlinks.
  out_yaml_f    : str   Output joint h5ads YAML.
  """
  print("building join h5ads YAML")
  pathlib.Path(h5ads_dir).mkdir(parents=True, exist_ok=True)

  joint_manifest = {}
  for pid, h5ads_yaml_f in zip(project_ids, h5ads_yaml_fs):
    with open(h5ads_yaml_f) as fh:
      h5ad_paths = yaml.safe_load(fh)

    for batch_key, h5ad_entry in h5ad_paths.items():
      src_path   = h5ad_entry if isinstance(h5ad_entry, str) else h5ad_entry['path']
      joint_key  = f"{pid}_{batch_key}"
      link_name  = pathlib.Path(h5ads_dir) / f"{joint_key}.h5ad"

      # create or update symlink
      if link_name.exists() or link_name.is_symlink():
        link_name.unlink()
      link_name.symlink_to(src_path)

      joint_manifest[joint_key] = {
        'path':       str(link_name),
        'project_id': pid
      }
      print(f"  {joint_key} -> {src_path}")

  with open(out_yaml_f, 'w') as fh:
    yaml.dump(joint_manifest, fh, default_flow_style=False)

  print(f"  saved to {out_yaml_f}")
  print("done!")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args():
  parser = argparse.ArgumentParser(description='scprocess join utilities')
  sub    = parser.add_subparsers(dest='cmd')

  # --- select_joint_hvgs ---
  p1 = sub.add_parser('select_joint_hvgs')
  p1.add_argument('--var_stats_fs', nargs='+', required=True,
    help='Per-project standardized variance stats CSV.gz files')
  p1.add_argument('--project_ids', nargs='+', required=True,
    help='Project IDs (parallel to --var_stats_fs)')
  p1.add_argument('--n_hvgs', type=int, required=True)
  p1.add_argument('--out_f', required=True)

  # --- build_joint_matrix ---
  p2 = sub.add_parser('build_joint_matrix')
  p2.add_argument('--joint_hvgs_f',     required=True)
  p2.add_argument('--h5ads_yaml_fs',    nargs='+', required=True,
    help='Per-project h5ads_clean_paths YAML files')
  p2.add_argument('--project_ids',      nargs='+', required=True)
  p2.add_argument('--integrated_dt_fs', nargs='+', required=True,
    help='Per-project integrated_dt CSV.gz files')
  p2.add_argument('--sample_meta_fs',   nargs='+', required=True,
    help='Per-project sample metadata CSV files')
  p2.add_argument('--metadata_vars',    default='',
    help='Space-separated metadata variable names')
  p2.add_argument('--out_h5_f',         required=True)
  p2.add_argument('--out_coldata_f',    required=True)
  p2.add_argument('--out_sample_meta_f', required=True)

  # --- build_join_h5ads_yaml ---
  p3 = sub.add_parser('build_join_h5ads_yaml')
  p3.add_argument('--h5ads_yaml_fs', nargs='+', required=True)
  p3.add_argument('--project_ids',   nargs='+', required=True)
  p3.add_argument('--h5ads_dir',     required=True)
  p3.add_argument('--out_yaml_f',    required=True)

  return parser.parse_args()


if __name__ == '__main__':
  args = _parse_args()

  if args.cmd == 'select_joint_hvgs':
    select_joint_hvgs(
      var_stats_fs = args.var_stats_fs,
      project_ids  = args.project_ids,
      n_hvgs       = args.n_hvgs,
      out_f        = args.out_f
    )

  elif args.cmd == 'build_joint_matrix':
    build_joint_matrix(
      joint_hvgs_f       = args.joint_hvgs_f,
      h5ads_yaml_fs      = args.h5ads_yaml_fs,
      project_ids        = args.project_ids,
      integrated_dt_fs   = args.integrated_dt_fs,
      sample_meta_fs     = args.sample_meta_fs,
      metadata_vars_str  = args.metadata_vars,
      out_h5_f           = args.out_h5_f,
      out_coldata_f      = args.out_coldata_f,
      out_sample_meta_f  = args.out_sample_meta_f
    )

  elif args.cmd == 'build_join_h5ads_yaml':
    build_join_h5ads_yaml(
      h5ads_yaml_fs = args.h5ads_yaml_fs,
      project_ids   = args.project_ids,
      h5ads_dir     = args.h5ads_dir,
      out_yaml_f    = args.out_yaml_f
    )

  else:
    print(f"unknown command: {args.cmd}", file=sys.stderr)
    sys.exit(1)
