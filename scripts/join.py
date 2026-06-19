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
  """Select top HVGs by mean rank across projects."""
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
  """Return a boolean filter expression for non-doublet, demultiplexed cells.
  If is_dbl/in_dbl_cl columns are absent (e.g. zoom integrated_dt where
  doublets were already removed upstream), all cells are considered clean.
  Cells with null/empty sample_id are excluded (undemultiplexed HTO cells).
  """
  filt = pl.lit(True)
  if 'is_dbl' in int_dt.columns and 'in_dbl_cl' in int_dt.columns:
    filt = filt & (pl.col('is_dbl') == False) & (pl.col('in_dbl_cl') == False)
  if 'sample_id' in int_dt.columns:
    filt = filt & pl.col('sample_id').is_not_null() & (pl.col('sample_id').cast(pl.Utf8) != 'None') & (pl.col('sample_id').cast(pl.Utf8) != '')
  return filt


def _check_sample_id_uniqueness(project_ids, integrated_dt_fs):
  """Raise ValueError if any sample_id (among clean cells) appears in >1 project."""
  seen = {}
  for pid, int_f in zip(project_ids, integrated_dt_fs):
    int_dt  = pl.read_csv(int_f)
    samples = int_dt.filter(_ok_cells_filter(int_dt))['sample_id'].unique().to_list()
    for s in samples:
      if s in seen and seen[s] != pid:
        raise ValueError(
          f"sample_id '{s}' appears in both project '{seen[s]}' and '{pid}'. "
          f"Use unique sample IDs across projects."
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


def _build_project_coldata(int_dt, pid):
  """Filter to clean cells and prefix sample_id with project ID."""
  proj_cells = int_dt.filter(_ok_cells_filter(int_dt)).select(
    ['cell_id', 'sample_id'] +
    [c for c in int_dt.columns if c not in ['cell_id', 'sample_id', 'project_id']])
  proj_cells = proj_cells.with_columns([
    pl.concat_str([pl.lit(f"{pid}_"), pl.col('sample_id').cast(pl.Utf8)]).alias('sample_id'),
    pl.lit(pid).alias('project_id')
  ])
  return proj_cells


def build_joint_matrix(joint_hvgs_f, h5ads_yaml_fs, project_ids, integrated_dt_fs,
                       out_h5_f):
  """Assemble a joint HVG count matrix from per-project h5ads."""
  import numpy as np
  import h5py
  from scipy.sparse import hstack

  print("building joint count matrix")

  hvg_df   = pl.read_csv(joint_hvgs_f)
  hvg_list = hvg_df['gene_id'].to_list()
  print(f"  joint HVGs: {len(hvg_list)}")

  _check_sample_id_uniqueness(project_ids, integrated_dt_fs)

  all_mats     = []
  all_barcodes = []

  for pid, h5ads_yaml_f, int_f in zip(project_ids, h5ads_yaml_fs, integrated_dt_fs):
    print(f"  processing project: {pid}")
    with open(h5ads_yaml_f) as fh:
      h5ad_paths = yaml.safe_load(fh)
    int_dt   = pl.read_csv(int_f)
    ok_cells = set(int_dt.filter(_ok_cells_filter(int_dt))['cell_id'].to_list())
    print(f"    clean cells: {len(ok_cells)}")

    for batch_key, h5ad_entry in h5ad_paths.items():
      h5ad_path = h5ad_entry if isinstance(h5ad_entry, str) else h5ad_entry['path']
      result = _load_batch_hvg_matrix(h5ad_path, ok_cells, hvg_list, pid, batch_key)
      if result is None:
        continue
      csc_sub, kept_bcs = result
      all_mats.append(csc_sub)
      all_barcodes.extend(kept_bcs)

  if len(all_mats) == 0:
    raise ValueError("No cells loaded from any project")

  print("  concatenating matrices")
  joint_csc = hstack(all_mats, format='csc').tocsc() if len(all_mats) > 1 else all_mats[0]
  print(f"  joint matrix shape: {joint_csc.shape} (genes x cells)")

  if joint_csc.shape[1] != len(all_barcodes):
    raise ValueError(
      f"Matrix has {joint_csc.shape[1]} columns but {len(all_barcodes)} barcodes")

  if len(all_barcodes) != len(set(all_barcodes)):
    from collections import Counter
    counts = Counter(all_barcodes)
    dups = [bc for bc, n in counts.items() if n > 1]
    raise ValueError(f"Duplicate cell IDs: {dups[:10]}")

  print(f"  saving matrix to {out_h5_f}")
  pathlib.Path(out_h5_f).parent.mkdir(parents=True, exist_ok=True)
  with h5py.File(out_h5_f, 'w') as f:
    f.create_dataset('matrix/data',          data=joint_csc.data)
    f.create_dataset('matrix/indices',       data=joint_csc.indices)
    f.create_dataset('matrix/indptr',        data=joint_csc.indptr)
    f.create_dataset('matrix/shape',         data=joint_csc.shape)
    f.create_dataset('matrix/features/name', data=np.array(hvg_list, dtype='S'))
    f.create_dataset('matrix/features/id',   data=np.array(hvg_list, dtype='S'))
    f.create_dataset('matrix/barcodes',      data=np.array(all_barcodes, dtype='S'))

  print("done!")


def build_joint_coldata(h5_f, project_ids, integrated_dt_fs, sample_meta_fs,
                        out_coldata_f, out_sample_meta_f):
  """Build joint coldata and sample metadata from matrix barcodes."""
  import numpy as np
  import h5py

  print("building joint coldata")

  # read barcodes from HDF5 to define cell order
  with h5py.File(h5_f, 'r') as f:
    all_barcodes = [b.decode('utf-8') for b in f['matrix/barcodes'][:]]
  print(f"  {len(all_barcodes)} cells in matrix")

  # build per-project coldata and sample metadata
  all_coldata_dfs = []
  all_smeta_dfs   = []

  for pid, int_f, smeta_f in zip(project_ids, integrated_dt_fs, sample_meta_fs):
    print(f"  processing project: {pid}")
    int_dt   = pl.read_csv(int_f)
    smeta_dt = pl.read_csv(smeta_f)

    coldata_df = _build_project_coldata(int_dt, pid)
    all_coldata_dfs.append(coldata_df)

    smeta_df = smeta_dt.with_columns([
      pl.concat_str([pl.lit(f"{pid}_"), pl.col('sample_id').cast(pl.Utf8)]).alias('sample_id'),
      pl.lit(pid).alias('project_id')
    ])
    if 'bad_sample_id' not in smeta_df.columns:
      smeta_df = smeta_df.with_columns(pl.lit(False).alias('bad_sample_id'))
    all_smeta_dfs.append(smeta_df)

  # concatenate and reorder to match matrix barcodes
  print("  assembling coldata")
  coldata_df = pl.concat(all_coldata_dfs, how='diagonal')
  fixed_cols = ['cell_id', 'sample_id', 'project_id']
  other_cols = [c for c in coldata_df.columns if c not in fixed_cols]
  coldata_df = coldata_df.select(fixed_cols + other_cols)

  bc_order_df = pl.DataFrame({
    'cell_id': all_barcodes,
    '_order':  range(len(all_barcodes)),
  })
  coldata_df = coldata_df.join(bc_order_df, on='cell_id').sort('_order').drop('_order')
  with gzip.open(out_coldata_f, 'wb') as fh:
    coldata_df.write_csv(fh)

  print("  saving sample metadata")
  smeta_df = _smart_concat(all_smeta_dfs)
  smeta_df.write_csv(out_sample_meta_f)

  print("done!")


def _smart_concat(dfs):
  """Concat DataFrames with mixed schemas, promoting Int to Float on conflict."""
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
  for schema in truth_schemas:
    for col, dtype in schema.items():
      if col in master_schema:
        prev_dtype = master_schema[col]
        if prev_dtype != dtype:
          # Resolve Int vs Float conflict: promote to Float
          if (prev_dtype.is_integer() and dtype.is_float()):
            master_schema[col] = dtype
          elif (prev_dtype.is_float() and dtype.is_integer()):
            continue # Keep the Float already in master_schema
          else:
            # For String vs Int, etc., you might still want an error
            raise ValueError(f"Hard conflict for '{col}': {prev_dtype} vs {dtype}")
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
  """Create symlinks and a joint h5ads YAML manifest."""
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
  p2.add_argument('--out_h5_f',         required=True)

  # --- build_joint_coldata ---
  p2b = sub.add_parser('build_joint_coldata')
  p2b.add_argument('--h5_f',             required=True,
    help='Joint HVG count matrix HDF5 (for barcodes)')
  p2b.add_argument('--project_ids',      nargs='+', required=True)
  p2b.add_argument('--integrated_dt_fs', nargs='+', required=True,
    help='Per-project integrated_dt CSV.gz files')
  p2b.add_argument('--sample_meta_fs',   nargs='+', required=True,
    help='Per-project sample metadata CSV files')
  p2b.add_argument('--out_coldata_f',    required=True)
  p2b.add_argument('--out_sample_meta_f', required=True)

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
      out_h5_f           = args.out_h5_f
    )

  elif args.cmd == 'build_joint_coldata':
    build_joint_coldata(
      h5_f               = args.h5_f,
      project_ids        = args.project_ids,
      integrated_dt_fs   = args.integrated_dt_fs,
      sample_meta_fs     = args.sample_meta_fs,
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
