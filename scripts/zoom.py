import polars as pl

# zoom function: get list of all mean var files for zooms
def get_zoom_raw_mean_var_files(zoom_name, zoom_dir, ZOOM_PARAMS, FULL_TAG, DATE_STAMP):
  group_names = ZOOM_PARAMS[zoom_name]['hvg']['hvg_group_names']
  num_chunks  = ZOOM_PARAMS[zoom_name]['hvg']['hvg_num_chunks']

  return [
    f'{zoom_dir}/{zoom_name}/tmp_mean_var_{group}_group_chunk_{chunk}_{FULL_TAG}_{zoom_name}_{DATE_STAMP}.csv.gz'
    for group in group_names
    for chunk in range(num_chunks)
  ]


# zoom function: get list of all mean var files for zooms
def get_zoom_std_var_stats_files(zoom_name, zoom_dir, ZOOM_PARAMS, FULL_TAG, DATE_STAMP, BATCHES):
  hvg_method = ZOOM_PARAMS[zoom_name]['hvg']['hvg_method']

  if hvg_method == "sample":
    return [
      f'{zoom_dir}/{zoom_name}/tmp_std_var_stats_{batch}_sample_{FULL_TAG}_{zoom_name}_{DATE_STAMP}.csv.gz'
      for batch in BATCHES
    ]
  else:
    group_names = ZOOM_PARAMS[zoom_name]['hvg']['hvg_group_names']
    num_chunks  = ZOOM_PARAMS[zoom_name]['hvg']['hvg_num_chunks']

    return [
      f'{zoom_dir}/{zoom_name}/tmp_std_var_stats_{group}_group_chunk_{chunk}_{FULL_TAG}_{zoom_name}_{DATE_STAMP}.csv.gz'
      for group in group_names
      for chunk in range(num_chunks)
    ]


def zoom_check_clusters(labels_f, labels_col, sel_labels_str, check_f):
  lbls_df     = pl.read_csv(labels_f)
  labels      = lbls_df[labels_col].unique().to_list()
  labels      = [cl for cl in labels if cl is not None]
  sel_labels  = sel_labels_str.split(',')
  missing_cls = set(sel_labels) - set(labels)
  if len(missing_cls) > 0:
    raise ValueError(
      f"the following labels were specified in the zoom params yaml but are not present in the file:\n"
      f"  {', '.join(missing_cls)}")
  with open(check_f, 'w') as f:
    f.write('ok\n')


# zoom function: make df with good / bad sample labels for a specific zoom
def extract_zoom_sample_statistics(qc_stats_f, labels_f, labels_col, sel_labels, batches, batch_var, min_n_sample, ambient_method):
  # load inputs
  qc_df       = pl.read_csv(qc_stats_f).drop("n_cells")
  lbls_df     = pl.read_csv(labels_f)

  # keep selected labels
  lbls_df     = lbls_df.filter( pl.col(labels_col).is_in(sel_labels) )

  # count the number of cells per sample
  zoom_stats  = lbls_df[batch_var].value_counts(name = "n_cells")

  # add empty samples
  empty_ss    = list(set(batches) - set(zoom_stats[batch_var]))
  empty_df    = pl.DataFrame({ batch_var: empty_ss, "n_cells": 0 }).cast(zoom_stats.schema)
  zoom_stats  = pl.concat([zoom_stats, empty_df])

  # identify samples that do not meet the minimum cell threshold
  zoom_stats  = zoom_stats.with_columns( (pl.col('n_cells') < min_n_sample).alias('bad_zoom_qc') )

  # merge new and existing sample stats
  batches_df  = qc_df.join(zoom_stats, on=batch_var, how='left')

  # update 'bad_sample' column
  bad_batch_col = f'bad_{batch_var}'
  if ambient_method == 'cellbender':
    batches_df  = batches_df.with_columns( (pl.col('bad_bender') | pl.col('bad_qc') | pl.col('bad_zoom_qc')).alias(bad_batch_col) )
  else:
    batches_df  = batches_df.with_columns( (pl.col('bad_qc') | pl.col('bad_zoom_qc')).alias(bad_batch_col) )

  # check that at least 2 good samples remain
  good_batches_count = batches_df.filter(pl.col(bad_batch_col) == False).shape[0]
  assert good_batches_count >= 2, \
    "Fewer than 2 samples available for this zoom."

  return batches_df


def filter_zoom_labels_by_qc(labels_f, qc_all_f, output_f, labels_col, sel_labels,
    qc_min_counts=None, qc_min_feats=None, qc_min_mito=None, qc_max_mito=None,
    qc_min_splice=None, qc_max_splice=None):
  import math
  import gzip

  thresholds = {
    'qc_min_counts': qc_min_counts, 'qc_min_feats': qc_min_feats,
    'qc_min_mito': qc_min_mito, 'qc_max_mito': qc_max_mito,
    'qc_min_splice': qc_min_splice, 'qc_max_splice': qc_max_splice
  }
  has_thresholds = any(v is not None for v in thresholds.values())

  lbls_df = pl.read_csv(labels_f)

  if not has_thresholds:
    with gzip.open(output_f, 'wb') as f:
      lbls_df.write_csv(f)
    print(f"No zoom QC thresholds specified. Wrote {lbls_df.shape[0]} rows unchanged.")
    return

  # get cell_ids matching selected labels
  sel_cell_ids = lbls_df.filter(pl.col(labels_col).cast(pl.Utf8).is_in(sel_labels))['cell_id'].to_list()

  # load qc_all (has log_counts, log_feats, logit_mito, logit_spliced, keep, cell_id)
  qc_df = pl.read_csv(qc_all_f)
  qc_df = qc_df.filter(pl.col('keep') & pl.col('cell_id').is_in(sel_cell_ids))

  # apply thresholds using pre-computed transformed columns (same as main QC)
  if qc_min_counts is not None:
    qc_df = qc_df.filter(pl.col('log_counts') >= math.log10(qc_min_counts))
  if qc_min_feats is not None:
    qc_df = qc_df.filter(pl.col('log_feats') >= math.log10(qc_min_feats))
  if qc_max_mito is not None and qc_max_mito < 1:
    qc_df = qc_df.filter(pl.col('logit_mito') < math.log(qc_max_mito / (1 - qc_max_mito)))
  if qc_min_mito is not None and qc_min_mito > 0:
    qc_df = qc_df.filter(pl.col('logit_mito') > math.log(qc_min_mito / (1 - qc_min_mito)))
  if qc_max_splice is not None and qc_max_splice < 1:
    qc_df = qc_df.filter(pl.col('logit_spliced') < math.log(qc_max_splice / (1 - qc_max_splice)))
  if qc_min_splice is not None and qc_min_splice > 0:
    qc_df = qc_df.filter(pl.col('logit_spliced') > math.log(qc_min_splice / (1 - qc_min_splice)))

  surviving_ids = set(qc_df['cell_id'].to_list())
  failed_ids = set(sel_cell_ids) - surviving_ids
  n_removed = len(failed_ids)
  print(f"Zoom QC filtering: {n_removed} / {len(sel_cell_ids)} cells removed ({len(surviving_ids)} remaining)")

  filtered_lbls_df = lbls_df.filter(~pl.col('cell_id').is_in(list(failed_ids)))
  with gzip.open(output_f, 'wb') as f:
    filtered_lbls_df.write_csv(f)
  print(f"Wrote filtered labels file with {filtered_lbls_df.shape[0]} rows to {output_f}")


# zoom function: specify some optional outputs for zoom (at the moment only FGSEA outputs)
def get_zoom_conditional_fgsea_files(ref_txome, zoom_dir, FULL_TAG, DATE_STAMP, do_gsea):
  valid_refs = ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020',
                'human_v1', 'mouse_v1', 'human_v2', 'mouse_v2']
  if do_gsea and (ref_txome in valid_refs):
    return {
      'fgsea_go_bp_f': f'{zoom_dir}/{{zoom_name}}/fgsea_{FULL_TAG}_{{zoom_name}}_{{mkr_sel_res}}_go_bp_{DATE_STAMP}.csv.gz', 
      'fgsea_go_cc_f': f'{zoom_dir}/{{zoom_name}}/fgsea_{FULL_TAG}_{{zoom_name}}_{{mkr_sel_res}}_go_cc_{DATE_STAMP}.csv.gz',
      'fgsea_go_mf_f': f'{zoom_dir}/{{zoom_name}}/fgsea_{FULL_TAG}_{{zoom_name}}_{{mkr_sel_res}}_go_mf_{DATE_STAMP}.csv.gz'
    }
  else:
    return {}


def get_zoom_conditional_xgboost_files(zoom_params, zoom_dir, zoom_name, code_dir):
  if 'train_xgboost' in zoom_params:
    ref_tag = zoom_params['train_xgboost']['ref_tag']
    return {
      'xgb_preds_f':   f'{zoom_dir}/{zoom_name}/train_xgboost/{ref_tag}_predictions.csv.gz',
      'xgb_imp_f':     f'{zoom_dir}/{zoom_name}/train_xgboost/{ref_tag}_gene_importance.csv',
      'xgb_pb_f':      f'{zoom_dir}/{zoom_name}/train_xgboost/{ref_tag}_pseudobulk.h5ad',
      'xgb_r_train_f': f'{code_dir}/train_xgboost.R',
    }
  else:
    return {}

