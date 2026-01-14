import pandas as pd
import polars as pl

# zoom function: get list of all mean var files for zooms
def get_zoom_raw_mean_var_files(zoom_name, zoom_dir, ZOOM_PARAMS, FULL_TAG, DATE_STAMP):
  group_names = ZOOM_PARAMS[zoom_name]['hvg']['hvg_group_names']
  num_chunks  = ZOOM_PARAMS[zoom_name]['hvg']['hvg_num_chunks']

  return [
    zoom_dir + f'/{zoom_name}/tmp_mean_var_{group}_group_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
    for group in group_names
    for chunk in range(num_chunks)
  ]


# zoom function: get list of all mean var files for zooms
def get_zoom_std_var_stats_files(zoom_name, zoom_dir, ZOOM_PARAMS, FULL_TAG, DATE_STAMP, BATCHES):
  hvg_method = ZOOM_PARAMS[zoom_name]['hvg']['hvg_method']

  if hvg_method == "sample":
    return [
      zoom_dir + f'/{zoom_name}/tmp_std_var_stats_{batch}_sample_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
      for batch in BATCHES
    ]
  else:
    group_names = ZOOM_PARAMS[zoom_name]['hvg']['hvg_group_names']
    num_chunks  = ZOOM_PARAMS[zoom_name]['hvg']['hvg_num_chunks']

    return [
      zoom_dir + f'/{zoom_name}/tmp_std_var_stats_{group}_group_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.csv.gz'
      for group in group_names
      for chunk in range(num_chunks)
    ]


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
  empty_df    = pl.DataFrame({ batch_var: empty_ss, "n_cells": 0 })
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


# zoom function: specify some optional outputs for zoom (at the moment only FGSEA outputs)
def get_zoom_conditional_fgsea_files(species, zoom_dir, FULL_TAG, DATE_STAMP, do_gsea):
  if do_gsea and (species in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020']):
    return {
      'fgsea_go_bp_f': zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_go_bp_' + DATE_STAMP + '.csv.gz', 
      'fgsea_go_cc_f': zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_go_cc_' + DATE_STAMP + '.csv.gz',
      'fgsea_go_mf_f': zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_go_mf_' + DATE_STAMP + '.csv.gz'
    }
  else:
    return {}

