import pandas as pd
import polars as pl

# zoom function: get list of all mean var files for zooms
def get_zoom_raw_mean_var_files(zoom_name, ZOOM_PARAMS, FULL_TAG, DATE_STAMP):
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
def extract_zoom_sample_statistics(qc_stats_f, BATCHES, BATCH_VAR, LABELS_F, LABELS_VAR, LABELS, MIN_N_SAMPLE, AMBIENT_METHOD):
 
    qc_df   = pl.read_csv(qc_stats_f).drop('n_cells')
    lbls_dt = pl.read_csv(LABELS_F)

    # keep selected labels
    lbls_dt = lbls_dt.filter(lbls_dt[LABELS_VAR].is_in(LABELS))

    # count the number of cells per sample
    zoom_sample_stats = (lbls_dt.group_by(BATCH_VAR).agg(
        pl.len().alias("n_cells")
    ))

    # add empty samples
    empty_ss = list(set(BATCHES) - set(zoom_sample_stats[BATCH_VAR].to_list()))
    empty_df = pl.DataFrame({BATCH_VAR: empty_ss, "n_cells": [0] * len(empty_ss)})
    zoom_sample_stats = pl.concat([zoom_sample_stats, empty_df])

    # identify samples that do not meet the minimum cell thresho
    zoom_sample_stats = zoom_sample_stats.with_columns(
        (zoom_sample_stats["n_cells"] < MIN_N_SAMPLE).alias("bad_zoom_qc")
    )

    # merge new and existing sample stats
    sample_df = qc_df.join(zoom_sample_stats, on=BATCH_VAR, how="left")

    # update 'bad_[batch_var]' column
    if AMBIENT_METHOD == 'cellbender':
        sample_df = sample_df.with_columns(
            (sample_df["bad_bender"] | sample_df["bad_qc"] | sample_df["bad_zoom_qc"]).alias(f'bad_{BATCH_VAR}')
        )
    else:
        sample_df = sample_df.with_columns(
            (sample_df["bad_qc"] | sample_df["bad_zoom_qc"]).alias(f'bad_{BATCH_VAR}')
        )

    # Check that at least 2 good samples remain
    good_smpls_count = (sample_df[f'bad_{BATCH_VAR}'] == False).sum()
    if good_smpls_count < 2:
        raise AssertionError("Fewer than 2 samples available for this zoom.")

    return sample_df


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


