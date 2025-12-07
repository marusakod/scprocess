import pandas as pd

# zoom function: get list of all mean var files for zooms
def get_zoom_raw_mean_var_files(zoom_name, ZOOM_PARAMS, FULL_TAG, DATE_STAMP):
  group_names = ZOOM_PARAMS[zoom_name]['hvg']['hvg_group_names']
  num_chunks  = ZOOM_PARAMS[zoom_name]['hvg']['hvg_num_chunks']

  return [
    zoom_dir + f'/{zoom_name}/tmp_mean_var_{group}_group_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
    for group in group_names
    for chunk in range(num_chunks)
  ]


# zoom function: get list of all mean var files for zooms
def get_zoom_std_var_stats_files(zoom_name, zoom_dir, ZOOM_PARAMS, FULL_TAG, DATE_STAMP, SAMPLES):
  hvg_method = ZOOM_PARAMS[zoom_name]['hvg']['hvg_method']

  if hvg_method == "sample":
    return [
      zoom_dir + f'/{zoom_name}/tmp_std_var_stats_{sample}_sample_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
      for sample in SAMPLES
    ]
  else:
    group_names = ZOOM_PARAMS[zoom_name]['hvg']['hvg_group_names']
    num_chunks  = ZOOM_PARAMS[zoom_name]['hvg']['hvg_num_chunks']

    return [
      zoom_dir + f'/{zoom_name}/tmp_std_var_stats_{group}_group_chunk_{chunk}_' + FULL_TAG + '_' + DATE_STAMP + '.txt.gz'
      for group in group_names
      for chunk in range(num_chunks)
    ]


# zoom function: make df with good / bad sample labels for a specific zoom
def extract_zoom_sample_statistics(qc_stats_f, SAMPLES, LABELS_F, LABELS_VAR, LABELS, MIN_N_SAMPLE, AMBIENT_METHOD):
  # load inputs
  qc_df     = pd.read_csv(qc_stats_f)
  qc_df     = qc_df.drop('n_cells', axis=1)
  lbls_dt   = pd.read_csv(LABELS_F, compression='gzip')

  # keep selected labels
  lbls_dt   = lbls_dt[ lbls_dt[LABELS_VAR].isin(LABELS) ]
  
  # count the number of cells per sample
  zoom_sample_stats = (
    lbls_dt.groupby('sample_id')
    .size()
    .reset_index(name='n_cells')
  )
  
  # add empty samples
  empty_ss  = list(set(SAMPLES) - set(zoom_sample_stats["sample_id"].tolist()))
  empty_df  = pd.DataFrame({ "sample_id": empty_ss, "n_cells": 0 })
  zoom_sample_stats = pd.concat([zoom_sample_stats, empty_df])

  # identify samples that do not meet the minimum cell threshold
  zoom_sample_stats['bad_zoom_qc'] = zoom_sample_stats['n_cells'] < MIN_N_SAMPLE
  
  # merge new and existing sample stats
  sample_df = qc_df.merge(zoom_sample_stats, on='sample_id',how='left')
  
  # update 'bad_sample' column
  if AMBIENT_METHOD == 'cellbender':
    sample_df['bad_sample'] = (
      sample_df['bad_bender'] | sample_df['bad_qc'] | sample_df['bad_zoom_qc']
    )
  else:
    sample_df['bad_sample'] = (
      sample_df['bad_qc'] | sample_df['bad_zoom_qc']
    )

  # check that at least 2 good samples remain
  good_smpls_count = (sample_df['bad_sample'] == False).sum()
  assert good_smpls_count >= 2, \
    "Fewer than 2 samples available for this zoom."
  
  return sample_df


# zoom function: specify some optional outputs for zoom (at the moment only FGSEA outputs)
def get_zoom_conditional_fgsea_files(species, zoom_dir, FULL_TAG, DATE_STAMP, do_gsea):
    if do_gsea and (species in ['human_2024', 'human_2020', 'mouse_2024', 'mouse_2020']):
        return {
            'fgsea_go_bp_f': zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_go_bp_' + DATE_STAMP + '.txt.gz', 
            'fgsea_go_cc_f': zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_go_cc_' + DATE_STAMP + '.txt.gz',
            'fgsea_go_mf_f': zoom_dir + '/{zoom_name}/fgsea_' + FULL_TAG  + '_{mkr_sel_res}_go_mf_' + DATE_STAMP + '.txt.gz'
        }
    else:
        return {}


