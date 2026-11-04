# The project index lists classifier training reports separately from zooms;
# a zoom training report must not be parsed as a zoom resolution report.
suppressPackageStartupMessages(library(data.table))
source('scripts/render_htmls.R')

docs_dir = tempfile(); dir.create(docs_dir)
for (f in c('mini_integration.html', 'mini_zoom_neurons_0.5.html', 'mini_zoom_neurons_smol_1.html',
  'mini_train_xgboost_brain.html', 'mini_zoom_neurons_smol_train_xgboost_neurons_test.html'))
  file.create(file.path(docs_dir, f))
config_f = file.path(docs_dir, 'config.yaml'); writeLines('project: {}', config_f)

sub_ls = get_sub_ls('index', proj_dir = docs_dir, your_name = 'a', affiliation = 'b', short_tag = 'mini',
  docs_dir = docs_dir, full_tag = 'mini_x', date_stamp = '2026-01-01', mkr_sel_res = '0.2',
  config_f = config_f, show_arv_uuids = FALSE)
stopifnot(sub_ls$zoom_title == '## Subclustering', sub_ls$train_title == '## Classifier training',
  grepl('neurons_0.5 ([link](mini_zoom_neurons_0.5.html))', sub_ls$zoom_links, fixed = TRUE),
  !grepl('train', sub_ls$zoom_links),
  grepl('- xgboost_brain ([link](mini_train_xgboost_brain.html))', sub_ls$train_links, fixed = TRUE),
  grepl('- neurons_smol zoom: xgboost_neurons_test ([link](mini_zoom_neurons_smol_train_xgboost_neurons_test.html))',
    sub_ls$train_links, fixed = TRUE))
message('SUCCESS: index lists zoom and classifier training reports separately')
