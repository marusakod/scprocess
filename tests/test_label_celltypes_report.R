# Labelling report helpers: installed label styles, Labels / Grouped labels
# tables and colours for many categories.
suppressPackageStartupMessages(library(data.table))
source('scripts/utils.R')
source('scripts/label_celltypes.R')

# A registry with several models: only the requested one may be used. A bare
# `model` inside data.table brackets would match every row and drop the styling.
scdata_dir = tempfile(); dir.create(file.path(scdata_dir, 'xgboost'), recursive = TRUE)
specs_f = file.path(scdata_dir, 'xgboost', 'm2_label_specs.csv')
coarse_f = file.path(scdata_dir, 'xgboost', 'm2_coarse_label_specs.csv')
fwrite(data.table(label = c('Neuron A', 'Astrocyte', 'Oligodendrocyte'),
  colour = c('#111111', '#222222', '#333333'), coarse_label = c('Neuron', 'Glia', 'Glia')), specs_f)
fwrite(data.table(coarse_label = c('Neuron', 'Glia'), coarse_colour = c('#444444', '#555555')), coarse_f)
fwrite(data.table(model = c('m1', 'm2', 'm3'), label_specs_f = c(NA, specs_f, NA),
  coarse_label_specs_f = c(NA, coarse_f, NA)), file.path(scdata_dir, 'xgboost', 'xgboost_models.csv'))
Sys.setenv(SCPROCESS_DATA_DIR = scdata_dir)

style = installed_label_style('m2')
stopifnot(identical(style$order, c('Neuron A', 'Astrocyte', 'Oligodendrocyte')),
  identical(unname(style$colours[c('Astrocyte', 'Neuron A')]), c('#222222', '#111111')),
  identical(style$coarse_order, c('Neuron', 'Glia')),
  identical(unname(style$mapping[['Oligodendrocyte']]), 'Glia'),
  is.null(installed_label_style('m1')$order))

# The grouped table holds only the grouped columns, under the standard names
guesses_dt = data.table(cell_id = c('c1', 'c2', 'c3'),
  predicted_label_naive = c('Neuron A', 'Astrocyte', 'Oligodendrocyte'),
  predicted_label_agg = c('Neuron A', 'Astrocyte', 'Astrocyte'),
  coarse_predicted_label_naive = c('Neuron', 'Glia', 'Glia'),
  coarse_predicted_label_agg = c('Neuron', 'Glia', 'Glia'))
views = get_label_views(guesses_dt, style)
stopifnot(identical(names(views), c('Labels', 'Grouped labels')),
  identical(names(views[['Grouped labels']]$labels_dt), c('cell_id', 'predicted_label_naive', 'predicted_label_agg')),
  setequal(views[['Grouped labels']]$labels_dt$predicted_label_naive, c('Neuron', 'Glia')))
grouped_table = calc_labels_table(views[['Grouped labels']]$labels_dt, label_order = style$coarse_order)
stopifnot(identical(as.character(grouped_table[[1]]), c('Neuron', 'Glia')),
  identical(names(get_label_views(guesses_dt[, 1:3], style)), 'Labels'))

# Colours repeat nice_cols once it runs out, and supplied colours win
many = paste0('cl', 1:(length(nice_cols) + 2))
cols = xgboost_colours(many, c(cl2 = '#000000'))
stopifnot(identical(unname(cols[length(nice_cols) + 1]), nice_cols[[1]]), cols[['cl2']] == '#000000')
# Headings: a tabset with one tab per view only when there are grouped labels
two = capture.output({ cat_label_views_heading('scprocess, m2', views, level = 3); cat_label_view_tab('Labels', views, level = 3) })
one = capture.output({ cat_label_views_heading('scprocess, m2', views['Labels'], level = 3); cat_label_view_tab('Labels', views['Labels'], level = 3) })
stopifnot(identical(two[nzchar(two)], c('### scprocess, m2 {.tabset}', '#### Labels')),
  identical(one[nzchar(one)], '### scprocess, m2'))
message('SUCCESS: installed label styles, grouped label tables, view headings and repeated palette')
