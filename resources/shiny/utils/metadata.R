# METADATA UTILITY FUNCTIONS ---------------------------------------------------

# Create cross-product combination columns in a metadata data.table.
# vars: character(2) of column names to combine.
# vars_lvls: named list of factor levels for each var (names must match vars).
add_metadata_var_combns <- function(meta_dt, vars, vars_lvls) {
  assert_that( all(vars %in% colnames(meta_dt)) )
  assert_that( length(vars) == 2 )
  assert_that( all(vars == names(vars_lvls)) )

  combn_cols     = meta_dt[, vars, with = FALSE] %>% as.list()
  combn_lvls     = apply(rev(expand.grid(rev(vars_lvls))), 1, paste, collapse = ' ')
  combn_lvls_rev = apply(rev(expand.grid(vars_lvls)),      1, paste, collapse = ' ')

  meta_dt[, paste(vars, collapse = '.')]      = do.call(paste, c(combn_cols,       sep = " ")) %>%
    factor(., levels = combn_lvls)
  meta_dt[, paste(rev(vars), collapse = '.')] = do.call(paste, c(rev(combn_cols), sep = " ")) %>%
    factor(., levels = combn_lvls_rev)

  return(meta_dt)
}

# Bin a continuous column into a categorical factor and append it to meta_dt.
add_metadata_var_cont <- function(meta_dt, var, breaks, labels = NULL) {
  assert_that( var %in% colnames(meta_dt) )
  assert_that( is.numeric(meta_dt[[var]]) )

  var_cat = paste0(var, '_cat')
  meta_dt[[var_cat]] = cut(meta_dt[[var]], breaks = breaks, labels = labels)
  meta_dt = meta_dt[is.na(get(var_cat)), (var_cat) := 'not known']

  return(meta_dt)
}

# Derive factor levels and colour palettes for all metadata variables.
# Uses user-supplied YAML palettes where available; falls back to nice_cols.
# Also reads an optional annotation.csv for cluster labels and colours.
# Returns list(vars_lvls, vars_pals, cluster_labels).
get_lvls_and_colours <- function(yaml_data, cluster_meta, sample_meta,
                                 metadata_vars, data_dir, nice_cols) {
  vars_lvls  = list()
  vars_pals  = list()
  yaml_pals  = yaml_data$metadata$palettes

  for (v in metadata_vars) {
    obs_vals = sample_meta[[ v ]]
    def_lvls = obs_vals %>% as.character %>% fct_infreq %>% levels

    if (v %in% names(yaml_pals)) {
      this_pal = yaml_pals[[ v ]]
      assert_that( 'values' %in% names(this_pal) )
      assert_that( all(def_lvls %in% this_pal$values) )

      this_lvls = this_pal$values
      if ('colours' %in% names(this_pal)) {
        assert_that( length(this_pal$values) == length(this_pal$colours) )
        this_pal  = this_pal$colours %>% setNames(this_lvls)
      } else {
        this_pal  = nice_cols[ seq_along(this_lvls) ] %>% setNames(this_lvls)
      }
    } else {
      this_lvls = def_lvls
      this_pal  = nice_cols[ seq_along(this_lvls) ] %>% setNames(this_lvls)
    }

    vars_lvls[[ v ]] = this_lvls
    vars_pals[[ v ]] = this_pal
  }

  annotation_f = file.path(data_dir, 'annotation.csv')
  if (file.exists(annotation_f)) {
    annot_dt       = fread(annotation_f)
    cluster_labels = annot_dt$cluster_name %>% setNames(annot_dt$cluster)

    vars_lvls[[ 'cluster' ]] = cluster_labels
    if ('colour' %in% names(annot_dt)) {
      vars_pals[[ 'cluster' ]] = annot_dt$colour %>% setNames(cluster_labels)
    } else {
      vars_pals[[ 'cluster' ]] = nice_cols[ seq_along(cluster_labels) ] %>% setNames(cluster_labels)
    }
  } else {
    cluster_labels = NULL
    cluster_lvls   = cluster_meta[, .(n_cells = sum(n_cells)), by = cluster] %>%
      .[ order(-n_cells) ] %>% .$cluster
    vars_lvls[[ 'cluster' ]] = cluster_lvls
    vars_pals[[ 'cluster' ]] = nice_cols[ seq_along(cluster_lvls) ] %>% setNames(cluster_lvls)
  }

  list(
    vars_lvls      = vars_lvls,
    vars_pals      = vars_pals,
    cluster_labels = cluster_labels
  )
}
