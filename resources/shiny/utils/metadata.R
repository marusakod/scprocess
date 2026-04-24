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
#
# Palette spec per variable (yaml_data$metadata$palettes[[v]]):
#   palette:  name  — call resolve_palette(name, n) to generate colours
#   colours:  [c1, c2, ...]  — explicit hex/named colours
#   values:   [l1, l2, ...]  — optional level ordering (defaults to freq order)
#
# yaml_data$metadata$cluster_palette: named palette for clusters (when no annotation.csv).
get_lvls_and_colours <- function(yaml_data, cluster_meta, sample_meta,
                                 metadata_vars, data_dir, nice_cols) {
  vars_lvls  = list()
  vars_pals  = list()
  yaml_pals  = yaml_data$metadata$palettes

  for (v in metadata_vars) {
    obs_vals = sample_meta[[ v ]]
    def_lvls = obs_vals %>% as.character %>% fct_infreq %>% levels

    if (v %in% names(yaml_pals)) {
      this_pal  = yaml_pals[[ v ]]
      this_lvls = if (!is.null(this_pal$values)) {
        assert_that( all(def_lvls %in% this_pal$values),
          msg = paste0("metadata_palettes$", v, "$values does not cover all observed levels: ",
                       paste(setdiff(def_lvls, this_pal$values), collapse = ", ")) )
        this_pal$values
      } else {
        def_lvls
      }

      if (!is.null(this_pal$colours)) {
        these_cols = unlist(this_pal$colours)  # list or vector -> character vector
        assert_that( length(this_lvls) == length(these_cols),
          msg = paste0("metadata_palettes$", v, ": length(colours) must equal number of levels (", length(this_lvls), ")") )
        this_pal = these_cols %>% setNames(this_lvls)
      } else if (!is.null(this_pal$palette)) {
        this_pal = resolve_palette(this_pal$palette, length(this_lvls)) %>% setNames(this_lvls)
      } else {
        this_pal = rep_len(nice_cols, length(this_lvls)) %>% setNames(this_lvls)
      }
    } else {
      this_lvls = def_lvls
      this_pal  = rep_len(nice_cols, length(this_lvls)) %>% setNames(this_lvls)
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
      n_cl     = length(cluster_labels)
      cl_pal   = yaml_data$metadata$cluster_palette
      cl_cols  = if (!is.null(cl_pal) && nchar(cl_pal) > 0)
        resolve_palette(cl_pal, n_cl)
      else
        rep_len(nice_cols, n_cl)
      vars_pals[[ 'cluster' ]] = cl_cols %>% setNames(cluster_labels)
    }
  } else {
    cluster_labels = NULL
    cluster_lvls   = cluster_meta[, .(n_cells = sum(n_cells)), by = cluster] %>%
      .[ order(-n_cells) ] %>% .$cluster
    n_cl    = length(cluster_lvls)
    cl_pal  = yaml_data$metadata$cluster_palette
    cl_cols = if (!is.null(cl_pal) && nchar(cl_pal) > 0)
      resolve_palette(cl_pal, n_cl)
    else
      nice_cols[ seq_len(n_cl) ]
    vars_lvls[[ 'cluster' ]] = cluster_lvls
    vars_pals[[ 'cluster' ]] = cl_cols %>% setNames(cluster_lvls)
  }

  list(
    vars_lvls      = vars_lvls,
    vars_pals      = vars_pals,
    cluster_labels = cluster_labels
  )
}
