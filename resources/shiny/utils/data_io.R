# DATA I/O FUNCTIONS -----------------------------------------------------------

# Locate and validate all required shiny input files in data_dir.
# Returns a named list: HDF5 entries are file paths; others are loaded data.tables.
get_all_input_fs <- function(data_dir, date_stamp, paga = 0) {

  f_dir = gsub('\\/$', '', data_dir)
  f_ls  = list.files(f_dir, recursive = FALSE, full.names = TRUE)

  all_fs = c(
    "count_h5_f"       = '-shiny_norm_count',
    "sample_pb_h5_f"   = '-shiny_sample_pb_count',
    "cluster_pb_h5_f"  = '-shiny_cluster_pb_count',
    "row_indx"         = '-shiny_row_indices',
    "cluster_markers"  = '-shiny_markers',
    "pb_hvgs"          = '-shiny_pb_hvgs',
    "fgsea"            = '-shiny_fgsea_res',
    "go_terms"         = '-shiny_go_terms',
    "sample_meta"      = '-shiny_sample_meta',
    "cluster_meta"     = '-shiny_cluster_meta',
    "cell_meta"        = '-shiny_cell_meta',
    "centroids"        = '-shiny_centroids',
    "repel_pos"        = '-shiny_repel_pos'
  )

  if (paga == 1) {
    all_fs = c(all_fs,
      'paga_mat' = '-shiny_paga_mat',
      'paga_pos' = '-shiny_paga_pos'
    )
  }

  f_matches = sapply(all_fs, function(f) any(grepl(f, f_ls)))
  if (!all(f_matches)) {
    missing = all_fs[f_matches == FALSE]
    stop("Missing input files containing: ", paste(missing, collapse = ', '),
         '; Shiny app cannot be made:(')
  }

  f_ls  = f_ls[ sapply(all_fs, function(f) grep(f, f_ls)) %>% unlist() ]
  names(f_ls) = names(all_fs)

  non_h5_fs = f_ls[!grepl('h5', names(f_ls))]
  non_h5_ls = sapply(non_h5_fs, fread, simplify = FALSE, USE.NAMES = TRUE)

  out_ls = c(
    f_ls[!names(f_ls) %in% names(non_h5_fs)],
    non_h5_ls
  )

  return(out_ls)
}

# Extract a single row from a BPCells matrix directory and return it as a
# plain numeric vector.
get_row_from_h5 <- function(f_name, row_index) {
  mat     = BPCells::open_matrix_dir(f_name)
  row_mat = mat[as.integer(row_index), ]
  as.numeric(as(as(row_mat, "dgCMatrix"), "matrix"))
}
