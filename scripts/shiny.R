suppressPackageStartupMessages({
  library("assertthat")
  library("data.table")
  library("magrittr")
  library("yaml")
  library("BPCells")
  library("MASS")
  library("strex")
})


#' Build Shiny app from scprocess outputs
#'
#' @param integration_f  Path to integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz
#' @param h5ads_yaml_f   Path to h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml (maps batch -> h5ad path)
#' @param sample_meta_f  Path to sample metadata CSV (config['project']['sample_metadata'])
#' @param mkrs_f         Path to pb_marker_genes_{FULL_TAG}_{res}_{DATE_STAMP}.csv.gz
#' @param pb_hvgs_f      Path to pb_hvgs_{FULL_TAG}_{res}_{DATE_STAMP}.csv.gz
#' @param fgsea_bp_f     Path to fgsea go_bp CSV.GZ
#' @param fgsea_cc_f     Path to fgsea go_cc CSV.GZ
#' @param fgsea_mf_f     Path to fgsea go_mf CSV.GZ
#' @param deploy_dir     Output directory (e.g. public/shiny)
#' @param scprocess_dir  Path to the scprocess installation directory
#' @param app_tag        Short tag for output file names (config['project']['short_tag'])
#' @param date_stamp     Date stamp (config['project']['date_stamp'])
#' @param mkr_sel_res    Leiden resolution used for marker genes (e.g. "0.5")
#' @param ref_txome      Reference transcriptome name (used to determine species)
#' @param metadata_vars  Comma-separated metadata variable names
#' @param app_title      Shiny app title
#' @param email          Contact email for app footer
#' @param keyword        Short keyword appearing in plot descriptions (e.g. "cells")
#' @param default_gene   Default gene symbol shown in Explore Genes tab
#' @param n_keep         Number of cells to retain in the downsampled UMAP (default 30000)
#' @param var_names      Comma-separated display names for metadata_vars (defaults to metadata_vars)
#' @param var_combns     JSON-encoded list of variable combination pairs (default "[]")
#' @param home_md_f      Optional path to a Markdown file used as the landing page (overrides placeholder)
#' @param annotation_csv_f Optional path to a CSV with columns cluster, cluster_name, colour defining
#'                        display names, colour, and order for clusters
#' @param n_cores        Number of data.table threads
make_shiny_app_scprocess <- function(
  integration_f, h5ads_yaml_f, sample_meta_f,
  mkrs_f, pb_hvgs_f, fgsea_bp_f, fgsea_cc_f, fgsea_mf_f,
  deploy_dir, scprocess_dir,
  app_tag, date_stamp, mkr_sel_res,
  ref_txome, metadata_vars,
  app_title = app_tag,
  email     = "",
  keyword   = "cells",
  default_gene = "",
  n_keep    = 30000L,
  var_names = metadata_vars,
  var_combns = "[]",
  home_md_f = "",
  annotation_csv_f = "",
  n_cores   = 1L
) {
  setDTthreads(n_cores)

  # ---- Resolve species from ref_txome name --------------------------------
  if (grepl("GRCh|hg|human", ref_txome, ignore.case = TRUE)) {
    species  <- "human"
    gene_re  <- "-ENSG"
  } else if (grepl("GRCm|mm|mouse", ref_txome, ignore.case = TRUE)) {
    species  <- "mouse"
    gene_re  <- "-ENSMUSG"
  } else {
    stop("Cannot determine species from ref_txome name '", ref_txome,
         "'. Name must contain GRCh/hg/human or GRCm/mm/mouse.")
  }
  message("Species: ", species)

  # ---- Parse comma-separated args ----------------------------------------
  metadata_vars <- strsplit(metadata_vars, ",")[[1]] %>% trimws()
  var_names     <- strsplit(var_names, ",")[[1]] %>% trimws()
  var_combns    <- jsonlite::fromJSON(var_combns)

  assert_that(length(metadata_vars) == length(var_names),
    msg = "metadata_vars and var_names must have the same length")

  n_keep        <- as.integer(n_keep)
  cluster_col   <- paste0("RNA_snn_res.", mkr_sel_res)

  # ---- Set up output directories ------------------------------------------
  deploy_dir <- gsub("/$", "", deploy_dir)
  assert_that(dir.exists(deploy_dir), msg = paste("deploy_dir does not exist:", deploy_dir))
  data_dir   <- file.path(deploy_dir, "data")
  dir.create(data_dir, showWarnings = FALSE)
  dir.create(file.path(deploy_dir, "www"), showWarnings = FALSE)

  # ---- Copy app files from resources/shiny/ --------------------------------
  message("Copying app files")
  app_src <- file.path(scprocess_dir, "resources/shiny")
  assert_that(dir.exists(app_src),
    msg = paste("resources/shiny not found at", app_src))

  for (f in c("app.R", "constants.R")) {
    file.copy(file.path(app_src, f), deploy_dir, overwrite = TRUE)
  }
  for (subdir in c("utils", "modules")) {
    src  <- file.path(app_src, subdir)
    dest <- file.path(deploy_dir, subdir)
    dir.create(dest, showWarnings = FALSE, recursive = TRUE)
    file.copy(list.files(src, full.names = TRUE), dest, overwrite = TRUE)
  }

  # ---- Copy or create home.md --------------------------------------------
  home_md_dest <- file.path(data_dir, "home.md")
  if (nchar(home_md_f) > 0 && file.exists(home_md_f)) {
    file.copy(home_md_f, home_md_dest, overwrite = TRUE)
    message(" copied user home.md")
  } else if (!file.exists(home_md_dest)) {
    writeLines(c(
      "## Welcome",
      "",
      "This app allows you to explore the single-cell RNA-seq data interactively.",
      "Use the tabs above to:",
      "",
      "- **Explore genes** — visualise expression of individual genes on the UMAP and across clusters/samples",
      "- **Explore clusters** — browse marker genes for each cluster",
      "- **Explore genesets** — view heatmaps for gene sets of interest",
      "- **Explore prevalence** — examine cluster composition across samples and conditions",
      "",
      "_Edit this file (`data/home.md`) to customise the welcome message._"
    ), con = home_md_dest)
    message(" created placeholder home.md")
  }

  # ---- Read integrated_dt, extract UMAP + cluster + cell metadata ----------
  message("Reading integration outputs")
  int_dt    <- fread(integration_f)

  assert_that(cluster_col %in% colnames(int_dt),
    msg = paste0("Cluster column '", cluster_col, "' not found in integrated_dt. ",
                 "Available columns: ", paste(colnames(int_dt), collapse = ", ")))
  assert_that("cell_id"  %in% colnames(int_dt))
  assert_that("sample_id" %in% colnames(int_dt))
  assert_that("UMAP1"    %in% colnames(int_dt))
  assert_that("UMAP2"    %in% colnames(int_dt))

  # Clean cells: not a doublet and not in a doublet-enriched cluster
  if ("is_dbl" %in% colnames(int_dt) && "in_dbl_cl" %in% colnames(int_dt)) {
    int_dt <- int_dt[is_dbl == FALSE & in_dbl_cl == FALSE]
  } else {
    # If doublet columns absent (single-batch run without doublets), keep all
    message(" doublet columns not found - keeping all cells")
  }

  umaps  <- int_dt[, .(cell_id, umap_1 = UMAP1, umap_2 = UMAP2)]
  clusts <- int_dt[, .(cell_id, cluster = get(cluster_col))]

  # ---- Join with sample metadata to get per-cell metadata -----------------
  message("Joining with sample metadata")
  sample_meta <- fread(sample_meta_f)
  assert_that("sample_id" %in% colnames(sample_meta),
    msg = "sample_meta_f must contain a 'sample_id' column")

  keep_meta_cols <- unique(c("sample_id", metadata_vars))
  keep_meta_cols <- intersect(keep_meta_cols, colnames(sample_meta))

  cell_meta <- merge(
    int_dt[, .(cell_id, sample_id)],
    sample_meta[, keep_meta_cols, with = FALSE],
    by = "sample_id", all.x = TRUE
  )

  # ---- Merge all metadata into a single table -----------------------------
  all_meta <- cell_meta %>%
    merge(umaps,  by = "cell_id") %>%
    merge(clusts, by = "cell_id") %>%
    .[cluster %in% unique(fread(mkrs_f)$cluster), ]  # keep only clusters with markers

  # ---- Downsample UMAP -----------------------------------------------------
  message("Downsampling UMAP to ", n_keep, " cells")
  keep_cells <- .subsample_umap(all_meta[, .(cell_id, UMAP_1 = umap_1, UMAP_2 = umap_2)],
                                to_keep = n_keep)
  all_meta[, keep_cell := cell_id %in% keep_cells]

  # Convert cluster to factor (sorted)
  all_meta$cluster <- factor(all_meta$cluster, levels = sort(unique(all_meta$cluster)))

  # ---- Read h5ad files and combine into a single BPCells matrix -----------
  message("Reading h5ad files")
  h5ad_paths  <- yaml::yaml.load_file(h5ads_yaml_f)  # list: batch_name -> path

  mat_list <- lapply(h5ad_paths, function(h5ad_f) {
    message(" reading ", basename(h5ad_f))
    BPCells::open_matrix_anndata_hdf5(h5ad_f)
  })
  counts_mat <- do.call(cbind, mat_list)
  rm(mat_list)

  # Normalise rowname separator and subset to cells in all_meta
  rownames(counts_mat) <- gsub("_", "-", rownames(counts_mat))

  shared_cells <- intersect(all_meta$cell_id, colnames(counts_mat))
  assert_that(length(shared_cells) > 0, msg = "No cell IDs overlap between integrated_dt and h5ad files")
  if (length(shared_cells) < nrow(all_meta))
    message(" WARNING: ", nrow(all_meta) - length(shared_cells),
            " cells in integrated_dt not found in h5ad files — they will be dropped")

  all_meta   <- all_meta[cell_id %in% shared_cells]
  counts_mat <- counts_mat[, all_meta$cell_id]

  # ---- Save cell metadata --------------------------------------------------
  out_fs <- .create_out_f_names(data_dir, app_tag, date_stamp)

  message("Saving cell metadata")
  if (!("n_cells" %in% colnames(all_meta))) all_meta[, n_cells := 1L]
  setnames(all_meta, c("umap_1", "umap_2"), c("UMAP_1", "UMAP_2"))
  fwrite(all_meta, out_fs["out_cell_meta_f"])

  # ---- Copy annotation.csv if supplied ------------------------------------
  if (nchar(annotation_csv_f) > 0 && file.exists(annotation_csv_f)) {
    annot <- fread(annotation_csv_f)
    assert_that(all(c("cluster", "cluster_name") %in% colnames(annot)),
      msg = "annotation_csv must contain columns 'cluster' and 'cluster_name'")
    annot[, cluster := as.character(cluster)]
    data_clusters  <- as.character(unique(all_meta$cluster))
    annot_clusters <- as.character(annot$cluster)
    missing_in_annot <- setdiff(data_clusters, annot_clusters)
    missing_in_data  <- setdiff(annot_clusters, data_clusters)
    if (length(missing_in_annot) > 0)
      message(" WARNING: clusters in data not in annotation_csv: ",
              paste(missing_in_annot, collapse = ", "))
    if (length(missing_in_data) > 0)
      message(" WARNING: clusters in annotation_csv not in data: ",
              paste(missing_in_data, collapse = ", "))
    fwrite(annot, file.path(data_dir, "annotation.csv"))
    message(" copied annotation.csv")
  }

  # ---- Normalise counts and write BPCells matrix (downsampled cells only) --
  message("Normalising and saving downsampled counts as BPCells")
  keep_mat    <- counts_mat[, all_meta$cell_id[all_meta$keep_cell]]
  norm_counts <- BPCells::multiply_cols(keep_mat, 1e4 / BPCells::colSums(keep_mat)) %>%
    log1p()

  .write_bpcells(norm_counts, out_fs["out_count_h5_f"])
  rm(keep_mat, norm_counts); gc(full = TRUE)

  # ---- Sample pseudobulks -------------------------------------------------
  message("Creating sample pseudobulks")
  sample_pb_sum    <- BPCells::pseudobulk_matrix(counts_mat, all_meta$sample_id, method = "sum")
  sample_pb_counts <- as(sample_pb_sum, "IterableMatrix") %>%
    BPCells::multiply_cols(1e6 / colSums(sample_pb_sum))
  sample_pb_meta   <- unique(all_meta[, keep_meta_cols, with = FALSE], by = "sample_id")
  sample_pb_meta[, pb_sample_id := sample_id]
  sample_pb_meta[, n_cells := all_meta[, .N, by = sample_id][match(sample_pb_meta$sample_id, sample_id), N]]

  .write_bpcells(sample_pb_counts, out_fs["out_sample_pb_h5_f"])
  fwrite(sample_pb_meta, out_fs["out_sample_meta_f"])

  # ---- Cluster pseudobulks ------------------------------------------------
  message("Creating cluster pseudobulks")
  pb_group          <- paste0(all_meta$sample_id, ".", all_meta$cluster)
  cluster_pb_sum    <- BPCells::pseudobulk_matrix(counts_mat, pb_group, method = "sum")
  cluster_pb_counts <- as(cluster_pb_sum, "IterableMatrix") %>%
    BPCells::multiply_cols(1e6 / colSums(cluster_pb_sum))

  # Build per-pseudobulk metadata from sample + cluster info
  pb_ids          <- colnames(cluster_pb_counts)
  cluster_pb_meta <- data.table(
    pb_sample_id = gsub("\\.", "_", pb_ids),
    sample_id    = sub("\\..*$", "", pb_ids),
    cluster      = sub("^[^.]+\\.", "", pb_ids)
  )
  cluster_pb_meta <- merge(
    cluster_pb_meta,
    sample_meta[, keep_meta_cols, with = FALSE],
    by = "sample_id", all.x = TRUE
  )
  pb_cell_counts  <- all_meta[, .N, by = .(pb_sample_id = gsub("\\.", "_", paste0(sample_id, ".", cluster)))]
  cluster_pb_meta <- merge(cluster_pb_meta, pb_cell_counts, by = "pb_sample_id", all.x = TRUE)
  setnames(cluster_pb_meta, "N", "n_cells")
  colnames(cluster_pb_counts) <- gsub("\\.", "_", colnames(cluster_pb_counts))

  .write_bpcells(cluster_pb_counts, out_fs["out_cluster_pb_h5_f"])
  fwrite(cluster_pb_meta, out_fs["out_cluster_meta_f"])

  # ---- Gene (row) metadata -----------------------------------------------
  message("Saving gene metadata")
  rowd <- data.table(
    gene_id = rownames(counts_mat),
    symbol  = strex::str_before_first(rownames(counts_mat), pattern = gene_re),
    ensembl = strex::str_after_last(rownames(counts_mat),  pattern = "-"),
    index   = seq_len(nrow(counts_mat))
  )
  fwrite(rowd, out_fs["out_row_indices_f"])

  # ---- TF annotation -------------------------------------------------------
  message("Loading transcription factor annotations")
  genesets_dir <- file.path(scprocess_dir, "resources/shiny/extdata/genesets")
  tfs_f        <- file.path(genesets_dir, paste0("transcription_factors_", species, ".txt.gz"))
  if (!file.exists(tfs_f)) {
    message(" WARNING: TF file not found at ", tfs_f, " — is_tf will be FALSE for all genes")
    tf_symbols <- character(0)
  } else {
    tfs_dt    <- fread(tfs_f)
    tf_symbols <- rowd[ensembl %in% tfs_dt$ensembl, symbol]
  }

  # ---- Markers -------------------------------------------------------------
  message("Processing cluster markers")
  markers <- fread(mkrs_f)
  assert_that(all(c("symbol", "cluster", "FDR", "logFC", "logcpm.sel") %in% colnames(markers)),
    msg = "markers file missing expected columns (symbol, cluster, FDR, logFC, logcpm.sel)")
  markers[, is.tf := symbol %in% tf_symbols]
  markers <- markers[, .(symbol, cluster, FDR, log2fc = logFC, CPM = round(exp(logcpm.sel) - 10), is.tf)] %>%
    setorder(cluster, FDR)
  fwrite(markers, out_fs["out_cluster_markers_f"])

  # ---- HVGs ----------------------------------------------------------------
  message("Processing HVGs")
  pb_hvgs <- fread(pb_hvgs_f)
  assert_that("vst_var" %in% colnames(pb_hvgs),
    msg = "pb_hvgs file must contain 'vst_var' column")
  pb_hvgs[, is.tf := symbol %in% tf_symbols]
  pb_hvgs <- pb_hvgs[, .(symbol, vst_var, is.tf)] %>% setorder(-vst_var)
  fwrite(pb_hvgs, out_fs["out_pb_hvgs_f"])

  # ---- GSEA ----------------------------------------------------------------
  message("Processing GSEA results")
  gsea_res <- list(
    go_bp = fread(fgsea_bp_f),
    go_cc = fread(fgsea_cc_f),
    go_mf = fread(fgsea_mf_f)
  )

  top_paths <- lapply(names(gsea_res), function(go_cat) {
    gsea_dt   <- gsea_res[[go_cat]]
    assert_that(all(gsea_dt$path_set == go_cat),
      msg = paste("path_set column in", go_cat, "file does not match expected value"))
    gsea_dt[
      main_path == TRUE
    ][,
      min_p  := min(padj, na.rm = TRUE), by = pathway
    ][,
      signif := ifelse(padj < 0.05, "significant", "not")
    ][
      min_p < 0.1 & NES > 0
    ] %>%
      setorder(cluster, padj) %>%
      .[, p_rank := seq_len(.N), by = cluster] %>%
      .[, .(pathway, padj, NES, size, p_rank, signif, cluster, go_category = path_set)]
  }) %>% rbindlist()

  fwrite(top_paths, out_fs["out_fgsea_f"])

  # ---- GO term gene lists --------------------------------------------------
  message("Building GO term gene lists")
  go_terms_f <- file.path(genesets_dir, paste0("genes_go_pathways_", species, ".txt.gz"))
  if (!file.exists(go_terms_f)) {
    message(" WARNING: GO terms file not found at ", go_terms_f, " — geneset exploration will be empty")
    go_genes <- data.table(pathway = character(), pathway_nice = character(),
                           path_short = character(), go_category = character(),
                           genes = character())
  } else {
    signif_markers <- markers[
      FDR <= 0.05 & CPM >= 10
    ][,
      max_fc := max(abs(log2fc)), by = symbol
    ][,
      .(symbol, max_fc)
    ] %>% unique()

    go_genes <- fread(go_terms_f) %>%
      setnames("gene", "symbol") %>%
      .[symbol %in% signif_markers$symbol] %>%
      merge(signif_markers, by = "symbol", allow.cartesian = TRUE) %>%
      setorder(-max_fc) %>%
      .[, row.idx := seq_len(.N), by = pathway] %>%
      .[row.idx <= 100] %>%
      .[, .(genes = paste(symbol, collapse = " ")),
        by = .(pathway, pathway_nice, path_short, go_category)]
  }
  fwrite(go_genes, out_fs["out_go_terms_f"])

  # ---- Write shinyconfig.yaml ----------------------------------------------
  message("Writing shinyconfig.yaml")
  shinyconfig <- list(
    date_stamp = date_stamp,
    data_dir   = "data",
    build = list(
      sample_col   = "sample_id",
      logo_f       = NULL,
      gsets_f      = NULL,
      include_paga = NULL
    ),
    app = list(
      app_title    = app_title,
      email        = email,
      keyword      = keyword,
      default_gene = if (nchar(default_gene) > 0) default_gene else rowd$symbol[1]
    ),
    metadata = list(
      vars      = as.list(setNames(metadata_vars, var_names)),
      var_names = as.list(var_names),
      var_combns = if (length(var_combns) > 0) var_combns else NULL
    )
  )
  yaml::write_yaml(shinyconfig, file.path(deploy_dir, "shinyconfig.yaml"))

  message("Done! Shiny app written to: ", deploy_dir)
}


# HELPER FUNCTIONS -------------------------------------------------------------

.create_out_f_names <- function(data_dir, app_tag, date_stamp) {
  fs <- c(
    out_count_h5_f        = paste0(app_tag, "-shiny_norm_count-",     date_stamp),
    out_sample_pb_h5_f    = paste0(app_tag, "-shiny_sample_pb_count-", date_stamp),
    out_cluster_pb_h5_f   = paste0(app_tag, "-shiny_cluster_pb_count-", date_stamp),
    out_row_indices_f     = paste0(app_tag, "-shiny_row_indices-",     date_stamp, ".txt.gz"),
    out_cluster_markers_f = paste0(app_tag, "-shiny_markers-",         date_stamp, ".txt.gz"),
    out_pb_hvgs_f         = paste0(app_tag, "-shiny_pb_hvgs-",         date_stamp, ".txt.gz"),
    out_fgsea_f           = paste0(app_tag, "-shiny_fgsea_res-",       date_stamp, ".txt.gz"),
    out_go_terms_f        = paste0(app_tag, "-shiny_go_terms-",        date_stamp, ".txt.gz"),
    out_sample_meta_f     = paste0(app_tag, "-shiny_sample_meta-",     date_stamp, ".txt.gz"),
    out_cluster_meta_f    = paste0(app_tag, "-shiny_cluster_meta-",    date_stamp, ".txt.gz"),
    out_cell_meta_f       = paste0(app_tag, "-shiny_cell_meta-",       date_stamp, ".txt.gz")
  )
  data_dir <- gsub("/$", "", data_dir)
  sapply(fs, function(f) file.path(data_dir, f))
}


.subsample_umap <- function(umap_dt, to_keep = 5e4) {
  umap_dt <- umap_dt[, density := .get_density(UMAP_1, UMAP_2, n = 500)] %>%
    .[, inv_dens := 1 / density] %>%
    .[, p_keep   := inv_dens / sum(inv_dens)]

  set.seed(20230308)
  if (to_keep >= nrow(umap_dt)) {
    keep_cells <- umap_dt$cell_id
  } else {
    keep_cells <- sample(umap_dt$cell_id, prob = umap_dt$p_keep, size = to_keep, replace = FALSE)
  }
  return(keep_cells)
}


.get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix   <- findInterval(x, dens$x)
  iy   <- findInterval(y, dens$y)
  ii   <- cbind(ix, iy)
  dens$z[ii]
}


.write_bpcells <- function(mx, dir_name) {
  if (dir.exists(dir_name)) unlink(dir_name, recursive = TRUE)
  if (!is(mx, "IterableMatrix")) mx <- as(mx, "IterableMatrix")
  BPCells::transpose_storage_order(mx, outdir = dir_name)
}
