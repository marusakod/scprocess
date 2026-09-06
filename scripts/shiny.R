suppressPackageStartupMessages({
  library("assertthat")
  library("data.table")
  library("jsonlite")
  library("magrittr")
  library("yaml")
  library("BPCells")
  library("MASS")
  library("strex")
})


.validate_shiny_annotation <- function(annotation_csv_f, data_clusters) {
  if (!file.exists(annotation_csv_f))
    stop("annotation_csv does not exist: ", annotation_csv_f)

  annot <- fread(annotation_csv_f)
  required_cols <- c("cluster", "cluster_name")
  missing_cols  <- setdiff(required_cols, colnames(annot))
  if (length(missing_cols) > 0L)
    stop("annotation_csv must contain columns 'cluster' and 'cluster_name'; missing: ",
         paste(missing_cols, collapse = ", "))

  annot[, cluster      := trimws(as.character(cluster))]
  annot[, cluster_name := trimws(as.character(cluster_name))]

  blank_clusters <- is.na(annot$cluster) | annot$cluster == ""
  if (any(blank_clusters))
    stop("annotation_csv contains missing or blank cluster values at rows: ",
         paste(which(blank_clusters), collapse = ", "))

  blank_names <- is.na(annot$cluster_name) | annot$cluster_name == ""
  if (any(blank_names))
    stop("annotation_csv contains missing or blank cluster_name values at rows: ",
         paste(which(blank_names), collapse = ", "))

  duplicate_clusters <- unique(annot$cluster[duplicated(annot$cluster)])
  if (length(duplicate_clusters) > 0L)
    stop("annotation_csv contains duplicated cluster values: ",
         paste(duplicate_clusters, collapse = ", "))

  duplicate_names <- unique(annot$cluster_name[duplicated(annot$cluster_name)])
  if (length(duplicate_names) > 0L)
    stop("annotation_csv contains duplicated cluster_name values: ",
         paste(duplicate_names, collapse = ", "))

  data_clusters  <- unique(trimws(as.character(data_clusters)))
  annot_clusters <- annot$cluster
  only_in_annot  <- sort(setdiff(annot_clusters, data_clusters))
  only_in_data   <- sort(setdiff(data_clusters, annot_clusters))

  if (length(only_in_annot) > 0L || length(only_in_data) > 0L) {
    mismatch_details <- c(
      if (length(only_in_annot) > 0L)
        paste0("Present only in annotation_csv: ", paste(only_in_annot, collapse = ", ")),
      if (length(only_in_data) > 0L)
        paste0("Present only in analysis data: ", paste(only_in_data, collapse = ", "))
    )
    stop(
      "annotation_csv does not match the clusters in the current analysis.\n",
      paste(mismatch_details, collapse = "\n"),
      "\nThe annotation may belong to a different clustering resolution or analysis run. ",
      "Update shiny.annotation_csv before rebuilding the app."
    )
  }

  annot
}


#' Build Shiny app from scprocess outputs
#'
#' @param integration_f  Path to integrated_dt_{FULL_TAG}_{DATE_STAMP}.csv.gz
#' @param h5ads_yaml_f   Path to h5ads_clean_paths_{FULL_TAG}_{DATE_STAMP}.yaml (maps batch -> h5ad path)
#' @param sample_meta_f  Path to sample metadata CSV (config['project']['sample_metadata'])
#' @param mkrs_f         Path to pb_marker_genes_{FULL_TAG}_{res}_{DATE_STAMP}.csv.gz
#' @param pb_hvgs_f      Path to pb_hvgs_{FULL_TAG}_{res}_{DATE_STAMP}.csv.gz
#' @param fgsea_bp_f     Path to fgsea go_bp CSV.GZ (pass "" to skip GSEA)
#' @param fgsea_cc_f     Path to fgsea go_cc CSV.GZ (pass "" to skip GSEA)
#' @param fgsea_mf_f     Path to fgsea go_mf CSV.GZ (pass "" to skip GSEA)
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
#' @param metadata_combns     JSON-encoded list of variable combination pairs (default "[]")
#' @param home_md_f      Optional path to a Markdown file used as the landing page (overrides placeholder)
#' @param annotation_csv_f Optional path to a CSV with columns cluster, cluster_name, colour defining
#'                        display names, colour, and order for clusters
#' @param cluster_palette  Optional palette name applied to cluster colours (when no annotation_csv).
#'                        Any name in VALID_PALETTE_NAMES (nice_cols, MetBrewer, RColorBrewer, ggsci).
#' @param metadata_palettes JSON-encoded object mapping metadata variable names to palette specs.
#'                        Each value may be: a string (palette name), an array (explicit hex colours),
#'                        or an object with optional keys 'palette', 'colours', and 'values'.
#' @param n_cores        Number of data.table threads
make_shiny_app_scprocess <- function(
  integration_f, h5ads_yaml_f, sample_meta_f,
  mkrs_f, pb_hvgs_f, fgsea_bp_f = "", fgsea_cc_f = "", fgsea_mf_f = "",
  deploy_dir, scprocess_dir,
  app_tag, date_stamp, mkr_sel_res,
  ref_txome, metadata_vars,
  app_title = app_tag,
  email     = "",
  keyword   = "cells",
  default_gene = "",
  n_keep    = 30000L,
  var_names = metadata_vars,
  metadata_combns = "[]",
  home_md_f = "",
  annotation_csv_f = "",
  cluster_palette = "",
  metadata_palettes = "{}",
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
  metadata_combns    <- jsonlite::fromJSON(metadata_combns)
  if (is.matrix(metadata_combns))
    metadata_combns  <- lapply(seq_len(nrow(metadata_combns)), function(i) metadata_combns[i, ])

  assert_that(length(metadata_vars) == length(var_names),
    msg = "metadata_vars and var_names must have the same length")

  n_keep        <- as.integer(n_keep)
  cluster_col   <- paste0("RNA_snn_res.", mkr_sel_res)

  # ---- Parse and validate palette arguments ----------------------------------
  # Load valid palette names from the shared resource file (mirrors scprocess_utils.py)
  .pal_file       <- file.path(scprocess_dir, "resources", "valid_palettes.json")
  .pal_groups     <- jsonlite::read_json(.pal_file)
  .valid_pal_names <- unlist(.pal_groups[!startsWith(names(.pal_groups), "_")], use.names = FALSE)

  .check_palette_name <- function(name, context) {
    if (!name %in% .valid_pal_names)
      stop("Unknown palette '", name, "' in ", context, ". ",
           "See the scprocess documentation for valid palette names.")
  }

  # Validate cluster_palette
  if (nchar(cluster_palette) > 0)
    .check_palette_name(cluster_palette, "cluster_palette")

  # Parse metadata_palettes from JSON and normalise into shinyconfig format.
  # Input shapes (per variable):
  #   "VanGogh1"             -> {palette: "VanGogh1"}
  #   ["#FF0000", "#0000FF"] -> {colours: [...]}
  #   {palette: "Klimt", values: [...]} -> {palette: "Klimt", values: [...]}
  #   {colours: [...], values: [...]}   -> {colours: [...], values: [...]}
  raw_pals <- jsonlite::fromJSON(metadata_palettes, simplifyVector = FALSE)

  normalised_pals <- lapply(names(raw_pals), function(v) {
    spec <- raw_pals[[ v ]]
    if (is.character(spec) && length(spec) == 1L) {
      # scalar string: palette name
      .check_palette_name(spec, paste0("metadata_palettes$", v))
      list(palette = spec)
    } else if (is.list(spec) && !is.null(names(spec))) {
      # named object: {palette:, colours:, values:}
      if (!is.null(spec$palette))
        .check_palette_name(spec$palette, paste0("metadata_palettes$", v, "$palette"))
      spec
    } else if (is.list(spec) && is.null(names(spec))) {
      # unnamed list from JSON array: explicit colour list
      list(colours = spec)
    } else {
      stop("metadata_palettes$", v, ": unexpected format. ",
           "Expected a palette name (string), colour list (array), ",
           "or object with 'palette'/'colours'/'values' keys.")
    }
  }) %>% setNames(names(raw_pals))

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

  # Copy valid_palettes.json so the deployed server does not need scprocess_dir
  file.copy(
    file.path(scprocess_dir, "resources", "valid_palettes.json"),
    file.path(data_dir, "valid_palettes.json"),
    overwrite = TRUE
  )
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
    merge(clusts, by = "cell_id")

  # Preserve all clusters from the integration, including rare clusters for
  # which marker-gene testing did not produce results.
  analysis_clusters <- sort(unique(as.character(all_meta$cluster)))
  all_meta$cluster  <- factor(all_meta$cluster, levels = analysis_clusters)

  # Validate cluster annotations against the full analysis before UMAP
  # subsampling. Annotation row order controls display order.
  annot <- NULL
  if (nchar(annotation_csv_f) > 0L) {
    annot <- .validate_shiny_annotation(
      annotation_csv_f,
      analysis_clusters
    )
  }

  # ---- Downsample UMAP -----------------------------------------------------
  message("Downsampling UMAP to ", n_keep, " cells")
  keep_cells <- .subsample_umap(
    all_meta[, .(cell_id, UMAP_1 = umap_1, UMAP_2 = umap_2, cluster)],
    to_keep = n_keep,
    stratify_by = "cluster"
  )
  all_meta[, keep_cell := cell_id %in% keep_cells]

  # Remove stale shiny data files from prior builds (e.g. .txt.gz from older
  # versions) only after inputs, including annotations, have been validated.
  stale <- list.files(data_dir, pattern = paste0(app_tag, "-shiny_"),
                      full.names = TRUE)
  if (length(stale) > 0L) {
    unlink(stale, recursive = TRUE)
    message("Removed ", length(stale), " stale data files from prior build")
  }
  annotation_dest <- file.path(data_dir, "annotation.csv")
  if (is.null(annot) && file.exists(annotation_dest))
    unlink(annotation_dest)

  # ---- Read h5ad files and combine into a single BPCells matrix -----------
  message("Reading h5ad files")
  h5ad_paths  <- yaml::yaml.load_file(h5ads_yaml_f)  # list: batch_name -> path (or {path, project_id})

  # Helper: load one h5ad; if entry is a join dict {path, project_id}, log the project
  # Cell IDs are kept as-is (the join pipeline already ensures uniqueness across projects)
  .load_one_h5ad_for_shiny <- function(h5ad_entry) {
    if (is.character(h5ad_entry)) {
      h5ad_f  <- h5ad_entry
    } else {
      h5ad_f  <- h5ad_entry[["path"]]
    }
    # skip empty files (excluded batches produce zero-byte h5ad placeholders)
    if (file.size(h5ad_f) == 0) {
      message(" skipping empty ", basename(h5ad_f))
      return(NULL)
    }
    message(" reading ", basename(h5ad_f))
    BPCells::open_matrix_anndata_hdf5(h5ad_f)
  }

  mat_list   <- Filter(Negate(is.null), lapply(h5ad_paths, .load_one_h5ad_for_shiny))
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

  # ---- Pre-compute cluster centroids from all cells (before downsampling) --
  message("Computing cluster centroids")
  centroids_dt <- all_meta[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = cluster] %>%
    setorder(cluster)
  fwrite(centroids_dt, out_fs["out_centroids_f"])

  # ---- Pre-compute repelled label positions from full UMAP density ----------
  message("Computing repelled label positions")
  repel_pos_dt <- .compute_repel_positions(all_meta, centroids_dt)
  fwrite(repel_pos_dt, out_fs["out_repel_pos_f"])

  # ---- Copy annotation.csv if supplied ------------------------------------
  if (!is.null(annot)) {
    fwrite(annot, annotation_dest)
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
  if ("gene_type" %in% colnames(markers)) {
    markers <- markers[, .(symbol, cluster, FDR, log2fc = logFC, CPM = round(exp(logcpm.sel) - 10), gene_type, is.tf)]
  } else {
    markers <- markers[, .(symbol, cluster, FDR, log2fc = logFC, CPM = round(exp(logcpm.sel) - 10), is.tf)]
  }
  markers <- markers %>% setorder(cluster, FDR)
  fwrite(markers, out_fs["out_cluster_markers_f"])

  # ---- HVGs ----------------------------------------------------------------
  message("Processing HVGs")
  pb_hvgs <- fread(pb_hvgs_f)
  variability_col <- intersect(c("logcpm_var", "vst_var"), colnames(pb_hvgs))[1]
  assert_that(!is.na(variability_col),
    msg = "pb_hvgs file must contain 'logcpm_var' or 'vst_var'")
  pb_hvgs[, is.tf := symbol %in% tf_symbols]
  pb_hvgs <- pb_hvgs[, .(symbol, variability = get(variability_col), is.tf)] %>%
    setorder(-variability)
  fwrite(pb_hvgs, out_fs["out_pb_hvgs_f"])

  # ---- GSEA ----------------------------------------------------------------
  do_gsea <- nchar(fgsea_bp_f) > 0 && nchar(fgsea_cc_f) > 0 && nchar(fgsea_mf_f) > 0

  if (do_gsea) {
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

  } else {
    message("Skipping GSEA (no fgsea files provided)")
    top_paths <- data.table(
      pathway     = character(), padj    = numeric(),
      NES         = numeric(),   size    = integer(),
      p_rank      = integer(),   signif  = character(),
      cluster     = character(), go_category = character()
    )
  }
  fwrite(top_paths, out_fs["out_fgsea_f"])

  # ---- GO term gene lists --------------------------------------------------
  go_terms_f <- file.path(genesets_dir, paste0("genes_go_pathways_", species, ".txt.gz"))
  if (!do_gsea) {
    message("Skipping GO term gene lists (no GSEA results)")
    go_genes <- data.table(
      pathway      = character(), pathway_nice = character(),
      path_short   = character(), go_category  = character(),
      genes        = character()
    )
  } else if (!file.exists(go_terms_f)) {
    message(" WARNING: GO terms file not found at ", go_terms_f, " — geneset exploration will be empty")
    go_genes <- data.table(
      pathway      = character(), pathway_nice = character(),
      path_short   = character(), go_category  = character(),
      genes        = character()
    )
  } else {
    message("Building GO term gene lists")
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
  meta_section <- list(
    vars      = as.list(metadata_vars),
    var_names = as.list(var_names),
    metadata_combns = if (length(metadata_combns) > 0) metadata_combns else NULL
  )
  if (nchar(cluster_palette) > 0)
    meta_section$cluster_palette <- cluster_palette
  if (length(normalised_pals) > 0)
    meta_section$palettes <- normalised_pals

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
    metadata = meta_section
  )
  yaml::write_yaml(shinyconfig, file.path(deploy_dir, "shinyconfig.yaml"))

  message("Done! Shiny app written to: ", deploy_dir)
}


#' Configure an existing scprocess Shiny deployment
#'
#' This is the lightweight counterpart to make_shiny_app_scprocess(). It copies
#' application resources and presentation-only inputs, validates annotations
#' and the default gene against already-built data, and rewrites
#' shinyconfig.yaml without rebuilding BPCells matrices or pseudobulks.
configure_shiny_app_scprocess <- function(
  deploy_dir, scprocess_dir,
  app_tag, date_stamp, metadata_vars,
  app_title = app_tag,
  email = "",
  keyword = "cells",
  default_gene = "",
  var_names = metadata_vars,
  metadata_combns = "[]",
  home_md_f = "",
  annotation_csv_f = "",
  cluster_palette = "",
  metadata_palettes = "{}"
) {
  deploy_dir <- gsub("/$", "", deploy_dir)
  data_dir   <- file.path(deploy_dir, "data")
  assert_that(dir.exists(data_dir),
    msg = paste("Shiny data directory does not exist:", data_dir))

  metadata_vars   <- strsplit(metadata_vars, ",")[[1]] %>% trimws()
  var_names       <- strsplit(var_names, ",")[[1]] %>% trimws()
  metadata_combns <- jsonlite::fromJSON(metadata_combns)
  if (is.matrix(metadata_combns))
    metadata_combns <- lapply(seq_len(nrow(metadata_combns)), function(i) metadata_combns[i, ])

  assert_that(length(metadata_vars) == length(var_names),
    msg = "metadata_vars and var_names must have the same length")

  pal_file       <- file.path(scprocess_dir, "resources", "valid_palettes.json")
  pal_groups     <- jsonlite::read_json(pal_file)
  valid_pal_names <- unlist(pal_groups[!startsWith(names(pal_groups), "_")], use.names = FALSE)

  check_palette_name <- function(name, context) {
    if (!name %in% valid_pal_names)
      stop("Unknown palette '", name, "' in ", context, ". ",
           "See the scprocess documentation for valid palette names.")
  }

  if (nchar(cluster_palette) > 0)
    check_palette_name(cluster_palette, "cluster_palette")

  raw_pals <- jsonlite::fromJSON(metadata_palettes, simplifyVector = FALSE)
  normalised_pals <- lapply(names(raw_pals), function(v) {
    spec <- raw_pals[[v]]
    if (is.character(spec) && length(spec) == 1L) {
      check_palette_name(spec, paste0("metadata_palettes$", v))
      list(palette = spec)
    } else if (is.list(spec) && !is.null(names(spec))) {
      if (!is.null(spec$palette))
        check_palette_name(spec$palette, paste0("metadata_palettes$", v, "$palette"))
      spec
    } else if (is.list(spec) && is.null(names(spec))) {
      list(colours = spec)
    } else {
      stop("metadata_palettes$", v, ": unexpected format. ",
           "Expected a palette name (string), colour list (array), ",
           "or object with 'palette'/'colours'/'values' keys.")
    }
  }) %>% setNames(names(raw_pals))

  out_fs       <- .create_out_f_names(data_dir, app_tag, date_stamp)
  cell_meta_f  <- unname(out_fs["out_cell_meta_f"])
  row_indices_f <- unname(out_fs["out_row_indices_f"])
  assert_that(file.exists(cell_meta_f),
    msg = paste("Built Shiny cell metadata does not exist:", cell_meta_f))
  assert_that(file.exists(row_indices_f),
    msg = paste("Built Shiny gene metadata does not exist:", row_indices_f))

  cell_meta <- fread(cell_meta_f, select = "cluster")
  analysis_clusters <- sort(unique(as.character(cell_meta$cluster)))

  annotation_dest <- file.path(data_dir, "annotation.csv")
  if (nchar(annotation_csv_f) > 0L) {
    annot <- .validate_shiny_annotation(annotation_csv_f, analysis_clusters)
    fwrite(annot, annotation_dest)
    message("Copied annotation.csv")
  } else if (file.exists(annotation_dest)) {
    unlink(annotation_dest)
    message("Removed annotation.csv")
  }

  rowd <- fread(row_indices_f, select = "symbol")
  if (nchar(default_gene) == 0L) {
    default_gene <- rowd$symbol[1]
  } else {
    assert_that(default_gene %in% rowd$symbol,
      msg = paste0("default_gene '", default_gene,
                   "' is not present in the built Shiny expression matrix"))
  }

  app_src <- file.path(scprocess_dir, "resources", "shiny")
  assert_that(dir.exists(app_src),
    msg = paste("resources/shiny not found at", app_src))
  for (f in c("app.R", "constants.R"))
    file.copy(file.path(app_src, f), deploy_dir, overwrite = TRUE)
  file.copy(pal_file, file.path(data_dir, "valid_palettes.json"), overwrite = TRUE)
  for (subdir in c("utils", "modules")) {
    src  <- file.path(app_src, subdir)
    dest <- file.path(deploy_dir, subdir)
    dir.create(dest, showWarnings = FALSE, recursive = TRUE)
    file.copy(list.files(src, full.names = TRUE), dest, overwrite = TRUE)
  }

  home_md_dest <- file.path(data_dir, "home.md")
  if (nchar(home_md_f) > 0L && file.exists(home_md_f)) {
    file.copy(home_md_f, home_md_dest, overwrite = TRUE)
    message("Copied user home.md")
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
    message("Created placeholder home.md")
  }

  meta_section <- list(
    vars = as.list(metadata_vars),
    var_names = as.list(var_names),
    metadata_combns = if (length(metadata_combns) > 0) metadata_combns else NULL
  )
  if (nchar(cluster_palette) > 0)
    meta_section$cluster_palette <- cluster_palette
  if (length(normalised_pals) > 0)
    meta_section$palettes <- normalised_pals

  shinyconfig <- list(
    date_stamp = date_stamp,
    data_dir = "data",
    build = list(
      sample_col = "sample_id",
      logo_f = NULL,
      gsets_f = NULL,
      include_paga = NULL
    ),
    app = list(
      app_title = app_title,
      email = email,
      keyword = keyword,
      default_gene = default_gene
    ),
    metadata = meta_section
  )
  yaml::write_yaml(shinyconfig, file.path(deploy_dir, "shinyconfig.yaml"))
  message("Configured Shiny app at: ", deploy_dir)
}


# HELPER FUNCTIONS -------------------------------------------------------------

.create_out_f_names <- function(data_dir, app_tag, date_stamp) {
  fs <- c(
    out_count_h5_f        = paste0(app_tag, "-shiny_norm_count-",     date_stamp),
    out_sample_pb_h5_f    = paste0(app_tag, "-shiny_sample_pb_count-", date_stamp),
    out_cluster_pb_h5_f   = paste0(app_tag, "-shiny_cluster_pb_count-", date_stamp),
    out_row_indices_f     = paste0(app_tag, "-shiny_row_indices-",     date_stamp, ".csv.gz"),
    out_cluster_markers_f = paste0(app_tag, "-shiny_markers-",         date_stamp, ".csv.gz"),
    out_pb_hvgs_f         = paste0(app_tag, "-shiny_pb_hvgs-",         date_stamp, ".csv.gz"),
    out_fgsea_f           = paste0(app_tag, "-shiny_fgsea_res-",       date_stamp, ".csv.gz"),
    out_go_terms_f        = paste0(app_tag, "-shiny_go_terms-",        date_stamp, ".csv.gz"),
    out_sample_meta_f     = paste0(app_tag, "-shiny_sample_meta-",     date_stamp, ".csv.gz"),
    out_cluster_meta_f    = paste0(app_tag, "-shiny_cluster_meta-",    date_stamp, ".csv.gz"),
    out_cell_meta_f       = paste0(app_tag, "-shiny_cell_meta-",       date_stamp, ".csv.gz"),
    out_centroids_f       = paste0(app_tag, "-shiny_centroids-",       date_stamp, ".csv.gz"),
    out_repel_pos_f       = paste0(app_tag, "-shiny_repel_pos-",       date_stamp, ".csv.gz")
  )
  data_dir <- gsub("/$", "", data_dir)
  sapply(fs, function(f) file.path(data_dir, f))
}


.subsample_umap <- function(umap_dt, to_keep = 5e4, stratify_by = NULL) {
  umap_dt <- umap_dt[, density := .get_density(UMAP_1, UMAP_2, n = 500)] %>%
    .[, inv_dens := 1 / density] %>%
    .[, p_keep   := inv_dens / sum(inv_dens)]

  set.seed(20230308)
  if (to_keep >= nrow(umap_dt)) {
    keep_cells <- umap_dt$cell_id
  } else if (!is.null(stratify_by)) {
    assert_that(stratify_by %in% colnames(umap_dt),
      msg = paste("stratify_by column not found:", stratify_by))
    group_idx <- split(seq_len(nrow(umap_dt)), umap_dt[[stratify_by]])
    assert_that(to_keep >= length(group_idx),
      msg = paste0("Cannot retain every ", stratify_by, " group when to_keep (",
                   to_keep, ") is smaller than the number of groups (",
                   length(group_idx), ")"))

    mandatory_idx <- vapply(group_idx, function(idx) {
      idx[sample.int(length(idx), size = 1L, prob = umap_dt$p_keep[idx])]
    }, integer(1))
    remaining_idx <- setdiff(seq_len(nrow(umap_dt)), mandatory_idx)
    n_remaining   <- to_keep - length(mandatory_idx)
    sampled_idx   <- if (n_remaining > 0L) {
      sample(remaining_idx, size = n_remaining,
             prob = umap_dt$p_keep[remaining_idx], replace = FALSE)
    } else {
      integer(0)
    }
    keep_cells <- umap_dt$cell_id[c(mandatory_idx, sampled_idx)]
  } else {
    keep_cells <- sample(umap_dt$cell_id, prob = umap_dt$p_keep, size = to_keep, replace = FALSE)
  }
  return(keep_cells)
}


# Compute non-overlapping label positions in low-density UMAP space.
# Uses a grid density map to find positions whose label footprint is clear of
# cells, then allocates those positions globally so labels cannot overlap.
# Returns a data.table with columns: cluster, repel_1, repel_2.
.compute_repel_positions <- function(all_meta, centroids_dt,
                                     grid_res = 200L, max_radius_frac = 0.5,
                                     r_clear = 6L, r_label = r_clear,
                                     label_gap = 1L, padding_frac = 0.05) {
  stopifnot(grid_res > 2L * r_clear,
            r_clear >= 0L, r_label >= 0L, label_gap >= 0L,
            padding_frac >= 0)

  raw_u1_min <- min(all_meta$UMAP_1); raw_u1_max <- max(all_meta$UMAP_1)
  raw_u2_min <- min(all_meta$UMAP_2); raw_u2_max <- max(all_meta$UMAP_2)
  u1_span <- raw_u1_max - raw_u1_min
  u2_span <- raw_u2_max - raw_u2_min
  if (u1_span == 0 || u2_span == 0)
    stop("UMAP coordinates must span both dimensions to place cluster labels")

  # A small empty margin gives crowded layouts somewhere to place labels
  # without falling back to positions on top of cells.
  u1_min <- raw_u1_min - u1_span * padding_frac
  u1_max <- raw_u1_max + u1_span * padding_frac
  u2_min <- raw_u2_min - u2_span * padding_frac
  u2_max <- raw_u2_max + u2_span * padding_frac

  x_breaks  <- seq(u1_min, u1_max, length.out = grid_res + 1L)
  y_breaks  <- seq(u2_min, u2_max, length.out = grid_res + 1L)
  x_centers <- (x_breaks[-1L] + x_breaks[-(grid_res + 1L)]) / 2
  y_centers <- (y_breaks[-1L] + y_breaks[-(grid_res + 1L)]) / 2

  # Bin all cells (vectorised)
  gx  <- pmax(1L, pmin(grid_res, findInterval(all_meta$UMAP_1, x_breaks, rightmost.closed = TRUE)))
  gy  <- pmax(1L, pmin(grid_res, findInterval(all_meta$UMAP_2, y_breaks, rightmost.closed = TRUE)))
  idx <- gx + (gy - 1L) * grid_res
  density_grid <- matrix(tabulate(idx, nbins = grid_res^2L),
                         nrow = grid_res, ncol = grid_res)

  # Build clearance grid: separable 2D box sum over (2*r_clear+1)^2 neighbourhood.
  # clearance_grid[i,j] == 0  iff  ALL cells within r_clear steps are empty —
  # i.e. the whole label footprint would sit in white space.
  # Boundary cells get Inf so they are never chosen as label positions.
  k <- rep(1, 2L * r_clear + 1L)
  row_sums <- matrix(Inf, nrow = grid_res, ncol = grid_res)
  for (i in seq_len(grid_res)) {
    v <- stats::filter(density_grid[i, ], k, sides = 2)
    row_sums[i, !is.na(v)] <- v[!is.na(v)]
  }
  clearance_grid <- matrix(Inf, nrow = grid_res, ncol = grid_res)
  for (j in seq_len(grid_res)) {
    v <- stats::filter(row_sums[, j], k, sides = 2)
    clearance_grid[!is.na(v), j] <- v[!is.na(v)]
  }

  make_candidates <- function(i, max_r) {
    cx  <- centroids_dt$UMAP_1[i]
    cy  <- centroids_dt$UMAP_2[i]
    cgx <- max(1L, min(grid_res, findInterval(cx, x_breaks, rightmost.closed = TRUE)))
    cgy <- max(1L, min(grid_res, findInterval(cy, y_breaks, rightmost.closed = TRUE)))

    seq_r <- seq(-max_r, max_r)
    candidates <- CJ(dx = seq_r, dy = seq_r)
    candidates[, `:=`(
      gx = cgx + dx,
      gy = cgy + dy,
      r = pmax(abs(dx), abs(dy)),
      dist2 = dx^2 + dy^2
    )]
    candidates <- candidates[
      gx > r_clear & gx <= grid_res - r_clear &
      gy > r_clear & gy <= grid_res - r_clear
    ]
    candidates[, density := clearance_grid[cbind(gx, gy)]]
    setorder(candidates, density, r, dist2, gx, gy)
    candidates
  }

  local_max_r <- max(1L, as.integer(round(grid_res * max_radius_frac)))
  local_candidates <- lapply(seq_len(nrow(centroids_dt)), make_candidates,
                             max_r = local_max_r)

  # Place clusters with the fewest empty local candidates first. This reduces
  # order sensitivity when several nearby centroids compete for the same gap.
  constraints <- rbindlist(lapply(seq_along(local_candidates), function(i) {
    candidates <- local_candidates[[i]]
    empty <- candidates[density == 0]
    data.table(
      idx = i,
      n_empty = nrow(empty),
      nearest_empty = if (nrow(empty) > 0L) min(empty$r) else Inf
    )
  }))
  setorder(constraints, n_empty, -nearest_empty, idx)

  min_label_sep <- 2L * r_label + label_gap
  placed <- data.table(idx = integer(), gx = integer(), gy = integer())

  is_clear_of_labels <- function(candidates) {
    if (nrow(placed) == 0L) return(rep(TRUE, nrow(candidates)))
    is_clear <- rep(TRUE, nrow(candidates))
    for (k in seq_len(nrow(placed))) {
      is_clear <- is_clear &
        pmax(abs(candidates$gx - placed$gx[k]),
             abs(candidates$gy - placed$gy[k])) >= min_label_sep
    }
    is_clear
  }

  for (i in constraints$idx) {
    candidates <- local_candidates[[i]]
    available <- which(is_clear_of_labels(candidates))

    # A crowded local neighbourhood may have no collision-free candidate.
    # Search the complete usable grid before declaring the layout impossible.
    if (length(available) == 0L) {
      candidates <- make_candidates(i, grid_res - 1L)
      available <- which(is_clear_of_labels(candidates))
    }
    if (length(available) == 0L) {
      stop(
        "Unable to place non-overlapping UMAP label for cluster ",
        as.character(centroids_dt$cluster[i]),
        ". Reduce the label footprint or increase grid_res."
      )
    }

    chosen <- candidates[available[1L]]
    placed <- rbind(
      placed,
      data.table(idx = i, gx = chosen$gx, gy = chosen$gy)
    )
  }

  repel_dt <- merge(
    data.table(idx = seq_len(nrow(centroids_dt))),
    placed,
    by = "idx",
    sort = FALSE
  )
  setorder(repel_dt, idx)
  repel_dt[, .(
    cluster = centroids_dt$cluster[idx],
    repel_1 = x_centers[gx],
    repel_2 = y_centers[gy]
  )]
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
