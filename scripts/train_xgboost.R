# Train an XGBoost classifier to predict cell type labels from scRNA-seq data.
# Workflow: load reference dataset -> Harmony batch correction -> HVG extraction -> XGBoost training -> diagnostic UMAP plots.
# The trained model and allowed-cluster table are consumed by the label_celltypes step of scprocess.

suppressPackageStartupMessages({
  library('magrittr')
  library('data.table')
  library('forcats')
  library('stringr')
  library('assertthat')
  library('BiocParallel')
  library('Seurat')
  library('harmony')
  library('scran')
  library('uwot')
  library('future')
  library('ggplot.multistats')
  library('xgboost')
  library('ggplot2')
})

options( future.globals.maxSize = 2^35 )

nice_cols   = c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6",
  "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A",
  "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4",
  "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3",
  "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"
)

yaml_f = 'test_config.yaml'

# Main entry point. Reads a YAML config, orchestrates the full training workflow,
# and saves the XGBoost classifier plus diagnostic plots to the output directory.
make_classifier <- function(yaml_f) {

  yaml_ls = .get_yaml_ls(yaml_f)

  # unpack config parameters
  samples_f   = yaml_ls$samples_f
  annots_f    = yaml_ls$annots_f
  xgboost_dir = yaml_ls$xgboost_dir
  ref_tag     = yaml_ls$ref_tag
  ref_date    = yaml_ls$ref_date

 # make output directory
  ref_dir     = sprintf('%s/%s_%s', yaml_ls$xgboost_dir, ref_tag, ref_date)
  if ( !dir.exists(ref_dir) )
    dir.create(ref_dir)
 # make sub directory for classifier
  xgboost_out_dir = file.path(ref_dir, 'xgboost')
  if ( !dir.exists(xgboost_out_dir) )
    dir.create(xgboost_out_dir)

  # define files to make
  meta_ref_f  = sprintf('%s/harmonized_metadata_%s_%s.csv',
                        ref_dir, ref_tag, ref_date)
  sce_all_f   = sprintf('%s/sce_input_%s_combined_%s.rds',
                        ref_dir, ref_tag, ref_date)
  keep_ref_f  = sprintf('%s/keep_dt_%s_combined_%s.txt.gz',
                        ref_dir, ref_tag, ref_date)
  cls_ref_f   = sprintf('%s/integrated_dt_%s_combined_%s.txt.gz',
                        ref_dir, ref_tag, ref_date)
  hvgs_ls_f   = sprintf('%s/harmony_hvgs_%s_combined_%s.txt.gz',
                        ref_dir, ref_tag, ref_date)
  hvgs_mat_f  = sprintf('%s/hvg_values_%s_%s.rds',
                        ref_dir, ref_tag, ref_date)
  sce_ref_f   = sprintf('%s/sce_clean_%s_combined_%s.rds',
                        ref_dir, ref_tag, ref_date)
  # xgboost output files
  xgb_f       = sprintf('%s/xgboost_obj_hvgs_%s_%s.rds',
                        xgboost_out_dir, ref_tag, ref_date)

  allowed_f   = sprintf('%s/allowed_cls_%s_%s.csv',
                        xgboost_out_dir, ref_tag, ref_date)
  # plot files
  pred_annots_pl_f           = sprintf('%s/predicted_annotations_%s_%s.png', ref_dir, ref_tag, ref_date)
  pred_true_annots_pl_f      = sprintf('%s/predicted_and_true_annotations_%s_%s.png', ref_dir, ref_tag, ref_date)
  ref_annots_pl_f            = sprintf('%s/umap_reference_annotations_%s_%s.png', ref_dir, ref_tag, ref_date)
  ref_annots_by_clust_pl_f   = sprintf('%s/umap_reference_annotations_by_cluster_%s_%s.png', ref_dir, ref_tag, ref_date)
  ref_clusts_pl_f            = sprintf('%s/umap_reference_clusters_%s_%s.png', ref_dir, ref_tag, ref_date)

  # get list of scprocess-generated input files
  paths_ls    = get_paths_ls(samples_f, split_var = yaml_ls$split_var, scprocess_v = 'old')

  # parameters for training xgboost
  n_cores   = yaml_ls$n_cores
  min_n_cl  = yaml_ls$min_n_cl
  n_train   = yaml_ls$n_train
  n_dims    = yaml_ls$n_dims
  n_hvgs    = yaml_ls$n_hvgs
  seed      = yaml_ls$seed
  sel_cl    = yaml_ls$sel_cl
  meta_vars = yaml_ls$split_var
  clust_var = 'annotation'


  setDTthreads(n_cores)

  # get table with original annotations
  annots_dt   = get_annots(annots_f)

  # make combined sce object
  prep_files_for_integrating(sce_all_f, meta_ref_f, keep_ref_f, paths_ls)

  # integrate this object (saves hvgs and hvg norm counts)
  if ( !all(file.exists(cls_ref_f, sce_ref_f, hvgs_ls_f, hvgs_mat_f) ) )
    integrate_ref_object( sce_all_f, keep_ref_f,
      cls_ref_f, hvgs_ls_f, hvgs_mat_f, sce_ref_f, n_hvgs, n_cores
    )

  # guess at cluster level
  guesses_dt  = assign_annotations_to_clusters(cls_ref_f, sel_cl, annots_dt, min_cl = 10)

  # load  metadata
  meta_ref    = meta_ref_f %>% fread %>%
    .[, c('sample_id', yaml_ls$split_var), with = FALSE]

  # train
  xgb_obj  = train_celltype_labeller(sce_ref_f, hvgs_mat_f, xgb_f, allowed_f,
    guesses_dt, meta_ref, clust_var = clust_var, use_all_samples = TRUE,
    meta_vars = meta_vars, min_n_cl = min_n_cl, n_train = n_train,
    n_dims = n_dims, n_hvgs = n_hvgs, seed = seed, n_cores = n_cores
    )

  # do diagnostic plots
  make_many_diagnostic_plots(ref_dir, ref_tag, ref_date,
    xgb_obj, annots_dt, allowed_f, cls_ref_f, hvgs_mat_f, hvgs_ls_f, sel_cl,
    pred_annots_pl_f, pred_true_annots_pl_f ,  ref_annots_pl_f, ref_annots_by_clust_pl_f,
    ref_clusts_pl_f
  )

}


# Load, validate, and return the YAML config as a named list.
# Required keys: annots_f, samples_f, xgboost_dir, ref_tag, ref_date.
# Optional keys are filled with defaults defined inside the function.
# (Originally adapted from scjoin)
.get_yaml_ls <- function(yaml_f) {
  # set up default list
  yaml_ls   = list(
    # required parameters
    ## input files
    annots_f        = NULL,
    samples_f       = NULL,

    ## output files
    xgboost_dir     = NULL,
    ref_tag         = NULL,
    ref_date        = NULL,

    # optional parameters
    split_var       = NULL,
    sel_cl          = 0.2,
    n_cores         = 16,
    min_n_cl        = 1000,
    n_train         = 500,
    n_dims          = 50,
    n_hvgs          = 3000,
    seed            = 20240311
  )

  # load yaml file
  yaml_in   = yaml::read_yaml(yaml_f)
  assert_that( length(setdiff(names(yaml_in), names(yaml_ls))) == 0 )

  # check that all required params are specified
  req_params = c('annots_f', 'samples_f', 'xgboost_dir', 'ref_tag', 'ref_date')
  for (r in req_params){
    assert_that(r %in% names(yaml_in), msg = paste0(r, " parameter missing from config file"))
  }

  # add parameters
  for (v in names(yaml_ls)) {
    if (v %in% names(yaml_in))
      yaml_ls[[ v ]]    = yaml_in[[ v ]]
  }

  # change sel_cl
  yaml_ls[['sel_cl']] = paste0('RNA_snn_res.', yaml_ls[['sel_cl']])

  # do some checks
  assert_that( file.exists(yaml_ls$annots_f) )
  assert_that( file.exists(yaml_ls$samples_f) )
  assert_that( dir.exists(yaml_ls$xgboost_dir) )

  return(yaml_ls)
}


# Build a named list of per-sample scprocess output file paths.
# Reads a samples CSV (must contain 'sample_id' and 'scprocess_config_f' columns),
# then resolves SCE, doublet, keep, and metadata paths for each sample.
get_paths_ls <- function(samples_f, split_var = NULL, scprocess_v = c('old', 'new')){

  scprocess_v = match.arg(scprocess_v)

  # load sample data
  samples_tmp = samples_f %>% fread

  # check if columns present
  smpl_cols = c('sample_id', 'scprocess_config_f')
  assert_that(all(smpl_cols %in% colnames(samples_tmp)))
  if(!is.null(split_var)) {
    assert_that(split_var %in% colnames(samples_tmp))
    smpl_cols = c(smpl_cols, split_var)
  }

  samples_dt  = samples_tmp[, smpl_cols, with = FALSE] %>% unique
  sample_ls   = samples_dt$sample_id

  scprocess_configs = samples_tmp$scprocess_config_f %>% unique
  assert_that( all(file.exists(scprocess_configs)) )
  if (scprocess_v == 'old') {
    scprocess_params_ls = scprocess_configs %>% lapply(FUN = .get_scprocess_params_old)

    # make output
    paths_ls    = list(
      sample_ls   = sample_ls,
      sce_paths   = unlist(lapply(scprocess_params_ls, function(x) x$sce_path)),
      dbl_paths   = unlist(lapply(scprocess_params_ls, function(x) x$dbl_path)),
      keep_paths  = unlist(lapply(scprocess_params_ls, function(x) x$keep_path)),
      meta_paths  = unlist(lapply(scprocess_params_ls, function(x) x$meta_f))
    )
  }else{
    scprocess_params_ls = samples_ls %>% lapply(FUN = .get_scprocess_params_new, samples_dt = samples_tmp)

    # make output
    paths_ls = list(
      sample_ls = sample_ls,
      sce_paths = unlist(lapply(scprocess_params_ls, function(x) x$sce_path)),
      meta_paths  = unlist(lapply(scprocess_params_ls, function(x) x$meta_f))
    )
  }
  return(paths_ls)
}


# Parse a single scprocess YAML config (old pipeline layout) and return
# paths to the doublet, keep, SCE, and sample metadata files.
# Asserts that all expected output files exist on disk.
.get_scprocess_params_old <- function(scprocess_config_f) {

  params = yaml::read_yaml(scprocess_config_f)

  # get project directory and check whether outputs are present
  proj_dir  = params$proj_dir
  short_tag = params$short_tag
  full_tag  = params$full_tag
  date      = params$date_stamp
  meta_f    = params$sample_metadata
  assert_that(dir.exists(proj_dir))
  assert_that(file.exists(meta_f))

  dbl_path  = sprintf('%s/output/%s04_doublet_id/scDblFinder_combined_outputs_%s_%s.txt.gz',
                      proj_dir, short_tag, full_tag, date)
  keep_path = sprintf('%s/output/%s05_qc/keep_dt_%s_%s.txt.gz',
                      proj_dir, short_tag, full_tag, date)
  sce_path  = sprintf('%s/output/%s06_integration/sce_clean_%s_%s.rds',
                      proj_dir, short_tag, full_tag, date)

  assert_that(
    file.exists(dbl_path),
    file.exists(keep_path),
    file.exists(sce_path)
  )

  params_ls = list(
    sce_path  = sce_path,
    dbl_path  = dbl_path,
    keep_path = keep_path,
    meta_f    = meta_f
  )

  return(params_ls)
}


# Parse a single scprocess YAML config (new pipeline layout) for a given sample
# and return paths to the SCE and sample metadata files.
.get_scprocess_params_new <- function(sel_s, samples_dt) {
  config_f = samples_dt[sample_id == sel_s, scprocess_config_f ]
  assert_that(length(config_f) == 1)
  params = yaml::read_yaml(config_f)

  # get project directory and check whether outputs are present
  proj_dir  = params$proj_dir
  short_tag = params$short_tag
  full_tag  = params$full_tag
  date      = params$date_stamp
  meta_f    = params$sample_metadata
  assert_that(dir.exists(proj_dir))
  assert_that(file.exists(meta_f))

  sce_path  = sprintf(
    '%s/output/%s_integration/sce_cells_clean_%s_%s_%s.rds',
    sel_s, proj_dir, short_tag, full_tag, date
    )

  assert_that(file.exists(sce_path))

  params_ls = list(
    sce_path  = sce_path,
    meta_f    = meta_f
  )

  return(params_ls)
}



# Load reference cell type annotations from a CSV file.
# Expects columns: sample_id, cell_id, annotation.
# Asserts that no annotation label contains spaces (XGBoost class names must be clean).
get_annots <- function(annots_f) {
  annots_dt = fread(annots_f)
  keep_cols = c('sample_id', 'cell_id', 'annotation')

  assert_that(all(keep_cols %in% colnames(annots_dt)))

  # check for spaces in annotations
  annots = annots_dt$annotation %>% unique()
  assert_that(
    all(!str_detect(annots, pattern = ' ')),
    msg = 'Some annotations contain spaces'
  )

  return(annots_dt)
}


# Prepare three cached files required before Harmony integration (each step is
# skipped if the output file already exists on disk):
#   1. sce_all_f  — merged SCE object across all samples
#   2. meta_ref_f — combined sample-level metadata table
#   3. keep_ref_f — table of cell IDs passing QC in the merged SCE
prep_files_for_integrating <- function(sce_all_f, meta_ref_f, keep_ref_f, paths_ls) {
  sample_ls   = paths_ls$sample_ls
  sce_paths   = paths_ls$sce_paths
  meta_paths  = paths_ls$meta_paths

  # combine sce files
  if (!file.exists(sce_all_f))
    combine_sce_objects(sce_all_f, sce_paths, sample_ls = sample_ls)

  # combine sample metadata
  if (!file.exists(meta_ref_f)) {
    dt_ls       = lapply(meta_paths, fread)
    col_ns      = lapply(dt_ls, colnames)
    assert_that( sapply( col_ns, function(ns) all(ns == col_ns[[1]]) ) %>% all )
    dt_tmp      = rbindlist(dt_ls) %>%
      .[ sample_id %in% sample_ls ]
    assert_that( all(sample_ls %in% unique(dt_tmp$sample_id)))
    fwrite(dt_tmp, file = meta_ref_f)
  }


  # save a txt file with good cells across all samples
  if (!file.exists(keep_ref_f)) {
    dt_tmp      = sce_all_f %>% readRDS %>% colData %>% as.data.table %>%
      .[, .(sample_id, cell_id)]
    fwrite(dt_tmp, file = keep_ref_f)
  }
}


# Assign the majority reference annotation to each Harmony cluster.
# Builds a (cluster x annotation) confusion table, selects the top annotation
# per cluster by proportion, then propagates that label back to individual cells.
# Returns a data.table with columns: sample_id, cell_id, cluster, annotation.
assign_annotations_to_clusters <- function(cls_ref_f, sel_cl, annots_dt, min_cl = 10) {
  # check new clusters
  cls_ref_dt  = cls_ref_f %>% fread %>%
    .[, .(sample_id, cell_id, cluster = get(sel_cl) %>% fct_infreq %>% fct_lump_min(min_cl)) ] %>%
    .[ cluster != "Other" ]

  confuse_dt  = merge(cls_ref_dt, annots_dt, by = 'cell_id', all.x = TRUE) %>%
    .[, .N, by = .(annotation, cluster) ] %>%
    .[ is.na(annotation), annotation := "missing" ] %>%
    dcast( cluster ~ annotation, value.var = 'N', fill = 0) %>%
    melt( id = c('cluster', 'missing'), var = 'annotation', val = 'N' ) %>%
    .[, N0          := N + 1 ] %>%
    .[, log_N       := log(N0) ] %>%
    .[, p_cluster   := N / sum(N), by = cluster ] %>%
    .[, log_p_cl1   := log(p_cluster) ] %>%
    .[, p_annot     := N / sum(N), by = annotation ] %>%
    .[, log_p_annot := log(p_annot) ]

  top_guesses = confuse_dt[ order(cluster, -p_cluster) ] %>%
    .[, .SD[1:min(.N, 1)], by = cluster ] %>%
    .[, .(cluster, missing, top_annotation = annotation,
          N_annot = N, pct_cl = round(p_cluster * 100, 1),
          log2_ratio = log2((missing + 1) / (N + 1)) %>% round(1) ) ] %>%
    .[ order(log2_ratio) ]
  print(top_guesses)

  # guess all cells with this
  guesses_dt  = cls_ref_dt %>%
    merge(top_guesses[, .(cluster, annotation = top_annotation)], by = 'cluster') %>%
    .[, .(sample_id, cell_id, cluster, annotation) ]

  return(guesses_dt)
}


# Load multiple SCE RDS files, optionally subset to sample_ls, validate that
# rownames and colData column names match across objects, cbind them, and
# save the merged SCE to sce_new_f.
combine_sce_objects <- function(sce_new_f, sce_fs, sample_ls = NULL) {
  # check everything makes sense
  assert_that( all(file.exists(sce_fs)) )
  # load up all sce files, restricting to samples as we go
  .load_sce <- function(sce_f) {
    message('  ', sce_f)
    sce         = sce_f %>% readRDS
    if (!is.null(sample_ls)) {
      keep_idx    = sce$sample_id %in% sample_ls
      sce         = sce[, keep_idx ]
    }
    return(sce)
  }
  message('loading up sce files')
  sce_ls      = lapply(sce_fs, .load_sce)

  # check that row names match
  all_row_ns  = lapply(sce_ls, rownames)
  assert_that( sapply( all_row_ns, function(rn) all(rn == all_row_ns[[1]]) ) %>% all )

  # check that coldata names match
  coldata_ns  = lapply(sce_ls, function(sce) sce %>% colData %>% names)
  assert_that( sapply( coldata_ns, function(ns) all(ns == coldata_ns[[1]]) ) %>% all )

  # join everything together
  message('joining them together')
  sce_joined  = do.call(cbind, sce_ls)

  if (!is.null(sample_ls)) {
    missing_ls  = setdiff( sample_ls, sce_joined$sample_id )
    if ( length(missing_ls) > 0 ) {
      warnings('the following samples were requested but not found in the sce objects:',
               paste(missing_ls, collapse = " "))
    }
  }

  message('saving joined object')
  saveRDS(sce_joined, file = sce_new_f, compress = FALSE)
}

# Run Harmony batch correction on the merged SCE, compute UMAP and Seurat clusters,
# and save four output files: cluster/UMAP table (cls_ref_f), HVG gene list
# (hvgs_ls_f), HVG normalised count matrix (hvgs_mat_f), and annotated SCE (sce_ref_f).
integrate_ref_object <- function(sce_all_f, keep_ref_f,
                                 cls_ref_f, hvgs_ls_f, hvgs_mat_f, sce_ref_f, n_hvgs, n_cores) {
  # define parameters
  exc_regex   = "^MT-"     # regex for genes to exclude (mitochondrial)
  n_dims      = 50         # number of PCA / Harmony dimensions to use
  theta       = 0.1        # Harmony regularisation (lower = less regularisation)
  res_ls      = c(0.1, 0.2, 0.5, 1, 2)  # Seurat clustering resolutions

  message('running Harmony')
  # load sce, restrict to ok cells
  message('  setting up cluster')
  plan("multicore", workers = n_cores)
  options( future.globals.maxSize = 2^35 )

  message('  loading relevant cell ids')
  keep_ids    = fread(keep_ref_f) %>% .$cell_id

  message('  loading sce')
  assert_that( length(keep_ids) > 0 )
  sce         = readRDS(sce_all_f) %>% .[, keep_ids ]
  assert_that("sample_id" %in% colnames(colData(sce)))

  # exclude genes if requested
  if (!is.null(exc_regex)) {
    exc_idx     = rownames(sce) %>% str_detect(exc_regex)
    if (sum(exc_idx) > 0) {
      exc_gs      = rowData(sce)$symbol[ exc_idx ]
      sprintf("    excluding %d genes: %s", sum(exc_idx), paste0(exc_gs, collapse = " ")) %>%
        message
      sce         = sce[ !exc_idx, ]
    }
  }

  # turn into seurat object
  message('  prepping seurat object')
  suppressWarnings({
    seu_obj     = Seurat::CreateSeuratObject(
      counts      = counts(sce),
      meta.data   = data.frame(colData(sce)),
      project     = "dummy"
    )
  })
  # rm(sce); gc()

  # run harmony including doublets
  message('  running Harmony')
  hmny_dt     = .run_harmony(seu_obj, n_hvgs, n_dims, hvgs_ls_f, hvgs_mat_f, theta = theta, res_ls)

  # save outputs
  fwrite(hmny_dt, file = cls_ref_f)


  # make sce file with only ok cells, and nice clusters
  message('  making nice clean sce')
  sce         = .annotate_sce_w_harmony(sce, hmny_dt)
  message('    saving')
  saveRDS(sce, file = sce_ref_f, compress = FALSE)
  message('done!')
}


# Generate five diagnostic UMAP plots after training:
#   1. pred_annots_pl_f          — naive per-cell XGBoost predictions
#   2. pred_true_annots_pl_f     — true vs. cluster-smoothed predictions side-by-side
#   3. ref_annots_pl_f           — reference UMAP coloured by raw annotation
#   4. ref_annots_by_clust_pl_f  — reference UMAP coloured by cluster-smoothed annotation
#   5. ref_clusts_pl_f           — reference UMAP coloured by Harmony cluster ID
make_many_diagnostic_plots <- function(ref_dir, ref_tag, ref_date, xgb_obj,
  annots_dt,allow_f, cls_ref_f, hvgs_mat_f, hvgs_ls_f, sel_cl,
  pred_annots_pl_f, pred_true_annots_pl_f ,  ref_annots_pl_f, ref_annots_by_clust_pl_f,
  ref_clusts_pl_f) {

  annot_lvls  = annots_dt$annotation %>% unique()
  cls_ref_dt  = cls_ref_f %>% fread %>%
    .[, .(cell_id, UMAP1, UMAP2)]

  # use classifier to make predictions
  preds_dt    = .apply_classifier(xgb_obj, allow_f,  hvgs_mat_f)
  hmny_dt     = .load_clusters(cls_ref_f)
  min_cl_prop = 0.5
  min_cl_size = 100
  guesses_dt  = .apply_labels_by_cluster(hmny_dt, preds_dt, min_cl_prop, min_cl_size)

  confuse_dt  = merge(
    annots_dt[, .(cell_id, annot_true = annotation)],
    guesses_dt[, .(cell_id, annot_pred = cl_pred_RNA_snn_res.1)],
    # preds_dt[, .(cell_id, annot_pred = cl_pred)],
    by = 'cell_id', all = TRUE) %>%
    .[, .N, by = .(annot_true, annot_pred) ] %>%
    .[ is.na(annot_true), annot_true := 'missing' ] %>%
    dcast( annot_pred ~ annot_true, value.var = 'N', fill = 0) %>%
    melt( id = 'annot_pred', var = 'annot_true', val = 'N' ) %>%
    .[, N0          := N + 1 ] %>%
    .[, log_N       := log(N0) ] %>%
    .[, p_pred      := N / sum(N), by = annot_pred ] %>%
    .[, log_p_cl1   := log(p_pred) ] %>%
    .[, p_true      := N / sum(N), by = annot_true ] %>%
    .[, log_p_true  := log(p_true) ]

  confuse_dt[ (N > 10) & (annot_true != 'missing') ] %>%
    .[ as.character(annot_pred) != as.character(annot_true) ] %>%
    .[ order(-p_true) ] %>% .[ (N >= 10) & (p_true > 0.1) ] %>%
    .[, .(
      annot_true  = abbreviate(annot_true, minlength = 30),
      annot_pred  = abbreviate(annot_pred, minlength = 30),
      N, pct_true = round(p_true * 100, 1) ) ]

  # make plot of predictions
  plot_dt     = cls_ref_dt %>%
    merge(preds_dt[, .(cell_id, cl_pred_naive, p_pred)],
          id = 'cell_id', all = TRUE)
  g = ggplot( plot_dt ) +
    aes( x = UMAP1, y = UMAP2, colour = cl_pred_naive ) +
    geom_point( size = 0.1 ) +
    scale_colour_manual( values = c(nice_cols, nice_cols),
                         guide = guide_legend(override.aes = list(size = 3)) ) +
    theme_bw() +
    theme( axis.text = element_blank(), legend.text = element_text( size = 6 ),
           panel.grid = element_blank() ) +
    labs(colour = "Predicted clusters")
    ggsave(pred_annots_pl_f, g, h = 7, w = 10)

  plot_dt     = cls_ref_dt %>%
    merge(guesses_dt[, .(cell_id, annot_pred = cl_pred_RNA_snn_res.1)],
          id = 'cell_id', all = TRUE) %>%
    merge(annots_dt[, .(cell_id, annot_true = annotation)],
          id = 'cell_id', all = TRUE) %>%
    melt( id = c('cell_id', 'UMAP1', 'UMAP2'), var = 'variable', val = 'annotation' ) %>%
    .[, variable    := variable %>% fct_relevel('annot_true') ] %>%
    .[, annotation  := annotation %>% factor( levels = c(annot_lvls, 'unknown') ) %>%
        fct_relevel('unknown', after = Inf) ]
  g = ggplot( plot_dt[ !is.na(annotation) ] ) +
    aes( x = UMAP1, y = UMAP2, colour = annotation ) +
    geom_point( size = 0.1 ) +
    scale_colour_manual( values = c(nice_cols[ seq(length(annot_lvls) - 1) ], 'grey', 'grey20'),
                         guide = guide_legend(override.aes = list(size = 3)) ) +
    facet_wrap( ~ variable ) +
    theme_bw() +
    theme( axis.text = element_blank(), aspect.ratio = 1,
           panel.grid = element_blank(),
           strip.background = element_rect( fill = 'white' ),
           legend.text = element_text( size = 6 ), legend.position = 'bottom' ) +
    labs(colour = "Annotation")
    ggsave(pred_true_annots_pl_f, g,  h = 7, w = 11)


  annot_cols  = nice_cols[ seq_along(annot_lvls) ] %>% setNames(annot_lvls)
  annot_cols[[ 'Neurons' ]] = 'grey'

  hmny_dt     = cls_ref_f %>% fread %>%
    .[, .(cell_id, UMAP1, UMAP2, cluster = get(sel_cl) %>% fct_lump_min(10))]
  plot_dt     = merge(hmny_dt, annots_dt, by = 'cell_id', all = TRUE)
  g = ggplot(plot_dt[ !is.na(annotation) ]) +
    aes( x = UMAP1, y = UMAP2, colour = annotation ) +
    geom_point( size = 0.1 ) +
    scale_colour_manual( values = annot_cols,
                         guide = guide_legend( override.aes = list(size = 3))) +
    theme_bw() +
    theme( aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank(),
           panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    labs( colour = 'annotation, raw' )
    ggsave(ref_annots_pl_f, g, height = 6, width = 8)

  guesses_dt  = assign_annotations_to_clusters(cls_ref_f, sel_cl, annots_dt, min_cl = 10)
  plot_dt     = merge(hmny_dt, guesses_dt, by = 'cell_id', all = TRUE)
  g = ggplot(plot_dt[ !is.na(annotation) ]) +
    aes( x = UMAP1, y = UMAP2, colour = annotation ) +
    geom_point( size = 0.1 ) +
    scale_colour_manual( values = annot_cols,
                         guide = guide_legend( override.aes = list(size = 3))
    ) +
    theme_bw() +
    theme( aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank(),
           panel.grid.major = element_blank(), panel.grid.minor = element_blank() ) +
    labs( colour = 'annotation, smoothed' )
    ggsave(ref_annots_by_clust_pl_f, g, height = 6, width = 8)

  g = ggplot(hmny_dt) +
    aes( x = UMAP1, y = UMAP2, colour = cluster ) +
    geom_point( size = 0.1 ) +
    scale_colour_manual( values = c(nice_cols, nice_cols),
                         guide = guide_legend( override.aes = list(size = 3))
    ) +
    theme_bw() +
    theme( aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank(),
           panel.grid.major = element_blank(), panel.grid.minor = element_blank() )
    ggsave(  ref_clusts_pl_f, g, height = 6, width = 9)

}

# Apply a trained XGBoost model to all cells in hvgs_mat_f.
# Returns a data.table with per-cell predicted label (cl_pred_naive), raw best-guess
# label (cl_pred_raw), and max class probability (p_pred).
# Cells with p_pred < min_pred are labelled "unknown".
.apply_classifier <- function(xgb_obj, allow_f, hvgs_mat_f, min_pred = 0.5) {
  # load up xgboost object
  allow_dt    = fread(allow_f)
  xgb_names   = allow_dt$cluster

  # load up HVG values
  hvgs_mat    = hvgs_mat_f %>% readRDS
  hvgs_mat    = as.matrix(hvgs_mat)

  # get probabilities for each predicted cluster
  probs_mat   = predict(xgb_obj, hvgs_mat, reshape = TRUE)
  assert_that(
    length(rownames(hvgs_mat)) == nrow(probs_mat)
  )
  probs_mat   = probs_mat %>%
    set_colnames( xgb_names ) %>%
    set_rownames( rownames(hvgs_mat) )

  # prediction for each cell
  preds_dt    = data.table(
    cell_id     = rownames(probs_mat),
    cl_pred_raw = colnames(probs_mat)[ apply(probs_mat, 1, which.max) ],
    p_pred      = apply(probs_mat, 1, max)
  ) %>% cbind(probs_mat) %>%
    .[, cl_pred_naive := (p_pred > min_pred) %>% ifelse(cl_pred_raw, "unknown") %>%
        factor(levels = c(xgb_names, "unknown"))  ]

  return(preds_dt)
}



# Full XGBoost training pipeline for cell type labelling.
# Steps: optionally select a representative cell subset (balanced by metadata) ->
#        extract HVG expression matrix -> train/validation split ->
#        train XGBoost with early stopping -> print validation confusion metrics.
# Saves the XGBoost model to xgb_f and the allowed cluster table to allowed_f.
train_celltype_labeller <- function(sce_f, hvgs_xgb_f, xgb_f, allowed_f,
   clusters_dt, meta_dt, clust_var = "cluster", use_all_samples = FALSE, meta_vars = NULL,
   min_n_cl = 200, n_train = 1000, n_dims = 50, sel_gs = NULL, n_hvgs = 2000,
   seed = 123, n_cores = 4) {


  if (use_all_samples) {
    message('  using all samples')
    samples_dt  = copy(clusters_dt)
  } else {
    message('  getting representative subset of cells for training')
    set.seed(seed)
    samples_dt  = .get_representative_subset(clusters_dt, meta_dt,
     clust_var = clust_var, meta_vars = meta_vars,
     n_per_type = min_n_cl, min_n = 10)

    message(sprintf('    %d cells chosen, split like so:', nrow(samples_dt)))
    samples_dt[, .N, by = 'cluster'] %>% .[ order(-N) ] %>% print
  }
  # now subset
  sel_ids     = samples_dt$cell_id
  sel_cls     = clusters_dt[ cell_id %in% sel_ids ] %>%
    .[, .(sample_id, cell_id, cluster = get(clust_var))]

  # get HVGs for these
  message('  calculate HVGs on these')
  hvgs_mat    = .get_hvgs_mat(hvgs_xgb_f, sce_f, sel_ids, what = "hvgs",
    sel_gs = sel_gs, n_hvgs = n_hvgs, n_dims = n_dims, n_cores = n_cores)
  hvg_gs      = hvgs_mat %>% colnames %>% setdiff(c("cluster", "cell_id"))

  # make data for xgboost
  message('  split into train/test')
  set.seed(seed)
  data_broad  = .load_train_test_data(sel_cls, hvgs_mat, min_n_cl, n_train)
  train_dt    = data_broad$train
  valid_dt    = data_broad$valid

  # run xgboost
  message('  run XGBoost')
  set.seed(seed)
  xgb_obj     = .run_boost_watchlist(train_dt, valid_dt, n_cores)
  # xgb_obj$cl_lvls = levels(train_dt$cluster)
  message('    saving')
  saveRDS(xgb_obj, file = xgb_f, compress = FALSE)
  allowed_dt  = data.table( cluster = levels(train_dt$cluster) )
  fwrite(allowed_dt, file = allowed_f)

  # predict on validation data
  message('  print some outputs to show performance on validation data')
  valid_all   = rbind(data_broad$valid, data_broad$valid_rest)
  pred_valid  = .get_pred_valid(xgb_obj, valid_all)
  conf_dt     = .calc_confuse_xgboost_dt(pred_valid)
  conf_tmp    = .calc_confuse_xgboost_dt(pred_valid[ p_pred > 0.5 ])
  conf_dt[ (cl_pred != cl_true) & (prop > 0.1) ] %>%
    .[, .(cl_true, cl_pred, N_true, pct_pred = round(prop * 100, 1))] %>%
    .[ order(-pct_pred) ] %>%
    print
  conf_tmp[ (cl_pred != cl_true) & (prop > 0.01) ] %>%
    .[, .(cl_true, cl_pred, N_true, pct_pred = round(prop * 100, 1))] %>%
    .[ order(-pct_pred) ] %>%
    print
  conf_tmp[ (cl_pred == cl_true) ] %>% .[ order(prop) ] %>%
    .[, .(cl_true, N_true, pct_true = round(prop * 100, 1))] %>% print

  return(xgb_obj)

  message('done.')
}



# Run the Seurat + Harmony integration pipeline on a Seurat object.
# Normalises counts, identifies HVGs, runs PCA, applies Harmony batch correction
# (if >1 sample), computes UMAP, and calls FindClusters at multiple resolutions.
# Saves the HVG gene-list table to hvgs_ls_f and the normalised HVG count matrix
# to hvgs_mat_f. Returns a data.table of cluster assignments, UMAP coordinates,
# and Harmony PCA embeddings.
.run_harmony <- function(seu_obj, n_hvgs, n_dims, hvgs_ls_f, hvgs_mat_f, theta = 0, res_ls) {
  message('    normalizing and finding HVGs')

  options(future.globals.maxSize = 500 * 1024^3)
  seu_obj   = seu_obj %>%
    NormalizeData( verbose = TRUE ) %>%
    FindVariableFeatures( nfeatures = n_hvgs, verbose = TRUE ) %>%
    ScaleData( verbose = TRUE )

  # check whether we have one or more values of batch
  n_samples       = seu_obj$sample_id %>% unique %>% length
  is_one_batch    = n_samples == 1
  this_reduction  = ifelse(is_one_batch, "pca", "harmony")
  if (is_one_batch) {
    stop("only one sample; doesn't make sense to do this")
    # run Seurat pipeline, plus clustering
    warning('    only one sample; not doing Harmony')
    message('    running UMAP')
    seu_obj   = seu_obj %>%
      RunPCA( npcs = n_dims, verbose = FALSE ) %>%
      RunUMAP( reduction = this_reduction, dims = 1:n_dims, verbose = FALSE  )
  } else {
    # run Seurat pipeline, plus clustering
    message('    running Harmony and UMAP')
    seu_obj   = seu_obj %>%
      RunPCA( npcs = n_dims, verbose = FALSE )
    seu_obj = seu_obj %>% RunHarmony( c("sample_id"), theta = theta, verbose = TRUE )
    seu_obj = seu_obj %>% RunUMAP( reduction = this_reduction, dims = 1:n_dims, verbose = FALSE  )
  }

  message('    finding clusters')
  seu_obj   = seu_obj %>%
    FindNeighbors( reduction = this_reduction, dims = 1:n_dims, verbose = FALSE  ) %>%
    FindClusters( resolution = res_ls, verbose = FALSE ) %>%
    identity()

  # switch cluster off
  plan("sequential")

  # saving datatable with hvg stats
  message('    recording highly variable genes')
  message('     saving table with highly variable gene statistics')
  hvgs_dt     = seu_obj@assays$RNA@meta.data %>%
    as.data.table %>%
    .[, gene_id := rownames(seu_obj) ] %>%
    .[, is_hvg  := vf_vst_counts_variable ]
  fwrite( hvgs_dt, file = hvgs_ls_f  )

  # saving normalized counts for highly variable genes
  message('     saving matrix with normalized counts for hvgs')
  hvgs_ls = hvgs_dt[ is_hvg == TRUE, gene_id ]
  hvgs_mat    = GetAssayData(seu_obj, layer = "data", assay = "RNA") %>%
    .[ hvgs_ls, ] %>% t
  saveRDS(hvgs_mat, file = hvgs_mat_f)

  message('    recording clusters')
  clusts_dt   = seu_obj[[]] %>% as.data.table %>% .[, integration := "harmony" ]
  cl_vs       = names(clusts_dt) %>% str_subset('RNA_snn_res\\.')
  clusts_dt   = clusts_dt[, c('integration', 'cell_id', 'sample_id', cl_vs), with = FALSE ]
  for (cl_v in cl_vs) {
    clusts_dt[[ cl_v ]] = clusts_dt[[ cl_v ]] %>% fct_infreq %>% as.integer %>%
      sprintf("cl%02d", .)
  }

  # get embeddings
  message('    extracting other outputs')
  hmny_pca_dt = Embeddings(seu_obj, reduction = this_reduction) %>%
    as.data.table( keep.rownames = "cell_id" ) %>%
    setnames(names(.), names(.) %>% str_replace("_(?=[0-9]$)", "_0") %>%
               str_replace("harmony", "hmny_pca"))
  umap_dt     = Embeddings(seu_obj, reduction = "umap") %>%
    as.data.table( keep.rownames = "cell_id" ) %>%
    setnames(names(.), names(.) %>% str_replace("umap_", "UMAP"))

  # join together
  hmny_dt     = clusts_dt %>% merge(umap_dt, by = "cell_id") %>%
    merge(hmny_pca_dt, by = "cell_id")

  return(hmny_dt)
}

# Add Harmony-derived UMAP coordinates and cluster assignments to a SCE object's
# colData, restricting the SCE to only the cells present in hmny_dt.
.annotate_sce_w_harmony <- function(sce, hmny_dt) {
  # restrict to just ok cells
  sce       = sce[ , hmny_dt$cell_id ]

  # get useful harmony variables
  hmny_vs   = c('UMAP1', 'UMAP2', str_subset(names(hmny_dt), "RNA_snn_res"))

  # add these to sce object
  for (v in hmny_vs) {
    if (str_detect(v, "RNA_snn_res")) {
      colData(sce)[[ v ]] = hmny_dt[[ v ]] %>% factor
    } else {
      colData(sce)[[ v ]] = hmny_dt[[ v ]]
    }
  }

  return(sce)
}


# Greedily select a subset of samples that (a) covers at least n_per_type cells
# per cluster and (b) balances metadata covariates (meta_vars) across the selection.
# Samples are added one at a time via .pick_next_sample until each cluster is covered.
.get_representative_subset <- function(cl_dt, meta_dt, clust_var = "cluster",
  meta_vars = NULL, n_per_type = 100, min_n = 10) {
  # initialize
  cl_tmp      = copy(cl_dt) %>%
    .[, .(sample_id, cell_id, cluster = get(clust_var))]
  sample_list = NULL
  cl_sub      = data.table(
    cell_id = character(0),
    sample_id = character(0),
    cluster = character(0)
    )

  # define metadata combinations we want to balance
  if (is.null(meta_vars)) {
    meta_track  = meta_dt[, .(sample_id, combn_var = 'dummy')]
  } else {
    meta_track  = meta_dt[, c('sample_id', meta_vars), with = FALSE] %>%
      .[, combn_var := do.call(paste, c(.SD, sep = "_")),
        .SDcols = meta_vars, by = meta_vars ]
  }

  props_all  = meta_track[, .N, by = combn_var] %>%
    .[, p_all := N / sum(N) ]

  # add samples one by one until we have at least n_per_type cells per cluster
  totals_dt   = cl_tmp[, .N, by = .(sample_id, cluster)] %>%
    .[ N > min_n ]
  types_list  = cl_tmp[, .N, .(cluster)] %>%
    .[ order(N) ] %>%
    use_series('cluster') %>% as.character

  # loop
  for (this_type in types_list) {
    n_type    = cl_sub[ cluster == this_type ] %>% nrow
    n_total   = totals_dt[ cluster == this_type ]$N %>% sum
    while (n_type < min(n_per_type, n_total)) {
      # which samples would help?
      sample_opts = cl_tmp[ cluster == this_type, .N, by = sample_id] %>%
        setorder(-N) %>% .[ N > min_n ] %>% use_series("sample_id")

      # pick one which improves metadata representation
      sel_sample  = .pick_next_sample(meta_track, props_all, sample_list, sample_opts)

      # add to list
      sample_list = c(sample_list, sel_sample)
      cl_sub      = rbind(cl_sub, cl_tmp[sample_id == sel_sample])
      cl_tmp      = cl_tmp[!(sample_id %in% sample_list)]

      # update count
      n_type      = cl_sub[ cluster == this_type ] %>% nrow
    }
  }

  # check worked
  ns_all      = totals_dt[, .(n_all = sum(.N)), by = cluster]
  ns_sub      = cl_sub[, .(n_sub = .N), by = cluster]
  check_dt    = merge(ns_all, ns_sub, by = 'cluster', all = TRUE) %>%
    .[, n_per_type  := n_per_type ] %>%
    .[, is_ok       := n_sub >= pmin(n_per_type, n_all) ]
  assert_that( all(check_dt$is_ok) )

  return(cl_sub)
}


# Choose the next sample to add to the training set to improve metadata balance.
# On the first call (sample_list is NULL) one sample is picked uniformly at random.
# Subsequently, the metadata combination (combn_var) most under-represented relative
# to its target proportion is identified, and a matching sample is drawn at random.
.pick_next_sample <- function(meta_track, props_all, sample_list, sample_opts) {
  # to start pick one at random
  if (is.null(sample_list))
    return(sample(sample_opts, 1))

  # otherwise calc current props
  props_now   = meta_track[ sample_id %in% sample_list ] %>%
    .[, .N, by = combn_var ] %>%
    .[, p_now   := N / sum(N) ]

  # combine with target props
  props_now   = merge(props_now, props_all, by = 'combn_var', all.y = TRUE) %>%
    .[ is.na(p_now), p_now := 0 ] %>%
    .[, p_delta := p_all - p_now ]

  # which is most out of whack, in the samples where we see this celltype?
  vars_ok     = meta_track[ sample_id %in% sample_opts ]$combn_var %>% unique
  sel_val     = props_now[ order(-p_delta) ] %>%
    .[ combn_var %in% vars_ok ] %>%
    use_series('combn_var') %>% .[[1]]

  # pick a sample at random from these
  sel_sample  = meta_track[ (combn_var == sel_val) & (sample_id %in% sample_opts) ]$sample_id %>%
    sample(1)

  return(sel_sample)
}


# Extract a cell x feature matrix for XGBoost training from an SCE file.
# what = "hvgs": returns normalised counts for the top n_hvgs variable genes (cells x genes).
# what = "pca" : returns the top n_dims Seurat PCA embeddings instead (cells x PCs).
# The result is cached to hvgs_xgb_f and reloaded on subsequent calls
# (unless overwrite = TRUE).
.get_hvgs_mat <- function(hvgs_xgb_f, sce_f, sel_ids, what = c("pca", "hvgs"),
                          sel_gs = NULL, n_hvgs = 2000, n_dims = 50, n_cores = 4, overwrite = FALSE) {
  if (file.exists(hvgs_xgb_f) & overwrite == FALSE) {
    message('  already done')
    hvgs_mat    = readRDS(hvgs_xgb_f)

    return(hvgs_mat)
  }
  # check inputs
  what        = match.arg(what)

  message('  subsetting to HVGs')
  # load sce, restrict to ok cells
  message('    setting up cluster')
  plan("multicore", workers = n_cores)
  options( future.globals.maxSize = 2^35 )

  message('    loading sce')
  sce_sel     = readRDS(sce_f) %>% .[, sel_ids ]

  # restrict to specified genes
  if (!is.null(sel_gs)) {
    ref_gs      = rowData(sce_sel)[[ names(sel_gs) ]]
    sel_gs_v    = unlist(sel_gs)
    assert_that( all(sel_gs_v %in% ref_gs) )
    sel_idx     = ref_gs %in% sel_gs_v
    sce_sel     = sce_sel[ sel_idx, ]
  }
  # get counts
  counts_mat  = counts(sce_sel)
  if (!is.null(sel_gs)) {
    rownames(counts_mat) = ref_gs[ sel_idx ]
  }

  # turn into seurat object
  message('    converting to Seurat object')
  seu_obj     = Seurat::CreateSeuratObject(
    counts      = counts_mat,
    meta.data   = data.frame(colData(sce_sel)),
    project     = "MS2"
  )
  rm(sce_sel); gc()

  # run Seurat pipeline, plus clustering
  message('    finding HVGs')
  seu_obj     = NormalizeData(seu_obj, verbose = FALSE ) %>%
    FindVariableFeatures( nfeatures = n_hvgs, verbose = FALSE )
  var_feats   = VariableFeatures(seu_obj)
  seu_obj     = seu_obj %>%
    ScaleData( verbose = FALSE ) %>%
    RunPCA( features = var_feats, verbose = FALSE, npcs = n_dims ) %>%
    identity()

  # switch cluster off
  plan("sequential")

  if (what == "pca") {
    # now use these genes for all cells
    message('    extract PCs, save')
    save_mat    = Embeddings(seu_obj, reduction = "pca")
  } else if (what == "hvgs") {
    # now use these genes for all cells
    message('    extract HVGs, save')
    save_mat    = GetAssayData(seu_obj, slot = "data", assay = "RNA") %>%
      .[ var_feats, ] %>% t
  }
  saveRDS(save_mat, file = hvgs_xgb_f)

  message('  done')

  return(save_mat)
}



# Join cluster labels with the HVG matrix and split into train / validation / test sets.
# Clusters smaller than min_cells are dropped.
# Train and validation each sample up to ceil(n_train / 2) cells per cluster;
# remaining labelled cells go to valid_rest; unlabelled cells go to test.
.load_train_test_data <- function(clusts_dt, hvgs_mat, min_cells, n_train) {
  # check inputs
  assert_that( all(clusts_dt$cell_id %in% rownames(hvgs_mat)) )
  clusts_na   = hvgs_mat %>% as.data.table( keep.rownames = 'cell_id' ) %>%
    merge(clusts_dt, by = "cell_id", all = TRUE) %>%
    .[, sample_id := NULL ]

  # which clusters are too small to bother with?
  ns_dt       = clusts_na[ !is.na(cluster) ] %>% .[, .N, by = cluster]
  keep_cl     = ns_dt[ N >= min_cells ]$cluster %>% as.character

  # get train data, balance samples
  train_dt    = clusts_na[ cluster %in% keep_cl ] %>%
    .[, .SD[ sample(ceiling(min(.N, n_train) / 2)) ], by = cluster ] %>%
    .[, cluster := cluster %>% fct_drop ]

  # get validation data
  valid_dt    = clusts_na[ cluster %in% keep_cl ] %>%
    .[ !(cell_id %in% train_dt$cell_id) ] %>%
    .[, .SD[ sample(ceiling(min(.N, n_train) / 2)) ], by = cluster ] %>%
    .[, cluster := cluster %>% fct_drop ]
  valid_rest  = clusts_na[ cluster %in% keep_cl ] %>%
    .[ !(cell_id %in% c(train_dt$cell_id, valid_dt$cell_id)) ] %>%
    .[, cluster := cluster %>% fct_drop ]

  # get test data
  test_dt     = clusts_na[ is.na(cluster) ]

  # some checks
  assert_that( all( levels(valid_dt$cluster) == levels(train_dt$cluster) ) )
  chk_dt_1    = clusts_na[ is.na(cluster) | (cluster %in% keep_cl) ]
  chk_dt_2    = rbind(train_dt, valid_dt, valid_rest, test_dt)
  assert_that( nrow(chk_dt_1) == nrow(chk_dt_2) )
  assert_that( all( sort(chk_dt_1$cell_id) == sort(chk_dt_2$cell_id) ) )

  return(list(
    train       = train_dt,
    valid       = valid_dt,
    valid_rest  = valid_rest,
    test        = test_dt
  ))
}

# Train a multiclass XGBoost classifier (multi:softprob objective) with early stopping.
# Uses the validation set as a watchlist; training stops after 5 rounds without
# improvement on the test evals. Returns the trained xgb.Booster object.
.run_boost_watchlist <- function(train_dt, valid_dt, n_cores) {
  # convert training data to expected format
  assert_that( !is.null(train_dt$cluster) )
  train_vars  = colnames(train_dt) %>% setdiff(c("cluster", "cell_id"))
  assert_that( length(train_vars) > 0 )
  train_mat   = train_dt[, c("cell_id", train_vars), with = FALSE] %>%
    as.matrix( rownames = "cell_id" )
  train_cl    = train_dt$cluster
  train_y     = train_cl %>% as.integer %>% `-`(1)
  # weights_v   = 1 / table(train_y) * 1000
  # weights_y   = weights_v[ train_y + 1 ]
  assert_that(
    min(train_y) == 0,
    max(train_y) + 1 == length(levels(train_dt$cluster))
  )

  # convert validation data to expected format
  valid_mat   = valid_dt[, c("cell_id", train_vars), with = FALSE] %>%
    as.matrix( rownames = "cell_id" )
  valid_cl    = valid_dt$cluster
  valid_y     = valid_cl %>% as.integer %>% `-`(1)

  # blah
  dtrain      = xgb.DMatrix( data = train_mat, label = train_y )
  dvalid      = xgb.DMatrix( data = valid_mat, label = valid_y )
  watchlist   = list( train = dtrain, test = dvalid )

  boost_obj = xgb.train(
    data = dtrain,
    evals = watchlist,
    params = list(objective = "multi:softprob",
                  num_class = max(train_y) + 1,
                  nthread = n_cores),
    nrounds = 100,
    early_stopping_rounds = 5,
    verbose = 2
    )

  return(boost_obj)
}


# Generate XGBoost predictions on a validation data.table.
# Returns a data.table with the true label (cl_true), predicted label (cl_pred),
# max class probability (p_pred), and the full per-class probability matrix.
.get_pred_valid <- function(boost_obj, valid_dt) {
  # get validation data
  train_vars  = colnames(valid_dt) %>% setdiff(c("cluster", "cell_id"))
  assert_that( length(train_vars) > 0 )
  valid_mat   = valid_dt[, c("cell_id", train_vars), with = FALSE] %>%
    as.matrix( rownames = "cell_id" )

  # get probabilities for each predicted cluster
  probs_mat   = predict(boost_obj, valid_mat, reshape = TRUE)
  assert_that(
    length(levels(valid_dt$cluster)) == ncol(probs_mat),
    length(rownames(valid_mat)) == nrow(probs_mat)
  )
  probs_mat   = probs_mat %>%
    set_colnames(levels(valid_dt$cluster)) %>%
    set_rownames(rownames(valid_mat))

  # prediction for each cell
  pred_valid  = data.table(
    cell_id     = rownames(probs_mat),
    cl_true     = valid_dt$cluster,
    cl_pred     = colnames(probs_mat)[ apply(probs_mat, 1, which.max) ],
    p_pred      = apply(probs_mat, 1, max)
  ) %>% cbind(probs_mat)

  return(pred_valid)
}


# Compute a confusion matrix summary from XGBoost prediction output.
# Returns per-(true, predicted) pair: cell count (N), per-true-class total (N_true),
# fraction of true-class cells (prop), Laplace-smoothed logit of that fraction
# (logit), and per-true-class Shannon entropy (H).
.calc_confuse_xgboost_dt <- function(pred_valid) {
  n_cats    = length(unique(pred_valid$cl_true))
  confuse_dt  = pred_valid[, .N, by=.(cl_true, cl_pred)] %>%
    .[, N_true  := sum(N), by = cl_true] %>%
    .[, prop    := N / N_true ] %>%
    .[, logit   := qlogis((N+1) / (N_true + n_cats))] %>%
    .[, H       := -sum(prop*log2(prop)), by = cl_true ]

  return(confuse_dt)
}

# Load the Harmony cluster output file, drop rows without UMAP coordinates
# (i.e. cells excluded from the clean UMAP), and return a data.table with
# sample_id, cell_id, UMAP1/2, and all RNA_snn_res.* cluster columns.
.load_clusters <- function(cls_f) {
  cls_dt      = cls_f %>% fread(na.strings = "") %>% .[ !is.na(UMAP1) ]
  cl_cols     = colnames(cls_dt) %>% str_subset("RNA_snn_res")
  cls_dt      = cls_dt %>%
    .[, c("sample_id", "cell_id", "UMAP1", "UMAP2", cl_cols), with = FALSE]

  return(cls_dt)
}



# Smooth per-cell XGBoost predictions by assigning cluster-level majority labels.
# For each (resolution, cluster) pair the dominant prediction is used as the
# cluster label if it exceeds min_cl_prop and the cluster has >= min_cl_size cells;
# clusters that are too small are excluded from smoothing.
# Returns a wide data.table (one row per cell) with cl_int, cl_pred, and prop_pred
# columns for each clustering resolution.
.apply_labels_by_cluster <- function(int_dt, preds_dt, min_cl_prop, min_cl_size) {
  # melt clusters
  non_cl_vars = c("sample_id", "cell_id", "UMAP1", "UMAP2")
  int_cls    = int_dt %>%
    melt.data.table( id = non_cl_vars, var = "res_int", val = "cl_int")

  # exclude tiny clusters
  int_ns     = int_cls[, .(N_cl = .N), by = .(res_int, cl_int) ]
  keep_cls    = int_ns[ N_cl >= min_cl_size ]
  if ( nrow(keep_cls) > 0 ) {
    message("  excluding some clusters bc they are tiny:")
    int_ns[ N_cl < min_cl_size ] %>% .[ order(res_int, cl_int) ] %>% print
    int_cls  = int_cls %>% merge(int_ns, by = c("res_int", "cl_int")) %>%
      .[, N_cl := NULL ]
  }

  # match these up to predictions, calculate proportions for each cluster
  match_dt    = preds_dt[, .(cell_id, cl_pred = cl_pred_naive)] %>%
    merge(int_cls, by = "cell_id") %>%
    .[, .N,                 by = .(res_int, cl_int, cl_pred)] %>%
    .[, prop := N / sum(N), by = .(res_int, cl_int) ] %>%
    setorder(res_int, cl_int, -prop)

  # take top prediction for each cluster
  match_lu    = match_dt[, .SD[1], by = .(res_int, cl_int)] %>%
    .[ (cl_pred != "unknown") & (prop > min_cl_prop) ]


  # add these results to original cluster labels
  guesses_dt  = match_lu[, .(res_int, cl_int, cl_pred, prop_pred = prop)] %>%
    merge(int_cls, by = c("res_int", "cl_int"), all.y = TRUE) %>%
    merge(preds_dt[, .(cell_id, cl_pred_raw, cl_pred_naive, p_pred)], by = "cell_id") %>%
    setcolorder( c(non_cl_vars, "cl_pred_raw", "cl_pred_naive", "p_pred") ) %>%
    dcast.data.table( sample_id + cell_id + UMAP1 + UMAP2 +
                        cl_pred_raw + cl_pred_naive + p_pred ~ res_int,
                      value.var = c("cl_int", "cl_pred", "prop_pred") )

  return(guesses_dt)
}






