suppressPackageStartupMessages({
  library("Matrix")
  library("SingleCellExperiment")
  library("BiocParallel")
  RhpcBLASctl::omp_set_num_threads(1L)
  library("fgsea")
})

# define pathways
gsea_regex  = "^(HALLMARK_|GOBP_|GOCC_|GOMF_|BIOCARTA_|REACTOME_|KEGG_)(.+)"

run_fgsea <- function(mkrs_f, fgsea_go_bp_f, fgsea_go_cc_f, fgsea_go_mf_f,
  species, gsea_dir, min_cpm_go, max_zero_p, gsea_cut, not_ok_re, n_cores){

  # check some inputs
  assert_that(
    is.character(mkrs_f),
    is.character(gsea_dir),
    is.character(not_ok_re), 
    is.character(fgsea_go_bp_f), 
    is.character(fgsea_go_cc_f),
    is.character(fgsea_go_mf_f), 
    is.character(not_ok_re)
  )
  assert_that(
    dir.exists(gsea_dir)
  )
  assert_that(
    is.numeric(max_zero_p),
    is.numeric(gsea_cut),
    is.numeric(n_cores)
  )
  setDTthreads(n_cores)
  
  # get markers
  mkrs_dt = fread(mkrs_f)

  message("  running FGSEA")
  fgsea_fs    = c(fgsea_go_bp_f, fgsea_go_cc_f, fgsea_go_mf_f)
  names(fgsea_fs) = str_extract(fgsea_fs, "(go_bp|go_cc|go_mf)")

  # get genesets, check they match expected files
  gsets_list  = get_fgsea_genesets(gsea_dir, species)
  assert_that( all(names(fgsea_fs) == names(gsets_list)) )

  # restrict slightly
  mkrs_tmp    = mkrs_dt[ str_detect(gene_type, not_ok_re, negate = TRUE) ] %>%
    .[ logcpm.sel > log(min_cpm_go + 1) ] %>%
    .[ (n_zero/n_cl < max_zero_p) | (logFC < 0) ]
  calc_fgsea_dt(gsets_list, fgsea_fs, mkrs_tmp, gsea_cut, gsea_var = "z_score",
    n_cores = n_cores)
}


get_fgsea_genesets <- function(gsea_dir, species) {
  if (species %in% c("human_2024","human_2020")) {
    gsets_list  = list(
      go_bp   = "c5.go.bp.v2023.1.Hs.symbols.gmt",
      go_cc   = "c5.go.cc.v2023.1.Hs.symbols.gmt",
      go_mf   = "c5.go.mf.v2023.1.Hs.symbols.gmt"
    ) %>% lapply(function(p) file.path(gsea_dir, p))
  } else if (species %in% c("mouse_2024", "mouse_2020")) {
    gsets_list  = list(
      go_bp   = "m5.go.bp.v2023.1.Mm.symbols.gmt",
      go_cc   = "m5.go.cc.v2023.1.Mm.symbols.gmt",
      go_mf   = "m5.go.mf.v2023.1.Mm.symbols.gmt"
    ) %>% lapply(function(p) file.path(gsea_dir, p))
  }
  assert_that( all(file.exists(unlist(gsets_list))) )

  path_names  = c(
    "GO BP",
    "GO CC",
    "GO MF"
  ) %>% setNames(names(gsets_list))
  return(gsets_list)
}

calc_fgsea_dt <- function(gsets_list, fgsea_fs, markers_dt, gsea_cut,
  gsea_var = c("z_score", "logFC"), n_cores = 4, overwrite = FALSE) {
  gsea_var    = match.arg(gsea_var)

  # split muscat results up
  message('running FGSEA')
  if (gsea_var == "z_score") {
    # fix any 0s
    min_fdr   = markers_dt[ FDR > 0 ]$FDR %>% min
    mkrs_tmp  = copy(markers_dt) %>%
      .[, .(cluster, symbol, logFC, log10_FDR = log10(FDR)) ] %>%
      .[ log10_FDR == -Inf, log10_FDR := log10(min_fdr) - 1 ]
    dt_list   = mkrs_tmp %>%
      .[, .(cluster, symbol, gsea_var = gsea_var,
        gsea_val = log10_FDR * sign(logFC) * -1, decreases = TRUE)] %>%
      split( by = 'cluster', drop = TRUE )
  } else if (gsea_var == "logFC") {
    dt_list   = markers_dt %>%
      .[, .(cluster, symbol, gsea_var = gsea_var,
        gsea_val = logFC, decreases = TRUE)] %>%
      split( by = 'cluster', drop = TRUE )
  }

  # set up cluster
  bpparam     = MulticoreParam(workers = n_cores, tasks = length(dt_list))
  on.exit(bpstop(bpparam))

  # run for each set of paths
  for (p in names(gsets_list)) {
    # check if already done
    fgsea_f   = fgsea_fs[[p]]

    # get pathways
    message('  ', p, appendLF = FALSE)
    paths_f   = gsets_list[[p]]
    pathways  = gmtPathways(paths_f)

    # set up parallel
    fgsea_dt  = bplapply(seq_along(dt_list), function(i) {
     
      # get labels
      dt      = dt_list[[i]]
      spec_dt = dt[1][, .(cluster, gsea_var)]
      tmp_dt  = .calc_enrichment(dt, pathways, gsea_cut) %>%
        cbind(spec_dt) %>% .[, path_set := p]
      tmp_dt  = tmp_dt[ !is.na(pval) ] %>%
        .[ !is.null(leadingEdge) ]

      return(tmp_dt)
      }, BPPARAM = bpparam) %>% rbindlist

    # save results
    fwrite(fgsea_dt, file = fgsea_f)
  }
  bpstop()
}

.calc_enrichment <- function(dt, pathways, gsea_cut, min_size = 5, max_size = 500) {
  # get ordering
  decreases   = unique(dt$decreases)
  assert_that( length(decreases) == 1 )

  # sort by specified variable
  ranks_vec = dt$gsea_val %>%
    setNames(dt$symbol) %>% sort( decreasing = decreases )

  # sort out any duplicated symbols, and miRNA names
  dupes       = names(ranks_vec) %>% table %>% .[ . > 1 ] %>% names
  ranks_vec   = ranks_vec[ !(names(ranks_vec) %in% dupes) ]

  # sort out any miRNA names
  mir_idx     = names(ranks_vec) %>% str_detect('^MIR.+HG$')
  if (sum(mir_idx) > 0) {
    mir_names   = names(ranks_vec)[mir_idx]
    names(ranks_vec)[mir_idx] = mir_names %>% str_replace_all('HG$', '')
  }

  # calculate fgsea
  fgsea_dt  = fgseaMultilevel(
    pathways  = pathways,
    stats     = ranks_vec,
    minSize   = min_size,
    maxSize   = max_size,
    nproc     = 1,
    scoreType = 'std'
    ) %>% .[, main_path := FALSE]

  # calculate collapsed pathways
  if (is.null(fgsea_dt) | nrow(fgsea_dt) == 0)
    return(fgsea_dt)
  fgsea_sig = fgsea_dt[ padj < gsea_cut ]

  # collapse pathways
  if (nrow(fgsea_sig) > 0) {
    paths_col   = collapsePathways(fgsea_sig[order(pval)], pathways, ranks_vec,
      pval.threshold = gsea_cut)
    main_paths  = paths_col$mainPathways
  } else {
    main_paths  = NULL
  }

  # add as label
  fgsea_dt[, main_path := pathway %in% main_paths ]

  return(fgsea_dt)
}

load_gsea_ls <- function(gsets_list, fgsea_pat) {
  gsea_ls   = lapply(names(gsets_list), function(p)
    sprintf(fgsea_pat, p) %>% fread) %>% setNames(names(gsets_list))
}


plot_gsea_dotplot <- function(gsea_dt, n_top_paths = 10, gsea_cut = 0.05,
  maxp_cut = 0.1, n_chars = 50, max_nes = Inf, min_log10_padj = -10,
  size_range = c(0, Inf), what = c('both', 'pos_only'), cl_order = NULL) {

  # check inputs
  what      = match.arg(what)
  if (!is.null(cl_order)) {
    assert_that( all(sort(unique(gsea_dt$cluster)) == sort(cl_order)) )
  }

  # which terms?
  top_dt    = copy(gsea_dt) %>%
    .[ main_path == TRUE ]
  if (what == 'pos_only') {
    top_dt    = top_dt[ NES > 0 ]
  }
  top_dt    = top_dt %>%
    .[, min_p := min(padj, na.rm = TRUE), by = pathway ] %>%
    .[ min_p < maxp_cut ] %>%
    .[ size %between% size_range ]
  if (nrow(top_dt) == 0)
    return(NULL)
  top_dt    = top_dt %>%
    .[ order(cluster, padj) ] %>%
    .[, p_rank := 1:.N, by = cluster ]
  top_paths   = top_dt[ p_rank <= n_top_paths ]$pathway %>% unique

  # put p-values into wide format
  plot_dt   = gsea_dt[ pathway %in% top_paths ] %>%
    .[, cluster  := factor(cluster) ] %>%
    .[ is.na(padj), padj := 1 ] %>%
    .[, log_p       := log10(padj) ] %>%
    .[, pathway     := str_match(pathway, gsea_regex)[, 3] %>%
      tolower %>% str_replace_all("_", " ")] %>%
    .[, path_short  := pathway %>% str_sub(1, n_chars) ] %>%
    .[, signif      := ifelse(padj < gsea_cut, 'significant', 'not') ]
  assert_that(
    !any(is.na(plot_dt$pathway)),
    !any(is.na(plot_dt$path_short))
  )

  if (!is.null(cl_order))
    plot_dt    = plot_dt %>% .[, cluster := factor(cluster, levels = cl_order) ]

  # put in nice order
  if (what == 'pos_only') {
    tmp_dt    = copy(plot_dt) %>% .[ NES > 0 ]
  } else if (what == 'both') {
    tmp_dt    = copy(plot_dt)
  }
  max_nes_dt  = tmp_dt[, .(
    min_p       = log10(.SD[min(padj) == padj]$padj[1]) * sign(.SD[min(padj) == padj]$NES[1]),
    cluster     = .SD[min(padj) == padj]$cluster[1]
  ), by = path_short] %>% setorder(-cluster, -min_p)
  assert_that( nrow(max_nes_dt) == length(unique(max_nes_dt$path_short)))
  plot_dt[, path_short := factor(path_short, levels = max_nes_dt$path_short) ]

  # what limits?
  res         = 0.5
  max_nes     = min(ceiling(max(abs(plot_dt$NES)) * res) / res, max_nes)
  plot_dt     = plot_dt %>%
    .[, nes_trunc         := sign(NES) * pmin(abs(NES), max_nes) ] %>%
    .[, log10_padj_trunc  := pmax(log10(padj), min_log10_padj) ]

  # do dotplot
  g = ggplot(plot_dt[ order(NES) ]) +
    aes(x = cluster, y = path_short, fill = nes_trunc, size = -log10_padj_trunc,
      alpha = signif ) +
    geom_point(shape = 21, colour = 'black')

  g = g + scale_fill_distiller( palette = "RdBu", limits = c(-max_nes, max_nes),
    breaks = pretty_breaks() )
  g = g + scale_alpha_manual(values = c(signif = 1, not = 0.5)) +
    scale_size( range = c(1, 6), breaks = pretty_breaks() ) +
    # facet_grid( . ~ cluster ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid  = element_blank()
      ) +
    labs(
      x     = 'Cluster',
      y     = 'Pathway',
      fill  = "Normalized\nenrichment\nscore",
      size  = "-log10(adjusted p)",
      alpha = sprintf('Significant?\n(at %dpc)', round(100*gsea_cut))
    )

  return(g)
}