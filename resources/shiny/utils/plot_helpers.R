# SHARED PLOT HELPER FUNCTIONS -------------------------------------------------

# Return a blank ggplot (used when no data has been selected yet).
blank_plot <- function() {
  ggplot(data.frame(x = 1, y = 1), aes(x = x, y = y)) +
    geom_blank() +
    theme_classic() +
    theme(
      axis.text  = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.line  = element_blank()
    )
}

# Shared theme for pseudobulk dot-plots (gene and composition tabs).
theme_scjoin_dots <- function() {
  theme_classic() +
    theme(
      axis.text.x    = element_text(angle = -45, hjust = 0, vjust = 0.5, size = FONT_TEXT),
      axis.text.y    = element_text(size = FONT_SMALL),
      axis.title     = element_text(size = FONT_AXIS),
      legend.text    = element_text(size = FONT_SMALL),
      legend.box     = 'horizontal',
      plot.subtitle  = element_text(color = 'grey60', size = FONT_TEXT, hjust = 1),
      legend.title   = element_text(size = FONT_AXIS),
      plot.title     = element_text(size = FONT_TITLE),
      strip.text     = element_text(size = FONT_AXIS)
    )
}

# Build a labelled FDR annotation data.table for a single gene symbol.
# Returns a data.table(cluster, label) with arrow + FDR-value labels.
make_fdr_anno_dt <- function(sel_symbol, cluster_markers) {

  fdr_anno_dt = copy(cluster_markers) %>%
    .[ symbol == sel_symbol ] %>%
    .[, FDR_chr := sprintf("%.1e", FDR)] %>%
    .[, label   := fcase(
      sign(log2fc) ==  1 & FDR <= 0.05, paste(sprintf('\u2191'), FDR_chr),
      sign(log2fc) == -1 & FDR <= 0.05, paste(sprintf('\u2193'), FDR_chr),
      FDR > 0.05,                        paste('\u2248', FDR_chr)
    )] %>%
    .[, .(cluster, label)] %>%
    setorder(cluster)

  return(fdr_anno_dt)
}

# Render a list of pickerInput widgets for subsetting, one per variable.
# prefix: string prepended to the variable value to form the inputId
#         (e.g. "pb_subset_" -> inputId "pb_subset_cluster")
# vars:   named character vector (name = display label, value = column name)
# cluster_meta: data.table with factor columns matching vars values
make_subset_picker_ui <- function(prefix, vars, cluster_meta, ns = identity) {
  lapply(seq_along(vars), function(ii) {
    tags$div(
      class = "custom-picker",
      pickerInput(
        inputId  = ns(paste0(prefix, vars[[ii]])),
        label    = names(vars)[[ii]],
        multiple = TRUE,
        choices  = levels(cluster_meta[[ vars[[ii]] ]]),
        options  = pickerOptions(actionsBox = TRUE),
        selected = NULL
      )
    )
  })
}
