# SETUP ------------------------------------------------------------------------

suppressPackageStartupMessages({
  library('data.table')
  library('tibble')
  library('strex')
  library('yaml')
  library('markdown')
  library('shiny')
  library('shinyWidgets')
  library('shinydashboardPlus')
  library('shinydashboard')
  library('shinyjs')
  library('shinyBS')
  library('DT')
  library('esquisse')
  library('assertthat')
  library('forcats')
  library('stringr')
  library('ggplot2')
  library('ggrepel')
  library('scattermore')
  library('ggbeeswarm')
  library('scales')
  library('RColorBrewer')
  library('MetBrewer')
  library('ggsci')
  library('viridis')
  library('patchwork')
  library('ComplexHeatmap')
  library('circlize')
  library('BPCells')
})

# Source helpers and modules
source('constants.R')
source('utils/colors.R')
source('utils/data_io.R')
source('utils/metadata.R')
source('utils/plot_helpers.R')
source('utils/plots_genes.R')
source('utils/plots_clusters.R')
source('utils/plots_metadata.R')
source('modules/explore_genes.R')
source('modules/explore_clusters.R')
source('modules/explore_genesets.R')
source('modules/explore_prevalence.R')

www_dir = 'www'

# CONFIGURATION ----------------------------------------------------------------

yaml_data     = yaml::yaml.load_file('shinyconfig.yaml')
data_dir      = yaml_data$data_dir

date_stamp    = yaml_data$date_stamp
app_title     = yaml_data$app$app_title
email         = yaml_data$app$email
keyword       = yaml_data$app$keyword
default_gene  = yaml_data$app$default_gene

include_paga  = yaml_data$build$include_paga
logo_f        = yaml_data$build$logo_f
gsets_f       = yaml_data$build$gsets_f
sample_col    = yaml_data$build$sample_col

metadata_vars        = unlist(yaml_data$metadata$vars)
names(metadata_vars) = unlist(yaml_data$metadata$var_names)
metadata_var_combns  = yaml_data$metadata$var_combns
subset_vars          = c('cluster', metadata_vars)
names(subset_vars)   = c('cluster', names(metadata_vars))

# optional data sources table
data_sources_f = file.path(data_dir, "data_sources.csv")
if (file.exists(data_sources_f)) ds_dt = fread(data_sources_f)

# logo path
if (!is.null(logo_f))
  logo_f = file.path(www_dir, strex::str_after_last(logo_f, '\\/'))

# predefined genesets
if (!is.null(gsets_f)) {
  gsets_f      = file.path(data_dir, strex::str_after_last(gsets_f, '\\/'))
  gsets_dt     = fread(gsets_f)
  predef_gsets = 1L
} else {
  gsets_dt     = NULL
  predef_gsets = 0L
}

# PAGA
if (!is.null(include_paga)) {
  paga = 1L
  paga_weight_cut = if (!is.null(include_paga$paga_weight_cut))
    as.numeric(include_paga$paga_weight_cut)
  else
    0
} else {
  paga            = 0L
  paga_weight_cut = 0
}

# DATA LOADING -----------------------------------------------------------------

fs_ls = get_all_input_fs(data_dir, date_stamp, paga)
list2env(fs_ls, envir = .GlobalEnv)

tmp_ls         = get_lvls_and_colours(yaml_data, cluster_meta, sample_meta,
                                      metadata_vars, data_dir, nice_cols)
vars_lvls      = tmp_ls$vars_lvls
vars_pals      = tmp_ls$vars_pals
cluster_labels = tmp_ls$cluster_labels

# apply factor levels to metadata tables
for (v in metadata_vars) {
  cluster_meta[[v]] = factor(cluster_meta[[v]], levels = vars_lvls[[v]])
  sample_meta[[v]]  = factor(sample_meta[[v]],  levels = vars_lvls[[v]])
  cell_meta[[v]]    = factor(cell_meta[[v]],     levels = vars_lvls[[v]])
}

# add variable combination columns
if (!is.null(metadata_var_combns)) {
  for (vars in metadata_var_combns) {
    cluster_meta = add_metadata_var_combns(cluster_meta, vars, vars_lvls[vars])
    sample_meta  = add_metadata_var_combns(sample_meta,  vars, vars_lvls[vars])
    cell_meta    = add_metadata_var_combns(cell_meta,    vars, vars_lvls[vars])

    combn_vars = c(paste(vars, collapse = '.'), paste(rev(vars), collapse = '.')) %>%
      setNames(c(
        paste(names(metadata_vars[match(vars,      metadata_vars)]), collapse = ' + '),
        paste(names(metadata_vars[match(rev(vars), metadata_vars)]), collapse = ' + ')
      ))
    metadata_vars = c(metadata_vars, combn_vars)

    new_lvls  = lapply(combn_vars, function(v) sample_meta[[v]] %>% levels) %>%
      setNames(combn_vars)
    new_pals  = lapply(new_lvls, function(l)
      rep(nice_cols, times = 20)[ seq_along(l) ] %>% setNames(l)) %>%
      setNames(combn_vars)
    vars_lvls = c(vars_lvls, new_lvls)
    vars_pals = c(vars_pals, new_pals)
  }
}

# cluster labels / factor levels
cluster_lvls = sort(unique(cluster_meta$cluster))
if (!is.null(cluster_labels)) {
  cluster_meta$cluster    = cluster_labels[cluster_meta$cluster]
  cell_meta$cluster       = cluster_labels[cell_meta$cluster]
  cluster_markers$cluster = cluster_labels[cluster_markers$cluster]
  fgsea$cluster           = cluster_labels[fgsea$cluster]
  centroids$cluster       = cluster_labels[centroids$cluster]
  repel_pos$cluster       = cluster_labels[repel_pos$cluster]
  if (paga == 1) {
    paga_mat$cluster = cluster_labels[paga_mat$cluster]
    setnames(paga_mat, cluster_lvls, cluster_labels[cluster_lvls])
    paga_pos$cluster = cluster_labels[paga_pos$cluster]
  }
  cluster_lvls = cluster_labels
}

cluster_meta$cluster    = factor(cluster_meta$cluster,    levels = cluster_lvls)
cell_meta$cluster       = factor(cell_meta$cluster,       levels = cluster_lvls)
cluster_markers$cluster = factor(cluster_markers$cluster, levels = cluster_lvls)
fgsea$cluster           = factor(fgsea$cluster,           levels = cluster_lvls)
centroids$cluster       = factor(centroids$cluster,       levels = cluster_lvls)
repel_pos$cluster       = factor(repel_pos$cluster,       levels = cluster_lvls)
if (paga == 1) {
  paga_mat$cluster = factor(paga_mat$cluster, levels = cluster_lvls)
  paga_pos$cluster = factor(paga_pos$cluster, levels = cluster_lvls)
}

# merge GO term descriptions into fgsea
fgsea = merge(fgsea, go_terms[, .(pathway, pathway_nice, path_short)], by = 'pathway')

# heatmap column width (used by explore_genesets module)
heatmap_width = (length(cluster_lvls) * 0.25 + 3) %>% ceiling %>% max(1) %>% min(10)

# convenience: reverse lookup for display names of metadata variables
var_choice_names = setNames(names(metadata_vars), metadata_vars) %>%
  c(., 'cluster' = 'cluster')

# USER INTERFACE ---------------------------------------------------------------

header <- dashboardHeader(
  title = app_title,
  tags$li(class = "dropdown",
    tags$style(".main-header {height: 55px}"),
    tags$style(".main-header .logo {height:55px}")),
  tags$li(img(src = logo_f, height = "45px"),
    style = "padding-top:5px; padding-bottom:5px; padding-right:10px",
    class = "dropdown")
)

sidebar <- dashboardSidebar(
  collapsed = TRUE,
  sidebarMenu(id = "tabs",
    menuItem("Home",               tabName = "description",        icon = icon("house",       verify_fa = FALSE)),
    menuItem("Explore genes",      tabName = "explore_gene",       icon = icon("dna",         verify_fa = FALSE)),
    menuItem("Explore clusters",   tabName = "explore_clusters",   icon = icon("sitemap",     verify_fa = FALSE)),
    menuItem("Explore genesets",   tabName = "explore_gsets",      icon = icon("layer-group", verify_fa = FALSE)),
    menuItem("Explore prevalence", tabName = "explore_prevalence", icon = icon("percent",     verify_fa = FALSE))
  )
)

body <- dashboardBody(
  useSweetAlert(),
  useShinyjs(),
  tags$style(HTML(".redBox { background-color: #ff0000; color: #ffffff; }")),
  tags$style(HTML("
    /* Header */
    .main-header .navbar,
    .main-header .navbar .sidebar-toggle { background-color: #724DB8 !important; }
    .main-header .logo                   { background-color: #724DB8 !important; border-bottom: 0 !important; }
    .main-header .logo:hover             { background-color: #543988 !important; }
    /* Sidebar */
    .main-sidebar, .left-side            { background-color: #B6A9DF !important; }
    .sidebar-menu > li > a               { color: #3D3756 !important; }
    .sidebar-menu > li.active > a        { background-color: rgba(61,55,86,0.12) !important;
                                           border-left-color: #C9A561 !important;
                                           color: #3D3756 !important; }
    .sidebar-menu > li > a:hover         { background-color: rgba(61,55,86,0.08) !important; }
    /* Body background */
    .content-wrapper, .right-side        { background-color: #F0F0F0 !important; }
    /* infoBox icon panels and link text */
    .info-box .info-box-icon.bg-purple   { background-color: #724DB8 !important; }
    .info-box .info-box-content a        { color: #724DB8 !important; }
    /* Buttons */
    .btn-primary                         { background-color: #724DB8 !important; border-color: #724DB8 !important; }
    .btn-primary:hover,
    .btn-primary:focus,
    .btn-primary:active                  { background-color: #543988 !important; border-color: #543988 !important; }
  ")),
  tags$script(HTML("
    var openTab = function(tabName) {
      $('a', $('.sidebar')).each(function() {
        if (this.getAttribute('data-value') == tabName) { this.click() };
      });
    }
  ")),
  tags$style(HTML("
    .custom-picker .bootstrap-select .dropdown-menu li a { font-size: 12px; font-weight: normal; }
    .custom-picker .bootstrap-select .btn               { font-size: 12px; font-weight: normal; }
  ")),

  tabItems(

    # HOME ---------------------------------------------------------------------
    tabItem(tabName = "description",
      div(h2(strong(app_title)), style = "text-align: center; color:#724DB8"),
      h3("What would you like to do?"),
      fluidRow(
        infoBox("", a("Explore genes",      onclick = "openTab('explore_gene')",       href="#"), icon = icon("dna",         verify_fa = FALSE), color = "purple", width = 3),
        infoBox("", a("Explore clusters",   onclick = "openTab('explore_clusters')",   href="#"), icon = icon("sitemap",     verify_fa = FALSE), color = "purple", width = 3),
        infoBox("", a("Explore genesets",   onclick = "openTab('explore_gsets')",      href="#"), icon = icon("layer-group", verify_fa = FALSE), color = "purple", width = 3),
        infoBox("", a("Explore prevalence", onclick = "openTab('explore_prevalence')", href="#"), icon = icon("percent",     verify_fa = FALSE), color = "purple", width = 3)
      ),
      tags$style(".markdown-content { font-size: 17px; text-align: left; }"),
      div(class = "markdown-content", includeMarkdown("data/home.md"))
    ),

    # EXPLORE GENES ------------------------------------------------------------
    tabItem(tabName = "explore_gene",       genesUI("genes")),

    # EXPLORE CLUSTERS ---------------------------------------------------------
    tabItem(tabName = "explore_clusters",   clustersUI("clusters")),

    # EXPLORE GENESETS ---------------------------------------------------------
    tabItem(tabName = "explore_gsets",      genesetsUI("genesets")),

    # EXPLORE PREVALENCE -------------------------------------------------------
    tabItem(tabName = "explore_prevalence", prevalenceUI("prevalence"))
  )
)

footer <- shinydashboardPlus::dashboardFooter(
  right = if (!is.null(email) && nchar(email) > 0)
    div(style = "font-size: 18px;",
      "Contact: ",
      a(href = paste0("mailto:", email), email, style = "margin-right: 5px;")
    ),
  left = div(style = "color: #3D3756; font-size: 18px;",
    "This app was generated with ",
    a(href = "https://marusakod.github.io/scprocess/", strong("scprocess"),
      style = "color: #9785D0;"), " :)"
  )
)

ui <- function() {
  if (!dir.exists(www_dir)) dir.create(www_dir)
  addResourcePath("www", www_dir)
  tagList(shinydashboardPlus::dashboardPage(
    skin = 'purple', header = header, sidebar = sidebar,
    body = body, footer = footer
  ))
}

# SERVER -----------------------------------------------------------------------

server <- function(input, output, session) {

  # UMAP data for cluster overview (rendered reactively in each module)
  cluster_umap_dt = copy(cell_meta) %>%
    .[keep_cell == TRUE, .(UMAP_1, UMAP_2, cluster)]

  # shared data passed to all modules
  shared <- list(
    cluster_meta     = cluster_meta,
    cell_meta        = cell_meta,
    sample_meta      = sample_meta,
    cluster_markers  = cluster_markers,
    fgsea            = fgsea,
    go_terms         = go_terms,
    pb_hvgs          = pb_hvgs,
    row_indx         = row_indx,
    gsets_dt         = gsets_dt,
    count_h5_f       = count_h5_f,
    cluster_pb_h5_f  = cluster_pb_h5_f,
    sample_pb_h5_f   = sample_pb_h5_f,
    keyword          = keyword,
    default_gene     = default_gene,
    sample_col       = sample_col,
    metadata_vars    = metadata_vars,
    subset_vars      = subset_vars,
    vars_pals        = vars_pals,
    vars_lvls        = vars_lvls,
    var_choice_names = var_choice_names,
    paga             = paga,
    paga_weight_cut  = paga_weight_cut,
    predef_gsets     = predef_gsets,
    heatmap_width    = heatmap_width,
    cluster_umap_dt  = cluster_umap_dt,
    centroids        = centroids,
    repel_pos        = repel_pos,
    current_tab      = reactive(input$tabs)
  )
  if (paga == 1) {
    shared$paga_mat = paga_mat
    shared$paga_pos = paga_pos
  }

  # modal dialog for cluster annotation comments
  if (!is.null(cluster_labels)) {
    observeEvent(input$comment_button, {
      showModal(modalDialog(
        title = paste0("If you don't agree with our cluster annotations, ",
                       "please write an e-mail to ", email,
                       " and tell us your rationale. Thank you for collaboration!"),
        footer = tagList(modalButton("Cancel"))
      ))
    })
  }

  # reset all inputs when switching tabs
  observeEvent(input$tabs, { shinyjs::reset() })

  # register modules
  genesServer("genes",           shared)
  clustersServer("clusters",     shared)
  genesetsServer("genesets",     shared)
  prevalenceServer("prevalence", shared)
}

shinyApp(ui = ui, server = server)
