# MODULE: Explore Genes --------------------------------------------------------

genesUI <- function(id) {
  ns = NS(id)
  tagList(
    column(width = 2,
      fluidRow(
        box(width = 0, title = NULL, headerBorder = FALSE,
          awesomeRadio(
            inputId  = ns("search_gene"),
            label    = div(style = "font-size: 14px", "Search gene by:"),
            choices  = c("Gene Symbol" = 1, "Ensembl ID" = 2),
            selected = 1,
            inline   = FALSE
          ),
          conditionalPanel(
            sprintf("input['%s'] == 1", ns("search_gene")),
            selectizeInput(
              inputId = ns("Entry_gene"),
              label   = div(style = "font-size: 14px", "Gene symbol"),
              multiple = FALSE, size = 20, choices = NULL,
              options = list(
                placeholder  = 'e.g. gene',
                onInitialize = I('function() { this.setValue("");}')
              )
            )
          ),
          conditionalPanel(
            sprintf("input['%s'] == 2", ns("search_gene")),
            selectizeInput(
              inputId = ns("Entry_ID"),
              label   = div(style = "font-size: 14px", "Ensembl gene ID"),
              multiple = FALSE, size = 20, choices = NULL,
              options = list(
                placeholder  = 'e.g. ENSG...',
                onInitialize = I('function() { this.setValue(""); }')
              )
            )
          ),
          palettePicker(
            inputId   = ns("UMAP_pal"),
            label     = "UMAP palette",
            choices   = lapply(palette_map, '[[', 'cols'),
            textColor = rep(rgb(0, 0, 0, alpha = 0), 5)
          ),
          prettySwitch(
            inputId = ns("gene_repel_labels"),
            label   = "Repel cluster labels",
            value   = FALSE
          ),
          radioButtons(
            inputId  = ns("pb_dots_choice"),
            label    = "Dots are:",
            choices  = c("clusters", "samples"),
            selected = "clusters"
          ),
          uiOutput(ns("pb_x_axis_ui")),
          uiOutput(ns("pb_fill_ui")),
          p(strong("Show advanced plotting options")),
          prettySwitch(inputId = ns("gene_show_advanced"), label = '', value = FALSE),
          conditionalPanel(
            sprintf("input['%s'] == true", ns("gene_show_advanced")),
            uiOutput(ns("pb_facet_ui")),
            prettySwitch(
              inputId = ns("gene_order_plot"),
              label   = "Order by mean expression",
              value   = FALSE
            ),
            p(strong("Subset data on...")),
            uiOutput(ns("gene_subset_picker"))
          )
        )
      )
    ),
    column(width = 10,
      fluidRow(
        column(width = 6, uiOutput(ns("get_gene_umap_box"))),
        column(width = 6, uiOutput(ns("get_cluster_overview_box"))),
        column(width = 12, align = "center", uiOutput(ns("get_gene_dotplot_box")))
      )
    )
  )
}

genesServer <- function(id, shared) {
  moduleServer(id, function(input, output, session) {
    ns = session$ns

    # populate server-side selectize choices
    .init_entry_gene <- function()
      updateSelectizeInput(session, "Entry_gene",
        server  = TRUE, selected = shared$default_gene,
        options = list(placeholder = sprintf('e.g. %s', shared$default_gene)),
        choices = shared$row_indx$symbol)

    .init_entry_gene()

    # re-sync selectize after the global shinyjs::reset() that fires on tab switch
    observeEvent(shared$current_tab(), {
      req(shared$current_tab() == "explore_gene")
      .init_entry_gene()
    }, ignoreInit = TRUE)

    updateSelectizeInput(session, "Entry_ID",
      server   = TRUE, selected = NULL,
      options  = list(placeholder = 'e.g. ENSG...',
                      onInitialize = I('function() { this.setValue("");}')),
      choices  = shared$row_indx$ensembl)

    # compute axis/facet choices from current dots type
    dot_axis_choices <- reactive({
      if (input$pb_dots_choice == "clusters") {
        list(
          x_choices   = c("cluster" = "cluster", shared$metadata_vars),
          fct_choices = c("none", shared$metadata_vars, "cluster"),
          x_sel       = "cluster"
        )
      } else {
        list(
          x_choices   = shared$metadata_vars,
          fct_choices = c("none", shared$metadata_vars),
          x_sel       = shared$metadata_vars[[1]]
        )
      }
    })

    output$pb_x_axis_ui <- renderUI({
      ch = dot_axis_choices()
      selectizeInput(ns("pb_x_axis"),
        label    = div(style = "font-size: 14px", "x-axis variable"),
        choices  = ch$x_choices, selected = ch$x_sel, multiple = FALSE, size = 20)
    })

    output$pb_fill_ui <- renderUI({
      ch = dot_axis_choices()
      selectizeInput(ns("pb_fill"),
        label    = div(style = "font-size: 14px", "Colour"),
        choices  = ch$x_choices, selected = ch$x_sel, multiple = FALSE, size = 20)
    })

    output$pb_facet_ui <- renderUI({
      ch = dot_axis_choices()
      selectizeInput(ns("pb_facet"),
        label    = div(style = "font-size: 12px", "Facet by"),
        choices  = ch$fct_choices, selected = "none", multiple = FALSE)
    })

    # reset gene search inputs when search-by radio changes
    observeEvent(input$search_gene, {
      shinyjs::reset(ns("Entry_gene"))
      shinyjs::reset(ns("Entry_ID"))
    })

    # track most recently entered gene (symbol or ensembl).
    # initialised with default_gene so plots render immediately on first visit.
    # guards against empty strings from shinyjs::reset() on tab switch.
    Entry_variable <- reactiveVal(shared$default_gene)
    observeEvent(input$Entry_gene, {
      if (nchar(input$Entry_gene) > 0) Entry_variable(input$Entry_gene)
    })
    observeEvent(input$Entry_ID, {
      if (nchar(input$Entry_ID) > 0) Entry_variable(input$Entry_ID)
    })

    # row index for selected gene
    gene_idx <- reactive({
      req(Entry_variable(), nchar(Entry_variable()) > 0)
      if (input$search_gene == 1)
        shared$row_indx[symbol  == Entry_variable(), index]
      else
        shared$row_indx[ensembl == Entry_variable(), index]
    })

    # UMAP data with gene expression column
    gene_umap_dt <- reactive({
      dt = copy(shared$cell_meta) %>% .[keep_cell == TRUE, .(UMAP_1, UMAP_2)]
      if (length(gene_idx()) == 0) return(dt)
      dt$gene_exp = get_row_from_h5(shared$count_h5_f, gene_idx())
      setorder(dt, gene_exp)
    })

    # pseudobulk dotplot data
    gene_dotplot_dt <- reactive({
      req(length(gene_idx()) > 0)
      if (input$pb_dots_choice == 'clusters') {
        meta          = copy(shared$cluster_meta)
        meta$gene_exp = get_row_from_h5(shared$cluster_pb_h5_f, gene_idx())
        meta          = meta[n_cells >= 10, ]
      } else {
        meta          = copy(shared$sample_meta)
        meta$gene_exp = get_row_from_h5(shared$sample_pb_h5_f, gene_idx())
      }
      meta
    })

    # FDR annotations for dotplot (only when x-axis is cluster)
    gene_dotplot_anno_dt <- reactive({
      if (is.null(input$pb_x_axis) || input$pb_x_axis != 'cluster') return(NULL)
      sym = if (input$search_gene == 2)
        shared$row_indx[ensembl == Entry_variable(), symbol]
      else
        Entry_variable()
      make_fdr_anno_dt(sym, shared$cluster_markers)
    })

    # subset picker UI
    output$gene_subset_picker <- renderUI({
      make_subset_picker_ui("pb_subset_", shared$subset_vars, shared$cluster_meta, ns)
    })

    # UMAP box (only shown when a gene is selected)
    output$get_gene_umap_box <- renderUI({
      if (length(gene_idx()) != 0) {
        box(width = 0, title = NULL, height = "692px", headerBorder = FALSE,
            plotOutput(ns("gene_umap"), height = "680px"))
      }
    })

    output$gene_umap <- renderPlot({
      if (length(gene_idx()) == 0) return(blank_plot())
      plot_gene_umap(gene_umap_dt(), gene = Entry_variable(), pal_name = input$UMAP_pal)
    })

    # dotplot box (only shown when a gene is selected)
    output$get_gene_dotplot_box <- renderUI({
      if (length(gene_idx()) != 0) {
        box(width = 0, title = NULL, height = "500px", headerBorder = FALSE,
            plotOutput(ns("gene_dotplot"), height = "480px"))
      }
    })

    output$gene_dotplot <- renderPlot({
      if (length(gene_idx()) == 0) return(blank_plot())

      use_clusters = input$pb_dots_choice == 'clusters'
      default_x    = if (use_clusters) "cluster" else shared$metadata_vars[[1]]
      x_axis       = if (!is.null(input$pb_x_axis) && nchar(input$pb_x_axis) > 0) input$pb_x_axis else default_x
      fill         = if (!is.null(input$pb_fill)   && nchar(input$pb_fill)   > 0) input$pb_fill   else default_x

      if (use_clusters) {
        note     = sprintf('Dots represent %s clusters in each sample', shared$keyword)
        anno_dt  = gene_dotplot_anno_dt()
        sel_vars = shared$subset_vars
      } else {
        note     = sprintf('Dots represent %s in each sample', shared$keyword)
        anno_dt  = NULL
        sel_vars = setdiff(shared$subset_vars, 'cluster')
      }

      facet_by  = if (!is.null(input$pb_facet) && nchar(input$pb_facet) > 0) input$pb_facet else "none"

      subset_ls = sel_vars %>%
        lapply(function(var) input[[ paste0("pb_subset_", var) ]]) %>%
        setNames(sel_vars)

      plot_gene_dotplot(
        gene_dotplot_dt(),
        subset    = subset_ls,
        x_axis    = x_axis,
        fill      = fill,
        fill_name = shared$var_choice_names[ fill ],
        fill_pal  = shared$vars_pals[[ fill ]],
        facet_by  = facet_by,
        gene      = Entry_variable(),
        anno_dt   = anno_dt,
        note      = note,
        order     = input$gene_order_plot
      )
    })

    # cluster overview UMAP (shared pre-computed plot)
    output$get_cluster_overview_box <- renderUI({
      box(width = 0, title = NULL, height = "692px", headerBorder = FALSE,
          plotOutput(ns("cluster_overview_umap"), height = "680px"))
    })

    output$cluster_overview_umap <- renderPlot({
      label_mode = if (isTRUE(input$gene_repel_labels)) "repel" else "centroid"
      plot_cluster_umap(
        shared$cluster_umap_dt,
        col_pal    = shared$vars_pals[['cluster']],
        centroids  = shared$centroids,
        repel_pos  = shared$repel_pos,
        label_mode = label_mode
      )
    })
  })
}
