# MODULE: Explore Clusters -----------------------------------------------------

clustersUI <- function(id) {
  ns = NS(id)
  tagList(
    column(width = 2,
      fluidRow(
        box(width = 0, title = NULL, headerBorder = FALSE,
          uiOutput(ns("cl_select_cluster_ui")),
          uiOutput(ns('cl_display_paga_sw')),
          awesomeRadio(
            inputId  = ns("cl_gene_type_filt"),
            label    = div(style = "font-size: 14px", "Gene type"),
            choices  = c("all" = "all", "protein coding" = "protein_coding",
                         "excluding lncRNAs" = "no_lncRNA"),
            selected = "no_lncRNA"
          ),
          awesomeRadio(
            inputId  = ns("cl_cut_fdr"),
            label    = div(style = "font-size: 14px", "FDR threshold"),
            choices  = c("0.1" = 0.1, "0.05" = 0.05, "0.01" = 0.01,
                         "0.001" = 0.001, "10-4" = 0.0001, "10-5" = 0.00001),
            selected = 0.05
          ),
          sliderInput(
            inputId = ns("cl_cut_log2fc"),
            label   = div(style = "font-size: 14px", "log2fc threshold"),
            value   = 0, min = 0, max = 5, step = 0.5, ticks = TRUE
          ),
          awesomeRadio(
            inputId  = ns("cl_cut_cpm"),
            label    = div(style = "font-size: 14px", "CPM threshold"),
            choices  = c("0" = 0, "10" = 10, "20" = 20, "50" = 50, "100" = 100),
            selected = 50
          ),
          div(style = "text-align:center",
            downloadBttn(
              outputId = ns("cl_download_table"),
              label    = tags$span(id = "down_sel_cl", "Download table"),
              style    = "material-flat", color = "primary", size = "sm",
              icon     = icon("download")
            )
          ),
          tags$br(),
          prettySwitch(
            inputId = ns("cl_order_plot"),
            label   = "Order clusters by mean expression",
            value   = FALSE
          )
        )
      )
    ),
    column(width = 10,
      fluidRow(
        column(width = 6,
          box(width = 0, title = NULL, height = "710px", headerBorder = FALSE,
              plotOutput(ns("cl_umap_overview"), height = "680px"))
        ),
        column(width = 6,
          box(width = 0, title = NULL, height = "710px", headerBorder = FALSE,
              dataTableOutput(ns("cl_table_mkrs")))
        ),
        column(width = 12, uiOutput(ns("cl_marker_dotplot_box")))
      )
    )
  )
}

clustersServer <- function(id, shared) {
  moduleServer(id, function(input, output, session) {
    ns = session$ns

    # render cluster selector with initial choices
    cluster_lvls = levels(shared$cluster_markers$cluster)
    output$cl_select_cluster_ui <- renderUI({
      selectizeInput(
        inputId  = ns("cl_select_cluster"),
        label    = div(style = "font-size: 14px", "Select cluster:"),
        choices  = cluster_lvls,
        selected = cluster_lvls[1],
        multiple = FALSE,
        size     = 20
      )
    })

    # PAGA toggle (only shown when paga data was loaded)
    output$cl_display_paga_sw <- renderUI({
      if (shared$paga == 1) {
        prettySwitch(inputId = ns('paga_switch'), label = 'PAGA view', value = FALSE)
      }
    })

    # UMAP or PAGA graph
    output$cl_umap_overview <- renderPlot({
      if (!is.null(input$paga_switch) && shared$paga == 1 && input$paga_switch) {
        plot_cluster_paga(
          pos_dt     = shared$paga_pos,
          mat_dt     = shared$paga_mat,
          weight_cut = shared$paga_weight_cut,
          col_pal    = shared$vars_pals[['cluster']]
        )
      } else {
        shared$cl_cluster_umap
      }
    })

    # filtered markers table for selected cluster
    selected_cluster_mkrs_dt <- reactive({
      sel_cluster = if (!is.null(input$cl_select_cluster) && nchar(input$cl_select_cluster) > 0) input$cl_select_cluster else cluster_lvls[1]
      mkrs = copy(shared$cluster_markers) %>%
        .[cluster == sel_cluster &
            FDR    <= as.numeric(input$cl_cut_fdr)  &
            CPM    >= as.numeric(input$cl_cut_cpm)  &
            log2fc >= input$cl_cut_log2fc] %>%
        .[, -'is.tf'] %>%
        setorder(-log2fc)

      if (input$cl_gene_type_filt == 'protein_coding') {
        mkrs = mkrs[gene_type == 'protein_coding']
      } else if (input$cl_gene_type_filt == 'no_lncRNA') {
        mkrs = mkrs[gene_type != 'lncRNA']
      }
      mkrs
    })

    output$cl_table_mkrs <- renderDataTable({
      DT::datatable(
        selected_cluster_mkrs_dt(),
        options   = list(pageLength = 15, scrollX = TRUE, dom = 'ftp'),
        selection = "single"
      ) %>% formatSignif(columns = "FDR", digits = 2)
    })

    # gene selected by clicking a row in the markers table
    sel_marker <- reactive({
      sel_idx = input$cl_table_mkrs_rows_selected
      if (is.null(sel_idx) || length(sel_idx) == 0) return(character(0))
      selected_cluster_mkrs_dt()[sel_idx, symbol]
    })

    # dotplot box appears only when a marker is selected
    output$cl_marker_dotplot_box <- renderUI({
      if (length(sel_marker()) != 0) {
        box(width = 0, title = NULL, height = "520px", headerBorder = FALSE,
            plotOutput(ns("cl_marker_dotplot"), height = "500px"))
      }
    })

    output$cl_marker_dotplot <- renderPlot({
      if (length(sel_marker()) == 0) return(blank_plot())

      sel_sym     = sel_marker()
      gene_idx    = shared$row_indx[symbol == sel_sym, index]

      marker_dt   = copy(shared$cluster_meta)
      marker_dt$gene_exp = get_row_from_h5(shared$cluster_pb_h5_f, gene_idx)
      marker_dt   = marker_dt[n_cells >= 10, ]

      anno_dt     = make_fdr_anno_dt(sel_sym, shared$cluster_markers)

      plot_gene_dotplot(
        marker_dt,
        subset    = NULL,
        x_axis    = 'cluster',
        fill      = 'cluster',
        fill_name = 'cluster',
        fill_pal  = shared$vars_pals[['cluster']],
        facet_by  = "none",
        gene      = sel_sym,
        anno_dt   = anno_dt,
        note      = sprintf('Dots represent %s clusters in each sample', shared$keyword),
        order     = input$cl_order_plot
      )
    })

    # download handler for filtered markers table
    output$cl_download_table <- downloadHandler(
      filename = function() {
        sel_cluster = if (!is.null(input$cl_select_cluster) && nchar(input$cl_select_cluster) > 0) input$cl_select_cluster else cluster_lvls[1]
        paste0(sprintf("%s_markers_", shared$keyword), sel_cluster, '.csv')
      },
      content = function(file) {
        write.csv(selected_cluster_mkrs_dt(), file, row.names = FALSE, quote = FALSE)
      }
    )
  })
}
