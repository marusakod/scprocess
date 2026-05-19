# MODULE: Explore Genesets -----------------------------------------------------

genesetsUI <- function(id) {
  ns = NS(id)
  tagList(
    fluidRow(
      column(width = 2,
        fluidRow(column(width = 12, uiOutput(ns("get_gset_sel")))),
        fluidRow(column(width = 12, uiOutput(ns("get_heat_params_box"))))
      ),
      uiOutput(ns("get_gsea_stats_box")),
      # heatmap_width is computed in app.R and passed via shared
      uiOutput(ns("get_gset_heat_col"))
    )
  )
}

genesetsServer <- function(id, shared) {
  moduleServer(id, function(input, output, session) {
    ns = session$ns

    # ---------------------------------------------------------------------------
    # Helpers: one function per geneset-selection mode
    # ---------------------------------------------------------------------------

    genes_from_predefined <- function() {
      req(input$sel_predef_gset)
      shared$gsets_dt[name == input$sel_predef_gset, symbol]
    }

    genes_from_tfs <- function() {
      req(input$sel_tfs_opt)
      if (input$sel_tfs_opt == "hv_tfs") {
        shared$pb_hvgs %>%
          .[is.tf == TRUE, ] %>%
          setorder(-vst_var) %>%
          .[, idx := 1:.N] %>%
          .[idx <= input$top_n_hv_tfs, symbol]
      } else {
        shared$cluster_markers %>%
          .[is.tf == TRUE & CPM >= input$CPM_thr_tfs & FDR <= 0.05 & log2fc > 0, ] %>%
          setorder(cluster, FDR) %>%
          .[, idx := 1:.N, by = cluster] %>%
          .[idx <= input$top_n_mkrs_tfs, symbol] %>%
          unique()
      }
    }

    genes_from_markers <- function() {
      req(input$gene_type_filt_mkrs, input$CPM_thr_mkrs, input$sel_mkrs_opt)
      sel_dt = switch(input$gene_type_filt_mkrs,
        protein_coding = copy(shared$cluster_markers)[gene_type == 'protein_coding'],
        no_lncRNA      = copy(shared$cluster_markers)[gene_type != 'lncRNA'],
        copy(shared$cluster_markers)
      )
      if (input$sel_mkrs_opt == 'all') {
        sel_dt %>%
          .[(FDR <= 0.05) & (log2fc > 0) & (CPM >= as.numeric(input$CPM_thr_mkrs))] %>%
          setorder(cluster, FDR) %>%
          .[, idx := 1:.N, by = cluster] %>%
          .[idx <= input$top_n_mkrs_all, symbol] %>%
          unique()
      } else {
        sel_dt %>%
          .[(FDR <= 0.05) & (log2fc > 0) &
              (CPM >= as.numeric(input$CPM_thr_mkrs)) &
              (cluster == input$sel_mkrs_opt)] %>%
          setorder(FDR) %>%
          .[, idx := 1:.N] %>%
          .[idx <= input$top_n_mkrs, symbol] %>%
          unique()
      }
    }

    genes_from_hvgs <- function() {
      req(input$top_n_hvgs)
      shared$pb_hvgs %>% setorder(-vst_var) %>% .[1:input$top_n_hvgs, symbol]
    }

    genes_from_go <- function() {
      req(input$sel_opt_for_go, input$go_cat)
      if (input$sel_opt_for_go == 'all') {
        req(input$go_paths)
        shared$go_terms[pathway == input$go_paths, genes] %>% strsplit(' ') %>% .[[1]]
      } else {
        sel_idx = input$gsea_stats_table_rows_selected
        if (length(sel_idx) == 0) return(NULL)
        nice_path = gsea_stats()[sel_idx, pathway]
        pat_genes = shared$go_terms[pathway_nice == nice_path, genes] %>% strsplit(' ') %>% .[[1]]
        stats = shared$cluster_markers %>%
          .[cluster == input$sel_opt_for_go & symbol %in% pat_genes] %>%
          setorder(FDR)
        if (nrow(stats) > 100) stats$symbol[1:100] else stats$symbol
      }
    }

    parse_custom_geneset <- function() {
      req(input$customGset)
      input$customGset %>%
        stringr::str_split(',') %>% unlist() %>% stringr::str_replace_all(' ', '')
    }

    # ---------------------------------------------------------------------------
    # Dynamic UI
    # ---------------------------------------------------------------------------

    gset_sel_opts <- reactive({
      opts = c(
        "Transcription factors"       = "tfs",
        "Top cluster markers"         = "mkrs",
        "Top highly variable genes"   = "hvgs",
        "Top GO terms for clusters"   = "go",
        "Menu selection"              = "menu",
        "Custom gene list"            = "custom"
      )
      if (shared$predef_gsets == 1)
        opts = c("Predefined genesets" = "pre", opts)
      opts
    })

    output$get_gset_sel <- renderUI({
      box(width = 0, title = NULL, headerBorder = FALSE,
        pickerInput(
          inputId = ns("gsets_in_opt"),
          label   = div(style = "font-size: 14px", "Geneset selection"),
          choices = gset_sel_opts(), multiple = FALSE, selected = NULL,
          options = list(title = "select one option")
        )
      )
    })

    output$get_heat_params_box <- renderUI({
      pre_sel_in = NULL
      if (shared$predef_gsets == 1) {
        gset_opts = shared$gsets_dt$name %>% unique()
        pre_sel_in = if (length(gset_opts) <= 5) {
          awesomeRadio(inputId = ns("sel_predef_gset"),
                       label   = div(style = "font-size: 14px", "Select geneset:"),
                       choices = gset_opts, selected = gset_opts[1])
        } else {
          selectizeInput(inputId = ns("sel_predef_gset"),
                         label   = div(style = "font-size: 14px", "Select geneset:"),
                         choices = gset_opts, selected = gset_opts[1],
                         multiple = FALSE, size = 20)
        }
      }

      mkrs_opts  = c('all', levels(shared$cluster_markers$cluster)) %>%
        setNames(c('all markers', paste0('significant in ', levels(shared$cluster_markers$cluster))))
      fgsea_opts = c('all', levels(shared$fgsea$cluster)) %>%
        setNames(c('all terms', paste0('significant in ', levels(shared$fgsea$cluster))))

      box(width = 0, title = NULL, headerBorder = FALSE,

        conditionalPanel(sprintf("input['%s'] == 'pre'",    ns("gsets_in_opt")), pre_sel_in),

        conditionalPanel(sprintf("input['%s'] == 'tfs'",    ns("gsets_in_opt")),
          awesomeRadio(inputId = ns("sel_tfs_opt"),
                       label   = div(style = "font-size: 14px", "I would like to see:"),
                       choices = c('highly variable TFs' = 'hv_tfs', 'top TF markers' = 'mkrs_tfs'),
                       selected = 'hv_tfs'),
          conditionalPanel(sprintf("input['%s'] == 'hv_tfs'", ns("sel_tfs_opt")),
            numericInput(inputId = ns('top_n_hv_tfs'),
                         label = div(style = "font-size: 14px", "Top n TFs"),
                         value = 30, min = 1, max = 500, step = 1, width = '70px')),
          conditionalPanel(sprintf("input['%s'] == 'mkrs_tfs'", ns("sel_tfs_opt")),
            awesomeRadio(inputId = ns("CPM_thr_tfs"),
                         label = div(style = "font-size: 14px", "CPM threshold"),
                         choices = c("0" = 0, "10" = 10, "20" = 20, "50" = 50, "100" = 100),
                         selected = 50),
            numericInput(inputId = ns('top_n_mkrs_tfs'),
                         label = div(style = "font-size: 14px", "Top n TFs"),
                         value = 5, min = 1, max = 500, step = 1, width = '70px'))
        ),

        conditionalPanel(sprintf("input['%s'] == 'mkrs'",   ns("gsets_in_opt")),
          selectizeInput(inputId = ns("sel_mkrs_opt"),
                         label   = div(style = "font-size: 14px", "Selection criteria"),
                         choices = mkrs_opts, selected = mkrs_opts[1], multiple = FALSE, size = 20),
          conditionalPanel(sprintf("input['%s'] == 'all'",  ns("sel_mkrs_opt")),
            numericInput(inputId = ns('top_n_mkrs_all'),
                         label = div(style = "font-size: 14px", "Top n markers per cluster"),
                         value = 5, min = 1, max = 100, step = 1)),
          conditionalPanel(sprintf("input['%s'] != 'all'",  ns("sel_mkrs_opt")),
            numericInput(inputId = ns('top_n_mkrs'),
                         label = div(style = "font-size: 14px", "Top n markers"),
                         value = 20, min = 1, max = 500, step = 1)),
          awesomeRadio(inputId = ns("gene_type_filt_mkrs"),
                       label = div(style = "font-size: 14px", "Gene type"),
                       choices = c("all" = "all", "protein coding" = "protein_coding",
                                   "excluding lncRNAs" = "no_lncRNA"),
                       selected = "no_lncRNA"),
          awesomeRadio(inputId = ns("CPM_thr_mkrs"),
                       label = div(style = "font-size: 14px", "CPM threshold"),
                       choices = c("0" = 0, "10" = 10, "20" = 20, "50" = 50, "100" = 100),
                       selected = 50)
        ),

        conditionalPanel(sprintf("input['%s'] == 'hvgs'",   ns("gsets_in_opt")),
          numericInput(inputId = ns('top_n_hvgs'),
                       label = div(style = "font-size: 14px", "Top n HVGs"),
                       value = 30, min = 1, max = 500, step = 1)),

        conditionalPanel(sprintf("input['%s'] == 'go'",     ns("gsets_in_opt")),
          selectizeInput(inputId = ns("sel_opt_for_go"),
                         label   = div(style = "font-size: 14px", "Selection criteria"),
                         choices = fgsea_opts, selected = fgsea_opts[1], multiple = FALSE, size = 20),
          radioGroupButtons(inputId = ns("go_cat"),
                            label   = div(style = "font-size: 14px", "GO category"),
                            choices = c("BP" = 'go_bp', "CC" = 'go_cc', "MF" = 'go_mf'),
                            selected = "go_bp",
                            checkIcon = list(yes = icon("ok", lib = "glyphicon"))),
          conditionalPanel(sprintf("input['%s'] == 'all'",  ns("sel_opt_for_go")),
            selectizeInput(inputId = ns('go_paths'),
                           label = div(style = "font-size: 14px", "Pathway"),
                           choices = NULL, multiple = FALSE, size = 20))
        ),

        conditionalPanel(sprintf("input['%s'] == 'menu'",   ns("gsets_in_opt")),
          selectizeInput(inputId = ns("Entry_gene_for_gsets"),
                         label = div(style = "font-size: 14px", "Gene symbol"),
                         multiple = TRUE, size = 20, choices = NULL,
                         options = list(placeholder = sprintf('e.g. %s', shared$default_gene),
                                        onInitialize = I('function() { this.setValue("");}')))
        ),

        conditionalPanel(sprintf("input['%s'] == 'custom'", ns("gsets_in_opt")),
          textAreaInput(inputId = ns("customGset"),
                        label   = div(style = "font-size: 14px", "Enter gene symbols:",
                                      tags$p("Please enter gene symbols separated by commas",
                                             style = "color: #AAB7B8; font-size: 12px; margin-top: 5px;")),
                        rows = 3)
        ),

        prettySwitch(inputId = ns('heat_clust_rows'), label = 'Cluster rows',    value = TRUE),
        prettySwitch(inputId = ns("heat_clust_cols"), label = 'Cluster columns', value = TRUE),
        tags$br(),
        div(style = "text-align:center",
            actionBttn(ns("new_heat"), "Update plot", icon("redo"),
                       size = "sm", style = "material-flat", block = FALSE, color = "primary")),
        tags$br(),
        radioGroupButtons(inputId = ns("heat_show_exp"), label = "Show:",
                          choices = c("logFCs", "CPMs"))
      )
    })

    # ---------------------------------------------------------------------------
    # Observers
    # ---------------------------------------------------------------------------

    # populate GO pathway selectize when cluster/category changes
    observeEvent(c(input$sel_opt_for_go, input$go_cat), {
      req(input$sel_opt_for_go, input$go_cat)
      filtered = shared$go_terms[go_category == input$go_cat, .(path_short, pathway)] %>%
        tibble::deframe()
      updateSelectizeInput(session, "go_paths", choices = filtered, selected = filtered[1])
    }, ignoreNULL = FALSE)

    # populate menu gene selectize when that option is chosen
    observe({
      req(input$gsets_in_opt == "menu")
      updateSelectizeInput(session, "Entry_gene_for_gsets",
        server   = TRUE, selected = NULL,
        options  = list(placeholder = sprintf('e.g. %s', shared$default_gene),
                        onInitialize = I('function() { this.setValue("");}')),
        choices  = shared$row_indx$symbol)
    })

    # ---------------------------------------------------------------------------
    # GSEA stats table (shown when a specific cluster's GO terms are selected)
    # ---------------------------------------------------------------------------

    gsea_stats <- reactive({
      cl = input$sel_opt_for_go
      req(cl, cl != 'all')
      copy(shared$fgsea) %>%
        .[cluster == cl & go_category == input$go_cat, .(pathway_nice, padj, NES, size)] %>%
        setnames('pathway_nice', 'pathway') %>%
        setorder(padj)
    })

    output$gsea_stats_table <- renderDataTable({
      DT::datatable(gsea_stats(),
                    options   = list(pageLength = 15, scrollX = TRUE, dom = 'ftp'),
                    selection = "single") %>%
        DT::formatRound(columns = c('padj', 'NES'), digits = 3)
    })

    output$get_gsea_stats_box <- renderUI({
      req(input$sel_opt_for_go, input$gsets_in_opt)
      if (input$sel_opt_for_go != 'all' && input$gsets_in_opt == 'go') {
        column(width = 3,
          box(width = 0, title = NULL, headerBorder = FALSE,
              dataTableOutput(ns("gsea_stats_table"))))
      }
    })

    # ---------------------------------------------------------------------------
    # Geneset reactive: resolve selected genes
    # ---------------------------------------------------------------------------

    gset_ls <- reactive({
      req(input$gsets_in_opt, nchar(input$gsets_in_opt) > 0)

      genes = switch(input$gsets_in_opt,
        pre    = genes_from_predefined(),
        tfs    = genes_from_tfs(),
        mkrs   = genes_from_markers(),
        hvgs   = genes_from_hvgs(),
        go     = genes_from_go(),
        menu   = input$Entry_gene_for_gsets,
        custom = parse_custom_geneset(),
        NULL
      )

      if (is.null(genes) || all(genes == '')) return(NULL)

      genes = intersect(genes, unique(shared$cluster_markers$symbol))
      if (length(genes) == 0) NULL else genes
    })

    # ---------------------------------------------------------------------------
    # Heatmap rendering
    # ---------------------------------------------------------------------------

    heat_height <- reactive({
      input$new_heat
      isolate({
        gs_size = length(gset_ls())
        if (gs_size < 2) return(NULL)
        sprintf('%spx', gs_size * 20 + 160)
      })
    })

    heat_box_height <- reactive({
      input$new_heat
      isolate({
        if (is.null(heat_height())) return('50px')
        sprintf('%spx', as.numeric(gsub("px", "", heat_height())) + 50)
      })
    })

    gset_heatmap_obj <- reactive({
      input$new_heat
      isolate({
        req(input$gsets_in_opt)
        hmap_gs = gset_ls()
        if (length(hmap_gs) > 1) {
          plot_heatmap(shared$cluster_markers, hmap_gs,
                       cluster_rows    = input$heat_clust_rows,
                       cluster_columns = input$heat_clust_cols,
                       show = 'log2fc')
        }
      })
    })

    output$gset_heatmap <- renderPlot({
      hmap_gs = gset_ls()
      hmap    = gset_heatmap_obj()
      req(!is.null(hmap), length(hmap_gs) > 1)
      if (input$heat_show_exp == 'CPMs') {
        plot_heatmap(shared$cluster_markers, hmap_gs,
                     cluster_rows    = input$heat_clust_rows,
                     cluster_columns = input$heat_clust_cols,
                     show = 'logcpm')
      } else {
        hmap
      }
    })

    output$get_gset_heat_box <- renderUI({
      hmap    = gset_heatmap_obj()
      hmap_gs = gset_ls()
      input$new_heat
      isolate({
        req(!is.null(hmap), length(hmap_gs) > 1)
        box(width = 0, title = NULL, headerBorder = FALSE,
            height = heat_box_height(),
            plotOutput(ns('gset_heatmap'), height = heat_height()))
      })
    })

    output$get_gset_heat_col <- renderUI({
      column(width = shared$heatmap_width, uiOutput(ns("get_gset_heat_box")))
    })

    observeEvent(input$new_heat, {
      hmap    = gset_heatmap_obj()
      hmap_gs = gset_ls()
      if (is.null(hmap) || length(hmap_gs) <= 1) {
        sendSweetAlert(session, title = "Error...",
                       text = "Please provide at least two valid genes", type = "error")
      }
    })
  })
}
