# MODULE: Explore Prevalence ---------------------------------------------------

prevalenceUI <- function(id) {
  ns = NS(id)
  tagList(
    fluidRow(
      column(width = 2,
        box(width = 0, title = NULL, headerBorder = FALSE,
          awesomeRadio(
            inputId  = ns("meta_which_plot"),
            label    = div(style = "font-size: 14px", "Composition plot"),
            choices  = c("barplot" = 1, "dotplot" = 2, "density plot" = 3),
            selected = 1,
            inline   = FALSE
          ),

          # --- Barplot controls ---
          conditionalPanel(sprintf("input['%s'] == 1", ns("meta_which_plot")),
            uiOutput(ns("meta_selected_var_ui")),
            uiOutput(ns("meta_order_barplot_ui"))
          ),

          # --- Dotplot controls ---
          conditionalPanel(sprintf("input['%s'] == 2", ns("meta_which_plot")),
            uiOutput(ns("meta_comp_x_var_ui")),
            uiOutput(ns("meta_comp_col_var_ui")),
            uiOutput(ns("meta_comp_sel_cls_ui")),
            p(strong("Show advanced plotting options")),
            prettySwitch(inputId = ns("meta_comp_show_advanced"), label = '', value = FALSE),
            conditionalPanel(sprintf("input['%s'] == true", ns("meta_comp_show_advanced")),
              uiOutput(ns("meta_comp_facet_ui")),
              numericInput(inputId = ns('meta_comp_ncols'),
                           label = div(style = "font-size: 14px", "Number of columns"),
                           value = 2, min = 1, max = 8, step = 1, width = '150px'),
              prettySwitch(inputId = ns("meta_comp_order"),
                           label = "Order by mean proportion", value = FALSE),
              p(strong("Subset data on...")),
              uiOutput(ns("meta_comp_subset_picker"))
            ),
            tags$br(),
            div(style = "text-align:center",
                actionBttn(ns("meta_update_plot"), "Update plot", icon("redo"),
                           size = "sm", style = "material-flat", block = FALSE, color = "primary")),
            bsPopover(ns("meta_update_plot"), title = "Update plot", placement = "bottom",
                      options = list(container = "body"),
                      content = "This will generate a new plot with selected clusters")
          ),

          # --- Density controls ---
          conditionalPanel(sprintf("input['%s'] == 3", ns("meta_which_plot")),
            uiOutput(ns("meta_density_var_ui")),
            p(strong("Show advanced plotting options")),
            prettySwitch(inputId = ns("meta_dens_show_advanced"), label = '', value = FALSE),
            conditionalPanel(sprintf("input['%s'] == true", ns("meta_dens_show_advanced")),
              numericInput(inputId = ns('meta_dens_ncols'),
                           label = div(style = "font-size: 14px", "Number of columns"),
                           value = NULL, min = 1, max = 8, step = 1, width = '150px'),
              p(strong("Subset data on...")),
              uiOutput(ns("meta_dens_subset_picker"))
            )
          )
        )
      ),
      uiOutput(ns("get_comp_box"))
    )
  )
}

prevalenceServer <- function(id, shared) {
  moduleServer(id, function(input, output, session) {
    ns = session$ns

    # render selectize/picker inputs with correct initial choices
    meta_and_cluster = c(shared$metadata_vars, "cluster" = "cluster")

    output$meta_selected_var_ui <- renderUI({
      selectizeInput(ns("meta_selected_var"),
        label = div(style = "font-size: 14px", "Select variable:"),
        choices = shared$metadata_vars, selected = shared$metadata_vars[[1]],
        multiple = FALSE, size = 20)
    })
    output$meta_order_barplot_ui <- renderUI({
      selectizeInput(ns("meta_order_barplot"),
        label = div(style = "font-size: 14px", "Order rows by:"),
        choices = c('selected variable', levels(shared$cell_meta$cluster)),
        selected = 'selected variable', multiple = FALSE, size = 20)
    })
    output$meta_comp_x_var_ui <- renderUI({
      selectizeInput(ns("meta_comp_x_var"),
        label = div(style = "font-size: 14px", "x-axis variable"),
        choices = meta_and_cluster, selected = shared$metadata_vars[[1]],
        multiple = FALSE, size = 20)
    })
    output$meta_comp_col_var_ui <- renderUI({
      selectizeInput(ns("meta_comp_col_var"),
        label = div(style = "font-size: 14px", "Colour"),
        choices = meta_and_cluster, selected = shared$metadata_vars[[1]],
        multiple = FALSE, size = 20)
    })
    output$meta_comp_sel_cls_ui <- renderUI({
      pickerInput(ns("meta_comp_sel_cls"),
        label    = 'Select clusters:',
        multiple = TRUE,
        options  = list(`actions-box` = TRUE),
        choices  = levels(shared$cell_meta$cluster),
        selected = levels(shared$cell_meta$cluster)[[1]])
    })
    output$meta_comp_facet_ui <- renderUI({
      selectizeInput(ns("meta_comp_facet"),
        label = div(style = "font-size: 12px", "Facet by"),
        choices = c("cluster", shared$metadata_vars), selected = "cluster",
        multiple = FALSE)
    })
    output$meta_density_var_ui <- renderUI({
      selectizeInput(ns("meta_density_var"),
        label = div(style = "font-size: 14px", "x-axis variable"),
        choices = shared$metadata_vars, selected = shared$metadata_vars[[1]],
        multiple = FALSE, size = 20)
    })

    # subset picker UIs
    output$meta_comp_subset_picker <- renderUI({
      make_subset_picker_ui("meta_comp_subset_", shared$metadata_vars, shared$cluster_meta, ns)
    })

    output$meta_dens_subset_picker <- renderUI({
      make_subset_picker_ui("meta_dens_subset_", shared$subset_vars, shared$cluster_meta, ns)
    })

    # ---------------------------------------------------------------------------
    # Plot output
    # ---------------------------------------------------------------------------

    output$meta_comp_plot <- renderPlot({
      default_var  = shared$metadata_vars[[1]]
      sel_var      = if (!is.null(input$meta_selected_var)  && nchar(input$meta_selected_var)  > 0) input$meta_selected_var  else default_var
      order_bp     = if (!is.null(input$meta_order_barplot) && nchar(input$meta_order_barplot) > 0) input$meta_order_barplot else 'selected variable'
      comp_x_var   = if (!is.null(input$meta_comp_x_var)   && nchar(input$meta_comp_x_var)   > 0) input$meta_comp_x_var    else default_var
      comp_col_var = if (!is.null(input$meta_comp_col_var)  && nchar(input$meta_comp_col_var)  > 0) input$meta_comp_col_var  else default_var
      dens_var     = if (!is.null(input$meta_density_var)   && nchar(input$meta_density_var)   > 0) input$meta_density_var   else default_var

      if (input$meta_which_plot == 1) {
        # barplot
        cl_ord = ifelse(order_bp != 'selected variable',
          gsub('by: ', '', order_bp), 'no')
        plot_metadata_barplot(
          meta_dt  = shared$cell_meta,
          meta_var = sel_var,
          col_pal  = shared$vars_pals[['cluster']],
          cl_ord   = cl_ord
        )

      } else if (input$meta_which_plot == 2) {
        # dotplot
        input$meta_update_plot
        isolate({
          dot_pre  = shared$cell_meta %>%
            .[, c('cell_id', shared$sample_col, 'cluster', 'n_cells',
                  unname(shared$metadata_vars)), with = FALSE]
          props_dt = make_composition_dt(dot_pre, sample_col = shared$sample_col) %>%
            .[ cluster %in% input$meta_comp_sel_cls ]

          subset_ls = shared$metadata_vars %>%
            lapply(function(var) input[[ paste0("meta_comp_subset_", var) ]]) %>%
            setNames(shared$metadata_vars)
          for (nn in names(subset_ls)) {
            if (!is.null(subset_ls[[ nn ]])) {
              props_dt = props_dt[ props_dt[[ nn ]] %in% subset_ls[[ nn ]] ]
            }
          }

          if (nrow(props_dt) == 0) return(blank_plot())

          g = plot_metadata_dotplot(
            props_dt  = props_dt,
            x_var     = comp_x_var,
            facet_by  = if (!is.null(input$meta_comp_facet) && nchar(input$meta_comp_facet) > 0) input$meta_comp_facet else "cluster",
            leg_title = shared$var_choice_names[ comp_col_var ],
            col_var   = comp_col_var,
            col_pal   = shared$vars_pals[[ comp_col_var ]],
            order     = input$meta_comp_order,
            n_cols    = input$meta_comp_ncols
          )

          n_panels = length(unique(shared$cell_meta[[ comp_col_var ]]))
          if (n_panels > 15) {
            g = g + guides(fill = guide_legend(
              nrow = ceiling(n_panels / 5), byrow = TRUE,
              override.aes = list(size = 5), title.position = 'top'))
          } else {
            g = g + guides(fill = guide_legend(override.aes = list(size = 5)))
          }
          g
        })

      } else {
        # density plot
        subset_ls = shared$subset_vars %>%
          lapply(function(var) input[[ paste0("meta_dens_subset_", var) ]]) %>%
          setNames(shared$subset_vars)

        meta_tmp = copy(shared$cell_meta)
        for (nn in names(subset_ls)) {
          if (!is.null(subset_ls[[ nn ]]))
            meta_tmp = meta_tmp[ meta_tmp[[ nn ]] %in% subset_ls[[ nn ]] ]
        }

        x_range = c(min(shared$cell_meta$UMAP_1), max(shared$cell_meta$UMAP_1)) * 1.05
        y_range = c(min(shared$cell_meta$UMAP_2), max(shared$cell_meta$UMAP_2)) * 1.05

        plot_metadata_density(meta_tmp, meta_var = dens_var,
                              x_range = x_range, y_range = y_range,
                              n_cols = input$meta_dens_ncols)
      }
    })

    # ---------------------------------------------------------------------------
    # Dynamic box sizing
    # ---------------------------------------------------------------------------

    comp_plot_height <- reactive({
      if (input$meta_which_plot == 1) return('500px')
      if (input$meta_which_plot == 3) return('500px')

      # dotplot
      input$meta_update_plot
      isolate({
        facet = if (!is.null(input$meta_comp_facet) && nchar(input$meta_comp_facet) > 0) input$meta_comp_facet else "cluster"
        facet_vals  = shared$cell_meta[[ facet ]] %>% unique %>% as.character
        subset_var  = if (facet == "cluster") "meta_comp_sel_cls"
                      else paste0("meta_comp_subset_", facet)
        subset_vals = input[[ subset_var ]]
        if (is.null(subset_vals) && facet == "cluster")
          subset_vals = levels(shared$cell_meta$cluster)[[1]]
        if (!is.null(subset_vals))
          facet_vals = intersect(facet_vals, subset_vals)
        n_rows = ceiling(length(facet_vals) / input$meta_comp_ncols)
        sprintf('%spx', 150 + 200 * n_rows)
      })
    })

    comp_plot_box_height <- reactive({
      height_pl = as.numeric(gsub("px", "", comp_plot_height()))
      if (input$meta_which_plot == 1) return('550px')
      sprintf('%spx', height_pl + 50)
    })

    comp_plot_column_width <- reactive({
      switch(as.character(input$meta_which_plot), "1" = 8L, "2" = 10L, 8L)
    })

    output$get_comp_box <- renderUI({
      column(
        width = comp_plot_column_width(),
        box(width = 0, title = NULL, headerBorder = FALSE,
            height = comp_plot_box_height(),
            plotOutput(ns("meta_comp_plot"), height = comp_plot_height()))
      )
    })
  })
}
