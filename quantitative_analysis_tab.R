# quantitative_analysis_tab.R --------------------------------------------------
# Runs PiProteline::quantitative_analysis() on the already-preprocessed data
# (never on the raw upload). The volcano and pairwise-MDS pickers are built
# from the actual result names, so this works for any number of groups (2..n),
# not just a fixed set. The gene column is always column 1 of data_norm, since
# preprocessing_data() reorders it there.

library(shiny)
library(DT)
library(shinyWidgets)
library(PiProteline)

quantitative_analysis_ui <- function() {
  tabPanel(
    "Quantitative Analysis",
    sidebarLayout(
      sidebarPanel(
        numericInput("significance_manova", "Significance threshold (MANOVA)",
                     value = 0.05, min = 0.001, max = 1, step = 0.001),
        actionBttn("run_analysis", "Run Quantitative Analysis", style = "fill", color = "primary"),
        helpText("Run 'Preprocess data' first."),
        tags$hr(),

        actionButton("show_manova", "Show MANOVA Results"),
        tags$hr(),

        h4("Volcano (pairwise)"),
        selectInput("qa_volcano_pick", "Contrast", choices = NULL),
        actionButton("show_volcano", "Show volcano"),
        tags$hr(),

        h4("MDS"),
        actionButton("show_mds", "Show global MDS"),
        tags$br(), tags$br(),
        selectInput("qa_mds_pick", "Pairwise MDS contrast", choices = NULL),
        actionButton("show_mds_pairwise", "Show pairwise MDS")
      ),
      mainPanel(
        conditionalPanel("input.show_manova > 0",
                         DTOutput("manova_table"),
                         downloadButton("download_manova", "Download MANOVA Results")),
        tags$hr(),

        h4(textOutput("qa_volcano_title")),
        plotOutput("qa_volcano"),
        downloadButton("download_volcano", "Download volcano"),
        tags$hr(),

        h4(textOutput("qa_mds_title")),
        plotOutput("qa_mds_global"),
        downloadButton("download_mds_global", "Download global MDS"),
        tags$hr(),

        h4(textOutput("qa_mds_pw_title")),
        plotOutput("qa_mds_pairwise"),
        downloadButton("download_mds_pairwise", "Download pairwise MDS")
      )
    )
  )
}

quantitative_analysis_server <- function(input, output, session, shared_data) {

  qa_result <- reactiveVal(NULL)

  observeEvent(input$run_analysis, {
    req(shared_data$preproc_data, shared_data$group_names)
    pp <- shared_data$preproc_data

    withProgress(message = "Running Quantitative Analysis", value = 0, {
      setProgress(0.4, detail = "Computing MANOVA, volcano, MDS ...")
      tryCatch({
        result <- PiProteline::quantitative_analysis(
          dataset             = pp$data_norm,
          names_of_groups     = shared_data$group_names,
          gene_column         = 1,                       # gene col = column 1 of data_norm
          data_grouped_full   = pp$data_grouped_full,
          significance_manova = input$significance_manova
        )

        # Compatibility shim for the pairwise field name across package versions.
        if (is.null(result$manova_pairw_results)) {
          if (!is.null(result$manova_pairw))            result$manova_pairw_results <- result$manova_pairw
          if (!is.null(result$manova_pairwise_results)) result$manova_pairw_results <- result$manova_pairwise_results
        }

        qa_result(result)
        shared_data$quant_results <- result   # share with the network tab

        # Populate the contrast pickers from the real results.
        updateSelectInput(session, "qa_volcano_pick",
                          choices = contrast_labels(result$volcano_plots,  shared_data$group_names))
        updateSelectInput(session, "qa_mds_pick",
                          choices = contrast_labels(result$mds_plot_pairw, shared_data$group_names))

        setProgress(1)
        showNotification("Quantitative analysis completed.", type = "message")
      }, error = function(e) {
        showNotification(paste("Error in quantitative analysis:", e$message), type = "error")
        qa_result(NULL)
      })
    })
  })

  output$manova_table <- renderDT({ req(qa_result()); qa_result()$manova_results })
  output$download_manova <- downloadHandler("manova_results.csv",
    function(file) write.csv(qa_result()$manova_results, file, row.names = FALSE))

  # ---- Volcano (dynamic) ----
  sel_volcano <- reactive({
    req(qa_result(), input$qa_volcano_pick)
    lst  <- qa_result()$volcano_plots
    labs <- contrast_labels(lst, shared_data$group_names)
    idx  <- match(input$qa_volcano_pick, labs)
    validate(need(!is.na(idx), "Selected contrast not found."))
    lst[[idx]]
  })
  output$qa_volcano_title <- renderText({
    if (isTruthy(input$qa_volcano_pick)) paste("Volcano:", input$qa_volcano_pick) else ""
  })
  output$qa_volcano <- renderPlot({ req(input$show_volcano); print(sel_volcano()) })
  output$download_volcano <- downloadHandler(
    filename = function() paste0("volcano_", gsub("[^A-Za-z0-9]+","_", input$qa_volcano_pick), ".png"),
    content  = function(file) { png(file, 1400, 900, res = 150); print(sel_volcano()); dev.off() }
  )

  # ---- Global MDS ----
  output$qa_mds_title <- renderText({ if (isTruthy(input$show_mds)) "Global MDS" else "" })
  output$qa_mds_global <- renderPlot({ req(input$show_mds, qa_result()); print(qa_result()$mds_plot) })
  output$download_mds_global <- downloadHandler("mds_global.png",
    function(file) { png(file, 1400, 900, res = 150); print(qa_result()$mds_plot); dev.off() })

  # ---- Pairwise MDS (dynamic) ----
  sel_mds_pw <- reactive({
    req(qa_result(), input$qa_mds_pick)
    lst  <- qa_result()$mds_plot_pairw
    labs <- contrast_labels(lst, shared_data$group_names)
    idx  <- match(input$qa_mds_pick, labs)
    validate(need(!is.na(idx), "Selected pairwise MDS not found."))
    lst[[idx]]
  })
  output$qa_mds_pw_title <- renderText({
    if (isTruthy(input$qa_mds_pick)) paste("Pairwise MDS:", input$qa_mds_pick) else ""
  })
  output$qa_mds_pairwise <- renderPlot({ req(input$show_mds_pairwise); print(sel_mds_pw()) })
  output$download_mds_pairwise <- downloadHandler(
    filename = function() paste0("mds_", gsub("[^A-Za-z0-9]+","_", input$qa_mds_pick), ".png"),
    content  = function(file) { png(file, 1400, 900, res = 150); print(sel_mds_pw()); dev.off() }
  )
}
