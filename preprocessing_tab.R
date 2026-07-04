# preprocessing_tab.R (v2) -----------------------------------------------------
# Changes vs v1:
#  * FIX #4: the user's chosen gene column is stored ONCE in `gene_column_raw`
#    and is never overwritten. v1 clobbered it with colnames(data_norm)[1] right
#    after preprocessing, so by the time other tabs read it, it was wrong.
#  * The chosen normalization (`norm_type`) is now stored in shared_data so the
#    rest of the app uses the SAME normalization the user selected here
#    (part of FIX #3 — no tab silently falls back to a different default).

library(shiny)
library(DT)
library(openxlsx)
library(shinyWidgets)
library(PiProteline)

preprocessing_ui <- function() {
  tabPanel(
    "Preprocessing",
    sidebarLayout(
      sidebarPanel(
        fileInput("dataset", "Upload dataset (.xlsx, .xls, .csv)",
                  accept = c(".xlsx", ".xls", ".csv")),
        uiOutput("sheet_ui"),          # shown only for Excel files with >1 sheet
        uiOutput("gene_column_ui"),

        textInput("group_names", "Group names (comma-separated)",
                  value = "", placeholder = "e.g. Control,Treated"),
        helpText("Type the group tags exactly as they appear inside your sample ",
                 "column names. Avoid the reserved words listed in the README ",
                 "(fc, specific, centrality, weighted, fdr, ...)."),

        selectInput("norm_type", "Normalization type",
                    choices = c("ln","Znorm","MinMax","Robust","UnitVector",
                                "TotSigNorm","MaxSigNorm","RowSigmaNorm"),
                    selected = "TotSigNorm"),
        helpText("Note: Znorm / Robust / MinMax can produce zero or negative ",
                 "values, which break the log2 fold change used downstream. ",
                 "TotSigNorm keeps values positive."),

        actionBttn("preprocess", "Preprocess data", style = "fill", color = "primary"),

        tags$hr(),
        actionButton("show_data_unique", "Show data_unique"),
        actionButton("show_data_norm",   "Show data_norm"),
        tags$hr(),
        actionBttn("desc_stats", "Descriptive statistics", style = "fill", color = "success"),
        actionButton("show_desc_stats_col", "Show DS_col"),
        actionButton("show_desc_stats_row", "Show DS_row")
      ),
      mainPanel(
        verbatimTextOutput("debug_output"),
        conditionalPanel("input.show_data_unique > 0",
                         h4("data_unique"), DTOutput("data_unique_table"),
                         downloadButton("download_data_unique","Download")),
        conditionalPanel("input.show_data_norm > 0",
                         h4("data_norm"), DTOutput("data_norm_table"),
                         downloadButton("download_data_norm","Download")),
        conditionalPanel("input.show_desc_stats_col > 0",
                         h4("Descriptive Stats - Columns"), DTOutput("desc_stats_col_table"),
                         downloadButton("download_desc_stats_col","Download")),
        conditionalPanel("input.show_desc_stats_row > 0",
                         h4("Descriptive Stats - Rows"), DTOutput("desc_stats_row_table"),
                         downloadButton("download_desc_stats_row","Download"))
      )
    )
  )
}

preprocessing_server <- function(input, output, session) {

  shared_data <- reactiveValues(
    dataset         = NULL,
    group_names     = NULL,
    gene_column_raw = NULL,   # the user's ORIGINAL choice; never overwritten
    norm_type       = NULL,   # the normalization the user picked
    preproc_data    = NULL,
    desc_stats      = NULL,
    quant_results   = NULL,   # filled by the quantitative tab (shared)
    pipelineResults = NULL    # filled by the network tab (shared)
  )

  output$debug_output <- renderText(paste("Temp dir:", tempdir()))

  # Hold the uploaded file so we can re-read it when the chosen sheet changes.
  upload <- reactiveValues(path = NULL, ext = NULL, sheets = NULL)

  # Step 1: a file is uploaded. For Excel, list the sheets and pick the most
  # likely DATA sheet by default (the one with the most columns) - this stops a
  # "Legend"/"Legenda" sheet from silently hijacking the upload. For CSV, read
  # it straight away (no sheets).
  observeEvent(input$dataset, {
    req(input$dataset)
    upload$ext  <- tolower(tools::file_ext(input$dataset$name))
    upload$path <- input$dataset$datapath

    if (upload$ext %in% c("xlsx", "xls")) {
      sheets <- tryCatch(openxlsx::getSheetNames(upload$path), error = function(e) character(0))
      upload$sheets <- sheets

      # Guess the data sheet: the one with the most columns in its header row.
      default_sheet <- sheets[1]
      if (length(sheets) > 1) {
        ncols <- vapply(sheets, function(s) {
          h <- tryCatch(openxlsx::read.xlsx(upload$path, sheet = s, rows = 1),
                        error = function(e) NULL)
          if (is.null(h)) 0L else ncol(h)
        }, integer(1))
        default_sheet <- sheets[which.max(ncols)]
      }

      output$sheet_ui <- renderUI({
        if (length(sheets) <= 1) return(NULL)   # nothing to choose
        selectInput("sheet", "Worksheet (the data sheet)",
                    choices = sheets, selected = default_sheet)
      })
      # If there is only one sheet, load it now; otherwise the sheet observer loads it.
      if (length(sheets) <= 1) {
        shared_data$dataset <- as.data.frame(openxlsx::read.xlsx(upload$path))
      }
    } else if (upload$ext == "csv") {
      output$sheet_ui <- renderUI(NULL)
      upload$sheets <- NULL
      shared_data$dataset <- as.data.frame(read.csv(upload$path, check.names = FALSE))
    } else {
      output$sheet_ui <- renderUI(NULL)
      showNotification("Unsupported file type", type = "error")
      shared_data$dataset <- NULL
    }
  })

  # Step 2: (re)load the chosen Excel sheet whenever the selection changes.
  observeEvent(input$sheet, {
    req(upload$path, input$sheet)
    shared_data$dataset <- tryCatch(
      as.data.frame(openxlsx::read.xlsx(upload$path, sheet = input$sheet)),
      error = function(e) { showNotification(paste("Could not read sheet:", e$message), type = "error"); NULL }
    )
  })

  output$gene_column_ui <- renderUI({
    req(shared_data$dataset)
    cols <- colnames(shared_data$dataset)
    default <- if ("GeneName" %in% cols) "GeneName" else cols[1]
    selectInput("gene_column", "Select gene column", choices = cols, selected = default)
  })

  observe({
    req(input$group_names)
    shared_data$group_names <- trimws(strsplit(input$group_names, ",")[[1]])
  })
  # Store the user's choice once. This is the value used for preprocessing.
  observe({ req(input$gene_column); shared_data$gene_column_raw <- input$gene_column })
  observe({ req(input$norm_type);   shared_data$norm_type       <- input$norm_type })

  observeEvent(input$preprocess, {
    req(shared_data$dataset, shared_data$group_names,
        shared_data$gene_column_raw, shared_data$norm_type)

    withProgress(message = "Preprocessing", value = 0, {
      tryCatch({
        res <- PiProteline::preprocessing_data(
          dataset         = shared_data$dataset,
          names_of_groups = shared_data$group_names,
          gene_column     = shared_data$gene_column_raw,  # user's ORIGINAL column
          norm_type       = shared_data$norm_type
        )
        shared_data$preproc_data <- res
        # NOTE: we deliberately do NOT overwrite gene_column_raw here.
        # After preprocessing, the gene column is column 1 of data_norm /
        # data_grouped (the backend reorders it), so downstream steps use 1 -
        # exactly as PiProteline::pipeline() does internally.
        showNotification("Preprocessing done. You can now run the other tabs.",
                         type = "message")
      }, error = function(e) {
        showNotification(paste("Preprocessing error:", e$message), type = "error")
        shared_data$preproc_data <- NULL
      })
    })
  })

  observeEvent(input$desc_stats, {
    req(shared_data$preproc_data)
    tryCatch({
      shared_data$desc_stats <- PiProteline::descriptive_statistics(
        shared_data$preproc_data$data_unique,
        shared_data$preproc_data$data_grouped_full
      )
      showNotification("Descriptive statistics computed", type = "message")
    }, error = function(e) {
      showNotification(paste("Desc. stats error:", e$message), type = "error")
      shared_data$desc_stats <- NULL
    })
  })

  output$data_unique_table   <- renderDT({ req(shared_data$preproc_data, input$show_data_unique); shared_data$preproc_data$data_unique })
  output$data_norm_table     <- renderDT({ req(shared_data$preproc_data, input$show_data_norm);   shared_data$preproc_data$data_norm })
  output$desc_stats_col_table <- renderDT({ req(shared_data$desc_stats, input$show_desc_stats_col); shared_data$desc_stats$DS_col })
  output$desc_stats_row_table <- renderDT({ req(shared_data$desc_stats, input$show_desc_stats_row); shared_data$desc_stats$DS_row })

  output$download_data_unique   <- downloadHandler("data_unique.csv", function(f) write.csv(shared_data$preproc_data$data_unique, f, row.names = FALSE))
  output$download_data_norm     <- downloadHandler("data_norm.csv",   function(f) write.csv(shared_data$preproc_data$data_norm,   f, row.names = FALSE))
  output$download_desc_stats_col <- downloadHandler("DS_col.csv", function(f) write.csv(shared_data$desc_stats$DS_col, f, row.names = FALSE))
  output$download_desc_stats_row <- downloadHandler("DS_row.csv", function(f) write.csv(shared_data$desc_stats$DS_row, f, row.names = FALSE))

  return(shared_data)
}
