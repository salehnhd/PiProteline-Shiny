# network_analysis_tab.R -------------------------------------------------------
# Runs network + functional analysis on the ALREADY-preprocessed data
# (shared_data$preproc_data), not on the raw upload. This matters because
# PiProteline::pipeline() re-runs preprocessing internally with its own
# defaults (TotSigNorm, gene_column = 1), which would silently ignore
# whatever normalization/gene-column the user picked in the Preprocessing
# tab. Instead this tab replicates pipeline()'s downstream steps directly, so
# nothing is recomputed and the user's choices carry through end to end. The
# gene column is column 1 of the preprocessed tables (PiProteline reorders it
# there), matching what pipeline() does internally.
# The Group dropdowns are populated from shared_data$group_names, so this
# works for any dataset with any number of groups.

library(shiny)
library(shinyWidgets)
library(CoPPIs)        # builds the interactome (see README for citation)
library(igraph)
library(dplyr)
library(PiProteline)
library(DT)

network_analysis_ui <- function() {
  tabPanel(
    "Network Analysis",
    sidebarLayout(
      sidebarPanel(

        helpText("Run 'Preprocess data' first. This step reuses that output - ",
                 "it does not re-run preprocessing."),
        numericInput("na_quantile", "Critical-node quantile threshold",
                     value = 0.75, min = 0.50, max = 0.99, step = 0.05),
        textInput("na_save_prefix", "Save results as (prefix)", value = "PiProteline_report"),
        actionBttn("na_run_pipeline", "Run network + functional analysis",
                   style = "fill", color = "primary"),
        helpText("The Functional Analysis tab will use these results."),
        tags$hr(),

        h4("Centralities"),
        selectInput("na_cent_group", "Group", choices = NULL),
        actionButton("na_show_unweighted_cent", "Unweighted_centralities"),
        actionButton("na_show_weighted_cent",   "Weighted_centralities"),
        tags$hr(),

        h4("Critical Nodes"),
        selectInput("na_cn_scope", "Scope", choices = c("NotSpecific","Specific","CentralitySpecific")),
        selectInput("na_cn_group", "Group", choices = NULL),
        actionButton("na_show_unweighted_cn", "Unweighted_criticalNodes"),
        actionButton("na_show_weighted_cn",   "Weighted_criticalNodes"),
        tags$hr(),

        h4("PPI graphs"),
        selectInput("na_ppi_group", "Group", choices = NULL),
        actionButton("na_plot_ppi_unweighted",   "Show PPI_unweighted"),
        actionButton("na_plot_ppi_correlations", "Show PPI_correlations")
      ),
      mainPanel(
        h4(textOutput("na_table_title")),
        DTOutput("na_table"),
        downloadButton("na_table_dl", "Download CSV"),
        tags$hr(),
        h4(textOutput("na_plot_title")),
        plotOutput("na_plot", height = 500),
        downloadButton("na_plot_dl", "Download Plot")
      )
    )
  )
}

network_analysis_server <- function(input, output, session, shared_data) {

  # Keep all three Group selectors in sync with the real group names.
  observe({
    req(shared_data$group_names)
    for (id in c("na_cent_group", "na_cn_group", "na_ppi_group")) {
      updateSelectInput(session, id, choices = shared_data$group_names)
    }
  })

  # ---- Run: consume the preprocessed object, reproduce pipeline()'s steps ----
  observeEvent(input$na_run_pipeline, {
    req(shared_data$preproc_data, shared_data$group_names)
    pp <- shared_data$preproc_data

    withProgress(message = "Running network + functional analysis", value = 0, {
      tryCatch({

        setProgress(0.15, detail = "Building interactome (CoPPIs)...")
        g_interactome <- CoPPIs::interactome.hs %>%
          CoPPIs::filter_interactome(scores_threshold = c(experimental = 150, database = 300)) %>%
          dplyr::select(3, 4) %>%
          igraph::graph_from_data_frame(directed = FALSE)

        # Reuse the quantitative results if the user already ran that tab;
        # otherwise compute them here on the SAME normalized data.
        setProgress(0.35, detail = "Quantitative analysis...")
        quantitativeAnalysis <- shared_data$quant_results
        if (is.null(quantitativeAnalysis)) {
          quantitativeAnalysis <- PiProteline::quantitative_analysis(
            dataset             = pp$data_norm,
            names_of_groups     = shared_data$group_names,
            gene_column         = 1,                    # gene col = column 1 of data_norm
            data_grouped_full   = pp$data_grouped_full,
            # `significance_manova` is defined in the Quantitative Analysis
            # tab's UI, not this tab's. It works because Shiny's `input` is
            # shared across the whole app (all tabs get the same `input`
            # object from app.R), but it means this tab depends on that one
            # still existing there.
            significance_manova = input$significance_manova %||% 0.05
          )
        }
        manova_pw <- quantitativeAnalysis$manova_pairw_results %||%
                     quantitativeAnalysis$manova_pairw %||%
                     quantitativeAnalysis$manova_pairwise_results

        setProgress(0.60, detail = "Network analysis...")
        set.seed(123)
        networkAnalysis <- PiProteline::network_analysis(
          data_grouped            = pp$data_grouped,
          data_grouped_even_dim   = pp$data_grouped_even_dim,
          g_interactome           = g_interactome,
          fun_list = c(Betweenness = igraph::betweenness,
                       Centroids   = PiProteline::centroids,
                       Bridging    = PiProteline::bridging_centrality),
          quantile_critical_nodes = input$na_quantile %||% 0.75,
          names_of_groups         = shared_data$group_names
        )

        setProgress(0.85, detail = "Functional analysis...")
        functionalAnalysis <- PiProteline::functional_analysis(
          dataset                 = pp$data_unique,
          manova_pairwise_results = manova_pw,
          unweighted_CN           = networkAnalysis$Unweighted_criticalNodes,
          weighted_CN             = networkAnalysis$Weighted_criticalNodes,
          names_of_groups         = shared_data$group_names,
          tax_ID                  = 9606,   # human; matches CoPPIs::interactome.hs
          categories              = c("Component","Function","Process","KEGG","RCTM","WikiPathways")
        )

        # Same top-level structure pipeline() returns, so the Functional tab works unchanged.
        res <- list(
          quantitativeAnalysis = quantitativeAnalysis,
          networkAnalysis      = networkAnalysis,
          functionalAnalysis   = functionalAnalysis
        )
        shared_data$pipelineResults <- res
        shared_data$quant_results   <- quantitativeAnalysis  # share back

        # Optional: save to disk like pipeline() does.
        if (nzchar(input$na_save_prefix %||% "")) {
          tryCatch(
            PiProteline::save_results(quantitativeAnalysis, networkAnalysis,
                                      functionalAnalysis,
                                      save_results_as = input$na_save_prefix),
            error = function(e) showNotification(paste("Save skipped:", e$message), type = "warning")
          )
        }

        setProgress(1, detail = "Done.")
        showNotification("Network + functional analysis finished. Functional Analysis tab is ready.",
                         type = "message")
      }, error = function(e) {
        showNotification(paste("Analysis error:", e$message), type = "error")
        shared_data$pipelineResults <- NULL
      })
    })
  })

  # ---- Result viewers (unchanged logic, now data-driven group choices) ----
  current_table <- reactiveVal(NULL)
  current_table_name <- reactiveVal("")

  observeEvent(input$na_show_unweighted_cent, {
    req(shared_data$pipelineResults); grp <- req(input$na_cent_group)
    tab <- get_path(shared_data$pipelineResults, c("networkAnalysis","Unweighted_centralities", grp))
    validate(need(!is.null(tab), "Unweighted_centralities table not found for this group."))
    current_table(as.data.frame(tab)); current_table_name(paste("Unweighted_centralities -", grp))
  })
  observeEvent(input$na_show_weighted_cent, {
    req(shared_data$pipelineResults); grp <- req(input$na_cent_group)
    tab <- get_path(shared_data$pipelineResults, c("networkAnalysis","Weighted_centralities", grp))
    validate(need(!is.null(tab), "Weighted_centralities table not found for this group."))
    current_table(as.data.frame(tab)); current_table_name(paste("Weighted_centralities -", grp))
  })

  output$na_table_title <- renderText({ req(current_table_name()); current_table_name() })
  output$na_table <- renderDT({ req(current_table()); current_table() },
                              options = list(pageLength = 10, scrollX = TRUE))
  output$na_table_dl <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_]+","_", current_table_name()), ".csv"),
    content  = function(file) write.csv(current_table(), file, row.names = FALSE)
  )

  show_cn <- function(weighted = FALSE) {
    req(shared_data$pipelineResults)
    scope <- req(input$na_cn_scope); grp <- req(input$na_cn_group)
    path <- if (!weighted) c("networkAnalysis","Unweighted_criticalNodes", scope, grp)
            else            c("networkAnalysis","Weighted_criticalNodes",   scope, grp)
    obj <- get_path(shared_data$pipelineResults, path)
    validate(need(!is.null(obj), "Critical nodes not found for this selection."))
    rn <- rownames(obj)
    validate(need(length(rn) > 0, "No critical nodes available."))
    current_table(data.frame(Node = rn, stringsAsFactors = FALSE))
    current_table_name(paste(ifelse(weighted,"Weighted","Unweighted"), "criticalNodes -", scope, "-", grp))
  }
  observeEvent(input$na_show_unweighted_cn, { show_cn(FALSE) })
  observeEvent(input$na_show_weighted_cn,   { show_cn(TRUE)  })

  output$na_plot_title <- renderText("")
  output$na_plot <- renderPlot({})

  plot_ppi <- function(kind = c("PPI_unweighted","PPI_correlations")) {
    req(shared_data$pipelineResults); kind <- match.arg(kind)
    grp <- req(input$na_ppi_group)
    gobj <- get_path(shared_data$pipelineResults, c("networkAnalysis", kind, grp))
    validate(need(!is.null(gobj), paste("Graph not found for", kind, "-", grp)))
    output$na_plot_title <- renderText(paste(kind, "-", grp))
    output$na_plot <- renderPlot({ plot(gobj) })
    output$na_plot_dl <- downloadHandler(
      filename = function() paste0(gsub("[^A-Za-z0-9_]+","_", paste(kind, grp)), ".png"),
      content  = function(file) { png(file, 1400, 900, res = 150); plot(gobj); dev.off() }
    )
  }
  observeEvent(input$na_plot_ppi_unweighted,   { plot_ppi("PPI_unweighted")   })
  observeEvent(input$na_plot_ppi_correlations, { plot_ppi("PPI_correlations") })
}
