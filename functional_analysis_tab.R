# =============================
# Functional Analysis Tab 
# =============================
library(shiny)
library(shinyWidgets)
library(CoPPIs)
library(igraph)
library(dplyr)
library(PiProteline)
library(ggplot2)
library(DT)
library(grid)

functional_analysis_ui <- function() {
  tabPanel(
    "Functional Analysis",
    sidebarLayout(
      sidebarPanel(
        
        tags$strong("Note:"), helpText("Run the pipeline in the Network Analysis tab first."),
        tags$hr(),
        
        
        uiOutput("fa_key_ui"),
        uiOutput("fa_cat_ui"),
        actionButton("fa_show", "Show selected MANOVA plot"),
        
        tags$hr(), h4("Tables"),
        actionButton("fa_show_spe", "Show singleProfileEnrichment"),
        tags$br(), tags$br(),
        
        uiOutput("ui_tbl_enrmanova"),
        actionButton("btn_tbl_enrmanova", "Show enr_manova"),
        tags$br(), tags$br(),
        
        uiOutput("ui_tbl_enrmanova_diff"),
        actionButton("btn_tbl_enrmanova_diff", "Show enr_manova.diff"),
        tags$br(), tags$br(),
        
        uiOutput("ui_tbl_enrNetw"),
        actionButton("btn_tbl_enrNetw", "Show enr_Netw"),
        tags$br(), tags$br(),
        
        uiOutput("ui_tbl_enrNetw_diff"),
        actionButton("btn_tbl_enrNetw_diff", "Show enr_Netw.diff"),
        tags$br(), tags$br(),
        
        uiOutput("ui_tbl_enr_single_Netw"),
        actionButton("btn_tbl_enr_single_Netw", "Show enr_single_Netw"),
        tags$br(), tags$br(),
        
        uiOutput("ui_tbl_enr_manovanetw"),
        actionButton("btn_tbl_enr_manovanetw", "Show enr_manovanetw"),
        tags$hr(),
        
        h4("Extra plots"),
        uiOutput("ui_plot_network"),
        uiOutput("ui_plot_network_cat"),
        actionButton("btn_plot_network", "Show network plot"),
        tags$br(), tags$br(),
        
        uiOutput("ui_plot_singleProfile"),
        actionButton("btn_plot_singleProfile", "Show singleProfile plot")
      ),
      mainPanel(
        verbatimTextOutput("fa_debug"),
        
        h4(textOutput("fa_title")),
        plotOutput("fa_plot_single", height = 500),
        downloadButton("fa_plot_single_dl", "Download PNG"),
        tags$hr(),
        
        h4(textOutput("fa_table_title")),
        DTOutput("fa_table"),
        downloadButton("fa_table_dl", "Download CSV"),
        tags$hr(),
        
        h4(textOutput("ttl_plot_network")),
        plotOutput("fa_plot_network", height = 500),
        downloadButton("fa_plot_network_dl", "Download Network Plot"),
        tags$hr(),
        
        h4(textOutput("ttl_plot_single")),
        plotOutput("fa_plot_singleProfile", height = 500),
        downloadButton("fa_plot_singleProfile_dl", "Download Single-Profile Plot")
      )
    )
  )
}

functional_analysis_server <- function(input, output, session, shared_data) {
  
  
  get_path <- function(x, path) {
    out <- x
    for (nm in path) {
      if (is.null(out) || is.null(out[[nm]])) return(NULL)
      out <- out[[nm]]
    }
    out
  }
  extract_drawable <- function(x) {
    if (is.null(x)) return(NULL)
    if (inherits(x, c("gg","ggplot"))) return(x)
    if (inherits(x, c("grob","gTree","gList","gtable"))) return(x)
    if (is.function(x)) return(x)
    if (is.list(x)) {
      if (!is.null(x$plot)) {
        out <- extract_drawable(x$plot); if (!is.null(out)) return(out)
      }
      for (nm in names(x)) {
        out <- extract_drawable(x[[nm]]); if (!is.null(out)) return(out)
      }
    }
    NULL
  }
  draw_any <- function(obj) {
    if (is.null(obj)) { plot.new(); text(0.5,0.5,"No plot", cex=1.2)
    } else if (inherits(obj, c("gg","ggplot"))) {
      print(obj)
    } else if (inherits(obj, c("grob","gTree","gList","gtable"))) {
      grid::grid.newpage(); grid::grid.draw(obj)
    } else if (is.function(obj)) {
      obj()
    } else {
      plot.new(); text(0.5,0.5,paste("Unsupported type:", class(obj)[1]), cex=1.1)
    }
  }
  canon_order <- c("Component","Function","KEGG","Process","Reactome","WikiPathways","All")
  normalize_category <- function(x) {
    lx <- tolower(x)
    dplyr::case_when(
      grepl("^comp", lx) ~ "Component",
      grepl("^func|^mf|go:mf|molecular", lx) ~ "Function",
      grepl("^proc|^bp|go:bp|process", lx) ~ "Process",
      grepl("kegg", lx) ~ "KEGG",
      grepl("react|rctm", lx) ~ "Reactome",
      grepl("wiki", lx) ~ "WikiPathways",
      TRUE ~ x
    )
  }
  list_categories <- function(node) {
    if (is.list(node)) {
      nms <- names(node); nms <- nms[!is.na(nms) & nzchar(nms)]
      if (length(nms)) {
        out <- unique(normalize_category(nms))
        out <- unique(c(canon_order[canon_order %in% out], setdiff(out, canon_order)))
        return(out)
      }
    }
    "All"
  }
  pick_category_node <- function(node, category) {
    if (identical(category, "All")) return(extract_drawable(node))
    if (!is.list(node)) return(extract_drawable(node))
    nms <- names(node); if (length(nms)) {
      norm <- normalize_category(nms); hit <- which(norm == category)
      if (length(hit)) return(extract_drawable(node[[hit[1]]]))
    }
    extract_drawable(node)
  }
  flatten_tables <- function(x, prefix = "") {
    out <- list()
    if (is.data.frame(x)) { nm <- if (nzchar(prefix)) prefix else "table"; out[[nm]] <- x; return(out) }
    if (is.list(x)) {
      nms <- names(x)
      for (i in seq_along(x)) {
        key <- if (!is.null(nms) && nzchar(nms[i])) nms[i] else paste0("idx", i)
        out <- c(out, flatten_tables(x[[i]], if (nzchar(prefix)) paste0(prefix, ".", key) else key))
      }
    }
    out
  }
  load_table_by_key <- function(container_path, key) {
    obj <- get_path(shared_data$pipelineResults, container_path)
    if (is.null(obj)) { path2 <- container_path; path2[1] <- ifelse(path2[1]=="functionalAnalysis","functional_analysis","functionalAnalysis"); obj <- get_path(shared_data$pipelineResults, path2) }
    flat <- flatten_tables(obj); if (length(flat) == 0) return(NULL)
    flat[[key]]
  }
  
  
  observeEvent(shared_data$pipelineResults, {
    res <- shared_data$pipelineResults
    if (is.null(res)) {
      showNotification("No pipeline results found. Run it in Network Analysis tab.", type = "warning")
      return()
    }
    
    
    mb <- get_path(res, c("functionalAnalysis","plot_enr","manova"))
    if (is.null(mb)) mb <- get_path(res, c("functional_analysis","plot_enr","manova"))
    keys <- if (!is.null(mb)) names(mb) else character(0)
    output$fa_key_ui <- renderUI({
      selectInput("fa_key", "MANOVA contrast", choices = keys, selected = if (length(keys)) keys[1] else NULL)
    })
    
    if (length(keys)) {
      cats <- list_categories(mb[[keys[1]]])
      output$fa_cat_ui <- renderUI({ selectInput("fa_cat", "Category", choices = cats, selected = cats[1]) })
    } else {
      output$fa_cat_ui <- renderUI(NULL)
    }
    
    
    enrmanova        <- get_path(res, c("functionalAnalysis","enrmanova","enr_manova"))
    if (is.null(enrmanova)) enrmanova <- get_path(res, c("functional_analysis","enrmanova","enr_manova"))
    output$ui_tbl_enrmanova <- renderUI({ selectInput("sel_tbl_enrmanova", "enr_manova key", choices = names(flatten_tables(enrmanova))) })
    
    enrmanova_diff   <- get_path(res, c("functionalAnalysis","enrmanova","enr_manova.diff"))
    if (is.null(enrmanova_diff)) enrmanova_diff <- get_path(res, c("functional_analysis","enrmanova","enr_manova.diff"))
    output$ui_tbl_enrmanova_diff <- renderUI({ selectInput("sel_tbl_enrmanova_diff", "enr_manova.diff key", choices = names(flatten_tables(enrmanova_diff))) })
    
    enrNetw <- get_path(res, c("functionalAnalysis","enrNetw","enr_Netw"))
    if (is.null(enrNetw)) enrNetw <- get_path(res, c("functional_analysis","enrNetw","enr_Netw"))
    output$ui_tbl_enrNetw <- renderUI({ selectInput("sel_tbl_enrNetw", "enr_Netw key", choices = names(flatten_tables(enrNetw))) })
    
    enrNetw_diff <- get_path(res, c("functionalAnalysis","enrNetw","enr_Netw.diff"))
    if (is.null(enrNetw_diff)) enrNetw_diff <- get_path(res, c("functional_analysis","enrNetw","enr_Netw.diff"))
    output$ui_tbl_enrNetw_diff <- renderUI({ selectInput("sel_tbl_enrNetw_diff", "enr_Netw.diff key", choices = names(flatten_tables(enrNetw_diff))) })
    
    enr_single_Netw <- get_path(res, c("functionalAnalysis","enrNetw","enr_single_Netw"))
    if (is.null(enr_single_Netw)) enr_single_Netw <- get_path(res, c("functional_analysis","enrNetw","enr_single_Netw"))
    output$ui_tbl_enr_single_Netw <- renderUI({ selectInput("sel_tbl_enr_single_Netw", "enr_single_Netw key", choices = names(enr_single_Netw)) })
    
    enr_manovanetw <- get_path(res, c("functionalAnalysis","enr_manovanetw"))
    if (is.null(enr_manovanetw)) enr_manovanetw <- get_path(res, c("functional_analysis","enr_manovanetw"))
    output$ui_tbl_enr_manovanetw <- renderUI({ selectInput("sel_tbl_enr_manovanetw", "enr_manovanetw key", choices = names(enr_manovanetw)) })
    
    
    plot_network <- get_path(res, c("functionalAnalysis","plot_enr","network"))
    if (is.null(plot_network)) plot_network <- get_path(res, c("functional_analysis","plot_enr","network"))
    output$ui_plot_network <- renderUI({
      selectInput("sel_plot_network", "Network plot key", choices = names(plot_network))
    })
    output$ui_plot_network_cat <- renderUI(NULL)
    
    plot_single <- get_path(res, c("functionalAnalysis","plot_enr","singleProfile"))
    if (is.null(plot_single)) plot_single <- get_path(res, c("functional_analysis","plot_enr","singleProfile"))
    output$ui_plot_singleProfile <- renderUI({
      selectInput("sel_plot_singleProfile", "Single-Profile plot key", choices = names(plot_single))
    })
  }, ignoreInit = FALSE)
  
  
  observeEvent(input$fa_key, {
    req(shared_data$pipelineResults, input$fa_key)
    res <- shared_data$pipelineResults
    mb  <- get_path(res, c("functionalAnalysis","plot_enr","manova"))
    if (is.null(mb)) mb <- get_path(res, c("functional_analysis","plot_enr","manova"))
    node <- mb[[input$fa_key]]
    cats <- list_categories(node)
    output$fa_cat_ui <- renderUI({ selectInput("fa_cat", "Category", choices = cats, selected = cats[1]) })
  }, ignoreInit = TRUE)
  
  
  output$fa_debug <- renderPrint({
    list(
      have_pipeline = !is.null(shared_data$pipelineResults),
      manova_key    = input$fa_key %||% NA,
      manova_cat    = input$fa_cat %||% NA
    )
  })
  
  
  output$fa_title <- renderText({
    if (!is.null(input$fa_key) && !is.null(input$fa_cat)) paste(input$fa_key, "—", input$fa_cat) else ""
  })
  observeEvent(input$fa_show, {
    req(shared_data$pipelineResults, input$fa_key, input$fa_cat)
    res <- shared_data$pipelineResults
    base <- get_path(res, c("functionalAnalysis","plot_enr","manova", input$fa_key))
    if (is.null(base)) base <- get_path(res, c("functional_analysis","plot_enr","manova", input$fa_key))
    validate(need(!is.null(base), "Selected contrast not found."))
    plt <- pick_category_node(base, input$fa_cat)
    validate(need(!is.null(plt), "No drawable plot for the selected contrast/category."))
    output$fa_plot_single <- renderPlot({ draw_any(plt) })
    key_part <- gsub("[^A-Za-z0-9_]+","_", input$fa_key)
    cat_part <- gsub("[^A-Za-z0-9_]+","_", input$fa_cat)
    output$fa_plot_single_dl <- downloadHandler(
      filename = function() paste0(key_part, "_", cat_part, ".png"),
      content  = function(file) { png(file, 1400, 900, res = 150); draw_any(plt); dev.off() }
    )
  })
  
  
  current_table <- reactiveVal(NULL)
  current_table_name <- reactiveVal("")
  
  observeEvent(input$fa_show_spe, {
    req(shared_data$pipelineResults)
    tab <- get_path(shared_data$pipelineResults, c("functionalAnalysis","singleProfileEnrichment"))
    if (is.null(tab)) tab <- get_path(shared_data$pipelineResults, c("functional_analysis","singleProfileEnrichment"))
    validate(need(!is.null(tab), "singleProfileEnrichment not found."))
    current_table(as.data.frame(tab)); current_table_name("singleProfileEnrichment")
  })
  observeEvent(input$btn_tbl_enrmanova, {
    req(input$sel_tbl_enrmanova)
    tab <- load_table_by_key(c("functionalAnalysis","enrmanova","enr_manova"), input$sel_tbl_enrmanova)
    validate(need(!is.null(tab), "Selected enr_manova table not found."))
    current_table(as.data.frame(tab)); current_table_name(paste0("enr_manova: ", input$sel_tbl_enrmanova))
  })
  observeEvent(input$btn_tbl_enrmanova_diff, {
    req(input$sel_tbl_enrmanova_diff)
    tab <- load_table_by_key(c("functionalAnalysis","enrmanova","enr_manova.diff"), input$sel_tbl_enrmanova_diff)
    validate(need(!is.null(tab), "Selected enr_manova.diff table not found."))
    current_table(as.data.frame(tab)); current_table_name(paste0("enr_manova.diff: ", input$sel_tbl_enrmanova_diff))
  })
  observeEvent(input$btn_tbl_enrNetw, {
    req(input$sel_tbl_enrNetw)
    tab <- load_table_by_key(c("functionalAnalysis","enrNetw","enr_Netw"), input$sel_tbl_enrNetw)
    validate(need(!is.null(tab), "Selected enr_Netw table not found."))
    current_table(as.data.frame(tab)); current_table_name(paste0("enr_Netw: ", input$sel_tbl_enrNetw))
  })
  observeEvent(input$btn_tbl_enrNetw_diff, {
    req(input$sel_tbl_enrNetw_diff)
    tab <- load_table_by_key(c("functionalAnalysis","enrNetw","enr_Netw.diff"), input$sel_tbl_enrNetw_diff)
    validate(need(!is.null(tab), "Selected enr_Netw.diff table not found."))
    current_table(as.data.frame(tab)); current_table_name(paste0("enr_Netw.diff: ", input$sel_tbl_enrNetw_diff))
  })
  observeEvent(input$btn_tbl_enr_single_Netw, {
    req(input$sel_tbl_enr_single_Netw)
    tab <- load_table_by_key(c("functionalAnalysis","enrNetw","enr_single_Netw"), input$sel_tbl_enr_single_Netw)
    validate(need(!is.null(tab), "Selected enr_single_Netw table not found."))
    current_table(as.data.frame(tab)); current_table_name(paste0("enr_single_Netw: ", input$sel_tbl_enr_single_Netw))
  })
  observeEvent(input$btn_tbl_enr_manovanetw, {
    req(input$sel_tbl_enr_manovanetw)
    tab <- load_table_by_key(c("functionalAnalysis","enr_manovanetw"), input$sel_tbl_enr_manovanetw)
    validate(need(!is.null(tab), "Selected enr_manovanetw table not found."))
    current_table(as.data.frame(tab)); current_table_name(paste0("enr_manovanetw: ", input$sel_tbl_enr_manovanetw))
  })
  output$fa_table_title <- renderText({ req(current_table_name()); current_table_name() })
  output$fa_table <- renderDT({ req(current_table()); current_table() }, options = list(pageLength = 10, scrollX = TRUE))
  output$fa_table_dl <- downloadHandler(
    filename = function() paste0(gsub("[^A-Za-z0-9_]+","_", current_table_name()), ".csv"),
    content  = function(file) write.csv(current_table(), file, row.names = FALSE)
  )
  
  
  observeEvent(input$sel_plot_network, {
    req(shared_data$pipelineResults, input$sel_plot_network)
    res  <- shared_data$pipelineResults
    base <- get_path(res, c("functionalAnalysis","plot_enr","network", input$sel_plot_network))
    if (is.null(base)) base <- get_path(res, c("functional_analysis","plot_enr","network", input$sel_plot_network))
    cats <- list_categories(base)
    output$ui_plot_network_cat <- renderUI({ selectInput("sel_plot_network_cat", "Network category", choices = cats, selected = cats[1]) })
  }, ignoreInit = TRUE)
  
  observeEvent(input$btn_plot_network, {
    req(shared_data$pipelineResults, input$sel_plot_network)
    res  <- shared_data$pipelineResults
    base <- get_path(res, c("functionalAnalysis","plot_enr","network", input$sel_plot_network))
    if (is.null(base)) base <- get_path(res, c("functional_analysis","plot_enr","network", input$sel_plot_network))
    validate(need(!is.null(base), "Selected network plot key not found."))
    node <- if (!is.null(input$sel_plot_network_cat)) pick_category_node(base, input$sel_plot_network_cat) else extract_drawable(base)
    validate(need(!is.null(node), "No drawable plot found for the selected network key/category."))
    ttl <- if (!is.null(input$sel_plot_network_cat)) paste("Network plot:", input$sel_plot_network, "—", input$sel_plot_network_cat) else paste("Network plot:", input$sel_plot_network)
    output$ttl_plot_network <- renderText(ttl)
    output$fa_plot_network <- renderPlot({ draw_any(node) })
    output$fa_plot_network_dl <- downloadHandler(
      filename = function() {
        key_part <- gsub("[^A-Za-z0-9_]+","_", input$sel_plot_network)
        cat_part <- if (!is.null(input$sel_plot_network_cat)) paste0("_", gsub("[^A-Za-z0-9_]+","_", input$sel_plot_network_cat)) else ""
        paste0("network_", key_part, cat_part, ".png")
      },
      content  = function(file) { png(file, 1400, 900, res = 150); draw_any(node); dev.off() }
    )
  })
  
  observeEvent(input$btn_plot_singleProfile, {
    req(shared_data$pipelineResults, input$sel_plot_singleProfile)
    res  <- shared_data$pipelineResults
    node <- get_path(res, c("functionalAnalysis","plot_enr","singleProfile", input$sel_plot_singleProfile))
    if (is.null(node)) node <- get_path(res, c("functional_analysis","plot_enr","singleProfile", input$sel_plot_singleProfile))
    validate(need(!is.null(node), "Selected singleProfile plot not found."))
    output$ttl_plot_single <- renderText(paste("Single-profile plot:", input$sel_plot_singleProfile))
    output$fa_plot_singleProfile <- renderPlot({ draw_any(extract_drawable(node)) })
    output$fa_plot_singleProfile_dl <- downloadHandler(
      filename = function() paste0("singleProfile_", gsub("[^A-Za-z0-9_]+","_", input$sel_plot_singleProfile), ".png"),
      content  = function(file) { png(file, 1400, 900, res = 150); draw_any(extract_drawable(node)); dev.off() }
    )
  })
}
