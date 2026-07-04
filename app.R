library(shiny)

# helpers.R must be sourced FIRST: it defines %||% and other utilities that the
# tab server functions rely on.
source("helpers.R")
source("preprocessing_tab.R")
source("quantitative_analysis_tab.R")
source("network_analysis_tab.R")
source("functional_analysis_tab.R")

ui <- fluidPage(
  titlePanel("PiProteline Shiny App"),
  tabsetPanel(
    preprocessing_ui(),
    quantitative_analysis_ui(),
    network_analysis_ui(),
    functional_analysis_ui()
  )
)

server <- function(input, output, session) {
  # The preprocessing tab owns the single shared_data store; every other tab
  # reads from it. This is what keeps the user's choices flowing through.
  shared_data <- preprocessing_server(input, output, session)

  quantitative_analysis_server(input, output, session, shared_data)
  network_analysis_server(input, output, session, shared_data)
  functional_analysis_server(input, output, session, shared_data)
}

shinyApp(ui = ui, server = server)
