# Set www directory for custom resources
addResourcePath("customwww", "www")

# Run the Shiny app
source("ui.R")
source("server.R")
shinyApp(ui = ui, server = server)


