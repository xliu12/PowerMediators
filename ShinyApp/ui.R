addResourcePath("customwww", "www")
ui <- fluidPage(
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
        tags$title("Power Analysis for Two Mediators")
    ),
    
    fluidRow(
        column(12, align = "center",
               h3("Monte Carlo Power Analysis for Interventional Indirect Effects with Two Mediators"),
               p("Developed by Jasmine Zeng, Xiao Liu")
        )
    ),
    
    tabsetPanel(
        tabPanel("Run Simulation",
                 
                 fluidRow(
                     column(3,
                            selectInput("coef_type", "Coefficient Input Type:",
                                        choices = c("Standardized Coefficient" = "standard", "Raw Coefficient" = "raw")),
                            radioButtons("objective", "Select Objective:",
                                         choices = c("Calculate Power" = "power", "Estimate Sample Size" = "samplesize"),
                                         selected = "power", inline = TRUE),
                            
                            uiOutput("model_panels")  # Dynamically render collapsible panels here
                     ),
                     
                     column(9,
                            actionButton("run", label = textOutput("run_label"), class = "btn-primary", width = "100%"),
                            br(),
                            verbatimTextOutput("check"),
                            DT::dataTableOutput("power_table"),
                            conditionalPanel(
                                condition = "input.objective == 'samplesize'",
                                plotOutput("power_plot", height = "500px"),
                                br(),
                                downloadButton("download_plot", "Download Power Curve Plot")
                            ),
                            br(),
                            downloadButton("download_table", "Download Results Table")
                     )
                 )
        ),
        
        tabPanel("Help",
                 div(style = "text-align: center;",
                     tags$p("Figure: Causal diagram with two correlated mediators."),
                     tags$img(src = "customwww/pic.png", style = "max-height:200px;")
                 ),
                 fluidRow(
                     column(8, offset = 2,
                            tags$h4("Instructions"),
                            tags$ol(
                                tags$li("Set simulation and model parameters using the collapsible sections."),
                                tags$li("Choose input method: Standardized coefficients with R-squared or raw coefficients."),
                                tags$li("Select your objective: calculate power or estimate required sample size."),
                                tags$li("Click 'Run' to perform the simulation and view/download results.")
                            ),
                            tags$br(),
                            tags$p("This Shiny app was designed for conducting power analysis and sample size estimation in two-mediator causal inference models. It supports both binary and continuous outcomes under different multiple testing correction methods."),
                            tags$h4("References")

                            
                     )
                 )
        )
    )
)

