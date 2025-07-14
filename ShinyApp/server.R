library(shiny)
library(shinyBS)

# Core tidyverse packages: includes dplyr, ggplot2, purr.
library(tidyverse)

# Source all R scripts inside the "function/" directory
function_files <- list.files("function", pattern = "\\.R$", full.names = TRUE)
sapply(function_files, source)

# Additional required packages
library(glue)
library(cubature)
library(DT)
library(mvtnorm)     # For multivariate normal functions
library(future)      # For asynchronous processing
library(promises)    # For working with futures
library(parallel)    # For parallel computing
library(shiny)

server <- function(input, output, session) {
    
    result <- reactiveVal(NULL)
    plot_result <- reactiveVal(NULL)
    
    output$run_label <- renderText({
        if (input$objective == "samplesize") "Estimate Sample Size" else "Calculate Power"
    })
    
    output$model_panels <- renderUI({
        bsCollapse(id = "collapse_panels", multiple = TRUE,
                   bsCollapsePanel("Simulation Settings",
                                   tagList(
                                       conditionalPanel(
                                           condition = "input.objective == 'power'",
                                           numericInput("n", "Sample Size", value = 50)
                                       ),
                                       conditionalPanel(
                                           condition = "input.objective == 'samplesize'",
                                           tagList(
                                               numericInput("TarPow", "Target Power", value = 0.8),
                                               numericInput("n", "Starting Sample Size", value = 50),
                                               numericInput("steps", "Sample Size Step", value = 20),
                                               numericInput("max_n", "Max Sample Size", value = 200),
                                               checkboxGroupInput("sig.adjust", "Significance Adjustment Method:",
                                                                  choices = c("no_adjust", "bonferroni", "modified_bon1", "modified_bon2"),
                                                                  selected = c("bonferroni")),
                                               checkboxGroupInput("mediation", "Mediation Effects:",
                                                                  choices = c("IIE_M1", "IIE_M2"), selected = c("IIE_M1")),
                                               checkboxGroupInput("power", "Power Type:",
                                                                  choices = c("familywise", "per-test"), selected = c("familywise")),
                                               checkboxGroupInput("effect", "Effects to Estimate:",
                                                                  choices = c("all", "IIE_M1(1,,1)", "IIE_M1(1,,0)", "IIE_M1(0,,1)",
                                                                              "IIE_M1(0,,0)", "IIE_M2(1,1,)", "IIE_M2(1,0,)", "IIE_M2(0,1,)", "IIE_M2(0,0,)"),
                                                                  selected = c("all"))
                                           )
                                       ),
                                       numericInput("nsims", "Number of Simulations", value = 1000),
                                       numericInput("seed", "Random Seed", value = 1234),
                                       numericInput("conf", "Significance Level (%)", value = 5, min = 1, max = 99),
                                       radioButtons("ci_type", "Confidence Interval Method:",
                                                    choices = c("Monte Carlo" = 'mc', "Bootstrap" = 'bootstrap'),
                                                    selected = 'mc', inline = TRUE),
                                       numericInput("n_rep", "Number of Replications", value = 1000, min = 1)
                                   )
                   ),
                   
                   bsCollapsePanel("Treatment & Covariates",
                                   tagList(
                                       numericInput("num_x", "Number of Covariates (X)", value = 2),
                                       numericInput("treat.prop", "Proportion of Treatment (A)", value = 0.5, min = 0, max = 1, step = 0.01),
                                       checkboxInput("treat.randomized", "Randomized Treatment Assignment", value = TRUE),
                                       checkboxInput("M1_binary", "Binary Mediator 1", value = FALSE),
                                       checkboxInput("M2_binary", "Binary Mediator 2", value = FALSE),
                                       checkboxInput("Y_binary", "Binary Outcome", value = FALSE)
                                   )
                   ),
                   
                   bsCollapsePanel("Exposure Model Parameters",
                                   if (input$coef_type == "raw") {
                                       numericInput("a_on_x", "Raw Coefficient: X → A", value = 0)
                                   } else {
                                       numericInput("R2.ax", "R-squared: X → A", value = 0)
                                   }
                   ),
                   
                   bsCollapsePanel("Mediators Model Parameters",
                                   tagList(
                                       if (input$coef_type == "raw") {
                                           tagList(
                                               numericInput("m1_on_a", "Raw Coefficient: A → M₁", value = 0.27),
                                               numericInput("m1_on_x", "Raw Coefficient: X → M₁", value = 0.1),
                                               numericInput("m2_on_a", "Raw Coefficient: A → M₂", value = 0.36),
                                               numericInput("m2_on_x", "Raw Coefficient: X → M₂", value = 0.1)
                                           )
                                       } else {
                                           tagList(
                                               numericInput("R2.m1x", "R-squared: X → M₁", value = 0.02),
                                               numericInput("R2.m2x", "R-squared: X → M₂", value = 0.02),
                                               numericInput("std.m1_on_a", "Standardized Coefficient: A → M₁", value = 0.27),
                                               numericInput("std.m2_on_a", "Standardized Coefficient: A → M₂", value = 0.36)
                                           )
                                       },
                                       numericInput("em_corr", "Residual Corr(M₁, M₂)", value = 0.1)
                                   )
                   ),
                   
                   bsCollapsePanel("Outcome Model Parameters",
                                   if (input$coef_type == "raw") {
                                       tagList(
                                           numericInput("y_on_x", "Raw Coefficient: X → Y", value = 0.1),
                                           numericInput("y_on_a", "Raw Coefficient: A → Y", value = 0.36),
                                           numericInput("y_on_m1", "Raw Coefficient: M₁ → Y", value = 0.36),
                                           numericInput("y_on_m2", "Raw Coefficient: M₂ → Y", value = 0.36),
                                           numericInput("y_on_am1_2way", "Interaction: A × M₁ → Y", value = 0.1),
                                           numericInput("y_on_am2_2way", "Interaction: A × M₂ → Y", value = 0.1),
                                           numericInput("y_on_m_2way", "Interaction: M₁ × M₂ → Y", value = 0.1),
                                           numericInput("y_on_am_3way", "Interaction: A × M₁ × M₂ → Y", value = 0.1)
                                       )
                                   } else {
                                       tagList(
                                           numericInput("R2.yx", "R-squared: X → Y", value = 0.02),
                                           numericInput("std.y_on_a", "Standardized Coefficient: A → Y", value = 0.03),
                                           numericInput("std.y_on_m1", "Standardized Coefficient: M₁ → Y", value = 0.25),
                                           numericInput("std.y_on_m2", "Standardized Coefficient: M₂ → Y", value = 0.15),
                                           numericInput("std.y_on_am1_2way", "Standardized Interaction: A × M₁ → Y", value = 0.1),
                                           numericInput("std.y_on_am2_2way", "Standardized Interaction: A × M₂ → Y", value = 0.1),
                                           numericInput("std.y_on_m_2way", "Standardized Interaction: M₁ × M₂ → Y", value = 0.1),
                                           numericInput("std.y_on_am_3way", "Standardized Interaction: A × M₁ × M₂ → Y", value = 0.1)
                                       )
                                   }
                   )
        )
    })
    observeEvent(input$run, {
        result(NULL)
        plot_result(NULL)
        
        progress <- Progress$new(session, min = 0, max = 1)
        on.exit(progress$close())
        progress$set(message = "Running simulation...", value = 0.1)
        
        args <- list(
            seed_user = input$seed,
            num_x = input$num_x,
            treat.prop = input$treat.prop,
            treat.randomized = input$treat.randomized,
            em_corr = input$em_corr,
            M_binary = c(input$M1_binary, input$M2_binary),
            Y_binary = input$Y_binary,
            sig.level = input$conf / 100,
            nsims = input$nsims,
            mc.cores = 1,
            n = input$n
        )
        
        if (input$ci_type == "mc") {
            args$n.draws <- input$n_rep
            args$nboot <- NULL
        } else if (input$ci_type == "bootstrap") {
            args$nboot <- input$n_rep
            args$n.draws <- NULL
        }
        
        if (input$coef_type == "standard") {
            args <- c(args, list(
                R2.ax = input$R2.ax,
                std.m_on_a = c(input$std.m1_on_a, input$std.m2_on_a),
                R2.mx = c(input$R2.m1x, input$R2.m2x),
                std.y_on_a = input$std.y_on_a,
                std.y_on_m = c(input$std.y_on_m1, input$std.y_on_m2),
                std.y_on_am_2way = c(input$std.y_on_am1_2way, input$std.y_on_am2_2way),
                std.y_on_m_2way = input$std.y_on_m_2way,
                std.y_on_am_3way = input$std.y_on_am_3way,
                R2.yx = input$R2.yx,
                a_on_x = NULL,
                m_on_a = NULL,
                m_on_x = NULL,
                y_on_a = NULL,
                y_on_m = NULL,
                y_on_am_2way = NULL,
                y_on_m_2way = NULL,
                y_on_am_3way = NULL,
                y_on_x = NULL
            ))
        } else {
            args <- c(args, list(
                a_on_x = input$a_on_x,
                m_on_a = c(input$m1_on_a, input$m2_on_a),
                m_on_x = c(input$m1_on_x, input$m2_on_x),
                y_on_a = input$y_on_a,
                y_on_m = c(input$y_on_m1, input$y_on_m2),
                y_on_am_2way = c(input$y_on_am1_2way, input$y_on_am2_2way),
                y_on_m_2way = input$y_on_m_2way,
                y_on_am_3way = input$y_on_am_3way,
                y_on_x = input$y_on_x,
                R2.ax = NULL,
                std.m_on_a = NULL,
                R2.mx = NULL,
                std.y_on_a =NULL,
                std.y_on_m = NULL,
                std.y_on_am_2way = NULL,
                std.y_on_m_2way = NULL,
                std.y_on_am_3way = NULL,
                R2.yx = NULL,
            ))
        }
        
        if (input$objective == "power") {
            tryCatch({
                cat("=== Arguments Passed to Function ===\n")
                print(args)
                res <- do.call(runPower, args)
                result(res)
            }, error = function(e) {
                showNotification(paste("Error in runPower:", e$message), type = "error")
            })
        } else {
            args <- c(args, list(
                steps = input$steps,
                TarPow = input$TarPow,
                max_n = input$max_n,
                sig.adjust = input$sig.adjust,
                mediation = input$mediation,
                effect = input$effect,
                power = input$power,
                plot = TRUE,
                verbose = TRUE
            ))
            tryCatch({
                cat("=== Arguments Passed to Function ===\n")
                print(args)
                res <- do.call(runSampsize, args)
                result(res$target_power_table)
                plot_result(res$plot)
                output_text <- capture.output(print(res)) 
                if (any(grepl("Target power not achieved", output_text))) {
                    showNotification("Target power not achieved. Max sample size used.", type = "warning", duration = 10)
                }
            }, error = function(e) {
                showNotification(paste("Error in runSampsize:", e$message), type = "error")
            })
        }
        
        progress$set(value = 1)
    })
    
    output$power_table <- DT::renderDataTable({
        req(result())
        DT::datatable(result(), options = list(scrollX = TRUE, scrollY = TRUE))
    })
    
    output$power_plot <- renderPlot({
        req(plot_result())
        print(plot_result())
    })
    
    output$download_table <- downloadHandler(
        filename = function() paste0("Power_Table_", Sys.Date(), ".csv"),
        content = function(file) {
            write.csv(result(), file, row.names = FALSE)
        }
    )
    
    output$download_plot <- downloadHandler(
        filename = function() paste0("Power_Curve_", Sys.Date(), ".png"),
        content = function(file) {
            ggsave(file, plot = plot_result(), width = 7, height = 5)
        }
    )
}


