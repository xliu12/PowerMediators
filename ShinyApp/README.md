# Two-Med-Causal-Power-Analysis

Shiny App: Power Analysis for Two-Mediator Causal Mediation Models

The online app can be accessed here: [https://jasminezqy.shinyapps.io/shinyapp/](https://jasminezqy.shinyapps.io/shinyapp/)  


The app can also be run locally using the following code in **R**:

```r
# Install required packages
install.packages(c(
  "shiny", "tidyverse", "mvtnorm", "glue", "cubature",
  "shinyBS", "DT", "future", "promises", "parallel"
))


# Run the app (if downloaded from GitHub)
shiny::runApp("path/to/downloaded/folder")

# Or run the app directly from the R console
shiny::runGitHub("Causal-Two-Mediator-power-Rshiny", "JasmineZqy")


