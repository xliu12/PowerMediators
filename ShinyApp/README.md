# 📈 Shiny App: Power Analysis for Two-Mediator Causal Mediation Models

This Shiny app provides **simulation-based power analysis** for causal mediation models with **two mediators**, supporting:

-   ✅ Familywise and per-test power estimation
-   ✅ Multiple testing corrections
-   ✅ Flexible model specification for mediator/outcome type
-   ✅ Informative plots and output tables

------------------------------------------------------------------------

## 🚀 How to Run the App

### 🔹 Option 1: Run Locally

#### Step 1: Install Required Packages

``` r
install.packages(c(
  "shiny", "tidyverse", "mvtnorm", "glue", "cubature",
  "shinyBS", "DT", "future", "promises", "parallel"
))
```

#### Step 2: Run from Local Folder

Navigate to the root folder (e.g., PowerMediators/) in your R console or RStudio, then run:


``` r
shiny::runApp("./PowerMediators/ShinyApp")
```
Alternatively, you can open the app.R file located in the ShinyApp/ folder in RStudio and run it.


### 🔹 Option 2: Run Directly from GitHub

``` r
shiny::runGitHub("PowerMediators", "xliu12", subdir = "ShinyApp")
```
