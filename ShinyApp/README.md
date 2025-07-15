# 📈 Shiny App: Power Analysis for Two-Mediator Causal Mediation Models

This Shiny app provides **simulation-based power analysis** for causal mediation models with **two mediators**, supporting:

-   ✅ Familywise and per-test power estimation\
-   ✅ Multiple testing corrections\
-   ✅ Flexible model specification for mediator/outcome type\
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

Make sure you're in the PowerMediators/ root folder, then run:

``` r
shiny::runApp("ShinyApp")
```

### 🔹 Option 2: Run Directly from GitHub

``` r
shiny::runGitHub("PowerMediators", "xliu12", subdir = "ShinyApp")
```
