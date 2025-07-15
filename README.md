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
packages <- c("shiny", "tidyverse", "mvtnorm",
"glue", "cubature", "shinyBS",
"DT", "future", "promises", "parallel")

# install a package if it has not been installed previously
lapply(packages, \(pkg) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)) 

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




# R Package: PowerMediators
Power Analysis for Causal Mediation Analysis with Multiple Mediators


# Install and load the package
```{r}
remotes::install_github("xliu12/PowerMediators", subdir = "PowerMediators")

```

# Example: Power Calculation

```{r}
library(mvtnorm)
library(tidyverse)
library(glue)
library(parallel)
library(calculus)

library(PowerMediators)

?PowerMediators::runPower

# other arguments use the default specification
power_res <- runPower(n = 100, nsims = 1000, mc.cores = 5)

power_res

```

The output includes a data.frame containing the statistical power for testing interventional indirect effects (IIE) via each mediator, including (1) the familywise ("FW") power of testing all four indirect effects for a mediator, and (2) the per-test ("PT") power of testing only one indirect effect for a mediator. The data.frame also contains the multiple testing adjustments (sig.adjust), Type I error rate (which is 0 if the true value is non-null), confidence interval coverage rate (cover_FW and cover_PT), the mediator (IIE_M1 for mediator M1 and IIE_M2 for mediator M2), and the indirect effects (IIEs; effect).


# Example: Sample Size Calculation

```{r}
library(mvtnorm)
library(tidyverse)
library(glue)
library(parallel)
library(calculus)

library(PowerMediators)

?PowerMediators::cal.sample_size

# other arguments use the default specification
result <- cal.sample_size(n.min = 200, 
                          n.max = 300, 
                          n.step = 20, 
                          power_target = 0.8,
                          FW_or_PT = c("FW", "PT"),
                          which_mediator = c("M1", "M2"),
                          which_effect = c("IIE_M1(1,,0)", "IIE_M2(1,1,)"),
                          mc.cores = 5, 
                          nsims = 1000) # may take some time to run 

result$res_n

result$power_FW_M1M2
result$power_FW_eachM

result$power_PT_someIIE
result$power_PT_allIIE
result$power_PT_eachIIE

```

The output includes a list containing the sample size calculation results (see `?PowerMediators::cal.sample_size`).



