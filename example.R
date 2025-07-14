library(mvtnorm)
library(tidyverse)
library(glue)
library(parallel)
library(calculus)

# download and unzip the package; and run
# devtools::load_all("PowerMediators")


# Alternatively: 
# # install from github
# remotes::install_github("xliu12/PowerMediators", subdir = "PowerMediators")
# library(PowerMediators)


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


?PowerMediators::runPower

# other arguments use the default specification
power_res <- runPower(n = 100, nsims = 1000, mc.cores = 5)

power_res
