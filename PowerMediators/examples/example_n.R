
\donttest{
    library(mvtnorm)
    library(tidyverse)
    library(glue)
    library(parallel)
    library(calculus)

    # other arguments use the default specification
    result <- cal.sample_size(n.min = 200,
                              n.max = 300,
                              n.step = 20,
                              power_target = 0.8,
                              FW_or_PT = c("FW", "PT"),
                              which_mediator = c("M1", "M2"),
                              which_effect = c("IIE_M1(1,,0)", "IIE_M2(1,1,)"),
                              mc.cores = 5,
                              nsims = 1000)

    result$res_n

    result$power_FW_M1M2
    result$power_FW_eachM

    result$power_PT_someIIE
    result$power_PT_allIIE
    result$power_PT_eachIIE
}




