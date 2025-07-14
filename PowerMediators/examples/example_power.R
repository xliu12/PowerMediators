
\donttest{
    library(mvtnorm)
    library(tidyverse)
    library(glue)
    library(parallel)
    library(calculus)

    # other arguments use the default specification
    power_res <- runPower(n = 100, nsims = 1000, mc.cores = 5)

    power_res
}




