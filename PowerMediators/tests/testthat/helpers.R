library(mvtnorm)
library(tidyverse)
library(glue)
library(parallel)


# with small number of simulations and bootstrap draws
# with user-specified parameters

# test_run1 <- runPower(
#   n = 200, 
#   num_x = 2,
#   treat.prop = 0.5,
#   treat.randomized = FALSE,
#   a_on_x = 0.1, # standardized coefficient in the treatment model
#   
#   m_on_a = c(0.3, 0.2), # standardized coefficients in the mediators' model
#   m_on_x = c(0.1, 0.1),
#   em_corr = 0,
#   
#   y_on_a = 0.1, # standardized coefficients in the outcome model
#   y_on_m = c(0, 0),
#   y_on_am_2way = c(0, 0.3),
#   y_on_m_2way = 0,
#   y_on_am_3way = 0.3,
#   y_on_x = 0.1,
#   
#   nboot = 10, # default is 1000 bootstrap draws for confidence intervals
#   sig.level = 0.05, # nominal Type I error rate for each mediator
#   
#   nsims = 2, # default is 1000 simulaiton replications
#   mc.cores = 1 # number of cores for parallel computation
# )
# 
# test_run1
