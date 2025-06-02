library(mvtnorm)
library(tidyverse)
library(glue)
library(parallel)
library(integral)



# download and unzip the package; and run
devtools::load_all("PowerMediators")

# Alternatively: 
# # install from github
# remotes::install_github("xliu12/PowerMediators", subdir = "PowerMediators")
# library(PowerMediators)

# with small number of simulations and bootstrap draws
# with user-specified parameters

test_run1 <- runPower(
  n = 200,
  num_x = 2,
  treat.prop = 0.5,
  treat.randomized = FALSE,
  
  R2.ax = 0.1, R2.mx = 0.1, R2.yx = 0.1, 
  std.m_on_a = c(0.36, 0.36),
  std.y_on_a = 0.14,
  std.y_on_m = c(0, 0),
  std.y_on_am_2way = c(0, 0.36),
  std.y_on_m_2way = 0,
  std.y_on_am_3way = 0.36,
  
  # a_on_x = sqrt(0.13), 
  # m_on_a = rep(0.14, 2), 
  # m_on_x = rep(sqrt(0.13), 2),
  em_corr = 0,
  M_binary = c(TRUE, FALSE),
  
  # y_on_a = 0.1, # standardized coefficient
  # y_on_m = rep(0.39, 2),
  # y_on_am_2way = rep(0.1, 2),
  # y_on_m_2way = rep(0.1, 2 * (2-1)/2),
  # y_on_am_3way = rep(0, 2 * (2-1)/2),
  # y_on_x = sqrt(0.02),
  Y_binary = FALSE,
  
  # for confidence interval and estimation
  nboot = NULL, # for bootstrap CI
  n.draws = 1000, # for Monte Carlo CI
  
  # for power
  sig.level = 0.05,
  # for simulation-based power calculation
  nsims = 2, mc.cores = 1
)

test_run1
