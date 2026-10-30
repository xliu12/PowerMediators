# The empirical example specification

EXAMPLE_SEED  <- 20261002L
EXAMPLE_SEEDS <- c(data  = EXAMPLE_SEED,
                   draws = EXAMPLE_SEED + 1000003L,
                   power = EXAMPLE_SEED + 2000003L)


# Variables -------------------------------------------------------------------

example_labels <- c(
  A  = "Family Foundations (1) vs control (0), randomized",
  M1 = "Low coparenting competition (fathers, child age 1; competition reverse-scored)",
  M2 = "Coparenting positivity (fathers, child age 1)",
  Y  = "Child adjustment (age 3; adjustment problems reverse-scored)",
  X1 = "Maternal education (pretest)",
  X2 = "Social desirability (pretest)"
)


# Inputs --------

example_inputs <- list(
  n                = 169,           
  num_x            = 2,
  treat.prop       = 0.5,
  treat.randomized = TRUE,
  R2.ax            = 0,
  R2.mx            = c(0.02, 0.02),
  R2.yx            = 0.02,
  std.m_on_a       = c(0.27, 0.36),
  std.y_on_a       = -0.03,
  std.y_on_m       = c(0.25, 0.15),
  std.y_on_am_2way = c(0.10, 0.10),
  std.y_on_m_2way  = 0.10,
  std.y_on_am_3way = 0.10,
  em_corr          = 0.20,
  M_binary         = c(FALSE, FALSE),
  Y_binary         = FALSE
)


# Decomposition pair ------------------------------

# TE = IDE(,0) + IIE_M1(1,,0) + IIE_M2(1,1,) + IIE_mu, the decomposition
# written in the manuscript (the IIE_mu definition after Eq. 2).
EXAMPLE_PAIR <- c("IIE_M1(1,,0)", "IIE_M2(1,1,)")

# sig.level per mediator family; `joint` applies when probing which mediator operates
EXAMPLE_ALPHA <- c(per_mediator = 0.05, joint = 0.025)

# Estimation on the one dataset 
example_estimation <- list(
  n.draws   = 1000,
  center    = TRUE,                 # mean-centre covariates (center_covariates())
  decisions = "p_MC",
  intervals = c(pointwise_95  = 0.95,   pointwise_975    = 0.975,
                bonferroni_9875 = 0.9875, bonferroni_99375 = 0.99375)
)

# Power and required n 
example_power <- list(
  nsims         = 1000,
  n.draws       = 1000,
  center        = TRUE,
  power_target  = 0.80,
  n_range       = c(n.min = 50, n.max = 1000, n.step = 10),
  sig.level     = unname(EXAMPLE_ALPHA),
  methods       = c("no_adjust", "bonferroni", "holm"),
  tools_methods = c("no_adjust", "bonferroni")
)

# Scenarios of R2_interaction_grid.R
EXAMPLE_GRID_VALUES <- c(-0.10, 0, 0.05, 0.10, 0.14)

example_scenarios <- list(
  maintext        = list(),
  published_a = list(std.m_on_a = c(0.17, 0.05)),
  baseline_R2 = list(R2.mx = c(0.30, 0.30), R2.yx = 0.30)
)


# gen.data()/runPower() arguments for a scenario
example_args <- function(scenario = "maintext") {
  stopifnot(scenario %in% names(example_scenarios))
  utils::modifyList(example_inputs, example_scenarios[[scenario]])
}
