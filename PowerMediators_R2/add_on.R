# Additional functions around the PowerMediators package version 0.1.0, added during the Round-2 revision of the paper; PowerMediators version 0.2.0 (the current version on github) adds Holm's procedure to runPower() and cal.sample_size(); either version of the package can run these additional functions.

# Packages ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mvtnorm)
  library(glue)
  library(tidyverse)
  library(PowerMediators)
})


# Data ----------------------------------------------------------------------

# Mean-centre every covariate (columns whose names start with "X"). Needed before estimating IIEs for continuous mediators and outcome, because cal.IIE_M1()/cal.IIE_M2() evaluate the mediator means at X = 0.
center_covariates <- function(data) {
  data %>% mutate(across(matches("^X"), \(x) x - mean(x)))
}

# Simulate one dataset with PowerMediators:::gen.data(). 
gen_data_r2 <- function(seed, ...) {
  set.seed(seed)
  g <- PowerMediators:::gen.data(...)
  list(
    data      = as_tibble(g$dat),
    true_vals = as_tibble(g$true_vals) %>% rename(true_val = true_vals),
    M_binary  = g$M_binary,
    Y_binary  = g$Y_binary,
    coefs     = list(m1_coefs = g$m1_coefs, m2_coefs = g$m2_coefs, y_coefs = g$y_coefs,
                     var_em = g$var_em, var_ey = g$var_ey,
                     std.m_on_x = g$std.m_on_x, std.y_on_x = g$std.y_on_x)
  )
}


# Monte Carlo draws, p-values and tests --------------------------------------

# n.draws x 8 tibble of Monte Carlo draws of the IIEs, from PowerMediators:::estimate_med2.reg(). 
mc_draws <- function(data, M_binary = c(FALSE, FALSE), Y_binary = FALSE,
                     n.draws = 5000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  PowerMediators:::estimate_med2.reg(as.data.frame(data), M_binary, Y_binary,
                                     n.draws = n.draws) %>%
    as_tibble()
}

# Two-sided Monte Carlo p-value 
p_mc <- function(draws) {
  draws <- as.matrix(draws)
  pmin(1, 2 * pmin(colMeans(draws <= 0, na.rm = TRUE),
                   colMeans(draws >= 0, na.rm = TRUE)))
}

# Bonferroni and Holm decisions from a named vector of p-values (one family).
holm_decisions <- function(p, sig.level = 0.05) {
  tibble(
    effect      = names(p) %||% seq_along(p),
    p           = unname(p),
    p_bonf      = p.adjust(p, "bonferroni"),
    p_holm      = p.adjust(p, "holm"),
    reject_bonf = p_bonf < sig.level,
    reject_holm = p_holm < sig.level
  )
}

# Per-test levels of the four non-step-down methods for one test in a family of k,
# using the same expressions as med2.reg().
alpha_pt <- function(sig.level, k, rbar) {
  c(no_adjust     = sig.level,
    bonferroni    = sig.level / k,
    modified_bon1 = sig.level / (k - (k - 1) * sqrt(rbar)),
    modified_bon2 = sig.level / (k^(1 - sqrt(rbar))))
}

# Tests for one family (the four IIEs of one mediator). Holm has 
# no simple simultaneous interval, so its CI columns are NA and its
# decision uses the Holm-adjusted p_MC. 
mc_tests_family <- function(draws_m, sig.level = 0.05) {
  draws_m <- as.matrix(draws_m)
  k <- ncol(draws_m)
  r <- cor(draws_m, use = "pairwise")
  rbar <- map_dbl(seq_len(k), \(i) mean(abs(r[i, -i])))
  p <- p_mc(draws_m)
  p_holm <- p.adjust(p, "holm")

  map(seq_len(k), \(i) {
    x <- draws_m[, i]
    a <- alpha_pt(sig.level, k, rbar[i])
    tibble(
      effect   = colnames(draws_m)[i],
      method   = c(names(a), "holm"),
      alpha_PT = c(unname(a), NA),
      SE       = sd(x, na.rm = TRUE),
      ci_low   = c(unname(quantile(x, probs = a / 2, na.rm = TRUE)), NA),
      ci_up    = c(unname(quantile(x, probs = 1 - a / 2, na.rm = TRUE)), NA),
      rbar     = rbar[i],
      p_mc     = p[i],
      p_adj    = c(p[i], min(1, p[i] * k), NA, NA, p_holm[i])
    )
  }) %>%
    list_rbind() %>%
    mutate(
      width  = ci_up - ci_low,
      if_sig = if_else(method == "holm", p_adj < sig.level, ci_low * ci_up > 0),
      sig_p  = p_adj < sig.level
    )
}

# Tests for both mediator families from a draws matrix (columns IIE_M1(...)
# and IIE_M2(...)). `est` is an optional named vector of point estimates.
mc_tests <- function(draws, est = NULL, sig.level = 0.05) {
  draws <- as.matrix(draws)
  out <- map(c("IIE_M1", "IIE_M2"), \(m) {
    cols <- str_starts(colnames(draws), fixed(m))
    mc_tests_family(draws[, cols, drop = FALSE], sig.level) %>%
      mutate(mediation = m, .before = 1)
  }) %>%
    list_rbind()
  if (!is.null(est)) {
    out <- out %>%
      left_join(tibble(effect = names(est), est = unname(est)), by = "effect") %>%
      relocate(est, .after = alpha_PT)
  }
  out
}

# Estimate and test the eight IIEs on one dataset: point estimates, Monte
# Carlo draws, the four non-step-down methods, Holm, p_MC and r-bar.

med2_reg_plus <- function(data, M_binary = c(FALSE, FALSE), Y_binary = FALSE,
                          sig.level = 0.05, n.draws = 5000, seed = NULL,
                          center = TRUE) {
  stopifnot(all(c("A", "M1", "M2", "Y") %in% names(data)))
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  if (center) data <- center_covariates(data)
  data <- as.data.frame(data)

  set.seed(seed)
  est <- unlist(PowerMediators:::estimate_med2.reg(data, M_binary, Y_binary, n.draws = NULL))
  draws <- PowerMediators:::estimate_med2.reg(data, M_binary, Y_binary, n.draws = n.draws)
  res <- mc_tests(draws, est, sig.level)

  res
}


# Interaction grids -------------------------------
# Grid of the four interaction coefficients (standardized). "joint" sets all
# four to each value; "one_at_a_time" varies one term while the others stay at `maintext`. 
interaction_grid <- function(values, type = c("joint", "one_at_a_time"), maintext = 0.10) {
  type <- match.arg(type)
  terms <- c("am1", "am2", "m1m2", "am1m2")
  grid <- if (type == "joint") {
    tibble(term = "all", value = values)
  } else {
    expand_grid(term = terms, value = values)
  }
  grid %>%
    mutate(
      coef = map2(term, value, \(t, v) set_names(if_else(terms == t | t == "all", v, maintext), terms)),
      std.y_on_am_2way = map(coef, \(b) unname(b[c("am1", "am2")])),
      std.y_on_m_2way  = map_dbl(coef, "m1m2"),
      std.y_on_am_3way = map_dbl(coef, "am1m2")
    ) %>%
    select(-coef) %>%
    mutate(grid_id = row_number(), .before = 1)
}

# Split a grid into the columns that match `fun`'s arguments (passed row by
# row) and the fixed arguments in `...` that the grid does not override.
grid_call <- function(grid, fun, fixed) {
  args <- grid %>% select(any_of(names(formals(fun))))
  fixed <- fixed[setdiff(names(fixed), names(args))]
  pmap(args, \(...) exec(fun, !!!fixed, ...))
}

# Power over a grid: calls the exported runPower() once per grid row
grid_power <- function(grid, ..., seed_user = 12) {
  res <- suppressMessages(
    grid_call(grid, PowerMediators::runPower, list(..., seed_user = seed_user))
  )
  grid %>%
    mutate(power = map(res, \(r) r %>% ungroup() %>% as_tibble() %>%
                         select(-any_of(names(grid))))) %>%
    unnest(power)
}

# Required sample size over a grid: calls the exported cal.sample_size() once per grid row 
grid_sample_size <- function(grid, ..., seed_user = 12) {
  targets <- c("power_FW_M1M2", "power_FW_eachM", "power_PT_someIIE", "power_PT_allIIE")
  res <- suppressMessages(
    grid_call(grid, PowerMediators::cal.sample_size, list(..., seed_user = seed_user))
  )
  grid %>%
    mutate(
      required_n = map(res, \(s) {
        compact(s[targets]) %>%
          imap(\(d, name) d %>% ungroup() %>% as_tibble() %>% mutate(target = name, .before = 1)) %>%
          list_rbind()
      }),
      res_n = map(res, \(s) as_tibble(ungroup(s$res_n)))
    )
}



