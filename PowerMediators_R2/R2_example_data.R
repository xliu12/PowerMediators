# The empirical example data
#
# The one simulated dataset of the main-text example, and the
# estimates, intervals and decisions for its eight interventional indirect
# effects (IIEs). The design is in R2_example_spec.R; the
# estimator and the tests come from add_on.R.

source("add_on.R")
source("R2_example_spec.R")


# Settings 

ALPHAS  <- unname(EXAMPLE_ALPHA)        # .05 per mediator; .025 per mediator when testing which mediators operate
N_DRAWS <- example_estimation$n.draws   


# Data ------------------------------------------------------------------------

# The example dataset
example_dataset <- function() {
  args <- example_args("maintext")
  g <- do.call(gen_data_r2, c(list(seed = EXAMPLE_SEEDS[["data"]]), args))
  g$args <- args
  g
}


# Estimation ------------------------------------------------------------------

example_estimates <- function(g) {
  map(ALPHAS, \(a) {
    med2_reg_plus(g$data, g$M_binary, g$Y_binary, sig.level = a, n.draws = N_DRAWS,
                  seed = EXAMPLE_SEEDS[["draws"]], center = example_estimation$center) %>%
      mutate(sig.level = a, .before = 1)
  }) %>%
    list_rbind() %>%
    left_join(g$true_vals, by = "effect") %>%
    relocate(true_val, .after = effect)
}


# Tables ----------------------------------------------------------------------

# Table A: one row per IIE. Intervals: pointwise 95% and 97.5% (no adjustment
# at .05 and .025) and Bonferroni 98.75% and 99.375% (at .05 and .025).
# Decisions are p_MC-based. `bonf_ci_excl0_*` is the interval-based Bonferroni
# decision of med2.reg(); it can differ from the p_MC-based one only at the
# one-draw boundary.
table_a <- function(res) {
  maintext <- res %>%
    filter(sig.level == ALPHAS[1], method == "no_adjust") %>%
    transmute(mediation, effect, in_pair = effect %in% EXAMPLE_PAIR, true_val, est, SE, p_mc, rbar)
  v <- function(a, m, col) {
    d <- res %>% filter(sig.level == a, method == m)
    set_names(d[[col]], d$effect)[maintext$effect]
  }
  a1 <- ALPHAS[1]
  a2 <- ALPHAS[2]
  maintext %>%
    mutate(
      holm_p                    = v(a1, "holm", "p_adj"),
      ci95_low                  = v(a1, "no_adjust", "ci_low"),
      ci95_up                   = v(a1, "no_adjust", "ci_up"),
      ci975_low                 = v(a2, "no_adjust", "ci_low"),
      ci975_up                  = v(a2, "no_adjust", "ci_up"),
      bonf9875_low              = v(a1, "bonferroni", "ci_low"),
      bonf9875_up               = v(a1, "bonferroni", "ci_up"),
      bonf99375_low             = v(a2, "bonferroni", "ci_low"),
      bonf99375_up              = v(a2, "bonferroni", "ci_up"),
      reject_none_05            = v(a1, "no_adjust", "sig_p"),
      reject_none_025           = v(a2, "no_adjust", "sig_p"),
      reject_bonf_p_05          = v(a1, "bonferroni", "sig_p"),
      reject_bonf_p_025         = v(a2, "bonferroni", "sig_p"),
      bonf_ci_excl0_05          = v(a1, "bonferroni", "if_sig"),
      bonf_ci_excl0_025         = v(a2, "bonferroni", "if_sig"),
      reject_holm_05            = v(a1, "holm", "sig_p"),
      reject_holm_025           = v(a2, "holm", "sig_p")
    ) %>%
    relocate(rbar, .after = last_col())
}

# Existence of mediation per mediator: the global null is rejected if any of
# the four IIEs is rejected under Bonferroni (p_MC < sig.level / 4), which is
# Holm's decision on the same p-values.
existence_table <- function(res) {
  res %>%
    filter(method %in% c("bonferroni", "holm")) %>%
    group_by(sig.level, mediation) %>%
    summarise(
      alpha_per_test_bonferroni = unique(alpha_PT[method == "bonferroni"]),
      min_p_mc                  = min(p_mc),
      reject_global_bonferroni  = any(sig_p[method == "bonferroni"]),
      reject_global_holm        = any(sig_p[method == "holm"]),
      n_rejected_bonferroni     = sum(sig_p[method == "bonferroni"]),
      n_rejected_holm           = sum(sig_p[method == "holm"]),
      iies_rejected_holm        = paste(effect[method == "holm" & sig_p], collapse = "; "),
      .groups = "drop"
    )
}

# The decomposition pair, one IIE per mediator: pointwise intervals at
# sig.level (95%) and at sig.level / 2 (97.5%, the level for testing which mediators operate).
pair_table <- function(res) {
  res %>%
    filter(effect %in% EXAMPLE_PAIR, method == "no_adjust") %>%
    group_by(sig.level) %>%
    mutate(interval = paste0(formatC(100 * (1 - alpha_PT), format = "fg"), "%"),
           both_rejected = all(sig_p)) %>%
    ungroup() %>%
    select(sig.level, interval, mediation, effect, true_val, est, SE, ci_low, ci_up,
           p_mc, reject = sig_p, both_rejected)
}



# Run -------------------------------------------------------------------------

g         <- example_dataset()
res       <- example_estimates(g)
tab       <- table_a(res)
existence <- existence_table(res)
pair      <- pair_table(res)

print(tab, n = Inf, width = Inf)         # the eight IIEs 
print(existence, n = Inf, width = Inf)   # the global test of each mediator
print(pair, n = Inf, width = Inf)        # the pre-specified pair
