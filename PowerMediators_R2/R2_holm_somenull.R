# The five multiple testing methods in the two conditions of the main-text Type I error table in which some IIEs via M1 are non-zero (its three-null and two-null blocks): three null IIEs (y_on_am1m2 = .36) and two null IIEs (y_on_m1m2 = .36), under the mediator-outcome null. All five methods are evaluated on the same replications: the four adjustment methods by their intervals, as in the package, and Holm's procedure by the Monte Carlo p-values. The per-test rates under no adjustment are those of the supplement's per-test table for these two conditions.
source("add_on.R")


# Settings 
NSIMS       <- 1000
N_DRAWS     <- 1000
MC_CORES    <- 1
MASTER_SEED <- 20261017
SIG_LEVEL   <- 0.05
N_VALUES    <- c(100, 300, 500)
CENTER_X    <- FALSE  # X has mean 0 by design


# Design ----------------------------------------------------------------------

some_null_configs <- tribble(
  ~config, ~null_type,                               ~m1_on_a, ~y_on_m1, ~y_on_am1, ~y_on_m1m2, ~y_on_am1m2,
  "N3",    "mediator-outcome null, three null IIEs", 0.36,     0,        0,         0,          0.36,
  "N2",    "mediator-outcome null, two null IIEs",   0.36,     0,        0,         0.36,       0
)

fixed_design <- tibble(   # as in R2_simulation.R
  m2_on_a = 0.36, y_on_a = 0.14, y_on_m2 = 0, y_on_am2 = 0.36, em_corr = 0.1,
  M1_binary = FALSE, M2_binary = FALSE, Y_binary = FALSE
) %>%
  mutate(R2 = if_else(Y_binary, 0.1, 0.3))  # as in simulation.R

condition <- some_null_configs %>%
  expand_grid(fixed_design, n = N_VALUES, n.draws = N_DRAWS) %>%
  mutate(n_reps = NSIMS) %>%
  mutate(cond_id = row_number(), .before = 1)

set.seed(MASTER_SEED)
rep_seeds <- expand_grid(n = N_VALUES, rep = seq_len(NSIMS)) %>%
  mutate(iseed = sample.int(1e8, n()))


# One replication and one cell ---------------------------------------------------

# Simulate one dataset and test the eight IIEs with every method, as in R2_simulation.R
one_rep <- function(iseed, cond) {
  g <- gen_data_r2(
    seed = iseed,
    n = cond$n,
    R2.ax = 0, R2.mx = cond$R2, R2.yx = cond$R2,
    std.m_on_a = c(cond$m1_on_a, cond$m2_on_a),
    std.y_on_a = cond$y_on_a,
    std.y_on_m = c(cond$y_on_m1, cond$y_on_m2),
    std.y_on_am_2way = c(cond$y_on_am1, cond$y_on_am2),
    std.y_on_m_2way = cond$y_on_m1m2,
    std.y_on_am_3way = cond$y_on_am1m2,
    em_corr = cond$em_corr,
    M_binary = c(cond$M1_binary, cond$M2_binary),
    Y_binary = cond$Y_binary
  )
  med2_reg_plus(g$data, g$M_binary, g$Y_binary,
                sig.level = SIG_LEVEL, n.draws = cond$n.draws,
                seed = iseed + 1000003L, center = CENTER_X) %>%
    left_join(g$true_vals, by = "effect")
}

# All replications of one condition
run_cell <- function(cond) {
  seeds <- rep_seeds %>% filter(n == cond$n)
  started <- Sys.time()
  reps <- parallel::mclapply(seq_len(nrow(seeds)), \(r) {
    tryCatch(
      one_rep(seeds$iseed[r], cond) %>%
        mutate(rep = seeds$rep[r], iseed = seeds$iseed[r], .before = 1),
      error = \(e) tibble(rep = seeds$rep[r], iseed = seeds$iseed[r],
                          error = conditionMessage(e))
    )
  }, mc.cores = MC_CORES, mc.preschedule = TRUE) %>%
    list_rbind()
  if (!"error" %in% names(reps)) reps <- reps %>% mutate(error = NA_character_)
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  message(glue("cell {cond$cond_id}/{nrow(condition)}: {cond$config}, n = {cond$n}, ",
               "{cond$n.draws} draws, {nrow(seeds)} replications in {round(seconds, 1)} s"))
  list(reps = reps %>% mutate(cond_id = cond$cond_id, .before = 1), seconds = seconds)
}


# Summaries ---------------------------------------------------------------------

# Family-level metrics per mediator and method. Holm has no intervals.
summarise_family <- function(reps_ok) {
  reps_ok %>%
    group_by(cond_id, config, null_type, n, n.draws, mediation, method, rep) %>%
    summarise(
      n_null         = sum(null),
      reject_null    = any(if_sig & null),
      reject_nonnull = any(if_sig & !null),
      all_covered    = all(covered),
      width          = mean(width),
      rbar           = mean(rbar),
      .groups = "drop"
    ) %>%
    group_by(cond_id, config, null_type, n, n.draws, mediation, method) %>%
    summarise(
      nreps    = n(),
      n_null   = max(n_null),
      alpha_FW = if (n_null > 0) mean(reject_null) else NA_real_,
      power_FW = if (n_null < 4) mean(reject_nonnull) else NA_real_,
      CR_FW    = mean(all_covered),
      width    = mean(width),
      rbar     = mean(rbar),
      .groups = "drop"
    )
}

# Per-IIE metrics
summarise_pertest <- function(reps_ok) {
  reps_ok %>%
    group_by(cond_id, config, null_type, n, n.draws, mediation, effect, method) %>%
    summarise(
      nreps    = n(),
      true_val = mean(true_val),
      null     = all(null),
      reject   = mean(if_sig),
      CR       = mean(covered),
      width    = mean(width),
      alpha_PT = mean(alpha_PT),
      rbar     = mean(rbar),
      .groups = "drop"
    ) %>%
    mutate(TypeI_PT = if_else(null, reject, NA_real_),
           power_PT = if_else(null, NA_real_, reject)) %>%
    select(-reject)
}


# Run -----------------------------------------------------------------------------

runs <- map(seq_len(nrow(condition)), \(i) run_cell(condition[i, ]))

reps <- runs %>% map("reps") %>% list_rbind()

reps_ok <- reps %>%
  filter(is.na(error)) %>%
  select(-error) %>%
  left_join(condition %>% select(cond_id, config, null_type, n, n.draws), by = "cond_id") %>%
  mutate(null = true_val == 0,
         covered = ci_low < true_val & ci_up > true_val) %>%
  select(cond_id, config, null_type, n, n.draws, rep, iseed, mediation, effect,
         method, alpha_PT, est, SE, ci_low, ci_up, width, if_sig, sig_p, rbar,
         p_mc, p_adj, true_val, null, covered)

family  <- summarise_family(reps_ok)
pertest <- summarise_pertest(reps_ok)

# The three-null and two-null blocks of Table 4 are the IIE_M1 rows: alpha_FW, Power_FW, CR and width of the four adjustment methods, and alpha_FW and Power_FW of method "holm" (CR and Width do not apply: Holm's procedure gives no intervals). The Bonferroni adjustment decides by its intervals and Holm's procedure by the p-values, which can differ at the boundary of an interval.
print(family %>%
        select(config, n, mediation, method, n_null, alpha_FW, power_FW, CR_FW, width) %>%
        arrange(mediation, config, n, method),
      n = Inf, width = Inf)

# The per-test rates: the Type I error rate of each null IIE (the rows of the supplement's per-test table under no adjustment), the power of each non-null IIE, coverage and width
print(pertest %>%
        select(config, n, mediation, effect, method, null, TypeI_PT, power_PT, CR, width) %>%
        arrange(mediation, config, n, effect, method),
      n = Inf, width = Inf)
