# Power at n = 169 and the required sample size for the main-text example
# over the interaction grid

source("add_on.R")
source("R2_example_spec.R")


# Settings --------------------------------------------------------------------

NSIMS     <- 1000
N_DRAWS   <- example_power$n.draws    
MC_CORES  <- max(1, parallel::detectCores() - 1)
if (.Platform$OS.type == "windows") MC_CORES <- 1  # mclapply() forks
RUN_REQN  <- TRUE                     
CACHE_DIR <- "R2_grid_cache"

MASTER_SEED  <- EXAMPLE_SEEDS[["power"]]
NSIMS_MAX    <- example_power$nsims          
ALPHAS       <- example_power$sig.level      # .05, .025
POWER_TARGET <- example_power$power_target   # .80
N_RANGE      <- example_power$n_range        # n.min 50, n.max 1000, n.step 10
N_EXAMPLE    <- example_inputs$n             # 169
CENTER_X     <- example_power$center         # TRUE
N_FINE       <- seq(N_RANGE[["n.min"]], N_RANGE[["n.max"]], by = N_RANGE[["n.step"]])
N_COARSE     <- seq(N_RANGE[["n.min"]], N_RANGE[["n.max"]], by = 50)
N_ALL        <- sort(unique(c(N_EXAMPLE, N_FINE)))
METHODS      <- c(no_adjust = "none", bonferroni = "bonf", holm = "holm")
stopifnot(NSIMS <= NSIMS_MAX, all(N_COARSE %in% N_FINE))

# The series whose required n is searched (each at every level in ALPHAS).
TRACKED <- tribble(
  ~target,                        ~column,       ~method,
  "T1_FW_both_mediators",         "FW_min",      "bonf",
  "T2_FW_M1",                     "FW_M1",       "bonf",
  "T3_min_individual_power_8",    "PT_min_all8", "none",
  "T3_min_individual_power_8",    "PT_min_all8", "holm",
  "T3_min_individual_power_8",    "PT_min_all8", "bonf",
  "T4_pair",                      "pair_min",    "none"
) %>%
  expand_grid(alpha = ALPHAS)


# Scenarios -------------------------------------------------------------------

# The labelled rows (joint grid, one-at-a-time grid, named scenarios) 
build_scenarios <- function() {
  maintext_args <- example_args("maintext")
  grid <- bind_rows(
    interaction_grid(EXAMPLE_GRID_VALUES, "joint", maintext = 0.10) %>% mutate(grid_type = "joint"),
    interaction_grid(EXAMPLE_GRID_VALUES, "one_at_a_time", maintext = 0.10) %>% mutate(grid_type = "one_at_a_time")
  ) %>%
    mutate(
      args = pmap(list(std.y_on_am_2way, std.y_on_m_2way, std.y_on_am_3way), \(am, mm, amm) {
        utils::modifyList(maintext_args, list(std.y_on_am_2way = am, std.y_on_m_2way = mm, std.y_on_am_3way = amm))
      }),
      label = if_else(grid_type == "joint", paste0("joint c = ", value), paste0(term, " = ", value)),
      candidate = if_else(grid_type == "joint", paste0("joint_", value), paste0(term, "_", value))
    ) %>%
    select(grid_type, term, value, label, candidate, args)
  named <- tibble(
    grid_type = "scenario", term = names(example_scenarios), value = NA_real_,
    label = names(example_scenarios), candidate = names(example_scenarios),
    args = map(names(example_scenarios), example_args)
  )
  labelled <- bind_rows(named, grid) %>%
    mutate(
      row_id = row_number(),
      key = map_chr(args, \(a) paste(names(a), map_chr(a, \(v) paste(format(v, digits = 15), collapse = ",")),
                                     sep = "=", collapse = ";")),
      priority = case_when(grid_type == "scenario" ~ 1L, grid_type == "joint" ~ 2L, TRUE ~ 3L)
    )
  labelled %>%
    arrange(key, priority, row_id) %>%
    group_by(key) %>%
    summarise(
      scenario_id   = first(candidate),
      label         = first(label),
      grid_type     = first(grid_type),
      term          = first(term),
      value         = first(value),
      members       = paste(candidate, collapse = " | "),
      n_members     = n(),
      in_tableB     = any(grid_type %in% c("scenario", "joint")),
      in_supplement = any(grid_type == "one_at_a_time"),
      args          = list(first(args)),
      .groups = "drop"
    ) %>%
    mutate(
      am1   = map_dbl(args, \(a) a$std.y_on_am_2way[1]),
      am2   = map_dbl(args, \(a) a$std.y_on_am_2way[2]),
      m1m2  = map_dbl(args, "std.y_on_m_2way"),
      am1m2 = map_dbl(args, "std.y_on_am_3way"),
      a1    = map_dbl(args, \(a) a$std.m_on_a[1]),
      a2    = map_dbl(args, \(a) a$std.m_on_a[2]),
      R2    = map_dbl(args, "R2.yx"),
      order_group = case_when(scenario_id == "maintext" ~ 0, grid_type == "joint" ~ 1,
                              scenario_id == "published_a" ~ 2, scenario_id == "baseline_R2" ~ 3,
                              TRUE ~ 4 + match(term, c("am1", "am2", "m1m2", "am1m2")))
    ) %>%
    arrange(order_group, value) %>%
    mutate(order = row_number(), .before = 1) %>%
    select(-order_group)
}

set.seed(MASTER_SEED)
rep_seeds <- expand_grid(n = N_ALL, rep = seq_len(NSIMS_MAX)) %>%
  mutate(iseed = sample.int(1e8, n()))


# One replication and one cell ---------------------------------------------------

# Simulate one dataset and test the eight IIEs
one_rep <- function(iseed, n, args) {
  g <- do.call(gen_data_r2, c(list(seed = iseed), utils::modifyList(args, list(n = n))))
  suppressMessages(
    med2_reg_plus(g$data, g$M_binary, g$Y_binary, sig.level = ALPHAS[1], n.draws = N_DRAWS,
                  seed = iseed + 1000003L, center = CENTER_X)
  ) %>%
    filter(method %in% names(METHODS)) %>%
    mutate(method = METHODS[method]) %>%
    select(mediation, effect, method, est, p_mc, p_adj) %>%
    pivot_wider(names_from = method, values_from = p_adj, names_prefix = "p_") %>%
    left_join(g$true_vals, by = "effect") %>%
    relocate(true_val, .after = effect)
}

# All replications of one cell 
run_cell <- function(scenario_id, args, n) {
  seeds <- rep_seeds %>% filter(n == !!n, rep <= NSIMS)
  started <- Sys.time()
  reps <- parallel::mclapply(seq_len(nrow(seeds)), \(r) {
    tryCatch(
      one_rep(seeds$iseed[r], n, args) %>%
        mutate(rep = seeds$rep[r], iseed = seeds$iseed[r], .before = 1),
      error = \(e) tibble(rep = seeds$rep[r], iseed = seeds$iseed[r], error = conditionMessage(e))
    )
  }, mc.cores = MC_CORES, mc.preschedule = TRUE)
  reps <- imap(reps, \(x, r) {
    if (is.data.frame(x)) x else tibble(rep = seeds$rep[r], iseed = seeds$iseed[r],
                                        error = paste("worker failed:", paste(as.character(x), collapse = " ")))
  }) %>%
    list_rbind()
  if (!"error" %in% names(reps)) reps <- reps %>% mutate(error = NA_character_)
  seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  ok     <- reps %>% filter(is.na(error)) %>% select(-error)
  failed <- reps %>% filter(!is.na(error)) %>% distinct(rep, iseed, error)
  summ   <- summarise_cell(ok)
  message(glue("  {scenario_id}, n = {n}: {n_distinct(ok$rep)} replications in {round(seconds, 1)} s",
               "{if (nrow(failed) > 0) glue(', {nrow(failed)} failed') else ''}"))
  list(scenario_id = scenario_id, n = n, seconds = seconds,
       nreps = n_distinct(ok$rep), n_failed = nrow(failed), failed = failed,
       by_iie = summ$by_iie, by_target = summ$by_target)
}

# Summaries of one cell from its per-replication rows
summarise_cell <- function(ok) {
  dec <- ok %>%
    select(rep, mediation, effect, true_val, est, p_none, p_bonf, p_holm) %>%
    pivot_longer(c(p_none, p_bonf, p_holm), names_to = "method", names_prefix = "p_", values_to = "p") %>%
    expand_grid(alpha = ALPHAS) %>%
    mutate(reject = p < alpha)
  by_iie <- dec %>%
    group_by(alpha, method, mediation, effect) %>%
    summarise(nreps = n(), true_val_mean = mean(true_val), true_val_sd = sd(true_val),
              true_val_min = min(true_val), true_val_max = max(true_val),
              est_mean = mean(est), power = mean(reject), .groups = "drop") 
  per_rep <- dec %>%
    group_by(alpha, method, rep) %>%
    summarise(any_M1 = any(reject[mediation == "IIE_M1"]),
              any_M2 = any(reject[mediation == "IIE_M2"]),
              all_8  = all(reject),
              pair_1 = reject[effect == EXAMPLE_PAIR[1]],
              pair_2 = reject[effect == EXAMPLE_PAIR[2]],
              .groups = "drop")
  min_iie <- by_iie %>%
    group_by(alpha, method) %>%
    summarise(PT_min_all8 = min(power), PT_min_effect = effect[which.min(power)], .groups = "drop")
  target_cols <- c("FW_M1", "FW_M2", "FW_min", "FW_joint", "PT_min_all8", "complete_all8",
                   "pair_M1", "pair_M2", "pair_min", "pair_both")
  by_target <- per_rep %>%
    group_by(alpha, method) %>%
    summarise(nreps = n(), FW_M1 = mean(any_M1), FW_M2 = mean(any_M2), FW_joint = mean(any_M1 & any_M2),
              complete_all8 = mean(all_8), pair_M1 = mean(pair_1), pair_M2 = mean(pair_2),
              pair_both = mean(pair_1 & pair_2), .groups = "drop") %>%
    mutate(FW_min = pmin(FW_M1, FW_M2), pair_min = pmin(pair_M1, pair_M2)) %>%
    left_join(min_iie, by = c("alpha", "method")) %>%
    select(alpha, method, nreps, all_of(target_cols), PT_min_effect)
  list(by_iie = by_iie, by_target = by_target)
}


# Cache 
cache_path <- function(scenario_id) file.path(CACHE_DIR, paste0(scenario_id, ".rds"))
design_key <- function() list(nsims = as.numeric(NSIMS), n_draws = as.numeric(N_DRAWS),
                              master_seed = as.numeric(MASTER_SEED),
                              center = CENTER_X, alphas = as.numeric(ALPHAS))

load_cache <- function(scenario_id) {
  p <- cache_path(scenario_id)
  if (!file.exists(p)) return(list(scenario_id = scenario_id, key = design_key(), cells = list()))
  cache <- readRDS(p)
  if (!identical(cache$key, design_key())) {
    stop(glue("{p} was written under other settings (nsims, n.draws, seed or centring); ",
              "delete it or use another CACHE_DIR."), call. = FALSE)
  }
  cache
}

save_cache <- function(cache) {
  p <- cache_path(cache$scenario_id)
  dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(p, ".tmp")
  saveRDS(cache, tmp)
  file.rename(tmp, p)
  invisible(cache)
}

# The cell (scenario x n) from the cache, or run and cached.
get_cell <- function(cache, scenario, n) {
  key <- as.character(n)
  if (is.null(cache$cells[[key]])) {
    cache$cells[[key]] <- run_cell(scenario$scenario_id, scenario$args[[1]], n)
    save_cache(cache)
  }
  cache
}

cell_targets <- function(cache) {
  map(cache$cells, \(cell) cell$by_target %>% mutate(scenario_id = cell$scenario_id, n = cell$n, .before = 1)) %>%
    list_rbind()
}


# Required-n search -------------------------------------------------------------

# Long table of the tracked series over the evaluated n (n = 169 excluded).
tracked_long <- function(cache) {
  bt <- cell_targets(cache)
  if (nrow(bt) == 0) return(TRACKED %>% mutate(n = NA_real_, power = NA_real_))
  cols <- unique(TRACKED$column)
  long <- bt %>%
    filter(n != N_EXAMPLE) %>%
    select(n, alpha, method, all_of(cols)) %>%
    pivot_longer(all_of(cols), names_to = "column", values_to = "power")
  TRACKED %>%
    left_join(long %>% select(n, alpha, method, column, power), by = c("column", "method", "alpha"),
              relationship = "many-to-many")
}

all_tracked_reached <- function(cache) {   
  tracked_long(cache) %>%
    filter(is.na(n) | n %in% N_COARSE) %>%
    group_by(target, column, method, alpha) %>%
    summarise(reached = any(power >= POWER_TARGET, na.rm = TRUE), .groups = "drop") %>%
    pull(reached) %>%
    all()
}

# First crossing on the coarse grid per tracked series (NA when never reached).
first_crossing <- function(cache) {
  tracked_long(cache) %>%
    filter(n %in% N_COARSE, power >= POWER_TARGET) %>%
    group_by(target, column, method, alpha) %>%
    summarise(n0 = min(n), .groups = "drop")
}

search_n <- function(scenario, cache) {
  for (n in N_COARSE) {
    if (all_tracked_reached(cache)) break
    cache <- get_cell(cache, scenario, n)
  }
  n0 <- unique(first_crossing(cache)$n0)
  refine <- sort(unique(unlist(map(n0, \(x) seq(x - 40, x - 10, by = 10)))))
  refine <- refine[refine >= N_RANGE[["n.min"]]]
  for (n in refine) cache <- get_cell(cache, scenario, n)
  cache
}

# Required n per tracked series
required_n_table <- function(cache) {
  tracked_long(cache) %>%
    group_by(target, column, method, alpha) %>%
    group_modify(\(d, k) {
      d <- d %>% filter(!is.na(n)) %>% arrange(n)
      none <- tibble(n_required = NA_real_, n_label = NA_character_, power_at_n = NA_real_,
                     n_last_fail = NA_real_, n_fail_above = NA_real_,
                     n_evaluated = nrow(d), n_max_evaluated = if (nrow(d)) max(d$n) else NA_real_,
                     status = "not searched")
      if (nrow(d) == 0) return(none)
      coarse <- d %>% filter(n %in% N_COARSE, power >= POWER_TARGET)
      if (nrow(coarse) == 0) {
        return(none %>% mutate(n_label = ">1000", status = "not reached by n = 1000",
                               n_last_fail = max(d$n[d$power < POWER_TARGET])))
      }
      n0 <- min(coarse$n)
      window <- d %>% filter(n > n0 - 50, n <= n0, power >= POWER_TARGET)
      nr <- min(window$n)
      below <- d$n[d$n < nr & d$power < POWER_TARGET]
      above <- d$n[d$n > nr & d$power < POWER_TARGET]
      tibble(n_required = nr, n_label = as.character(nr),
             power_at_n = d$power[d$n == nr],
             n_last_fail = if (length(below)) max(below) else NA_real_,
             n_fail_above = if (length(above)) min(above) else NA_real_,
             n_evaluated = nrow(d), n_max_evaluated = max(d$n), status = "reached")
    }) %>%
    ungroup()
}



# One row per scenario, showing power at n = 169 for the targets and the required n.
scenario_table <- function(scenarios, caches, required) {
  lvl <- \(a) if_else(a == ALPHAS[1], "05", "025")
  at169 <- map(caches, cell_targets) %>% list_rbind() %>%
    filter(n == N_EXAMPLE) %>%
    mutate(lvl = lvl(alpha))
  pick <- function(col, m) {
    at169 %>% filter(method == m) %>%
      transmute(scenario_id, name = paste0(col, "_", m, "_", lvl), value = .data[[col]]) %>%
      pivot_wider(names_from = name, values_from = value)
  }
  power_cols <- list(pick("FW_M1", "bonf"), pick("FW_M2", "bonf"), pick("FW_joint", "bonf"),
                     pick("PT_min_all8", "holm"), pick("PT_min_all8", "bonf"), pick("PT_min_all8", "none"),
                     pick("complete_all8", "holm"), pick("complete_all8", "bonf"),
                     pick("pair_M1", "none"), pick("pair_M2", "none"), pick("pair_both", "none")) %>%
    reduce(full_join, by = "scenario_id")
  min_effect <- at169 %>% filter(method == "holm") %>%
    transmute(scenario_id, name = paste0("PT_min_effect_holm_", lvl), PT_min_effect) %>%
    pivot_wider(names_from = name, values_from = PT_min_effect)
  n_cols <- required %>%
    mutate(short = c(T1_FW_both_mediators = "T1", T2_FW_M1 = "T2", T3_min_individual_power_8 = "T3", T4_pair = "T4")[target],
           name = paste0("n_", short, "_", method, "_", lvl(alpha))) %>%
    select(scenario_id, name, n_label) %>%
    pivot_wider(names_from = name, values_from = n_label)
  scenarios %>%
    select(order, scenario_id, label, in_tableB, in_supplement, term, value, am1, am2, m1m2, am1m2, a1, a2, R2) %>%
    left_join(power_cols, by = "scenario_id") %>%
    left_join(min_effect, by = "scenario_id") %>%
    left_join(n_cols, by = "scenario_id") %>%
    arrange(order)
}


# Run -------------------------------------------------------------------------

scenarios <- build_scenarios()
selected  <- scenarios %>%
  arrange(desc(in_tableB), order)

caches <- map(set_names(selected$scenario_id), load_cache)

for (i in seq_len(nrow(selected))) {
  s <- selected[i, ]
  caches[[s$scenario_id]] <- get_cell(caches[[s$scenario_id]], s, N_EXAMPLE)
}
if (RUN_REQN) {
  for (i in seq_len(nrow(selected))) {
    s <- selected[i, ]
    caches[[s$scenario_id]] <- search_n(s, caches[[s$scenario_id]])
  }
}

required <- map(caches, \(c) required_n_table(c) %>% mutate(scenario_id = c$scenario_id, .before = 1)) %>% list_rbind()
timing   <- map(caches, \(c) map(c$cells, \(cell) tibble(scenario_id = cell$scenario_id, n = cell$n, nreps = cell$nreps,
                                                          n_failed = cell$n_failed, seconds = cell$seconds)) %>%
                  list_rbind()) %>% list_rbind()
tab      <- scenario_table(selected, caches, required)

# Power at n = 169 by target (Table S20 of the paper) and the required n (Table S21)
print(tab %>% select(order, scenario_id, label, starts_with("FW_M"), starts_with("PT_min_all8"), starts_with("pair_M")),
      n = Inf, width = Inf)
print(tab %>% select(order, scenario_id, label, starts_with("n_")), n = Inf, width = Inf)
message(glue("Failed replications: {sum(timing$n_failed)}. ",
             "Total simulation time {round(sum(timing$seconds) / 60, 1)} min on {MC_CORES} core(s)."))
