#' Power Analysis for Causal Mediation Analysis
#'
#' Implement a simulation-based method for statistical power calculation in causal mediation analysis with two correlated mediators.
#'
#' @param n Sample size.
#'
#' @return A \code{data.frame} containing the Type I error rates and statistical power for testing interventional indirect effects via each mediator. Specifically, the function outputs (1) the familywise ("FW") Type I error rate and power of testing all four indirect effects for a mediator (i.e., the probability of finding significance for at least one of the indirect effects for a mediator, when that indirect effect were actually null and non-null, respectively), and (2) the per-test ("PT") Type I error rate and power of testing only one indirect effect for a mediator.
#'
#'
#' @export
#'
runPower <- function(seed_user = 12,
                     n = 100,
                     num_x = 2,
                     treat.prop = 0.5,
                     treat.randomized = FALSE,
                     
                     R2.ax = NULL, R2.mx = NULL, R2.yx = NULL, 
                     std.m_on_a = NULL,
                     std.y_on_a = NULL,
                     std.y_on_m = NULL,
                     std.y_on_am_2way = NULL,
                     std.y_on_m_2way = NULL,
                     std.y_on_am_3way = NULL,
                     
                     a_on_x = sqrt(0.13), 
                     m_on_a = rep(0.14, 2), 
                     m_on_x = rep(sqrt(0.13), 2),
                     em_corr = 0,
                     # R2x.m = rep(0.13, 2)
                     M_binary = rep(FALSE, 2),
                     
                     # y_formula = "Y ~ A + M1 + M2 + X + A:M1 + A:M2 + M1:M2 + A:M1:M2",
                     y_on_a = 0.1, # standardized coefficient
                     y_on_m = rep(0.39, 2),
                     y_on_am_2way = rep(0.1, 2),
                     y_on_m_2way = rep(0.1, 2 * (2-1)/2),
                     y_on_am_3way = rep(0, 2 * (2-1)/2),
                     y_on_x = sqrt(0.02),
                     Y_binary = FALSE,
                     
                     # for confidence interval and estimation
                     nboot = NULL, # for bootstrap CI
                     n.draws = 1000, # for Monte Carlo CI
                     nreps.binYM = 1000, # for imputation-based estimation of IIEs,
                     
                     # for power
                     sig.level = 0.05,
                     # for simulation-based power calculation
                     nsims = 1000, mc.cores = 1
){

  set.seed(seed_user)
  datseeds <- sample(1:1e6, nsims)

  res_reps <- parallel::mclapply(datseeds,
                       one.dataset,
                       # arguments for one.dataset()
                       n = n,
                       num_x = num_x,
                       treat.prop = treat.prop,
                       treat.randomized = treat.randomized,
                       
                       R2.ax = R2.ax, R2.mx = R2.mx, R2.yx = R2.yx, 
                       std.m_on_a = std.m_on_a,
                       std.y_on_a = std.y_on_a,
                       std.y_on_m = std.y_on_m,
                       std.y_on_am_2way = std.y_on_am_2way,
                       std.y_on_m_2way = std.y_on_m_2way,
                       std.y_on_am_3way = std.y_on_am_3way,
                       
                       a_on_x = a_on_x, 
                       m_on_a = m_on_a, 
                       m_on_x = m_on_x,
                       em_corr = em_corr,
                       M_binary = M_binary,
                       
                       y_on_a = y_on_a, # standardized coefficient
                       y_on_m = y_on_m,
                       y_on_am_2way = y_on_am_2way,
                       y_on_m_2way = y_on_m_2way,
                       y_on_am_3way = y_on_am_3way,
                       y_on_x = y_on_x,
                       Y_binary = Y_binary,
                       
                       # for confidence interval and estimation
                       nboot = nboot, # for bootstrap CI
                       n.draws = n.draws, # for Monte Carlo CI
                       nreps.binYM = nreps.binYM, # for imputation-based estimation of IIEs,
                       
                       # for power
                       sig.level = sig.level,

                       # arguments for mclapply
                       mc.preschedule = TRUE, mc.cores = mc.cores)


  resdf <- do.call(rbind, res_reps)


  # power ----

  result_FW <- resdf %>%
    group_by(n, effect) %>%
    mutate(true_val = mean(true_vals)) %>%
    group_by(n, effect, sig.adjust) %>%
    mutate(rep = 1:n(),
           mediation = substr(effect, 1, 6),
           nreps = n()) %>%
    group_by(n, mediation, sig.adjust
             , rep) %>%
    mutate(
      TypeIerror_FW = sum((true_val == 0) & (if_sig == TRUE))>=1,
      power_FW = sum((true_val != 0) & (if_sig == TRUE))>=1,
      cover_FW = sum((ci_low < (true_val)) & (ci_up > (true_val)))==4
    ) %>%
    group_by(n, mediation, sig.adjust) %>%
    summarise(
      across(c(contains("power"), contains("TypeIerror"), contains("cover")),
             ~mean(., na.rm = TRUE))
    )


  result_PT <- resdf %>%
    group_by(n, effect) %>%
    mutate(true_val = mean(true_vals)) %>%
    group_by(n, effect, sig.adjust) %>%
    mutate(rep = 1:n(),
           nreps = n()) %>%
    mutate(
      TypeIerror_PT = (true_val == 0) & (if_sig == TRUE),
      power_PT = (true_val != 0) & (if_sig == TRUE),
      cover_PT = (ci_low < (true_val)) & (ci_up > (true_val))
    ) %>%
    summarise(
      across(c(contains("power"), contains("TypeIerror"), contains("cover"), contains("true_val"), nreps),
             ~mean(., na.rm = TRUE))
    ) %>%
    mutate(mediation = substr(effect, 1, 6))

  res_power <- full_join(result_FW, result_PT, by = c("n", "mediation", "sig.adjust"))

  res_power
}




one.dataset <- function(iseed = 123,
                        n = 100,
                        num_x = 2,
                        treat.prop = 0.5,
                        treat.randomized = FALSE,
                        
                        R2.ax = NULL, R2.mx = NULL, R2.yx = NULL, 
                        std.m_on_a = NULL,
                        std.y_on_a = NULL,
                        std.y_on_m = NULL,
                        std.y_on_am_2way = NULL,
                        std.y_on_m_2way = NULL,
                        std.y_on_am_3way = NULL,
                        
                        a_on_x = sqrt(0.13), 
                        m_on_a = rep(0.14, 2), 
                        m_on_x = rep(sqrt(0.13), 2),
                        em_corr = 0,
                        M_binary = rep(FALSE, 2),
                        
                        # y_formula = "Y ~ A + M1 + M2 + X + A:M1 + A:M2 + M1:M2 + A:M1:M2",
                        y_on_a = 0.1, # standardized coefficient
                        y_on_m = rep(0.39, 2),
                        y_on_am_2way = rep(0.1, 2),
                        y_on_m_2way = rep(0.1, 2 * (2-1)/2),
                        y_on_am_3way = rep(0, 2 * (2-1)/2),
                        y_on_x = sqrt(0.02),
                        Y_binary = FALSE,
                        
                        # for confidence interval and estimation
                        nboot = NULL, # for bootstrap CI
                        n.draws = 1000, # for Monte Carlo CI
                        nreps.binYM = 1000, # for imputation-based estimation of IIEs,
                        
                        # for power
                        sig.level = 0.05
){
  
  gen_data <- gen.data(
    iseed = iseed,
    n = n,
    num_x = num_x,
    treat.prop = treat.prop,
    treat.randomized = treat.randomized,
    
    R2.ax = R2.ax, R2.mx = R2.mx, R2.yx = R2.yx, 
    std.m_on_a = std.m_on_a,
    std.y_on_a = std.y_on_a,
    std.y_on_m = std.y_on_m,
    std.y_on_am_2way = std.y_on_am_2way,
    std.y_on_m_2way = std.y_on_m_2way,
    std.y_on_am_3way = std.y_on_am_3way,
    
    a_on_x = a_on_x, 
    m_on_a = m_on_a,
    m_on_x = m_on_x,
    em_corr = em_corr,
    M_binary = M_binary,

    y_on_a = y_on_a, 
    y_on_m = y_on_m,
    y_on_am_2way = y_on_am_2way,
    y_on_m_2way = y_on_m_2way,
    y_on_am_3way = y_on_am_3way,
    y_on_x = y_on_x,
    Y_binary = Y_binary,
    nreps.binYM = nreps.binYM
  )


  dat <- gen_data$dat

  true_vals <- gen_data$true_vals

  res <- med2.reg(
    data = dat, M_binary, Y_binary, 
    sig.IIE = sig.level,
    sig.adjust = c("no_adjust", "bonferroni", "modified_bon1", "modified_bon2"),
    nboot = nboot, 
    n.draws = n.draws, # for Monte Carlo CI
    nreps.binYM = nreps.binYM # for imputation-based estimation of IIEs
  )

  res1 <- data.frame(n,
                     full_join(res, true_vals, by = "effect"), sig.level,
                     row.names = NULL)


  return(res1)
}

