#' Sample size planning to achieve a given power 
#' 
#' The function uses each of a sequence of sample sizes to conduct the power analysis for each mediator. The sample size suggestion is obtained as the smallest sample size with which a user-specified objective for a power analysis result is reached. 
#' 
#' @param n.min [\code{integer(1)}]\cr
#' An integer specifying the starting value of a sequence of sample sizes. Power will be calculated for each sample size of this sequence.
#' @param n.max [\code{integer(1)}]\cr
#' An integer specifying the end value of the sequence of sample sizes. 
#' @param by [\code{integer(1)}]\cr
#' An integer specifying the increment of the sequence of sample sizes.
#' @param power_target [\code{numeric(1)}]\cr
#' The targeted power level. Default is 0.8.
#' @param FW_or_PT [\code{character(1)}]\cr
#' The type of power to calculate. Options are \code{"FW"} for the familywise power of testing all indirect effects via a mediator, or \code{"PT"} for the per-test power of testing one of the indirect effects. Default is \code{FW_or_PT = "FW"}.
#' @param which_mediator [\code{character}]\cr
#' When \code{FW_or_PT = "FW"} is specified, \code{which_mediator} specifies the mediator(s) for which the familywise power should be used in the sample size determination. Options are \code{"M1"} or \code{"M2"} or \code{c("M1", "M2")}. 
#' If \code{which_mediator = c("M1", "M2")} and \code{FW_or_PT = "FW"} (the default), the sample size is determined such that both the familywise power for mediator M1 and the familywise power for mediator M2 are higher than or equal to \code{power_target}. 
#' @param which_effect [\code{character}]\cr
#' When \code{FW_or_PT = "PT"} is specified, \code{which_effect} specifies the indirect effect(s) for which the per-test power should be used in the sample size determination. Options are \code{"all"} or one or some of the following interventional indirect effects (IIEs): \code{"IIE_M1(1,,1)"}, \code{"IIE_M1(1,,0)"}, \code{"IIE_M1(0,,1)"}, \code{"IIE_M1(0,,0)"}, \code{"IIE_M2(1,1,)"}, \code{"IIE_M2(1,0,)"}, \code{"IIE_M2(0,1,)"}, \code{"IIE_M2(0,0,)"}. See the accompanying referrence manuscript for definitions of these effects. 
#' If \code{which_effect = "all"} and \code{FW_or_PT = "PT"}, the sample size is determined such that the per-test power for any of these indirect effects is higher than or equal to \code{power_target}. 
#'
#' @export



cal.sample_size <- function(seed_user = 12,
                     n.min = 200,
                     n.max = 300,
                     n.step = 20, 
                     power_target = 0.8,
                     FW_or_PT = "FW",
                     which_mediator = c("M1", "M2"),
                     which_effect = c("all"), 
                     multiple_testing_adjustment = c("no_adjust", "bonferroni", "modified_bon1", "modified_bon1"),
                     
                     num_x = 2,
                     treat.prop = 0.5,
                     treat.randomized = FALSE,
                     
                     R2.ax = 0, R2.mx = 0.1, R2.yx = 0.1, 
                     std.m_on_a = c(0.36, 0.36),
                     std.y_on_a = 0.14,
                     std.y_on_m = c(0, 0),
                     std.y_on_am_2way = c(0, 0.36),
                     std.y_on_m_2way = 0,
                     std.y_on_am_3way = 0.36,
                     
                     a_on_x = NULL, 
                     m_on_a = NULL, 
                     m_on_x = NULL,
                     
                     em_corr = 0.1,
                     M_binary = c(FALSE, FALSE),
                     
                     # y_formula = "Y ~ A + M1 + M2 + X + A:M1 + A:M2 + M1:M2 + A:M1:M2",
                     y_on_a = NULL, # standardized coefficient
                     y_on_m = NULL,
                     y_on_am_2way = NULL,
                     y_on_m_2way = NULL,
                     y_on_am_3way = NULL,
                     y_on_x = NULL,
                     
                     Y_binary = FALSE,
                     
                     # for confidence interval and estimation
                     nboot = NULL, # for bootstrap CI
                     n.draws = 1000, # for Monte Carlo CI
                     
                     # for power
                     sig.level = 0.05,
                     # for simulation-based power calculation
                     nsims = 1000, mc.cores = 1
) {
  n_seq <- seq(n.min, n.max, by = n.step)
  
  res_n_list <- lapply(n_seq, \(n) {
    runPower(
      seed_user = seed_user,
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
      
      # y_formula = "Y ~ A + M1 + M2 + X + A:M1 + A:M2 + M1:M2 + A:M1:M2",
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
      
      # for power
      sig.level = sig.level,
      # for simulation-based power calculation
      nsims = nsims, mc.cores = mc.cores
    )
  })
  
  res_n <- do.call(rbind, res_n_list)
  
  power_FW_M1M2 <- res_n %>% 
    group_by(n, sig.adjust) %>% 
    summarise(min_power_FW = min(power_FW), achieve_target = min(power_FW) >= power_target) %>% 
    group_by(sig.adjust, achieve_target) %>% 
    summarise(n_required = min(n), power = min(min_power_FW)) %>% 
    filter(achieve_target == TRUE)
  
  power_FW_eachM <- res_n %>% 
    group_by(n, sig.adjust, mediation) %>% 
    summarise(min_power_FW = min(power_FW), achieve_target = min(power_FW) >= power_target) %>% 
    group_by(sig.adjust, achieve_target, mediation) %>% 
    summarise(n_required = min(n), power = min(min_power_FW)) %>% 
    filter(achieve_target == TRUE)
  
  power_PT_eachIIE <- res_n %>% 
    group_by(n, sig.adjust, mediation, effect) %>% 
    summarise(min_power_PT = min(power_PT), achieve_target = min(power_PT) >= power_target) %>% 
    group_by(sig.adjust, achieve_target, mediation, effect) %>% 
    summarise(n_required = min(n), power = min(min_power_PT)) %>% 
    filter(achieve_target == TRUE)
  
  power_PT_allIIE <- res_n %>% 
    group_by(n, sig.adjust) %>% 
    summarise(min_power_PT = min(power_PT), achieve_target = min(power_PT) >= power_target) %>% 
    group_by(sig.adjust, achieve_target) %>% 
    summarise(n_required = min(n), power = min(min_power_PT)) %>% 
    filter(achieve_target == TRUE)
  
  if (which_effect[1] == "all") {
    which_effect <- unique(res_n$effect)
  }
  if (!is.null(which_effect)) {
    power_PT_someIIE <- res_n %>% 
      filter(effect %in% which_effect) %>% 
      group_by(n, sig.adjust) %>% 
      summarise(min_power_PT = min(power_PT), achieve_target = min(power_PT) >= power_target) %>% 
      group_by(sig.adjust, achieve_target) %>% 
      summarise(n_required = min(n), power = min(min_power_PT)) %>% 
      filter(achieve_target == TRUE)
  } else {
    power_PT_someIIE <- NULL
  }
  
  
  
  list(res_n = res_n,
       power_FW_M1M2 = power_FW_M1M2,
       power_FW_eachM = power_FW_eachM,
       power_PT_eachIIE = power_PT_eachIIE, 
       power_PT_allIIE = power_PT_allIIE, 
       power_PT_someIIE = power_PT_someIIE)
  
}