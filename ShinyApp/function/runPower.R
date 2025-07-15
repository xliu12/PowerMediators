#' Power analysis for causal mediation analysis with two mediators
#'
#' Simulation-based method for statistical power analysis of testing interventional indirect effects for each mediator in a causal mediation analysis with two correlated mediators. Both familywise power and per-test power can be calculated. Currently support binary treatments, binary or continuous mediators, and binary or continuous outcomes.
#'
#' @param n [\code{integer(1)}]\cr Sample size to be used for power calculation. Default is 100.
#' @param num_x [\code{integer(1)}]\cr The number of baseline covariates to be included in the models for power analysis. Default is \code{num_x = 2}.
#' @param treat.prop [\code{numeric(1)}]\cr The proportion of the sample allocated to the treatment condition. Default is \code{treat.prop = 0.5}, indicating equal allocation between the treatment and control conditions.
#' @param treat.randomized [\code{logic(1)}]\cr Whether the treatment variable is randomized. Default is \code{treat.randomized = FALSE}, allowing the treatment variable to depend on the baseline covariates. If \code{treat.randomized = TRUE} is specified, user-specification of \code{R2.ax} and \code{a_on_x} will be ignored and \code{R2.ax = 0} and \code{a_on_x = 0} will be used, indicating the treatment variable is independent of the covariates.
#' @param R2.ax [\code{numeric(1)}]\cr The proportion of variance in the treatment's latent continuous score explained by the covariates. If \code{treat.randomized = TRUE} is specified, \code{R2.ax = 0} will be used.
#' @param R2.mx [\code{numeric(2)}]\cr The proportion of variance in each of the two mediators (or the mediator's latent continuous score) that is explained by the covariates. Default is \code{R2.mx = c(0.1, 0.1)}.
#' @param R2.yx [\code{numeric(1)}]\cr The proportion of variance in the outcome (or the outcome's latent continuous score) that is explained by the covariates. Default is \code{R2.yx = 0.1}.
#' @param std.m_on_a [\code{numeric(2)}]\cr The standardized coefficients in the models regressing each mediator (or the mediator's latent continuous score) on the treatment variable, adjusting for the covariates. Default is \code{std.m_on_a = c(0.3, 0.3)}.
#' @param std.y_on_a [\code{numeric(1)}]\cr The standardized coefficient for the treatment in the outcome model that includes the treatment, mediators, interactions among the treatment and mediators, and the covariates. Default is \code{std.y_on_a = 0.1}.
#' @param std.y_on_m [\code{numeric(2)}]\cr The standardized coefficients for the two mediators in the outcome model that includes the treatment, mediators, interactions among the treatment and mediators, and the covariates. Default is \code{std.y_on_m = c(0.3, 0.3)}.
#' @param std.y_on_am_2way [\code{numeric(2)}]\cr The standardized coefficients for the two-way interactions of the treatment with each of the two mediators, in the outcome model that includes the treatment, mediators, interactions among the treatment and mediators, and the covariates. Default is \code{std.y_on_am_2way = c(0.1, 0.1)}.
#' @param std.y_on_m_2way [\code{numeric(1)}]\cr The standardized coefficient for the two-way interaction between the two mediators, in the outcome model that includes the treatment, mediators, interactions among the treatment and mediators, and the covariates. Default is \code{std.y_on_m_2way = 0.1}.
#' @param std.y_on_am_3way [\code{numeric(1)}]\cr The standardized coefficient for the three-way interaction of the treatment and two mediators, in the outcome model that includes the treatment, mediators, interactions among the treatment and mediators, and the covariates. Default is \code{std.y_on_am_3way = 0.1}.
#' @param em_corr [\code{numeric}]\cr The correlation between mediators' residuals. Default is \code{em_corr = 0}. If the mediators' residual correlation can differ between the treatment and control conditions, a user may specify a vector of 2 values, where the first and second elements correspond to the residual correlations under the treatment and control conditions, respectively.
#' @param M_binary [\code{logic(2)}]\cr Whether the mediators are binary. Options include (1) \code{M_binary = c(FALSE, FALSE)} for two continuous mediators (the default),
#' (2) \code{M_binary = c(TRUE, FALSE)} for a binary mediator (labeled with "M1") and a continuous mediator (labeled with "M2"), and
#' (3) \code{M_binary = c(TRUE, TRUE)} for two binary mediators.
#' @param Y_binary [\code{logic(1)}]\cr Whether the outcome is binary. Default is \code{Y_binary = FALSE}, indicating the outcome variable is continuous.
#' @param a_on_x [\code{numeric(1)}]\cr Optional. The unstandardized regression coefficient corresponding to \code{R2.ax = 0}.
#' @param m_on_a [\code{numeric(2)}]\cr Optional. The unstandardized coefficients for the treatment in the models for each mediator.
#' @param m_on_x [\code{numeric(2)}]\cr Optional. The unstandardized coefficients for each covariate in the models for each mediator. Within a mediator's model, the covariates have equal coefficients.
#' @param y_on_a [\code{numeric(1)}]\cr Optional. The unstandardized coefficient for the treatment in the outcome model.
#' @param y_on_x [\code{numeric(1)}]\cr Optional. The unstandardized coefficient for each covariate in the outcome model. The covariates have equal coefficients.
#' @param y_on_m [\code{numeric(2)}]\cr Optional. The unstandardized coefficients for the two mediators in the outcome model.
#' @param y_on_am_2way [\code{numeric(2)}]\cr Optional. The unstandardized coefficients for the two-way interactions of the treatment with each of the two mediators in the outcome model.
#' @param y_on_m_2way [\code{numeric(1)}]\cr Optional. The unstandardized coefficient for the two-way interaction between the two mediators in the outcome model.
#' @param y_on_am_3way [\code{numeric(1)}]\cr Optional. The unstandardized coefficient for the three-way interaction of the treatment and two mediators in the outcome model.
#' @param n.draws [\code{integer(1)}]\cr The number of random draws for obtaining the Monte Carlo confidence intervals for the indirect effects. Default is \code{n.draws = 1000}.
#' @param nboot [\code{integer(1)}]\cr Optional. The number of bootstrap samples if a user would like to use bootstrapping to obtain confidence intervals for the indirect effects.
#' @param sig.level [\code{numeric(1)}]\cr The significance level to be used as (i) the familywise Type I error rate in multiple testing adjustments and (ii) the per-test Type I error rate with no adjustment. Default is \code{sig.level = 0.05}. Currently implement the following multiple testing adjustments: the Bonferroni adjustment, and two modified Bonferroni adjustments that incorporate correlations among the indirect effect estimators.
#' @param nsims [\code{integer(1)}]\cr The number of replications to be used in the simulation-based power analysis. Default is \code{nsims = 1000}.
#' @param mc.cores [\code{integer(1)}]\cr The number of cores to use. Default is \code{mc.cores = 1}.
#' @param seed_user [code{integer}]\cr Seed to be used by the \code{set.seed()} for offsetting the random number generator. Default is to leave the random number generator alone.
#'
#'
#'
#' @return A \code{data.frame} containing the statistical power for testing interventional indirect effects (IIE) via each mediator, including (1) the familywise ("FW") power of testing all four indirect effects for a mediator, and (2) the per-test ("PT") power of testing only one indirect effect for a mediator. The \code{data.frame} also contains the multiple testing adjustments (\code{sig.adjust}), Type I error rate (which is 0 if the true value is non-null), confidence interval coverage rate (\code{cover_FW} and \code{cover_PT}), the mediator (\code{IIE_M1} for mediator M1 and \code{IIE_M2} for mediator M2), and the indirect effects (IIEs; \code{effect}).
#'
#'
#' @export
#'
#' @import dplyr
#' @importFrom purrr map
#'
#'
#'
#' @example examples/example_power.R

runPower <- function(seed_user = NULL,
                     n = 100,
                     num_x = 2,
                     treat.prop = 0.5,
                     treat.randomized = FALSE,
                     M_binary = c(FALSE, FALSE),
                     Y_binary = FALSE,

                     R2.ax = 0, R2.mx = c(0.1, 0.1), R2.yx = 0.1,

                     std.m_on_a = c(0.3, 0.3),
                     em_corr = 0,

                     std.y_on_a = 0.1,
                     std.y_on_m = c(0.3, 0.3),
                     std.y_on_am_2way = c(0.1, 0.1),
                     std.y_on_m_2way = 0.1,
                     std.y_on_am_3way = 0.1,

                     a_on_x = 0.13,
                     m_on_a = c(0.3, 0.3),
                     m_on_x = c(0.1, 0.1),
                     # R2x.m = rep(0.13, 2)

                     # y_formula = "Y ~ A + M1 + M2 + X + A:M1 + A:M2 + M1:M2 + A:M1:M2",
                     y_on_a = 0.1, # standardized coefficient
                     y_on_m = c(0.3, 0.3),
                     y_on_am_2way = c(0.1, 0.1),
                     y_on_m_2way = 0.1,
                     y_on_am_3way = 0.1,
                     y_on_x = 0.1,

                     # for confidence interval and estimation
                     nboot = NULL, # for bootstrap CI
                     n.draws = 1000, # for Monte Carlo CI

                     # for power
                     sig.level = 0.05,
                     # for simulation-based power calculation
                     nsims = 1000, mc.cores = 1
){

  if (!is.null(seed_user)) {
    set.seed(seed_user)
  }

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
    Y_binary = Y_binary
  )


  dat <- gen_data$dat

  true_vals <- gen_data$true_vals

  system.time(
    res <- med2.reg(
      data = dat, M_binary, Y_binary,
      sig.IIE = sig.level,
      sig.adjust = c("no_adjust", "bonferroni", "modified_bon1", "modified_bon2"),
      nboot = nboot,
      n.draws = n.draws # for Monte Carlo CI

    )
  )
  res1 <- data.frame(n,
                     full_join(res, true_vals, by = "effect"), sig.level,
                     row.names = NULL)


  return(res1)
}




