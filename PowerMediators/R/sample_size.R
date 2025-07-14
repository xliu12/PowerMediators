#' Sample size planning for causal mediation analysis with two mediators
#'
#' The function uses each of a sequence of sample sizes to conduct the power analysis for each mediator. The sample size result is obtained as the smallest sample size with which a user-specified target for the power analysis is reached.
#'
#' @param n.min [\code{integer(1)}]\cr
#' An integer specifying the starting value of a sequence of sample sizes. Power will be calculated for each sample size of this sequence.
#' @param n.max [\code{integer(1)}]\cr
#' An integer specifying the end value of the sequence of sample sizes.
#' @param n.step [\code{integer(1)}]\cr
#' An integer specifying the increment of the sequence of sample sizes.
#' @param power_target [\code{numeric(1)}]\cr
#' The targeted power level. Default is 0.8.
#' @param FW_or_PT [\code{character}]\cr
#' The type(s) of power to calculate. Options are \code{"FW"} for the familywise power of testing all indirect effects via a mediator, or \code{"PT"} for the per-test power of testing one of the indirect effects. Default is \code{FW_or_PT = c("FW", "PT")}.
#' @param which_mediator [\code{character}]\cr
#' When \code{FW_or_PT = "FW"} is specified, \code{which_mediator} specifies the mediator(s) for which the familywise power should be used in the sample size determination. Options are \code{"M1"} or \code{"M2"} or \code{c("M1", "M2")}.
#' If \code{which_mediator = c("M1", "M2")} and \code{FW_or_PT = "FW"} (the default), the sample size is determined such that both the familywise power for mediator M1 and the familywise power for mediator M2 are higher than or equal to \code{power_target}.
#' @param which_effect [\code{character}]\cr
#' When \code{FW_or_PT = "PT"} is specified, \code{which_effect} specifies the indirect effect(s) for which the per-test power should be used in the sample size determination. Options are \code{"all"} or one or some of the following interventional indirect effects (IIEs): \code{"IIE_M1(1,,1)"}, \code{"IIE_M1(1,,0)"}, \code{"IIE_M1(0,,1)"}, \code{"IIE_M1(0,,0)"}, \code{"IIE_M2(1,1,)"}, \code{"IIE_M2(1,0,)"}, \code{"IIE_M2(0,1,)"}, \code{"IIE_M2(0,0,)"}. See the accompanying referrence manuscript for definitions of these effects.
#' If \code{which_effect = "all"} and \code{FW_or_PT = "PT"}, the sample size is determined such that the per-test power for any of these indirect effects is higher than or equal to \code{power_target}.
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
#' @param multiple_testing_adjustment [\code{character}]\cr
#' Method(s) to be used for multiple testing adjustment.
#' Default is \code{multiple_testing_adjustment = c("no_adjust", "bonferroni", "modified_bon1", "modified_bon2")}.
#' Options include (1) \code{"no_adjust"}, which will use \code{sig.level} as the per-test Type I error rate;
#' (2)  \code{"bonferroni"}, which will use the Bonferroni adjustment with \code{sig.level} as the familywise Type I error rate for each mediator, leading to the per-test Type I error rate of \code{sig.level}/4;
#' (3)  \code{"modified_bon1"}, which will use the modified Bonferroni-1 adjustment with \code{sig.level} as the familywise Type I error rate for each mediator (\eqn{\alpha_{FW}}), leading to the per-test Type I error rate of \eqn{ \alpha_{PT,j} = \frac{\alpha_{FW}}{k -(k-1)\sqrt{|\bar r_j|}}};
#' (4)  \code{"modified_bon2"}, which will use the modified Bonferroni-2 adjustment with \code{sig.level} as the familywise Type I error rate for each mediator (\eqn{\alpha_{FW}}), leading to the per-test Type I error rate of \eqn{  \alpha_{PT ,j}= \frac{\alpha_{FW}}{k^{1-\sqrt{|\bar r_j|}}}}.
#' For the two modified Bonferroni adjustments,  \eqn{j=1,2,3,4} index the four indirect effects in the family for each mediator, \eqn{k=4} is the number of tests, and  \eqn{|\bar r_j|} is the average absolute correlation between the estimator of the j-th effect and the estimators of the other three effects in the family. These two adjustments are extended from Smith, C. E., & Cribbie, R. A. (2013). Multiplicity control in structural equation modeling: Incorporating parameter dependencies Structural Equation Modeling: A Multidisciplinary Journal, 20(1), 79–85.
#'
#'
#'
#'
#' @return A \code{list} containing the sample size calculation results.
#' \item{res_n}{A data frame of power calculation results for the sequence of sample sizes.}
#' \item{power_FW_M1M2}{A data frame containing the minimum sample size (\code{n_required}) with which both familywise power levels for mediators M1 and M2 reach at least 0.80. The data frame also contains the corresponding power level. If \code{FW_or_PT = "PT"} is specified, this will be \code{NULL}.}
#' \item{power_FW_eachM}{A data frame containing the minimum sample size (\code{n_required}) with which the familywise power level for each mediator specified in \code{which_mediator} reaches at least 0.80. The data frame also contains the corresponding power level. If \code{FW_or_PT = "PT"} is specified, this will be \code{NULL}.}
#' \item{power_PT_someIIE}{A data frame containing the minimum sample size (\code{n_required}) with which the per-test power levels reach at least 0.80 for the indirect effects specified in \code{which_effect}. The data frame also contains the corresponding power level. If \code{FW_or_PT = "FW"} is specified, this will be \code{NULL}.}
#' \item{power_PT_allIIE}{A data frame containing the minimum sample size (\code{n_required}) with which the per-test power levels reach at least 0.80 for all the indirect effects for mediators M1 and M2. The data frame also contains the corresponding power level. If \code{FW_or_PT = "FW"} is specified, this will be \code{NULL}.}
#' \item{power_PT_eachIIE}{A data frame containing,  for each indirect effect, the minimum sample size (\code{n_required}) with which the per-test power level reach at least 0.80. The data frame also contains the corresponding power level. If \code{FW_or_PT = "FW"} is specified, this will be \code{NULL}.}
#' \item{call}{The matched call.}
#'
#'
#' @export
#'
#' @example examples/example_n.R



cal.sample_size <- function(seed_user = 12,
                            n.min = 200,
                            n.max = 300,
                            n.step = 20,
                            power_target = 0.8,
                            FW_or_PT = c("FW", "PT"),
                            which_mediator = c("M1", "M2"),
                            which_effect = c("all"),

                            M_binary = c(FALSE, FALSE),
                            Y_binary = FALSE,

                            num_x = 2,
                            treat.prop = 0.5,
                            treat.randomized = FALSE,

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
                            multiple_testing_adjustment = c("no_adjust", "bonferroni", "modified_bon1", "modified_bon2"),

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

  res_n <- do.call(rbind.data.frame, res_n_list)

  res_n <- res_n %>% ungroup() %>%
    filter(sig.adjust %in% multiple_testing_adjustment)

  power_FW_M1M2 <- res_n %>%
    group_by(n, sig.adjust) %>%
    summarise(min_power_FW = min(power_FW), achieve_target = min(power_FW) >= power_target) %>%
    group_by(sig.adjust, achieve_target) %>%
    summarise(n_required = min(n), power = min(min_power_FW)) %>%
    filter(achieve_target == TRUE)

  if (!("FW" %in% FW_or_PT)) {
    power_FW_M1M2 <- NULL
  }


  power_FW_eachM <- res_n %>%
    group_by(n, sig.adjust, mediation) %>%
    summarise(min_power_FW = min(power_FW), achieve_target = min(power_FW) >= power_target) %>%
    group_by(sig.adjust, achieve_target, mediation) %>%
    summarise(n_required = min(n), power = min(min_power_FW)) %>%
    filter(achieve_target == TRUE)

  if (!("FW" %in% FW_or_PT)) {
    power_FW_eachM <- NULL
  } else {
    if (length(which_mediator) == 1) {
      power_FW_eachM <- power_FW_eachM %>% filter(mediation %in% glue("IIE_{which_mediator}"))
    }
  }

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

  if (!("PT" %in% FW_or_PT)) {
    power_PT_eachIIE <- NULL
    power_PT_allIIE <- NULL
    power_PT_someIIE <- NULL
  }

  n_result <- list(res_n = res_n,
       power_FW_M1M2 = power_FW_M1M2,
       power_FW_eachM = power_FW_eachM,
       power_PT_someIIE = power_PT_someIIE,
       power_PT_allIIE = power_PT_allIIE,
       power_PT_eachIIE = power_PT_eachIIE,
       call = match.call())


  # out <- mget(ls(envir = environment()))

  # return(out)
  n_result
}
