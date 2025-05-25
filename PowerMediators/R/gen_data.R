library(mvtnorm)
library(tidyverse)
library(glue)


gen.data <- function(iseed = 123,
  n = 100,
  num_x = 2,
  treat.prop = 0.5,
  treat.randomized = FALSE,
  a_on_x = sqrt(0.13), # standardized coefficient
  
  num_m = 2,
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
  Y_binary = FALSE
) {
  dat <- data.frame(id = 1:n)
  # Covariates --------------
  if (num_x == 0) {
    R2x <- 0
    num_x <- 1
  }
  X <- rmvnorm(n, sigma = diag(1/num_x, nrow = num_x))
  colnames(X) <- glue("X{1:num_x}")
  X_rowsum <- rowSums(X)
  # Treatment ----------------
  # a_on_x = sqrt(0.13)
  gen_a <- list(a_intercept = qnorm(treat.prop), a_on_x = a_on_x)
  a_continuous <- gen_a[["a_intercept"]] + gen_a[["a_on_x"]] * X_rowsum + rnorm(n, sd = sqrt(1 - gen_a[["a_on_x"]]^2))
  A <- 1*(a_continuous > 0)
  A_std <- as.numeric(scale(A))
  # Mediators -----------------
  # m_on_a = c(0, 0.4)
  # em_corr = 0.2
  
  var_em <- map_dbl(1:num_m, \(i) 1 - var(m_on_a[i] * A_std + m_on_x[i] * X_rowsum))
  
  if (any(var_em <= 0)) {
    stop("The residual of a mediator model has a negative variance in at least one replication. Please respecify the mediator model parameters.")
    var_em[var_em < 0] <- 0
  }
  
  if (length(em_corr)==1) {
    em_corr <- c(trt = em_corr, ctrl = em_corr)
  }
  
  Corr_em <- diag(1, nrow = num_m)
  
  em_list <- map(c("trt", "ctrl"), \(a="trt") {
    Corr_em[Corr_em == 0] <- em_corr[a]
    sigma_m <- Corr_em * (sqrt(var_em) %*% t(sqrt(var_em)))
    em <- rmvnorm(n, sigma = sigma_m)
  } )
  names(em_list) <- c("trt", "ctrl")
  
  m_intercept <- rep(0, num_m)
  gen_m <- list(m_intercept = rep(0, num_m), m_on_a = m_on_a, m_on_x = m_on_x)
  
  M <- map(1:num_m, \(i=1) {
    m_intercept[i] + m_on_a[i] * A_std + m_on_x[i] * X_rowsum + 
      em_list[["trt"]][, i] * A + em_list[["ctrl"]][, i] * (1 - A)
  }) %>% reduce(cbind)
  
  colnames(M) <- glue("M{1:num_m}")
  if (any(M_binary)) {
    M[, M_binary] <- 1*(M[, M_binary] > 0)
  }
  
  # Outcome ----------------
  y_formula <- "Y ~ A + M1 + M2 + X + A:M1 + A:M2 + M1:M2 + A:M1:M2"
  y_coefs = rep()
  
  wmat <- model.matrix(as.formula(gsub("Y", "", y_formula)), data = data.frame(A, M, X=X_rowsum))
  wmat[1, ]
  wmat_std <- scale(wmat[, -1])
  y_coefs <- numeric(length = ncol(wmat_std))
  names(y_coefs) <- colnames(wmat_std)
  y_coefs["A"] <- y_on_a
  y_coefs[glue("M{1:num_m}")] <- y_on_m
  y_coefs["X"] <- y_on_x
  
  y_coefs[glue("A:M{1:num_m}")] <- y_on_am_2way
  y_coefs[str_count(names(y_coefs), "M") == 2 & str_count(names(y_coefs), "A") == 0] <- y_on_m_2way
  y_coefs[str_count(names(y_coefs), "M") == 2 & str_count(names(y_coefs), "A") == 1] <- y_on_am_3way
  
  var_ey <- 1 - var(wmat_std %*% y_coefs)
  
  if (any(var_ey <= 0)) {
    stop("The residual of the outcome model has a negative variance in at least one replication. Please respecify the outcome model parameters.")
    var_ey <- 0
  }
  
  Y <- wmat_std %*% y_coefs + rnorm(n, sd = sqrt(var_ey))
  if (Y_binary) {
    Y <- 1*(Y>0)
  }
  
  # true values ------------------
  m1_coefs <- c(gen_m$m_intercept[1], gen_m$m_on_a[1], rep(gen_m$m_on_x[1], num_x))
  names(m1_coefs) <- c("(Intercept)", "A", glue("X{1:num_x}"))
  m2_coefs <- c(gen_m$m_intercept[2], gen_m$m_on_a[2], rep(gen_m$m_on_x[2], num_x))
  names(m2_coefs) <- c("(Intercept)", "A", glue("X{1:num_x}"))
  
  ## IIE_M1 --------
  a0a2 <- expand.grid(a0=c(0,1), a2=c(0,1))
  true_IIE_M1 <- map_dbl(1:nrow(a0a2), \(i) {
    cal.IIE_M1(a0 = a0a2$a0[i], a2 = a0a2$a2[i], y_coefs, m1_coefs, m2_coefs)
  })
  names(true_IIE_M1) <- glue("IIE_M1({a0a2$a0},,{a0a2$a2})")
  
  ## IIE_M2 --------
  a0a1 <- expand.grid(a0=c(0,1), a1=c(0,1))
  true_IIE_M2 <- map_dbl(1:nrow(a0a1), \(i) {
    cal.IIE_M2(a0 = a0a1$a0[i], a1 = a0a1$a1[i], y_coefs, m1_coefs, m2_coefs)
  })
  names(true_IIE_M2) <- glue("IIE_M2({a0a1$a0},{a0a1$a1},)")
  
  true_vals <- c(true_IIE_M1, true_IIE_M2)
  true_vals <- data.frame(true_vals) %>% rownames_to_column(var = "effect")
  
  # out -----
  dat <- data.frame(id = 1:n, A, M, Y, X)
  
  out <- mget(ls(envir = environment()))
  
  return(out)
}

condition_all <- data.frame(
  expand.grid(
    n = c(200),
    m1_on_a = c(0, 0.3),
    y_on_am1 = c(0, 0.1),
    y_on_m1m2 = c(0, 0.1),
    y_on_am1m2 = c(0, 0.1)
  )) 

condition <- condition_all %>% 
  filter(y_on_am1 != 0, y_on_m1m2 != 0, y_on_am1m2 != 0)

OneData <- function(iseed = 123, cond = 1){
  num_m <- 2
  gen_data <- gen.data(
    iseed = iseed,
    n = condition$n[cond],
    num_x = 2,
    treat.prop = 0.5,
    treat.randomized = FALSE,
    a_on_x = sqrt(0.13), # standardized coefficient
    
    num_m = 2,
    m_on_a = c(condition$m1_on_a[cond], 0.3),
    m_on_x = rep(sqrt(0.13), num_m),
    em_corr = 0.2,
    # R2x.m = rep(0.13, num_m)
    M_binary = rep(FALSE, num_m),
    
    # y_formula = "Y ~ A + M1 + M2 + X + A:M1 + A:M2 + M1:M2 + A:M1:M2",
    y_on_a = 0.1, # standardized coefficient
    y_on_m = rep(0.3, num_m),
    y_on_am_2way = c(condition$y_on_am1[cond], 0.1),
    y_on_m_2way = condition$y_on_m1m2[cond],
    y_on_am_3way = condition$y_on_am1m2[cond],
    y_on_x = sqrt(0.02),
    Y_binary = FALSE
  )
  
  
  dat <- gen_data$dat
  
  true_vals <- gen_data$true_vals
  
  res <- med2.reg(
    data = dat,
    sig.IIE = 0.05,
    sig.adjust = c("no_adjust", "bonferroni", "modified_bon1", "modified_bon2"),
    nboot = 1000
  )
  
  res1 <- data.frame(condition[cond, ],
                     full_join(res, true_vals, by = "effect"), 
                    row.names = NULL)
  
  
  return(res1)
}

