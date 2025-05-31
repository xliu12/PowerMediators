library(mvtnorm)
library(tidyverse)
library(glue)


gen.data <- function(iseed = 123,
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
  nreps.binYM = 1000
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
  if (!is.null(R2.ax)) {
    a_on_x <- sqrt(R2.ax)
  }
  gen_a <- list(a_intercept = qnorm(treat.prop), a_on_x = a_on_x)
  a_continuous <- gen_a[["a_intercept"]] + gen_a[["a_on_x"]] * X_rowsum + rnorm(n, sd = sqrt(1 - gen_a[["a_on_x"]]^2))
  A <- 1*(a_continuous > 0)
  
  # cor(a_continuous, X_rowsum)^2
  # Mediators -----------------
  num_m <- 2
  
  # m_on_a = c(0, 0.4)
  # em_corr = 0.2
  
  if(!is.null(R2.mx)) {
    std.m_on_x <- sqrt(R2.mx) - std.m_on_a * cor(A,X_rowsum)
    m_on_x <- std.m_on_x
  }
  if (!is.null(std.m_on_a)) {
    m_on_a <- std.m_on_a * 1 / sd(A)
  }
  
  var_em <- map_dbl(1:num_m, \(i) 1 - var(m_on_a[i] * A + m_on_x[i] * X_rowsum))
  
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
    m_intercept[i] + m_on_a[i] * A + m_on_x[i] * X_rowsum + 
      # em_list[["ctrl"]][, i]
      em_list[["trt"]][, i] * A + em_list[["ctrl"]][, i] * (1 - A)
  }) %>% reduce(cbind)
  
  
  colnames(M) <- glue("M{1:num_m}")
  if (any(M_binary)) {
    M[, M_binary] <- 1*(M[, M_binary] > 0)
  }
  
  # Outcome ----------------
  y_formula <- "Y ~ A + M1 + M2 + X + A:M1 + A:M2 + M1:M2 + A:M1:M2"
  
  wmat <- model.matrix(as.formula(gsub("Y", "", y_formula)), data = data.frame(A, M, X=X_rowsum))
  
  if (!is.null(R2.yx)) {
    std.y_on_x <- sqrt(R2.yx) - std.y_on_a * sqrt(R2.ax) - sum(std.y_on_m * sqrt(R2.mx)) - sum(std.y_on_am_2way * as.numeric(cor(A*M, X_rowsum))) - std.y_on_m_2way * cor(M[,1]*M[,2], X_rowsum) - std.y_on_am_3way * cor(A*M[,1]*M[,2], X_rowsum)
    y_on_x  <- std.y_on_x
  }
  if (!is.null(std.y_on_a)) {
    y_on_a <- std.y_on_a * 1 / sd(A)
    y_on_m <- std.y_on_m * 1 / apply(M, 2, sd)
    y_on_am_2way <- std.y_on_am_2way * 1 / apply(A*M, 2, sd)
    y_on_m_2way <- std.y_on_m_2way * 1 / sd(M[,1]*M[,2])
    y_on_am_3way <- std.y_on_am_3way * 1 / sd(A*M[,1]*M[,2])
  }
  
  y_coefs <- numeric(length = ncol(wmat))
  names(y_coefs) <- colnames(wmat)
  y_coefs["A"] <- y_on_a
  y_coefs[glue("M{1:num_m}")] <- y_on_m
  y_coefs["X"] <- y_on_x
  
  y_coefs[glue("A:M{1:num_m}")] <- y_on_am_2way
  y_coefs[str_count(names(y_coefs), "M") == 2 & str_count(names(y_coefs), "A") == 0] <- y_on_m_2way
  y_coefs[str_count(names(y_coefs), "M") == 2 & str_count(names(y_coefs), "A") == 1] <- y_on_am_3way
  
  var_ey <- 1 - var(wmat %*% y_coefs)
  
  if (any(var_ey <= 0)) {
    stop("The residual of the outcome model has a negative variance in at least one replication. Please respecify the outcome model parameters.")
    var_ey <- 0
  }
  
  Y <- wmat %*% y_coefs + rnorm(n, sd = sqrt(var_ey))
  if (Y_binary) {
    Y <- 1*(Y>0)
  }
  
  dat <- data.frame(id = 1:n, A, M, Y, X)
  data <- dat
  
  # true values ------------------
  m1_coefs <- c(gen_m$m_intercept[1], gen_m$m_on_a[1], gen_m$m_on_x[1])
  m1_coefs <- as.data.frame(t(m1_coefs))
  names(m1_coefs) <- c("(Intercept)", "A", "X")
  m2_coefs <- c(gen_m$m_intercept[2], gen_m$m_on_a[2], gen_m$m_on_x[2])
  m2_coefs <- as.data.frame(t(m2_coefs))
  names(m2_coefs) <- c("(Intercept)", "A", "X")
  
  y_coefs <- as.data.frame(t(y_coefs))

  if (sum(M_binary, Y_binary) == 0) {
    # gaussian Y, M -------------------
    # IIE_M1 
    a0a2 <- expand.grid(a0=c(0,1), a2=c(0,1))
    IIE_M1 <- map(1:nrow(a0a2), \(i) {
      cal.IIE_M1(a0 = a0a2$a0[i], a2 = a0a2$a2[i], y_coefs, m1_coefs, m2_coefs)
    }) %>% reduce(bind_cols)
    # names(IIE_M1) <- glue("IIE_M1({a0a2$a0},,{a0a2$a2})")
    
    # IIE_M2 
    a0a1 <- expand.grid(a0=c(0,1), a1=c(0,1))
    IIE_M2 <- map(1:nrow(a0a1), \(i) {
      cal.IIE_M2(a0 = a0a1$a0[i], a1 = a0a1$a1[i], y_coefs, m1_coefs, m2_coefs)
    }) %>% reduce(bind_cols)
    # names(IIE_M2) <- glue("IIE_M2({a0a1$a0},{a0a1$a1},)")
    
    true_vals <- cbind(IIE_M1, IIE_M2)[1, ] %>% unlist 
  } else {
    # binary Y or M -----------
    
    m1_formula <- as.formula(paste("M1 ~ A + X"))
    if (M_binary[1]) {
      m1_fit <- glm(m1_formula, data = data.frame(A, M, X=X_rowsum), family = binomial(link = "probit"))
    } else {
      m1_fit <- lm(m1_formula, data = data.frame(A, M, X=X_rowsum))
    }
    
    m2_formula <- as.formula(paste("M2 ~ A + X"))
    if (M_binary[2]) {
      m2_fit <- glm(m2_formula, data = data.frame(A, M, X=X_rowsum), family = binomial(link = "probit"))
    } else {
      m2_fit <- lm(m2_formula, data = data.frame(A, M, X=X_rowsum))
    }
    
    y_formula <- as.formula(paste("Y ~ A + M1 + M2 + X + A:M1 + A:M2 + M1:M2 + A:M1:M2"))
    if (Y_binary) {
      y_fit <- glm(y_formula, data = data.frame(A, M, Y, X=X_rowsum), family = binomial(link = "probit"))
    } else {
      y_fit <- lm(y_formula, data = data.frame(A, M, Y, X=X_rowsum))
    }
    true_vals <- purrr::map(1:nrow(y_coefs), \(i=1) {
      y_fit$coefficients <- y_coefs[i, ] %>% unlist
      m1_fit$coefficients <- m1_coefs[i, ] %>% unlist
      m2_fit$coefficients <- m2_coefs[i, ] %>% unlist
      
      IIEs_tmp <- bYM_cal.IIE(y_fit, m1_fit, m2_fit, data = data.frame(A, M, Y, X=X_rowsum),
                              M_binary,
                              Y_binary, nreps = nreps.binYM)
    }) %>% purrr::reduce(bind_rows)
    
  }
  
  # true_vals <- c(true_IIE_M1, true_IIE_M2)
  true_vals <- data.frame(true_vals) %>% rownames_to_column(var = "effect")
  
  # out -----
  
  out <- mget(ls(envir = environment()))
  
  return(out)
}



