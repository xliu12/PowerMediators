

estimate_med2.reg <- function(
    data = dat,
    M_binary,
    Y_binary,
    n.draws = 1000 
) {
  Xnames <- grep("^X", names(data), value = TRUE)
  
  # fit models ------------
  m1_formula <- as.formula(paste("M1 ~ A +", paste(Xnames, collapse = " + ")))
  if (M_binary[1]) {
    m1_fit <- glm(m1_formula, data = data, family = binomial(link = "probit"))
  } else {
    m1_fit <- lm(m1_formula, data = data)
  }
  
  m2_formula <- as.formula(paste("M2 ~ A +", paste(Xnames, collapse = " + ")))
  if (M_binary[2]) {
    m2_fit <- glm(m2_formula, data = data, family = binomial(link = "probit"))
  } else {
    m2_fit <- lm(m2_formula, data = data)
  }
  
  y_formula <- as.formula(paste("Y ~ A + M1 + M2 + A:M1 + A:M2 + M1:M2 + A:M1:M2 +", paste(Xnames, collapse = " + ")))
  if (Y_binary) {
    y_fit <- glm(y_formula, data = data, family = binomial(link = "probit"))
  } else {
    y_fit <- lm(y_formula, data = data)
  }
  
  # coefficient estimates
  y_coefs <- as.data.frame(t(coef(y_fit)))
  m1_coefs <- as.data.frame(t(coef(m1_fit)))
  m2_coefs <- as.data.frame(t(coef(m2_fit)))
  
  if (!is.null(n.draws)) {
    # draws of coefficient estimates for imputation-based estimator and MCCI
    y_coefs <- mvtnorm::rmvnorm(n.draws, mean = coef(y_fit), sigma = vcov(y_fit)) %>% as.data.frame()
    m1_coefs <- mvtnorm::rmvnorm(n.draws, mean = coef(m1_fit), sigma = vcov(m1_fit)) %>% as.data.frame()
    m2_coefs <- mvtnorm::rmvnorm(n.draws, mean = coef(m2_fit), sigma = vcov(m2_fit)) %>% as.data.frame()
  }
  
  
  if (sum(M_binary, Y_binary) == 0) {
    # gaussian Y, M -------------------
    # IIE_M1
    a0a2 <- expand.grid(a0=c(0,1), a2=c(0,1))
    est_IIE_M1 <- map(1:nrow(a0a2), \(i) {
      cal.IIE_M1(a0 = a0a2$a0[i], a2 = a0a2$a2[i], y_coefs, m1_coefs, m2_coefs)
    }) %>% reduce(bind_cols)
    
    # IIE_M2
    a0a1 <- expand.grid(a0=c(0,1), a1=c(0,1))
    est_IIE_M2 <- map(1:nrow(a0a1), \(i) {
      cal.IIE_M2(a0 = a0a1$a0[i], a1 = a0a1$a1[i], y_coefs, m1_coefs, m2_coefs)
    }) %>% reduce(bind_cols)
    
    IIEs <- cbind(est_IIE_M1, est_IIE_M2)
    
  } else {
    # binary Y or M -----------
    IIEs <- purrr::map(1:nrow(y_coefs), \(i=1) {
      y_fit$coefficients <- y_coefs[i, ] %>% unlist
      m1_fit$coefficients <- m1_coefs[i, ] %>% unlist
      m2_fit$coefficients <- m2_coefs[i, ] %>% unlist
      
      IIEs_tmp <- bYM_cal.IIE(y_coefs[i, ], m1_coefs[i, ], m2_coefs[i, ],
                              y_fit, m1_fit, m2_fit, data,
                              M_binary,
                              Y_binary)
    }) %>% purrr::reduce(bind_rows)
    
  }
  
  IIEs
}


med2.reg <- function(
    data = dat, M_binary, Y_binary,
    sig.IIE = 0.05,
    sig.adjust = c("no_adjust", "bonferroni", "modified_bon1", "modified_bon2"),
    nboot = NULL, # for bootstrap CI
    n.draws = 1000 # for Monte Carlo CI
) {
  est_IIE <- estimate_med2.reg(data, M_binary, Y_binary, n.draws = NULL)
  if (!is.null(nrow(est_IIE))) {
    est_IIE <- est_IIE %>% unlist
  }
  
  if (!is.null(nboot)) {
    boot_est <- t(replicate(nboot, {
      df_boot <- data[sample(1:nrow(data), replace = TRUE), ]
      est_IIE <- estimate_med2.reg(data = df_boot, M_binary, Y_binary, n.draws = NULL)
    })) %>% as.data.frame()
  } else {
    boot_est <- estimate_med2.reg(data, M_binary, Y_binary, n.draws)
  }
  
  
  # average absolute correlation r_IIE_M1 ----
  
  boot_est_IIE_M1 <- boot_est[, grep("IIE_M1", names(boot_est))]
  k <- ncol(boot_est_IIE_M1) # number of tests
  
  r_IIE_M1 <- map_dbl(1:k, \(i=1) {
    mean(abs(cor(boot_est_IIE_M1)[i, -i]))
  })
  names(r_IIE_M1) <- colnames(boot_est_IIE_M1)
  # modified_bon1
  # \alpha_\PT = \frac{\alpha_\FW}{k -(k-1)\sqrt{|\bar r_j|}}
  
  # modified_bon2
  # \alpha_\PT = \frac{\alpha_\FW}{k^{1-\sqrt{|\bar r_j|}}}
  
  
  sig.PT_IIE_M1 <- map(1:k, \(i=1) {
    sig.PT <- case_when(sig.adjust == "no_adjust" ~ sig.IIE,
                        sig.adjust == "bonferroni" ~ sig.IIE/k,
                        sig.adjust == "modified_bon1" ~ sig.IIE/(k-(k-1)*sqrt(r_IIE_M1[i])),
                        sig.adjust == "modified_bon2" ~ sig.IIE/(k^(1-sqrt(r_IIE_M1[i])))
    )
    
    names(sig.PT) <- sig.adjust
    sig.PT
  })
  
  names(sig.PT_IIE_M1) <- names(r_IIE_M1)
  
  test_IIE_M1 <- map_df(1:k, \(i=1) {
    ci_low <- quantile(boot_est_IIE_M1[[i]], probs = sig.PT_IIE_M1[[i]]/2)
    ci_up <- quantile(boot_est_IIE_M1[[i]], probs = 1-sig.PT_IIE_M1[[i]]/2)
    SE <- sd(boot_est_IIE_M1[[i]])
    
    if_sig <- (ci_low*ci_up>0)
    
    data.frame(effect = colnames(boot_est_IIE_M1)[i], SE = SE, ci_low = ci_low, ci_up = ci_up, if_sig = if_sig, sig.adjust = sig.adjust, row.names = NULL)
  })
  
  test_IIE_M1 <- test_IIE_M1 %>%
    group_by(sig.adjust) %>%
    mutate(any_sig = mean(if_sig)>0 )
  
  
  # average absolute correlation r_IIE_M2 ----
  
  boot_est_IIE_M2 <- boot_est[, grep("IIE_M2", names(boot_est))]
  k <- ncol(boot_est_IIE_M2) # number of tests
  
  r_IIE_M2 <- map_dbl(1:k, \(i=1) {
    mean(abs(cor(boot_est_IIE_M2)[i, -i]))
  })
  names(r_IIE_M2) <- colnames(boot_est_IIE_M2)
  # modified_bon1
  # \alpha_\PT = \frac{\alpha_\FW}{k -(k-1)\sqrt{|\bar r_j|}}
  
  # modified_bon2
  # \alpha_\PT = \frac{\alpha_\FW}{k^{1-\sqrt{|\bar r_j|}}}
  
  
  sig.PT_IIE_M2 <- map(1:k, \(i=1) {
    sig.PT <- case_when(sig.adjust == "no_adjust" ~ sig.IIE,
                        sig.adjust == "bonferroni" ~ sig.IIE/k,
                        sig.adjust == "modified_bon1" ~ sig.IIE/(k-(k-1)*sqrt(r_IIE_M2[i])),
                        sig.adjust == "modified_bon2" ~ sig.IIE/(k^(1-sqrt(r_IIE_M2[i])))
    )
    
    names(sig.PT) <- sig.adjust
    sig.PT
  })
  
  names(sig.PT_IIE_M2) <- names(r_IIE_M2)
  
  test_IIE_M2 <- map_df(1:k, \(i=1) {
    ci_low <- quantile(boot_est_IIE_M2[[i]], probs = sig.PT_IIE_M2[[i]]/2)
    ci_up <- quantile(boot_est_IIE_M2[[i]], probs = 1-sig.PT_IIE_M2[[i]]/2)
    SE <- sd(boot_est_IIE_M2[[i]])
    
    if_sig <- (ci_low*ci_up>0)
    
    data.frame(effect = colnames(boot_est_IIE_M2)[i], SE = SE, ci_low = ci_low, ci_up = ci_up, if_sig = if_sig, sig.adjust = sig.adjust, row.names = NULL)
  })
  
  test_IIE_M2 <- test_IIE_M2 %>%
    group_by(sig.adjust) %>%
    mutate(any_sig = mean(if_sig)>0 )
  
  
  test_IIE <- rbind(test_IIE_M1, test_IIE_M2)
  
  res <- full_join(
    data.frame(effect = names(est_IIE), est = est_IIE), test_IIE, by = "effect")
  
  res
  
}

# gaussian Y,M ---------------------

cal.IIE_M1 <- function(a0=1, a2=1, y_coefs, m1_coefs, m2_coefs,
                       x_means = rep(0, length(grep("X", names(y_coefs), value = TRUE)))
) {
  names(x_means) <- grep("X", names(m2_coefs), value = TRUE)
  
  meanM2_a2 <- m2_coefs[["(Intercept)"]] + m2_coefs[["A"]] * a2 + rowSums(m2_coefs[grep("X", names(m2_coefs), value = TRUE)] * t(replicate(nrow(m2_coefs), x_means)))
  
  est.IIE_M1 <- m1_coefs[["A"]] * (y_coefs[["M1"]] + y_coefs[["A:M1"]] * a0 + y_coefs[["M1:M2"]] * meanM2_a2+ y_coefs[["A:M1:M2"]] * a0 * meanM2_a2)
  
  est.IIE_M1 <- data.frame(est.IIE_M1)
  names(est.IIE_M1) <- glue("IIE_M1({a0},,{a2})")
  est.IIE_M1
}



cal.IIE_M2 <- function(a0=1, a1=0, y_coefs, m1_coefs, m2_coefs,
                       x_means = rep(0, length(grep("X", names(y_coefs), value = TRUE)))
) {
  meanM1_a1 <- m1_coefs[["(Intercept)"]] + m1_coefs[["A"]] * a1 + rowSums(m1_coefs[grep("X", names(m1_coefs), value = TRUE)] * t(replicate(nrow(m1_coefs), x_means)))
  
  est.IIE_M2 <- m2_coefs[["A"]] * (y_coefs[["M2"]] + y_coefs[["A:M2"]] * a0 + y_coefs[["M1:M2"]] * meanM1_a1 + y_coefs[["A:M1:M2"]] * a0 * meanM1_a1)
  
  est.IIE_M2 <- data.frame(est.IIE_M2)
  names(est.IIE_M2) <- glue("IIE_M2({a0},{a1},)")
  est.IIE_M2
}

# outcome regression-imputation -----
EY.a0a1a2 <- function(a0=1, a1=0, a2=1,
                      y_coefs, m1_coefs, m2_coefs,
                      y_fit, m1_fit, m2_fit, data,
                      M_binary = c(FALSE, FALSE),
                      Y_binary = TRUE) {
  
  if (Y_binary==TRUE & (M_binary[1]==TRUE) & (M_binary[2]==TRUE)) {
    pM1a1 <- predict(m1_fit, newdata = mutate(data, A=a1), type = "response")
    pM2a2 <- predict(m2_fit, newdata = mutate(data, A=a2), type = "response")
    
    EYx <- map2(.x = c(0,0,1,1), .y = c(0,1,0,1), .f = \(.x, .y) {
      predict(y_fit, newdata = mutate(data, A=a0, M1=.x, M2=.y), type = "response") * (pM1a1 * .x + (1-pM1a1) * (1-.x)) * (pM2a2 * .y + (1-pM2a2) * (1-.y))
    }) %>% reduce(.f = `+`)
    
    EYa0a1a2 <- mean(EYx)
  }
  
  if (Y_binary==TRUE & (M_binary[1]==TRUE) & (M_binary[2]==FALSE)) {
    EYa0a1a2 <- E.YbM1bM2g(a0, a1, a2, y_coefs, m1_coefs, m2_coefs, data)
  }
  
  if (Y_binary==TRUE & (M_binary[1]==FALSE) & (M_binary[2]==FALSE)) {
    EYa0a1a2 <- calE.YbM1gM2g(a0, a1, a2, y_coefs, m1_coefs, m2_coefs, data)
  }
  
  if (Y_binary==FALSE) {
    M1a1 <- predict(m1_fit, newdata = mutate(data, A=a1), type = "response")
    M2a2 <- predict(m2_fit, newdata = mutate(data, A=a2), type = "response")
    EYx <- predict(y_fit, newdata = mutate(data, A=a0, M1=M1a1, M2=M2a2), type = "response")
    EYa0a1a2 <- mean(EYx)
  }
  
  EYa0a1a2
}

bYM_cal.IIE <- function(y_coefs, m1_coefs, m2_coefs, 
                        y_fit, m1_fit, m2_fit, data,
                        M_binary = c(FALSE, FALSE),
                        Y_binary = TRUE) {
  
  a_vals <- expand.grid(a0=c(0,1), a1=c(0,1), a2=c(0,1))
  
  EYa0a1a2 <- sapply(1:nrow(a_vals), \(i) {
    EY.a0a1a2(a0 = a_vals$a0[i], a1 = a_vals$a1[i], a2 = a_vals$a2[i],
              y_coefs, m1_coefs, m2_coefs,
              y_fit, m1_fit, m2_fit, data,
              M_binary,
              Y_binary)
  })
  names(EYa0a1a2) <- glue("Y({a_vals$a0},M1({a_vals$a1}),M2({a_vals$a2}))")
  
  # contrast outcomes for IIEs
  IIE_M1 <- map2(.x=c(0,0,1,1), .y=c(0,1,0,1), .f = \(.x, .y) {
    IIE_M1_tmp <- EYa0a1a2[[glue("Y({.x},M1(1),M2({.y}))")]] - EYa0a1a2[[glue("Y({.x},M1(0),M2({.y}))")]]
    names(IIE_M1_tmp) <- glue("IIE_M1({.x},,{.y})")
    IIE_M1_tmp
  }) %>% unlist
  
  IIE_M2 <- map2(.x=c(0,0,1,1), .y=c(0,1,0,1), .f = \(.x, .y) {
    IIE_M2_tmp <- EYa0a1a2[[glue("Y({.x},M1({.y}),M2(1))")]] - EYa0a1a2[[glue("Y({.x},M1({.y}),M2(0))")]]
    names(IIE_M2_tmp) <- glue("IIE_M2({.x},{.y},)")
    IIE_M2_tmp
  }) %>% unlist
  
  c(IIE_M1, IIE_M2)
}


# NOT USED -------------
intEY.a0a1a2 <- function(a0=1, a1=0, a2=1,
                         y_fit, m1_fit, m2_fit, data,
                         M_binary = c(FALSE, FALSE),
                         Y_binary = TRUE, nreps = 1000) {
  
  if (Y_binary==TRUE & (M_binary[1]==TRUE) & (M_binary[2]==TRUE)) {
    pM1a1 <- predict(m1_fit, newdata = mutate(data, A=a1), type = "response")
    pM2a2 <- predict(m2_fit, newdata = mutate(data, A=a2), type = "response")
    
    EYx <- map2(.x = c(0,0,1,1), .y = c(0,1,0,1), .f = \(.x, .y) {
      predict(y_fit, newdata = mutate(data, A=a0, M1=.x, M2=.y), type = "response") * (pM1a1 * .x + (1-pM1a1) * (1-.x)) * (pM2a2 * .y + (1-pM2a2) * (1-.y))
    }) %>% reduce(.f = `+`)
    
  }
  
  if (Y_binary==TRUE & (M_binary[1]==TRUE) & (M_binary[2]==FALSE)) {
    pM1a1 <- predict(m1_fit, newdata = mutate(data, A=a1), type = "response")
    
    EYx <- replicate(nreps, {
      
      M2a2 <- predict(m2_fit, newdata = mutate(data, A=a2)) + rnorm(nrow(data), sd = sigma(m2_fit))
      
      EY <- map(.x = c(0,1), .f = \(.x) {
        predict(y_fit, newdata = mutate(data, A=a0, M1=.x, M2=M2a2), type = "response") * (pM1a1 * .x + (1-pM1a1) * (1-.x))
      }) %>% reduce(.f = `+`)
      
      EY
    })
  }
  
  if (Y_binary==TRUE & (M_binary[1]==FALSE) & (M_binary[2]==TRUE)) {
    pM2a2 <- predict(m2_fit, newdata = mutate(data, A=a2), type = "response")
    
    EYx <- replicate(nreps, {
      
      M1a1 <- predict(m1_fit, newdata = mutate(data, A=a1)) + rnorm(nrow(data), sd = sigma(m1_fit))
      
      EY <- map(.x = c(0,1), .f = \(.x) {
        predict(y_fit, newdata = mutate(data, A=a0, M1=M1a1, M2=.x), type = "response") * (pM2a2 * .x + (1-pM2a2) * (1-.x))
      }) %>% reduce(.f = `+`)
      
      EY
    })
  }
  
  if (Y_binary==TRUE & (M_binary[1]==FALSE) & (M_binary[2]==FALSE)) {
    EYx <- replicate(nreps, {
      M1a1 <- predict(m1_fit, newdata = mutate(data, A=a1)) + rnorm(nrow(data), sd = sigma(m1_fit))
      M2a2 <- predict(m2_fit, newdata = mutate(data, A=a2)) + rnorm(nrow(data), sd = sigma(m2_fit))
      EY <- predict(y_fit, newdata = mutate(data, A=a0, M1=M1a1, M2=M2a2), type = "response")
      
      EY
    })
  }
  
  if (Y_binary==FALSE) {
    M1a1 <- predict(m1_fit, newdata = mutate(data, A=a1), type = "response")
    M2a2 <- predict(m2_fit, newdata = mutate(data, A=a2), type = "response")
    EYx <- predict(y_fit, newdata = mutate(data, A=a0, M1=M1a1, M2=M2a2), type = "response")
    EYa0a1a2 <- mean(EYx)
  }
  
  EYa0a1a2
}

