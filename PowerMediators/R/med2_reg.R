

estimate_med2.reg <- function(
    data = dat
) {
  Xnames <- grep("X", names(data), value = TRUE)
  
  m1_formula <- paste("M1 ~ A +", paste(Xnames, collapse = " + "))
  m1_fit <- summary(lm(m1_formula, data = data))
  
  m2_formula <- paste("M2 ~ A +", paste(Xnames, collapse = " + "))
  m2_fit <- summary(lm(m2_formula, data = data))
  
  y_formula <- paste("Y ~ A + M1 + M2 + A:M1 + A:M2 + M1:M2 + A:M1:M2 +", paste(Xnames, collapse = " + "))
  y_fit <- summary(lm(y_formula, data = data))
  
  
  y_coefs <- y_fit$coefficients[, 1]
  m1_coefs <- m1_fit$coefficients[, 1]
  m2_coefs <- m2_fit$coefficients[, 1]
  
  # IIE_M1 --------
  a0a2 <- expand.grid(a0=c(0,1), a2=c(0,1))
  est_IIE_M1 <- map_dbl(1:nrow(a0a2), \(i) {
    cal.IIE_M1(a0 = a0a2$a0[i], a2 = a0a2$a2[i], y_coefs, m1_coefs, m2_coefs)
  })
  names(est_IIE_M1) <- glue("IIE_M1({a0a2$a0},,{a0a2$a2})")
  
  # IIE_M2 --------
  a0a1 <- expand.grid(a0=c(0,1), a1=c(0,1))
  est_IIE_M2 <- map_dbl(1:nrow(a0a1), \(i) {
    cal.IIE_M2(a0 = a0a1$a0[i], a1 = a0a1$a1[i], y_coefs, m1_coefs, m2_coefs)
  })
  names(est_IIE_M2) <- glue("IIE_M2({a0a1$a0},{a0a1$a1},)")
  
  c(est_IIE_M1, est_IIE_M2)
}


med2.reg <- function(
    data = dat,
    sig.IIE = 0.05,
    sig.adjust = c("no_adjust", "bonferroni", "modified_bon1", "modified_bon2"),
    nboot = 1000
) {
  est_IIE <- estimate_med2.reg(data = data)
  
  
  boot_est <- replicate(nboot, {
    df_boot <- data[sample(1:nrow(data), replace = TRUE), ]
    est_IIE <- estimate_med2.reg(data = df_boot)
  })
  
  # average absolute correlation r_IIE_M1 ----
  
  boot_est_IIE_M1 <- t(boot_est[grep("IIE_M1", rownames(boot_est)), ])
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
    ci_low <- quantile(boot_est_IIE_M1[, i], probs = sig.PT_IIE_M1[[i]]/2)
    ci_up <- quantile(boot_est_IIE_M1[, i], probs = 1-sig.PT_IIE_M1[[i]]/2)
    
    if_sig <- (ci_low*ci_up>0)
    
    data.frame(effect = colnames(boot_est_IIE_M1)[i], ci_low = ci_low, ci_up = ci_up, if_sig = if_sig, sig.adjust = sig.adjust, row.names = NULL)
    }) 
  
  test_IIE_M1 <- test_IIE_M1 %>% 
    group_by(sig.adjust) %>% 
    mutate(any_sig = mean(if_sig)>0 )
  
  
  # average absolute correlation r_IIE_M2 ----
  
  boot_est_IIE_M2 <- t(boot_est[grep("IIE_M2", rownames(boot_est)), ])
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
    ci_low <- quantile(boot_est_IIE_M2[, i], probs = sig.PT_IIE_M2[[i]]/2)
    ci_up <- quantile(boot_est_IIE_M2[, i], probs = 1-sig.PT_IIE_M2[[i]]/2)
    
    if_sig <- (ci_low*ci_up>0)
    
    data.frame(effect = colnames(boot_est_IIE_M2)[i], ci_low = ci_low, ci_up = ci_up, if_sig = if_sig, sig.adjust = sig.adjust, row.names = NULL)
  }) 
  
  test_IIE_M2 <- test_IIE_M2 %>% 
    group_by(sig.adjust) %>% 
    mutate(any_sig = mean(if_sig)>0 )
  
  
  test_IIE <- rbind(test_IIE_M1, test_IIE_M2)
  
  res <- full_join(
    data.frame(effect = names(est_IIE), est = est_IIE), test_IIE, by = "effect")
  
  res
  
}

cal.IIE_M1 <- function(a0=1, a2=1, y_coefs, m1_coefs, m2_coefs, 
                       x_means = rep(0, length(grep("X", names(y_coefs), value = TRUE)))
                       ) {
  meanM2_a2 <- m2_coefs["(Intercept)"] + m2_coefs["A"] * a2 + sum(m2_coefs[grep("X", names(m2_coefs), value = TRUE)] * x_means)
  
  est.IIE_M1 <- m1_coefs["A"] * (y_coefs["M1"] + y_coefs["A:M1"] * a0 + y_coefs["M1:M2"] * meanM2_a2+ y_coefs["A:M1:M2"] * a0 * meanM2_a2)
  
  est.IIE_M1
}

cal.IIE_M2 <- function(a0=1, a1=0, y_coefs, m1_coefs, m2_coefs, 
                       x_means = rep(0, length(grep("X", names(y_coefs), value = TRUE)))
) {
  meanM1_a1 <- m1_coefs["(Intercept)"] + m1_coefs["A"] * a1 + sum(m1_coefs[grep("X", names(m1_coefs), value = TRUE)] * x_means)
  
  est.IIE_M2 <- m2_coefs["A"] * (y_coefs["M2"] + y_coefs["A:M2"] * a0 + y_coefs["M1:M2"] * meanM1_a1 + y_coefs["A:M1:M2"] * a0 * meanM1_a1)
  
  est.IIE_M2
}
  