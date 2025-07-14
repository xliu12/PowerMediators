
devtools::load_all("PowerMediators")


# simulations ------------------
condition_all <- data.frame(
  expand.grid(
    n = c(100, 300, 500),
    m1_on_a = c(0.36),
    m2_on_a = c(0.36),
    y_on_m1 = c(0),
    y_on_m2 = c(0),
    y_on_am1 = c(0),
    y_on_am2 = c(0.36),
    y_on_m1m2 = c(0, 0.36),
    y_on_am1m2 = c(0, 0.36),
    em_corr = c(0.1),
    M1_binary = c(FALSE, TRUE),
    M2_binary = c(FALSE, TRUE),
    Y_binary = c(FALSE)
  )) %>% 
  # for TypeI binary Y
  bind_rows(data.frame(
    expand.grid(
      n = c(100, 300, 500),
      m1_on_a = c(0.4),
      m2_on_a = c(0.4),
      y_on_m1 = c(0),
      y_on_m2 = c(0),
      y_on_am1 = c(0),
      y_on_am2 = c(0.1),
      y_on_m1m2 = c(0, 0.59),
      y_on_am1m2 = c(0, 0.59),
      em_corr = c(0.1),
      M1_binary = c(FALSE, TRUE),
      M2_binary = c(FALSE, TRUE),
      Y_binary = c(TRUE)
    ))) %>% 
  # for power
  bind_rows(data.frame(expand.grid(
    n = c(100, 300, 500),
    m1_on_a = c(0.14, 0.27, 0.36),
    m2_on_a = c(0.27),
    y_on_m1 = c(0.14, 0.27, 0.36),
    y_on_m2 = c(0.27),
    y_on_am1 = c(0.1),
    y_on_am2 = c(0.1),
    y_on_m1m2 = c(0.1),
    y_on_am1m2 = c(0.1),
    em_corr = c(0.1),
    M1_binary = c(FALSE, TRUE),
    M2_binary = c(FALSE, TRUE),
    Y_binary = c(FALSE)
  ))) %>% 
  # for power
  bind_rows(data.frame(expand.grid(
    n = c(100, 300, 500),
    m1_on_a = c(0.14, 0.27, 0.36),
    m2_on_a = c(0.27),
    y_on_m1 = c(0.14, 0.27, 0.36),
    y_on_m2 = c(0.27),
    y_on_am1 = c(0.27),
    y_on_am2 = c(0.27),
    y_on_m1m2 = c(-0.27),
    y_on_am1m2 = c(-0.27),
    em_corr = c(0.1),
    M1_binary = c(FALSE, TRUE),
    M2_binary = c(FALSE, TRUE),
    Y_binary = c(TRUE)
  )))

condition <- condition_all %>%
  filter(
    # y_on_m1 == 0, ((y_on_m1m2 == 0 & y_on_am1m2 > 0) | (y_on_m1m2 > 0 & y_on_am1m2 == 0) | (y_on_m1m2 == 0 & y_on_am1m2 == 0)),
    y_on_m1 >= 0.14, y_on_m1 == m1_on_a,
    # y_on_m1m2 == 0 & y_on_am1m2 == 0,
    # y_on_m1m2 == 0 & y_on_am1m2 > 0,
    # y_on_m1m2 > 0, y_on_am1m2 == 0,
    # em_corr %in% c(0, 0.3),
    n %in% c(100),
    Y_binary %in% c(TRUE),
    ((!M1_binary) & (!M2_binary)) |
      ((M1_binary) & (!M2_binary)) |
      (M1_binary) & (M2_binary)
  )

iseed <- 119130 # datseeds[1]
cond <- 1

OneData <- function(iseed = 123, cond = 1){
  
  if (condition$Y_binary[cond]) {
    R2.mx = 0.1
    R2.yx = 0.1
  } else {
    R2.mx = 0.3
    R2.yx = 0.3
  }
  
  gen_data <- gen.data(
    iseed = iseed,
    n = condition$n[cond],
    R2.ax = 0, R2.mx = R2.mx, R2.yx = R2.yx, 
    std.m_on_a = c(condition$m1_on_a[cond], condition$m2_on_a[cond]),
    std.y_on_a = 0.14,
    std.y_on_m = c(condition$y_on_m1[cond], condition$y_on_m2[cond]),
    std.y_on_am_2way = c(condition$y_on_am1[cond], condition$y_on_am2[cond]),
    std.y_on_m_2way = condition$y_on_m1m2[cond],
    std.y_on_am_3way = condition$y_on_am1m2[cond],
    em_corr = condition$em_corr[cond],
    M_binary = c(condition$M1_binary[cond], condition$M2_binary[cond]),
    Y_binary = condition$Y_binary[cond]
  )
  
  
  dat <- gen_data$dat
  
  true_vals <- gen_data$true_vals
  
  res <- med2.reg(
    data = dat, 
    M_binary = gen_data$M_binary, 
    Y_binary = gen_data$Y_binary, 
    sig.IIE = 0.05,
    sig.adjust = c("no_adjust", "bonferroni", "modified_bon1", "modified_bon2"),
    n.draws = 1000
  )
  
  res1 <- data.frame(condition[cond, ],
                     full_join(res, true_vals, by = "effect"), 
                     row.names = NULL)
  
  
  return(res1)
}
