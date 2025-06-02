
# M1,M2 gaussian, Y binary ---------------
E.YbM1gM2g <- function(em2_vec=0, 
                        a0=1, a1=0, a2=1,
                      y_coefs, m1_coefs, m2_coefs, data) {
  n <- nrow(data)
  y_coefs <- unlist(y_coefs)
  m1_coefs <- unlist(m1_coefs)
  m2_coefs <- unlist(m2_coefs)
  # em2 <- rep(e.m2, n)
  X <- grep("^X", names(y_coefs), value = TRUE)
  Xmat <- as.matrix(data[, X, drop = FALSE])
  
  # calculate sigma_em
  hat.em1 <- data$M1 - m1_coefs[[1]] - as.matrix(data[, names(m1_coefs[-1])]) %*% m1_coefs[-1] #apply(data[, names(m1_coefs[-1])], 1, \(x) sum(x*m1_coefs[-1]))
  sigma_em1 <- sd(hat.em1) * sqrt(n-1)  / sqrt(n-length(m1_coefs))
  
  hat.em2 <- data$M2 - m2_coefs[[1]] - as.matrix(data[, names(m2_coefs[-1])]) %*% m2_coefs[-1] 
  
  sigma_em2 <- sd(hat.em2) * sqrt(n-1)  / sqrt(n-length(m2_coefs))
  
  # mean M1a1 given x
  M1a1x <- m1_coefs[[1]] + m1_coefs[["A"]]*a1 + Xmat %*% m1_coefs[X] 
  
  EYem2 <- numeric(length = length(em2_vec))
  for (i in 1:length(em2_vec)) {
    em2 <- em2_vec[i]
    
    M2a2xe <- em2 + m2_coefs[[1]] + m2_coefs[["A"]]*a2 + Xmat %*% m2_coefs[X]  
    # numerator
    g_num <- y_coefs[[1]] + y_coefs[["A"]]*a0 + (y_coefs[["M1"]] + y_coefs[["A:M1"]]*a0)*M1a1x + 
      (y_coefs[["M2"]] + y_coefs[["A:M2"]]*a0)*M2a2xe + (y_coefs[["M1:M2"]] + y_coefs[["A:M1:M2"]]*a0)*M1a1x*M2a2xe + Xmat %*% y_coefs[X]  
    
    g_den <- (y_coefs[["M1"]] + y_coefs[["A:M1"]]*a0 + (y_coefs[["M1:M2"]] + y_coefs[["A:M1:M2"]]*a0)*M2a2xe)^2
    
    EYem2[i] <- dnorm(em2/sigma_em2) * mean( pnorm(g_num / sqrt(g_den * sigma_em1^2 + 1)) )
  }
  
  return(EYem2)
}

calE.YbM1gM2g <- function(a0=1, a1=0, a2=1,
                          y_coefs, m1_coefs, m2_coefs, data) {
  
  value <- calculus::integral(
    E.YbM1gM2g, bounds = list(em2_vec = c(-3.3,3.3)), 
    params = list(a0, a1, a2,
                  y_coefs = y_coefs, m1_coefs = m1_coefs, m2_coefs = m2_coefs, data = data)
    , vectorize = TRUE)
  
  value$value
  
}


# M1 binary, M2 gaussian, Y binary ----------------
# label the binary mediator as M1
# label the continuous mediator as M2
E.YbM1bM2g <- function(a0=1, a1=0, a2=1,
                       y_coefs, m1_coefs, m2_coefs, data) {
  n <- nrow(data)
  y_coefs <- unlist(y_coefs)
  m1_coefs <- unlist(m1_coefs)
  m2_coefs <- unlist(m2_coefs)
  # em2 <- rep(e.m2, n)
  X <- grep("^X", names(y_coefs), value = TRUE)
  Xmat <- as.matrix(data[, X, drop = FALSE])
  # calculate sigma_em
  hat.em2 <- data$M2 - m2_coefs[[1]] - as.matrix(data[, names(m2_coefs[-1])]) %*% m2_coefs[-1] #apply(data[, names(m2_coefs[-1])], 1, \(x) sum(x*m2_coefs[-1]))
  sigma_em2 <- sd(hat.em2) * sqrt(n-1)  / sqrt(n-length(m2_coefs))
  
  # link M1a1 given x
  M1a1x <- m1_coefs[[1]] + m1_coefs[["A"]]*a1 + Xmat %*% m1_coefs[X]
    #apply(data[, X, drop=FALSE], 1, \(x) sum(x*m1_coefs[X]))
  
  
  # mean M2a2 given X
  M2a2x <- m2_coefs[[1]] + m2_coefs[["A"]]*a2 + Xmat %*% m2_coefs[X]
  #apply(data[, X, drop=FALSE], 1, \(x) sum(x*m2_coefs[X]))
  
  # with M1 | a1=1 
  g_num1 <- y_coefs[[1]] + y_coefs[["A"]]*a0 + (y_coefs[["M1"]] + y_coefs[["A:M1"]]*a0)*1 + 
    (y_coefs[["M2"]] + y_coefs[["A:M2"]]*a0)*M2a2x + (y_coefs[["M1:M2"]] + y_coefs[["A:M1:M2"]]*a0)*1*M2a2x + Xmat %*% y_coefs[X]
  #apply(data[, X, drop=FALSE], 1, \(x) sum(x*y_coefs[X]))
  
  g_den1 <- (y_coefs[["M2"]] + y_coefs[["A:M2"]]*a0 + y_coefs[["M1:M2"]] + y_coefs[["A:M1:M2"]]*a0)^2
  
  EY1 <- pnorm(M1a1x) * pnorm(g_num1 / sqrt(g_den1 * sigma_em2^2 + 1))
  
  # with M1 | a1=0 
  g_num0 <- y_coefs[[1]] + y_coefs[["A"]]*a0 + 
    (y_coefs[["M2"]] + y_coefs[["A:M2"]]*a0)*M2a2x + Xmat %*% y_coefs[X]
  #apply(data[, X, drop=FALSE], 1, \(x) sum(x*y_coefs[X]))
  
  g_den0 <- (y_coefs[["M2"]] + y_coefs[["A:M2"]]*a0)^2
  
  EY0 <- (1 - pnorm(M1a1x)) * pnorm(g_num0 / sqrt(g_den0 * sigma_em2^2 + 1))
  
  EY <- mean(EY1 + EY0)
  
  EY
}
