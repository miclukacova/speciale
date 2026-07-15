get_OLS_risk <- function(tau = 0.5, B = 10^4, n = 500) {

  set.seed(53794)

  # Definition of vectors for the collection of results
  riskss <- c()                                                                # This is a vector for every single calculation of the conditional probability
  ms <- c(seq(0.04, 1,by =0.02), seq(1, 10, by = 1), seq(10,100, by = 10))      # This is the sequence standard deviations
  risk_est <- c()

  for(j in seq_along(ms)){
    for(i in 1:B){
      X1 <- rnorm(n + 1, -10, 2)
      X2 <- rnorm(n + 1, -1, 1)
      X <- cbind(X1, X2)

      sd_X_exp <- sqrt((t(X[n + 1,]) %*% solve(t(X[1:n,]) %*% X[1:n,]) %*% X[n + 1,] + 1)) * ms[j]

      riskss[i] <- 2 * pnorm(-tau, 0, sd_X_exp)
    }
    risk_est[j] <- mean(riskss)
  }

  # We find gamma^*
  m_star <- ms[min(which(risk_est >= tau))]

  return(list(m_star = m_star, risk_est = data.frame(ms = ms, risk = risk_est)))
}
