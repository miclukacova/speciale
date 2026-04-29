get_OLS_risk <- function(tau = 0.5, B = 10^4, N = 1000, n = 500) {

  set.seed(53794)

  # Definition of vectors for the collection of results
  riskss <- c()                                                                # This is a vector for every single calculation of the conditional probability
  vars <- c(seq(0.05,1, by =0.01), seq(1,10, by = 1), seq(10,100, by = 10))    # This is the sequence standard deviations
  risk_est <- c()

  for(j in seq_along(vars)){
    for(i in 1:B){
      X1 <- rnorm(N+1,-10,2)
      X2 <- rnorm(N+1,-1,1)
      X <- cbind(X1[1:N],X2[1:N])

      riskss[i] <- 2* pnorm(-0.5, 0, t(X[n+1,]) %*% solve(t(X) %*% X) %*% X[n+1,] + vars[j])
    }
    risk_est[j] <- mean(riskss)
  }

  # We find gamma^*
  gamma_star <- vars[min(which(risk_est >= tau))]

  return(list(gamma_star = gamma_star, risk_est = data.frame(gamma = vars, risk = risk_est)))
}
