get_cv_test_data <- function(taus = seq(0.85,1.2, by = 0.005), B = 10^4, N = 500, K = 5, loss, alg = lm_alg) {
  #-------------------------------------------------------------------------------
  ### Rina Barber paper - LM example
  #-------------------------------------------------------------------------------
  # The idea is to employ the black box algorithm test, and check the type 1 error (when tau <= 1)
  # and the power (when tau > 1) of the test

  set.seed(6389)
  rej_prob <- c()
  i <- 1
  test_res <- c()

  # Estimator of risk
  risk_est <- function(f_hat, X, Y) {
    Y_hat <- f_hat(X)
    1/(N/K) * sum(loss(Y_hat, Y))
    mean(loss(Y_hat, Y))
  }

  for(tau in taus){
    for(b in 1:B){
      # Data
      #-------------------------------------------------------------------------------
      X1 <- rnorm(N,-10,2)
      X2 <- rnorm(N,-1,1)
      X <- cbind(X1,X2)
      Y <- X %*% beta + rnorm(N)

      # Black box test
      #-------------------------------------------------------------------------------

      # Fit models
      n <- N * (1 - 1 / K)
      parti <- sample(1:N, N)

      hat_f_k <- list()
      for(k in 1:K){
        indices <- ((k -1) * N/K + 1):(k * N/K)
        hat_f_k[[k]] <- alg(X[-indices,], Y[-indices])
      }

      risk_CV <- c()

      for(k in 1:K){
        risk_CV[k] <- risk_est(hat_f_k[[k]], X[((k -1) * N/K + 1):(k * N/K),], Y[((k -1) * N/K + 1):(k * N/K)])
      }

      test_res[b] <- mean(risk_CV) < tau
    }

    rej_prob[i] <- mean(test_res)
    i <- i+1
    print(i)
  }

  return(data.frame(taus = taus, rej_prob = rej_prob))
}
