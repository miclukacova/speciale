get_rina_bb_test_data <- function(alpha = 0.05, taus = seq(0.3,1, by = 0.01),
                                  B = 10000, N = 5000, n = 500, loss, alg){
  #-------------------------------------------------------------------------------
  ### Rina Barber paper - Binomial example
  #-------------------------------------------------------------------------------
  # The idea is to employ the binomial black box algorithm test

  set.seed(6389)

  rej_prob <- c()
  i <- 1
  test_res <- c()

  for(tau in taus){
    k_star <- 0
    a_star <- -1
    while(a_star < 0 | a_star > 1){
      a_star <- (alpha - pbinom(k_star -1, floor(N/(n+1)), prob = tau)) / dbinom(k_star, floor(N/(n+1)), prob = tau)
      k_star <- k_star + 1
      if(k_star > 100) stop("No solution")
    }

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
      K <- floor(N/(n+1))
      parti <- sample(1:N, N)

      hat_f_k <- list()
      for(k in 1:K){
        indices <- ((k -1) * (n+1) + 1):(k*(n+1) -1)
        hat_f_k[[k]] <- alg(X[indices,], Y[indices])
      }

      loss_k <- c()

      for(k in 1:K){
        y_hat <- hat_f_k[[k]](X[k*(n+1),])
        loss_k[k] <- loss(y_hat, Y[k*(n+1)])
      }

      S <- sum(loss_k)
      test_res[b] <- (S < k_star) + (S == k_star) * (runif(1) <= a_star)
    }

    rej_prob[i] <- mean(test_res)
    i <- i+1
  }

  return(data.frame(tau = taus, rej_prob = rej_prob))
}
