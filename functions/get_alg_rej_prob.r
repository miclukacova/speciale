get_alg_rej_prob <- function(alpha = 0.05, gammas = seq(0.6,0.9, by = 0.05),
                             tau = 0.5, B = 10^4, N = 3000, n = 500, loss, alg) {

  set.seed(6389)

  rej_prob <- c()
  i <- 1
  test_res <- c()

  for(gamma in gammas){
    k_star <- -1
    a_star <- -1
    while(a_star < 0 | a_star > 1){
      k_star <- k_star + 1
      a_star <- (alpha - pbinom(k_star -1, floor(N/(n+1)), prob = tau))/dbinom(k_star, floor(N/(n+1)), prob = tau)
    }

    for(b in 1:B){
      # Data
      #-------------------------------------------------------------------------------
      X1 <- rnorm(N,-10,2)
      X2 <- rnorm(N,-1,1)
      X <- cbind(X1,X2)
      Y <- X %*% beta + rnorm(N, sd = gamma)

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
    print(i)
  }

  return(data.frame(gamma = gammas, rej_prob = rej_prob))
}
