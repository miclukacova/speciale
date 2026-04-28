#-------------------------------------------------------------------------------
### Rina Barber paper - l_q stability
#-------------------------------------------------------------------------------
# The idea is to try and estimate the l_q stability of the OLS algorithm

get_l_q_stability_data <- function(beta = c(0.3,0.9,-2), B = 500, q = 2) {
  set.seed(6389)

  ns <- c(50, 60, 70, 80, 90, seq(100, 500, by = 50))
  i <- 1
  beta_q <- c()

  for(N in ns){
    loss_diff <- matrix(nrow=N, ncol = B)
    for(b in 1:B){
      X1 <- rnorm(N,-10,2)
      X2 <- rnorm(N,-1,1)
      X3 <- rbinom(N, 1, 0.4)
      X <- cbind(X1,X2,X3)
      Y <- X %*% beta + rnorm(N)

      X_n1 <- cbind(rnorm(1,-10,2),rnorm(1,-1,1),rbinom(1, 1, 0.4))
      Y_n1 <- X_n1 %*% beta + rnorm(1)

      Y_hat_n <- X_n1 %*% lm_alg(X,Y)

      for(j in 1:N){
        Y_hat_j <- X_n1 %*% lm_alg(X[-j,],Y[-j])
        loss_diff[j,b] <- abs(loss(Y_hat_n, Y_n1) - loss(Y_hat_j, Y_n1))^q
      }
    }
    beta_q[i] <- mean(loss_diff)^(1/q)
    i <- i +1
    print(i)
  }

  return(data.frame(n = ns, beta_q = beta_q))
}
