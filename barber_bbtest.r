#-------------------------------------------------------------------------------
### Rina Barber paper - LM example
#-------------------------------------------------------------------------------
set.seed(7357)
B <- 1000
tau <- 0.99

# Estimate Risk
loss <- function(y_hat, y) (y_hat - y)^2

risk_est <- function(f_hat, X, Y) {
  Y_hat <- X %*% f_hat
  1/(N/K) * sum(loss(Y_hat, Y))
}

lm_alg <- function(X, Y) {
  solve(t(X) %*% X) %*% t(X) %*% Y
}

pow_lm <- c()
for(b in 1:B){
  # Data
  #-------------------------------------------------------------------------------
  N <- 300
  X1 <- rnorm(N,-10,2)
  X2 <- rnorm(N,-1,1)
  X3 <- rbinom(N, 1, 0.4)
  X <- cbind(X1,X2,X3)
  beta <- c(0.3,0.9,-2)
  Y <- X %*% beta + rnorm(N)

  # Black box test
  #-------------------------------------------------------------------------------

  # Fit models
  K <- 5
  n <- N * (1 - 1 / K)
  parti <- sample(1:N, N)

  hat_f_k <- list()
  for(k in 1:K){
    indices <- ((k -1) * N/K + 1):(k * N/K)
    hat_f_k[[k]] <- lm_alg(X[-indices,], Y[-indices])
  }

  risk_CV <- list()

  for(k in 1:K){
    risk_CV[[k]] <- risk_est(hat_f_k[[k]], X[((k -1) * N/K + 1):(k * N/K),], Y[((k -1) * N/K + 1):(k * N/K)])
  }

  test[b] <- mean(unlist(risk_CV)) < tau
}

mean(test)
