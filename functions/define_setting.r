# The loss function
loss_bin <- function(y_hat, y) abs(y-y_hat) > 0.5
loss_sq_error <- function(y_hat, y) (y_hat - y)^2

# The algorithm
lm_alg <- function(X, Y) {
  beta_hat <- solve(t(X) %*% X) %*% t(X) %*% Y
  f_hat <- function(X_n1) X_n1 %*% beta_hat
  return(f_hat)
}
# Our particular choice of betas
beta <- c(0.3,0.9)

# Binomial data sampling function
sample_binom <- function(n, true_mean) {
  X <- rbinom(n, 1, true_mean)
}

# Normal data sampling function
sample_patient_norm <- function(N, m) {
  X <- MASS::mvrnorm(n = N, mu = m, Sigma = Sigma)
  X
}

# Normal data sampling function
sample_patient_norm <- function(N, m, Sigma) {
  X <- MASS::mvrnorm(n = N, mu = m, Sigma = Sigma)
  X
}

# Contaminated normal data sampling function
sample_patient_contaminated <- function(
    N, m, Sigma,
    eps = 0.20,
    inflation = 25
) {

  X <- MASS::mvrnorm(N, m, Sigma)

  ind <- sample.int(N, ceiling(eps * N))

  X[ind, ] <- MASS::mvrnorm(
    length(ind),
    m,
    inflation * Sigma
  )

  X
}

# log normal sampler
sample_patient_lognormal <- function(N, m, Sigma) {

  p <- length(m)

  Z <- MASS::mvrnorm(
    n = N,
    mu = rep(0, p),
    Sigma = Sigma
  )

  Y <- exp(Z)

  lognormal_mean <- exp(0.5 * diag(Sigma))

  X <- sweep(Y, 2, lognormal_mean, "-")
  X <- sweep(X, 2, m, "+")

  X
}
