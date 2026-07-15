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


sample_binom <- function(n, true_mean) {
  X <- rbinom(n, 1, true_mean)
}
