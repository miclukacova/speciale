# The loss function
loss <- function(y_hat, y) abs(y-y_hat) > 0.5
# The algorithm
lm_alg <- function(X, Y) {
  solve(t(X) %*% X) %*% t(X) %*% Y
}
# Our particular choice of betas
beta <- c(0.3,0.9)

