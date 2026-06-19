# The universal inference process
UIE <- function(X, f0_estimator, f1_estimator, N) {
  # Storage
  Es <- vector(length = N)
  f1s <- vector(length = N)

  # Compute the e-values
  for(i in 1:N){
    f1s[i] <- f1_estimator(X = X[1:i,], i = i)
    Es[i] <- prod(f1s[1:i] / f0_estimator(X[1:i,], i))
  }
  Es
}

UIE_test <- function(X, f0_estimator, f1_estimator, N, alpha) {
  UIE_res <- UIE(X, f0_estimator, f1_estimator, N)
  test_res <- UIE_res > 1/ alpha
  #test_fut <- UIE_res < alpha / 2
  if(any(test_res)){         #| any(test_fut)
    ESS <- which((test_res) == 1)[1]           # + test_fut
  } else {
    ESS <- N
  }
  return(c(Reject = any(test_res[1:ESS]), ESS = ESS))
}


Example <- FALSE
if(Example){
  # Variance matrix
  Sigma <- matrix(c(1,0,0,1), ncol = 2)

  # Data sampling function
  sample_patient <- function(N, m) {
    X <- MASS::mvrnorm(n = N, mu = m, Sigma = Sigma)
    X
  }

  N <- 200
  X <- sample_patient(N, c(0.3, 0.3))
  plot(1:N, UIE(X = X,
      f0_estimator = f0_estimator,
      f1_estimator = f1_estimator,
      N = N))
}


