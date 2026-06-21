# The universal inference process
UIE <- function(X, f0_estimator, f1_estimator, N, burnin = NULL) {
  # Storage
  Es <- vector(length = N)
  f1s <- vector(length = N)

  if(!is.null(burnin)){
    # Compute the e-values
    for(i in burnin:N){
      f1s[i - (burnin - 1)] <- f1_estimator(X = X[1:i,], i = i)
      f0s <- f0_estimator(X[1:i,], i)[burnin:i]

      if(any(f0s == 0)) stop("UIE(f0) is 0")

      Es[i - (burnin - 1)] <- prod(f1s[1:(i - (burnin - 1))] / f0s)
    }
    # Adapt sizes
    f1s <- f1s[1:(burnin - 1)]
    Es <- c(rep(1, (burnin - 1)), Es[1:(N - (burnin - 1))])

  } else {
    # Compute the e-values
    for(i in 1:N){
      f1s[i] <- f1_estimator(X = X[1:i,], i = i)
      f0s <- f0_estimator(X[1:i,], i)

      if(any(f0s == 0)) stop("UIE(f0) is 0")

      Es[i] <- prod(f1s[1:i] / f0s)
    }
  }
  if(any(f1s == 0)) stop("UIE(f1) is 0")
  Es
}

UIE_test <- function(X, f0_estimator, f1_estimator, N, alpha, burnin = FALSE) {
  UIE_res <- UIE(X, f0_estimator, f1_estimator, N, burnin = burnin)
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


