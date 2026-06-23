# The universal inference process
UIE <- function(X, log_f0, log_f1, N, burnin = NULL) {
  # Storage
  logE <- vector(length = N)
  f1s <- vector(length = N)

  if(!is.null(burnin)){
    # Compute the e-values
    for(i in burnin:N){
      f1s[i] <- log_f1(X = X[1:i,], i = i)
      f0s <- log_f0(X[1:i,], i)[burnin:i]

      if(!any(is.finite(f0s))) stop("log(f0) is infinite")

      logE[i] <- sum(f1s[burnin:i] - f0s)
    }
    # Adapt sizes
    f1s <- f1s[burnin:N]
    logE <- c(rep(0, (burnin - 1)), logE[burnin:N])

  } else {
    # Compute the e-values
    for(i in 1:N){
      f1s[i] <- log_f1(X = X[1:i,], i = i)
      f0s <- log_f0(X[1:i,], i)

      if(!any(is.finite(f0s))) {
        print(i)
        stop("log(f0) is infinite")
      }

      logE[i] <- sum(f1s[1:i] - f0s)
    }
  }
  if(!any(is.finite(f1s))) stop("log(f1) is infinite")
  logE
}

UIE_test <- function(X, log_f0, log_f1, N, alpha, burnin = NULL) {
  logUIE <- UIE(X, log_f0, log_f1, N, burnin = burnin)
  test_res <- logUIE > log(1 / alpha)
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


