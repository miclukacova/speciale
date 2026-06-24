# The universal inference process
UIE <- function(X, log_f0, log_f1, N, Sigma, sigmaUnknown, m_init, burnin) {
  # Storage
  logE <- vector(length = N)
  f1s <- vector(length = N)

  # Estimates
  # Mean
  mT_est  <- cumsum(X[,1]) / 1:N
  mC_est  <- cumsum(X[,2]) / 1:N
  m0_est <- (cumsum(X) / (1:(2 * N)))[seq(2, 2 * N, by = 2)]

  # Variance covariance matrix
  if(sigmaUnknown){
    # Moments
    mT2_est <- cumsum(X[,1]^2) / 1:N
    mC2_est <- cumsum(X[,2]^2) / 1:N
    mTC_est <- cumsum(X[,1] * X[,2]) / 1:N

    # Variance covariance
    s2_1 <- mT2_est - (mT_est)^2
    s2_2 <- mC2_est - (mC_est)^2
    s2_12 <- mTC_est - mT_est * mC_est

    # Initialize
    sigma_est <- diag(2)
    sigma_est_1 <-  diag(2)

  } else {
    sigma_est <- Sigma
    sigma_est_1 <- Sigma
  }

  # Since the f1 process should be predictable we add an initial mT and mC value
  mT_est  <- c(m_init, mT_est)
  mC_est  <- c(m_init, mC_est)

  # Calculate the e-process
  for(i in burnin:N){

    if(sigmaUnknown & i > 2){
      # Update variance estimate
      sigma_est_1 <- sigma_est
      sigma_est <- matrix(c(s2_1[i], s2_12[i], s2_12[i], s2_2[i]), ncol = 2)
    }

    # Compute the predictable density estimate evaluated in the newest observation
    f1s[i] <- log_f1(X = X[i,,drop = FALSE],
                     mT_est = mT_est[i],
                     mC_est = mC_est[i],
                     sigma_est = sigma_est_1)

    f0s <- log_f0(X = X[1:i,,drop = FALSE],
                  mu_est = m0_est[i],
                  sigma_est = sigma_est)

    if(!any(is.finite(f0s))) {
      print(i)
      stop("log(f0) is infinite")
    }

    # Compute the log e-values
    logE[i] <- sum(f1s[1:i] - f0s)
  }

  if(!any(is.finite(f1s))) stop("log(f1) is infinite")
  logE
}

UIE_test <- function(X, log_f0, log_f1, N, Sigma, sigmaUnknown = FALSE, m_init = 0.3, burnin = 1) {
  logUIE <- UIE(X = X,
                log_f0 = log_f0,
                log_f1 = log_f1,
                N = N,
                Sigma = Sigma,
                sigmaUnknown = sigmaUnknown,
                m_init = m_init,
                burnin = burnin)

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


