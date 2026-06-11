######################################################
######## A general version of the SPRT process #######
######################################################

run_sprt_test <- function(N, X, f0, f1, beta, alpha) {

  gamma0 <- beta / (1 - alpha)
  gamma1 <- (1 - alpha) / beta

  L <- cumprod(f1(X) / f0(X))
  ss <- min(which(L >= gamma1 | gamma0 >= L), N)

  return(c(L[ss] >= gamma1, ss))
}

# Example run
Example = FALSE
if(Example){
  N <- 500
  m_0 <- 0.5
  m_1 <- 0.6
  f0 <- function(X) dbinom(X, 1,  m_0)
  f1 <- function(X) dbinom(X, 1,  m_1)
  f1_adap <- function(X) {
    # OBS: estimate has to be predictable
    N <- length(X)
    m_est <- cumsum(c(0.5, X[1:(N-1)])) / seq_len(N)
    dbinom(X, 1,  m_est)
  }
  m_true <- 0.6
  alpha <- 0.05
  X <- rbinom(N, 1, m_true)

  run_sprt_test(N,
                X = X,
                f0 = f0,
                f1 = f1_adap,
                beta = alpha,
                alpha = alpha)
}

