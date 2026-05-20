################################################################################
# Core functions
################################################################################

# Function for sampling patients
sample_patient <- function(p_c, p_t, n) {
  pat_c <- rbinom(n,1,p_c)
  pat_t <- rbinom(n,1,p_t)
  pat_t - pat_c
}

#-------------------------------------------------------------------------------
# Obrian/Pocock boundaries

alphas_soko <- function(p_c, p_t, n_looks, Nmax, B = 10^4, alpha = 0.05) {

  # Define quantities
  look_times <- round(seq(Nmax/n_looks, Nmax, length.out = n_looks))
  info_frac <- look_times/Nmax

  # Sample data
  xT_null_mat <- matrix(stats::rbinom(B * Nmax, 1, p_c),nrow = B)
  xC_null_mat <- matrix(stats::rbinom(B * Nmax, 1, p_c),nrow = B)
  xT_alt_mat <- matrix(stats::rbinom(B * Nmax, 1, p_t),nrow = B)
  xC_alt_mat <- matrix(stats::rbinom(B * Nmax, 1, p_c),nrow = B)

  # Sum number of sucesses in each stage for all B replications
  cumT_null <- t(apply(xT_null_mat, 1, cumsum))[, look_times, drop = FALSE]
  cumC_null <- t(apply(xC_null_mat, 1, cumsum))[, look_times, drop = FALSE]
  cumT_alt <- t(apply(xT_alt_mat, 1, cumsum))[, look_times, drop = FALSE]
  cumC_alt <- t(apply(xC_alt_mat, 1, cumsum))[, look_times, drop = FALSE]

  # Number of observations in each stage for all B replications
  nn_looks <- matrix(look_times, nrow = B, ncol = length(look_times), byrow = TRUE)

  # Estimated effect sizes divided by standard error
  z_null <- evalinger:::.z_matrix(cumT_null, cumC_null, nn_looks)
  z_alt <- evalinger:::.z_matrix(cumT_alt, cumC_alt, nn_looks)

  # We find the maximal effect size across stages for each of the B trials
  m_null <- apply(z_null * matrix(sqrt(info_frac), nrow = B, ncol = length(info_frac), byrow = TRUE), 1, max)
  # Take the quantile of the maximal effect sizes
  gs_c <- as.numeric(stats::quantile(m_null, probs = 1 - alpha, names = FALSE))
  obf_bounds <- gs_c/sqrt(info_frac)

  return(obf_bounds)
}

################################################################################
# Simulation Functions
################################################################################

#-------------------------------------------------------------------------------
# E-process
#-------------------------------------------------------------------------------

eprocess_run <- function(p_c, p_t, Nmax, alpha, lambda = NULL, alphas, n_looks, null = FALSE) {

  if(is.null(lambda)) {
    lambda_star <- (
      p_t * (1 - p_c) - (1 - p_t) * p_c
    ) / (
      p_t * (1 - p_c) + (1 - p_t) * p_c
    )
  }
  else lambda_star <- lambda

  if(null) X <- sample_patient(p_c, p_c, Nmax)
  else X <- sample_patient(p_c, p_t, Nmax)

  e_vals <- cumprod(1 + lambda_star * X)

  n_stop <- min(which(e_vals >= 1 / alpha), Nmax)

  reject <- any(e_vals >= 1 / alpha)

  tibble(
    n = n_stop,
    reject = reject
  )
}


#-------------------------------------------------------------------------------
# SPRT
#-------------------------------------------------------------------------------

sprt_run <- function(p_c, p_t, Nmax, alpha, lambda = NULL, alphas, n_looks, null = FALSE){

  if(null) X <- sample_patient(p_c, p_c, Nmax)
  else X <- sample_patient(p_c, p_t, Nmax)

  f_1 <- p_t * (1 - p_c) * (X == 1) +  (1 - p_t) * p_c * (X == -1) + (X == 0) * ((1 - p_t) * (1 - p_c) + p_t * p_c)
  f_0 <-  p_c * (1 - p_c) * (X == 1 | X == -1) + (p_c * p_c + (1 - p_c)^2)* (X == 0)

  ratio <- cumprod(f_1 / f_0)

  n_stop <- min(
    which(ratio <= alpha | ratio >= 1 / alpha),
    Nmax
  )

  reject <- any(ratio >= 1 / alpha)

  tibble(
    n = n_stop,
    reject = reject
  )
}

#-------------------------------------------------------------------------------
# Group Sequential
#-------------------------------------------------------------------------------

gs_run <- function(p_c, p_t, Nmax, alpha, lambda = NULL, alphas, n_looks, null = FALSE) {

  # Information times / looks
  look_times <- round(seq(Nmax / n_looks, Nmax, length.out = n_looks))

  # Incremental sample sizes at each stage
  n_k <- c(look_times[1], diff(look_times))

  # Null SD
  sd_0 <- sqrt(2 * p_c * (1 - p_c))

  # Simulate data
  if(null) X <- sample_patient(p_c, p_c, Nmax)
  else X <- sample_patient(p_c, p_t, Nmax)

  # For storage
  Z_k_star <- vector()

  # Compute cumulative statistics
  for(k in 1:n_looks) {

    # Observations up to look k
    X_k <- X[1:look_times[k]]

    # Cumulative means
    Xbar_k <- mean(X_k)

    # Equation (1.2)
    Z_k_star[k] <- Xbar_k / sd_0 * sqrt(look_times[k])
  }

  # Boundary crossing
  test_res <- Z_k_star >= alphas

  # Sample size
  n <-  look_times[min(which(test_res), n_looks)]

  tibble(
    n = n,
    reject = any(test_res),
  )

}

#-------------------------------------------------------------------------------
# Main simulation wrapper
#-------------------------------------------------------------------------------

run_scenario <- function(
    p_c,
    p_t,
    Nmax,
    alpha,
    B = 10000,
    alphas,
    lambda = NULL,
    n_looks
) {

  methods <- c("E-process", "SPRT", "GS")

  res <- vector("list", length(methods))
  i <- 1
  for(method in methods){

    if(method == "E-process") run_func <- eprocess_run
    else if(method == "SPRT") run_func <- sprt_run
    else run_func <- gs_run

    # Under alternative
    tmp <- replicate(
      B,
      run_func(p_c = p_c, p_t = p_t, Nmax = Nmax, alpha = alpha, alphas = alphas,
               n_looks = n_looks, lambda = lambda),
      simplify = FALSE
    )

    # Under nul
    tmp2 <- replicate(
      B,
      run_func(p_c = p_c, p_t = p_t, Nmax = Nmax, alpha = alpha, alphas = alphas,
               n_looks = n_looks, lambda = lambda, null = TRUE),
      simplify = FALSE
    )

    tmp <- bind_rows(tmp)
    tmp2 <- bind_rows(tmp2)

    res[[i]] <- tibble(
      method = method,
      mean_n_alt = mean(tmp$n),
      sd_n_alt = sd(tmp$n),
      power = mean(tmp$reject),
      mean_n_nul = mean(tmp2$n),
      sd_n_nul = sd(tmp2$n),
      typeIerr = mean(tmp2$reject)
    )
    i <- i +1
  }

  bind_rows(res)
}

