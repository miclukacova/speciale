######################################################
######### A general version of the GS test ###########
######################################################

## To do employ this test

gs_run <- function(Nmax, alphas, n_looks, X, m_0, side = 1, sigmaUnknown = TRUE) {

  # Information times / looks
  look_times <- round(seq(Nmax / n_looks, Nmax, length.out = n_looks))

  # Incremental sample sizes at each stage
  n_k <- c(look_times[1], diff(look_times))

  # For storage
  Z_k_star <- vector()

  # Compute cumulative statistics
  for(k in 1:n_looks) {

    # Observations up to look k
    X_k <- X[1:look_times[k]]

    # Cumulative means
    Xbar_k <- mean(X_k)

    # Cumulative variance estimate
    if(sigmaUnknown) sigma_hat_k <- sqrt(1 / (look_times[k] - 1) * sum((X_k - Xbar_k)^2))
    else sigma_hat_k <- sigmaUnknown
    if(sigma_hat_k == 0) warning("Variance estimate is 0")

    # Z_k_star
    Z_k_star[k] <- (Xbar_k - m_0) / sigma_hat_k * sqrt(look_times[k])
  }

  # If we do not know sigma, we need to correct with the t distribution
  if(sigmaUnknown) alphas <- qt(p = pnorm(alphas), df = look_times - 1)

  # Boundary crossing
  if(side == 2) test_res <- abs(Z_k_star) >= alphas
  else test_res <- Z_k_star >= alphas

  # Sample size
  n <-  look_times[min(which(test_res), n_looks)]

  c(Reject = any(test_res), ESS = n)

}

# Obrian/Pocock boundaries

alphas_soko <- function(p_c, n_looks, Nmax, B = 10^4, alpha = 0.05) {

  # Define quantities
  look_times <- round(seq(Nmax/n_looks, Nmax, length.out = n_looks))
  info_frac <- look_times/Nmax

  # Sample data
  xT_null_mat <- matrix(stats::rbinom(B * Nmax, 1, p_c),nrow = B)
  xC_null_mat <- matrix(stats::rbinom(B * Nmax, 1, p_c),nrow = B)

  # Sum number of sucesses in each stage for all B replications
  cumT_null <- t(apply(xT_null_mat, 1, cumsum))[, look_times, drop = FALSE]
  cumC_null <- t(apply(xC_null_mat, 1, cumsum))[, look_times, drop = FALSE]

  # Number of observations in each stage for all B replications
  nn_looks <- matrix(look_times, nrow = B, ncol = length(look_times), byrow = TRUE)

  # Estimated effect sizes divided by standard error
  z_null <- evalinger:::.z_matrix(cumT_null, cumC_null, nn_looks)

  # We find the maximal effect size across stages for each of the B trials
  m_null <- apply(z_null * matrix(sqrt(info_frac), nrow = B, ncol = length(info_frac), byrow = TRUE), 1, max)
  # Take the quantile of the maximal effect sizes
  gs_c <- as.numeric(stats::quantile(m_null, probs = 1 - alpha, names = FALSE))
  obf_bounds <- gs_c / sqrt(info_frac)

  return(obf_bounds)
}

