######################################################
######### A general version of the GS test ###########
######################################################

gs_run <- function(Nmax,
                   alphas,
                   n_looks,
                   X,
                   m_0,
                   side = 1,
                   sigmaUnknown = TRUE,
                   sigma = NULL) {

  # Information times / looks
  look_times <- round(seq(Nmax / n_looks, Nmax, length.out = n_looks))

  cs  <- cumsum(X)
  cs2 <- cumsum(X^2)

  Xbar_k <- cs[look_times] / look_times

  if (sigmaUnknown) {
    var_k <- (cs2[look_times] - cs[look_times]^2 / look_times) /
      (look_times - 1)
    sigma_hat_k <- sqrt(var_k)
  } else {
    sigma_hat_k <- rep(sigma, length(look_times))
  }

  # Z_k_star
  Z_k_star <- (Xbar_k - m_0) / sigma_hat_k * sqrt(look_times)

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
  look_times <- round(seq(Nmax / n_looks, Nmax, length.out = n_looks))
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

