get_e_vs_p_power <- function(
    alphas = c(0.001, 0.01, 0.05, 0.08),
    ns = seq(50, 250, by = 50),
    true_means = seq(0.1, 0.9, by = 0.05),
    n = 100,                                  # sample size
    B = 5000,
    alpha = 0.05,
    true_mean = 0.7,
    p = 1/2,                                  # The value we test for
    c = 3/4,                                  # Truncation level
    sample_fct
) {

  if(FALSE){
    source("~/Desktop/Uni/Speciale/speciale/functions/HCP_test.R")
    source("~/Desktop/Uni/Speciale/speciale/functions/define_setting.r")
  }

  #Function for getting rejectiong probabilities
  e_p_vals <- function(n, B, alpha, true_mean, p, var){
    e_reject <- numeric(B)
    p_reject <- numeric(B)

    for (b in 1:B) {
      #Simulating data - binomial distribution cause its easy
      X <- sample_fct(n, true_mean = true_mean)

      # Calculate the lambda's
      lambda_tilde <- calculate_lambda_tilde(X, alpha = alpha)
      lambda_plus <- pmin(abs(lambda_tilde), c / p)[n]
      lambda_minus <- pmin(abs(lambda_tilde), c / (1 - p))[n]

      #one-sided t-test for mean = 1/2
      if (sd(X) == 0) {
        p_reject[b] <- mean(X) != 0.5
      } else {
        p_val <- t.test(X, alternative = "two.sided", mu = p)$p.value
        p_reject[b] <- p_val <= alpha
      }

      #e-value from Ramdas article
      e_val <- max(1/2 * prod(1 + lambda_plus * (X - 0.5)), 1/2 * prod(1 - lambda_minus * (X - 0.5)))
      e_reject[b] <- e_val >= 1/alpha
    }

    data.frame(
      p_power = mean(p_reject),
      e_power = mean(e_reject)
    )
  }

  #Getting values for different variables
  alpha_data <- do.call(rbind, lapply(alphas, function(a) {
    cbind(
      variable = "alpha",
      x = a,
      e_p_vals(n = n, B = B, alpha = a, true_mean = true_mean, var = alpha, p = p)
    )
  }))

  n_data <- do.call(rbind, lapply(ns, function(n_val) {
    cbind(
      variable = "N",
      x = n_val,
      e_p_vals(n = n_val, B = B, alpha = alpha, true_mean = true_mean, var = n, p = p)
    )
  }))

  true_mean_data <- do.call(rbind, lapply(true_means, function(true_mean_val) {
    cbind(
      variable = "p",
      x = true_mean_val,
      e_p_vals(n = n, B = B, alpha = alpha, true_mean = true_mean_val, var = true_mean, p = p)
    )
  }))

  power_data <- rbind(alpha_data, n_data, true_mean_data)

  return(data.frame(power_data, row.names = NULL))
}

