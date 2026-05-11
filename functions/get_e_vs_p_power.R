

get_e_vs_p_power <- function(
    alphas = c(0.001, 0.01, 0.05),
    Bs = c(1000, 5000, 10000),
    ns = seq(50, 250, by = 50),
    ps = seq(0.1, 0.9, by = 0.1),
    n = 100,
    B = 5000,
    alpha = 0.05,
    p = 0.7,
    lambda = 0.5
) {

  #Function for getting rejectiong probabilities
  e_p_vals <- function(n, B, alpha = 0.05, p = 0.7, mu = 1/2, lambda = 0.5, var){
    e_reject <- numeric(B)
    p_reject <- numeric(B)

    for (b in 1:B) {
      #Simulating data - binomial distribution cause its easy
      X <- rbinom(n, size = 1, prob = p)

      #one-sided t-test for mean = 1/2
      p_val <- t.test(X, alternative = "two.sided", mu = 1/2)$p.value
      p_reject[b] <- p_val <= alpha

      #e-value from Ramdas article
      e_val <- 1/2 * prod(1 + lambda * (X - 0.5)) + 1/2 * prod(1 - lambda * (X - 0.5))
      e_reject[b] <- e_val >= 1/alpha
    }

    data.frame(
      p_value_power = mean(p_reject),
      e_value_power = mean(e_reject)
    )
  }


  #Getting values for different variables
  alpha_data <- do.call(rbind, lapply(alphas, function(a) {
    cbind(
      variable = "alpha",
      x = a,
      e_p_vals(n = n, B = B, alpha = a, p = p, lambda = lambda, var = alpha)
    )
  }))

  B_data <- do.call(rbind, lapply(Bs, function(B_val) {
    cbind(
      variable = "B",
      x = B_val,
      e_p_vals(n = n, B = B_val, alpha = alpha, p = p, lambda = lambda, var = B)
    )
  }))

  n_data <- do.call(rbind, lapply(ns, function(n_val) {
    cbind(
      variable = "n",
      x = n_val,
      e_p_vals(n = n_val, B = B, alpha = alpha, p = p, lambda = lambda, var = n)
    )
  }))

  p_data <- do.call(rbind, lapply(ps, function(p_val) {
    cbind(
      variable = "p",
      x = p_val,
      e_p_vals(n = n, B = B, alpha = alpha, p = p_val, lambda = lambda, var = p)
    )
  }))

  power_data <- rbind(alpha_data, B_data, n_data, p_data)

  data.frame(power_data, row.names = NULL)
}


