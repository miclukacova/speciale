######################################################
######## A general version of the HCP process ########
######################################################

# Calculating the predictable sequence
calculate_lambda_tilde <- function(X, alpha){
  n <- length(X)
  mu_hat <- (1 /2 + cumsum(X)) / (1:n +1)
  sigma_hat <- (1 / 4 + cumsum((X - mu_hat) ^ 2)) / (1:n + 1)

  sigma_hat_t_1 <- c(1 / 4, sigma_hat)[1:n]
  lambda_tilde <- sqrt(2 * log(2 / alpha) / (sigma_hat_t_1 * 1:n * log(1:n + 1)))
  lambda_tilde

}

# The hedged capital process
HCP <- function(m_0, c, X, theta, alpha) {

  # Calculate the lambda's
  lambda_tilde <- calculate_lambda_tilde(X, alpha = alpha)
  lambda_plus <- pmin(abs(lambda_tilde), c / m_0)
  lambda_minus <- pmin(abs(lambda_tilde), c / (1 - m_0))

  # Calculate the capital process
  K_plus <- cumprod(1 + lambda_plus * (X - m_0))
  K_minus <- cumprod(1 - lambda_minus * (X - m_0))

  return(pmax(theta * K_plus, (1 - theta) * K_minus))
}

# The HCP test
run_HCP_test <- function(m_0, c, X, theta, alpha) {
  HCP_res <- HCP(m_0, c, X, theta, alpha)
  test_res <- HCP_res > 1/ alpha
  test_fut <- HCP_res < alpha / 2
  if(any(test_res) | any(test_fut)){
    ESS <- which((test_res + test_fut) == 1)[1]
  } else {
    ESS <- length(HCP_res)
  }
  return(c(Reject = any(test_res[1:ESS]), ESS = ESS))
}

# Example run
Example = FALSE
if(Example){
  # Parameters
  m_0 <- 0.5
  c <- 1 / 2
  N <- 200
  theta <- 3 / 4
  alpha <- 0.05
  m_true <- 0.6

  # Test run
  X <- rbinom(N, 1, m_true)
  run_HCP_test(m_0 = m_0,
               c = c,
               X = X,
               theta = theta,
               alpha)
}

# Plot process
plot_process = FALSE
if(plot_process){
  # Parameters
  m_0 <- 0.5
  c <- 1 / 2
  N <- 200
  theta <- 3 / 4
  alpha <- 0.05
  m_1 <- 0.6

  simulate_path <- function(m_true) {
    X <- rbinom(N, 1, m_true)
    HCP(m_0 = m_0, c = c, X = X, theta = theta)
  }

  HCP_mat <- replicate(10, simulate_path(m_true = m_1))

  df_alt <- data.frame(
    time = rep(1:N, 10),
    path = rep(1:10, each = N),
    HCP = as.vector(HCP_mat),
    scenario = "Alternative (m_true ≠ m_0)"
  )

  HCP_mat <- replicate(10, simulate_path(m_true = m_0))
  df_null <- data.frame(
    time = rep(1:N, 10),
    path = rep(1:10, each = N),
    HCP = as.vector(HCP_mat),
    scenario = "Null (m_true = m_0)"
  )

  df <- rbind(df_alt, df_null)

  p1 <- ggplot(df, aes(time, HCP, group = path)) +
    geom_line(alpha = 0.4) +
    geom_hline(yintercept = 1 / alpha,
               colour = "red",
               linetype = 2) +
    geom_hline(yintercept = alpha,
               colour = "red",
               linetype = 2) +
    scale_y_log10() +
    theme_minimal() +
    facet_wrap(~ scenario, ncol = 2)
}


