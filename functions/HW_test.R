######################################################
######## A general version of the HW-test ###########
######################################################

# sample_data_null should return the bootstrapped test statistic

run_HW_test <- function(N,
                        X,
                        sample_data_null = NULL,
                        calc_q_n = NULL,
                        gamma,
                        quanti,
                        B = 1000,
                        return_Q_n = FALSE) {

  if(is.null(calc_q_n)) {

    Q_n <- vector(length = N)
    i <- 2
    while(i <= (N-1) & all(Q_n < gamma)) {
      boot_data <- numeric(B)
      boot_data <- numeric(B)
      for(b in seq_len(B))
        boot_data[b] <- sample_data_null(
          N = N - i,
          X = X[1:i]
        )
      Q_n[i] <- mean(boot_data >= quanti)
      i <- i + 1
    }
    Q_n[N] <- sample_data_null(N = N, X = X) > quanti

  } else {
    Q_n <- calc_q_n(X, N, quanti)
  }

  Reject = any(Q_n >= gamma)

  if(Reject) ESS <- which(Q_n >= gamma)[1]
  else ESS <- N

  if(return_Q_n) return(Q_n)

  return(c(Reject = Reject, ESS = ESS))
}

Example = FALSE
if(Example){
  # Parameter values
  m_0 <- 0.5
  m_1 <- 0.6
  m_true <- 0.6
  alpha <- 0.05
  gamma <- 0.9
  N <- 300

  # Data sampling function, quantile, function to calculate Q_n
  X <- rbinom(N, 1, m_true)
  sample_data_null <- function(N, X) {
    X_boot <- rbinom(N, 1, m_0)
    return(sum(X, X_boot))
  }
  z_ag <- qbinom(p = 1 - (alpha * gamma), size = N, prob = m_0)
  calc_q_n <- function(X, N, z_ag) {
    1 - pbinom(q = z_ag - cumsum(X), size = seq(N-1, 0), prob = m_0)
  }

  # Run of test
  run_HW_test(N = N,
                  X = X,
                  sample_data_null = sample_data_null,
                  #calc_q_n = calc_q_n,
                  gamma = gamma,
                  quanti = z_ag,
                  B = 1000)
}


# Plot process
plot_process = FALSE
if(plot_process){
  # Parameter values
  m_0 <- 0.5
  m_1 <- 0.6
  m_true <- 0.6
  alpha <- 0.05
  gamma <- 0.9
  N <- 200
  n <- 8

  # Data sampling function, quantile, function to calculate Q_n
  z_ag <- qbinom(p = 1 - (alpha * gamma), size = N, prob = m_0)
  calc_q_n <- function(X, N) {
    1 - pbinom(q = z_ag - cumsum(X), size = seq(N-1, 0), prob = m_0)
  }

  # Run of test
  simulate_path <- function(m_true) {
    X <- rbinom(N, 1, m_true)
    run_HW_test(N = N,
                    X = X,
                    return_Q_n = TRUE,
                    calc_q_n = calc_q_n,
                    gamma = gamma,
                    quanti = z_ag)
  }

  set.seed(4825)

  Holm_mat <- replicate(n, simulate_path(m_true = m_1))

  df_alt <- data.frame(
    time = rep(1:N, n),
    path = rep(1:n, each = N),
    Q_n = as.vector(Holm_mat),
    scenario = "Alternative (m_true ≠ m_0)"
  )

  Holm_mat <- replicate(n, simulate_path(m_true = m_0))
  df_null <- data.frame(
    time = rep(1:N, n),
    path = rep(1:n, each = N),
    Q_n = as.vector(Holm_mat),
    scenario = "Null (m_true = m_0)"
  )

  df <- rbind(df_alt, df_null)
  df <- df[df$Q_n > 0.005,]

  p1 <- ggplot(df, aes(time, Q_n, group = path)) +
    geom_line(alpha = 0.4) +
    geom_hline(yintercept = gamma,
               colour = "red",
               linetype = 2) +
    scale_y_log10() +
    theme_minimal() +
    facet_wrap(~ scenario, ncol = 2)
}
