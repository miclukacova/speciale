get_ass_log_opt_e_proc <- function(mu0,
                                   theta_nu,
                                   N) {

  set.seed(1)

  # Simulate under nu
  X <- rnorm(N, theta_nu, 1)

  # Plug-in estimator
  mu_hat <- cumsum(X) / seq_len(N)

  # Adaptive e-process
  logE <- cumsum(
    dnorm(X, mu_hat, log = TRUE) -
      dnorm(X, mu0, log = TRUE)
  )

  # Oracle likelihood ratio process
  logE_nu <- cumsum(
    dnorm(X, theta_nu, log = TRUE) -
      dnorm(X, mu0, log = TRUE)
  )

  # Asymptotic log-optimality quantity
  regret <- (logE - logE_nu) / seq_len(N)

  D <- tibble(
    t = 1:N,
    logE = logE,
    logE_nu = logE_nu,
    regret = regret
  )

  D_long <- D %>%
    select(t, logE, logE_nu) %>%
    pivot_longer(-t)

  p1 <- ggplot(D_long,
               aes(t, value, color = name)) +
    geom_line(linewidth = 1.1) +
    labs(
      x = "Time",
      y = expression(log(E[n])),
      color = NULL
    ) +
    scale_color_manual(
      values = c(
        "logE" = "steelblue",
        "logE_nu" = "firebrick"
      ),
      labels = c(
        "Adaptive",
        "Oracle"
      )
    ) +
    theme_minimal(base_size = 15)

  p2 <- ggplot(D, aes(t, regret)) +
    geom_line(linewidth = 1.1, color = "firebrick") +
    geom_hline(yintercept = 0,
               linetype = 2) +
    labs(
      x = "Time",
      y = expression(
        frac(log(E[n]) - log(E[n]^nu), n)
      )
    ) +
    theme_minimal(base_size = 15)

  return(list(p1 = p1, p2 = p2))
}

