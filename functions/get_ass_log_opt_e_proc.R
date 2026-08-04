get_ass_log_opt_e_proc <- function(mu0,
                                   theta_nu,
                                   N,
                                   num_var = 100) {

  set.seed(63628)

  # Simulate
  X <- matrix(rnorm(num_var * N, theta_nu, 1), ncol = num_var)

  # Estimate means
  mu_hat <- apply(X, 2, \(x) cumsum(x) / seq_len(N))

  # Compute log e-processes
  logE <- sapply(1:num_var, \(i)
                 cumsum(dnorm(X[, i], mu_hat[, i], log = TRUE) -
                          dnorm(X[, i], mu0, log = TRUE))
  )

  logE_nu <- sapply(1:num_var, \(i)
                    cumsum(dnorm(X[, i], theta_nu, log = TRUE) -
                             dnorm(X[, i], mu0, log = TRUE))
  )

  # Store everything
  D <- tibble(
    n = rep(1:N, num_var),
    path = rep(1:num_var, each = N),
    logE = c(logE),
    logE_nu = c(logE_nu),
    regret = c((logE - logE_nu) / seq_len(N))
  )

  # Long format for plotting
  D_long <- D %>%
    pivot_longer(
      c(logE, logE_nu),
      names_to = "process",
      values_to = "value"
    ) %>%
    mutate(type = "individual")

  # Add averages
  D_avg <- D_long %>%
    group_by(n, process) %>%
    summarise(value = mean(value), .groups = "drop") %>%
    mutate(
      path = NA,
      type = "average"
    )

  D_plot <- bind_rows(D_long, D_avg)

  # Plot log e-processes
  p1 <- ggplot(D_plot, aes(n, value, colour = process)) +
    geom_line(
      data = subset(D_plot, type == "individual"),
      aes(group = interaction(process, path)),
      alpha = 0.15,
      linewidth = 0.4
    ) +
    geom_line(
      data = subset(D_plot, type == "average"),
      aes(group = process),
      linewidth = 1.2
    ) +
    scale_colour_manual(
      values = c(logE = "steelblue", logE_nu = "firebrick"),
      labels = c(logE = "Adaptive", logE_nu = "Oracle")
    ) +
    labs(x = "n", y = "", colour = NULL) +
    theme_minimal(base_size = 15)


  # Regret plot
  D_regret <- D %>%
    group_by(n) %>%
    summarise(
      regret = mean(regret),
      .groups = "drop"
    ) %>%
    mutate(type = "average")

  p2 <- ggplot(D, aes(n, regret, group = path)) +
    geom_line(alpha = 0.15, linewidth = 0.4,
              colour = "firebrick") +
    geom_line(
      data = D_regret,
      aes(group = 1),
      linewidth = 1.1,
      colour = "black"
    ) +
    geom_hline(yintercept = 0, linetype = 2) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "n",
      y = expression(
        frac(log(Lambda[n]) - log(Lambda[n]^theta[1]), n)
      )
    ) +
    theme_minimal(base_size = 15)

  return(list(p1 = p1, p2 = p2))
}

