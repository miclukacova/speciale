get_comp_soko_eproc <- function(alpha,
                                B,
                                N,
                                p_c,
                                p_t,
                                p_t_ms,
                                c,
                                theta) {
  if(FALSE) {
    library(tidyverse)
    source("~/Desktop/Uni/Speciale/speciale/functions/HCP_test.R")
  }

  #------------------
  # Settings
  #------------------

  set.seed(92638)
  m_0    <- (p_c - p_c + 1) / 2
  lambda_star <- (p_t * (1 - p_c) - (1 - p_t) * p_c) / (p_t * (1 - p_c) + (1 - p_t) * p_c)
  lambda_star_ms <- (p_t_ms * (1 - p_c) - (1 - p_t_ms) * p_c) / (p_t_ms * (1 - p_c) + (1 - p_t_ms) * p_c)

  sample_patient <- function(N, p_t) {
    pat_c <- rbinom(N, 1, p_c)
    pat_t <- rbinom(N, 1, p_t)
    (pat_t - pat_c + 1) / 2
  }

  #--------------------
  # Sokolova E-process
  #--------------------

  eprocess_run <- function(X,
                           N,
                           alpha,
                           lambda,
                           returnPath = FALSE) {

    e_vals <- cumprod(1 + lambda * (2 * X - 1))

    if (returnPath) {
      return(tibble(n = seq_len(N), value = e_vals))
    }

    reject <- any(e_vals >= 1 / alpha)

    n_stop <- if (reject) {
      min(which(e_vals >= 1 / alpha))
    } else {
      N
    }

    tibble(n = n_stop, reject = reject)
  }

  #--------------------------------
  # Wrapper for a single simulation
  #--------------------------------

  run_one_sim <- function(p_c,
                          p_t,
                          N,
                          alpha,
                          m_0,
                          c,
                          theta,
                          lambda_star,
                          lambda_star_ms) {

    X <- sample_patient(N, p_t)

    e_out <- eprocess_run(X = X,
                          N = N,
                          alpha = alpha,
                          lambda = lambda_star)

    e_out_ms <- eprocess_run(X = X,
                             N = N,
                             alpha = alpha,
                             lambda = lambda_star_ms)

    hcp_out <- run_HCP_test(m_0 = m_0,
                            c = c,
                            X = X,
                            theta = theta,
                            alpha = alpha)

    tibble(method = c("E-process", "MS E-process", "HCP"),
           reject = c(e_out$reject, e_out_ms$reject, hcp_out[1]),
           n_stop = c(e_out$n, e_out_ms$n, hcp_out[2]))
  }

  #------------------------------------------------
  # Monte Carlo comparison
  #------------------------------------------------

  compare_tests <- function(p_c,
                            p_t,
                            N,
                            alpha,
                            m_0,
                            c,
                            theta,
                            B = 5000) {

    map_dfr(seq_len(B), ~run_one_sim(p_c = p_c,
                                     p_t = p_t,
                                     N = N,
                                     alpha = alpha,
                                     m_0 = m_0,
                                     c = c,
                                     theta = theta,
                                     lambda_star = lambda_star,
                                     lambda_star_ms = lambda_star_ms), .id = "sim")
  }


  #------------------------------------------------
  # Type I error
  #------------------------------------------------

  # Choose parameters satisfying H0

  null_results <- compare_tests(p_c = p_c,
                                p_t = p_c,
                                N = N,
                                alpha = alpha,
                                m_0 = m_0,
                                c = c,
                                theta = theta,
                                B = B)

  type1_summary <- null_results |>
    group_by(method) |>
    summarise(
      typeI = mean(reject),
      avg_stop = mean(n_stop)
    )

  #------------------------------------------------
  # Power curve
  #------------------------------------------------

  alt_grid <- tibble(p_t = seq(0.35, 0.7, by = 0.05))

  power_results <- map_dfr(alt_grid$p_t, function(p_t_alt) {

        out <- compare_tests(p_c = p_c,
                             p_t = p_t_alt,
                             N = N,
                             alpha = alpha,
                             m_0 = m_0,
                             c = c,
                             theta = theta,
                             B = B)

        out |> group_by(method) |>
          summarise(
            power = mean(reject),
            avg_stop = mean(n_stop),
            .groups = "drop") |>
          mutate(p_t = p_t_alt)})

  #------------------------------------------------
  # Plot power
  #------------------------------------------------

  plot1 <- ggplot(power_results, aes(p_t, power, colour = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_manual(values = c("E-process" = "steelblue",
                                   "HCP" = "firebrick",
                                   "MS E-process" = "darkgreen"),
                        labels = c("E-process" = "Soko(0.45)-test",
                                   "HCP" = "HCP-test",
                                   "MS E-process" = "Soko(0.6)-test"))+
    labs(x = expression(p[t]), y = "Power") +
    theme_bw()

  #------------------------------------------------
  # Plot expected stopping time
  #------------------------------------------------

  plot2 <- ggplot(power_results, aes(p_t, avg_stop, colour = method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_colour_manual(values = c("E-process" = "steelblue",
                                   "HCP" = "firebrick",
                                   "MS E-process" = "darkgreen"),
                        labels = c("E-process" = "Soko(0.45)-test",
                                   "HCP" = "HCP-test",
                                   "MS E-process" = "Soko(0.6)-test"))+
    labs(x = expression(p[t]), y = "Expected sample size") +
    theme_bw()


  #------------------------------------------------
  # Example path under alternative
  #------------------------------------------------

  set.seed(3978435)
  n_rep = 10
  sim_one <- function(rep_id, lambda_star) {
    X <- sample_patient(N, p_t)

    e_path <- eprocess_run(X, N, alpha, lambda_star, returnPath = TRUE) |>
        mutate(process = "E-process")

    hcp_path <- tibble(n = 1:N,
                       value = HCP(m_0 = m_0,
                                   c = c,
                                   X = X,
                                   theta = theta,
                                   alpha = alpha),
                       process = "HCP")

    bind_rows(e_path, hcp_path) |>
      mutate(rep = rep_id)
    }

    df <- bind_rows(lapply(1:n_rep, function(rep) sim_one(rep, lambda_star)))
    df_ms <- bind_rows(lapply(1:n_rep, function(rep) sim_one(rep, lambda_star_ms)))

    y_range <- range(c(df$value, df_ms$value), na.rm = TRUE)

    plot3 <- ggplot(df, aes(n, value, colour = process, group = interaction(process, rep))) +
      geom_line(alpha = 0.5) +
      geom_hline(yintercept = 1 / alpha, linetype = 2) +
      scale_y_log10(limits = y_range) +
      labs(x = "N", y = "Process value (log scale)")+
      scale_colour_manual(values = c("E-process" = "steelblue",
                                     "HCP" = "firebrick"),
                          labels = c("E-process" = "Soko(0.45)-test",
                                     "HCP" = "HCP"))+
      theme_bw()

    plot4 <- ggplot(df_ms, aes(n, value, colour = process, group = interaction(process, rep))) +
      geom_line(alpha = 0.5) +
      geom_hline(yintercept = 1 / alpha, linetype = 2) +
      scale_y_log10(limits = y_range) +
      labs(x = "N", y = "Process value (log scale)")+
      scale_colour_manual(values = c("E-process" = "darkgreen",
                                     "HCP" = "firebrick"),
                          labels = c("E-process" = "Soko(0.6)-test",
                                     "HCP" = "HCP"))+
      theme_bw()

    return(list(type1_summary = type1_summary,
                plot1 = plot1,
                plot2 = plot2,
                plot3 = plot3,
                plot4 = plot4))

}


