library(tidyverse)
source("~/Desktop/Uni/Speciale/speciale/functions/HCP_test.R")

############################################################
# Settings
############################################################

set.seed(92638)

alpha  <- 0.05
N      <- 200
p_c    <- 0.3
p_t    <- 0.45
m_0    <- (p_c - p_c + 1) / 2
c      <- 3/4
theta  <- 1
B      <- 5000
lambda_star <- (p_t * (1 - p_c) - (1 - p_t) * p_c) / (p_t * (1 - p_c) + (1 - p_t) * p_c)

sample_patient <- function(N, p_t) {
  pat_c <- rbinom(N, 1, p_c)
  pat_t <- rbinom(N, 1, p_t)
  (pat_t - pat_c + 1) / 2
}

############################################################
# E-process
############################################################

eprocess_run <- function(X,
                         N,
                         alpha,
                         lambda,
                         returnPath = FALSE) {

  e_vals <- cumprod(1 + lambda * (2 * X - 1))

  if (returnPath) {
    return(
      tibble(
        n = seq_len(N),
        value = e_vals
      )
    )
  }

  reject <- any(e_vals >= 1 / alpha)

  n_stop <- if (reject) {
    min(which(e_vals >= 1 / alpha))
  } else {
    N
  }

  tibble(
    n = n_stop,
    reject = reject
  )
}

############################################################
# Wrapper for a single simulation
############################################################

run_one_sim <- function(p_c,
                        p_t,
                        N,
                        alpha,
                        m_0,
                        c,
                        theta,
                        lambda_star) {

  X <- sample_patient(N, p_t)

  e_out <- eprocess_run(
    X = X,
    N = N,
    alpha = alpha,
    lambda = lambda_star
  )

  hcp_out <- run_HCP_test(
    m_0 = m_0,
    c = c,
    X = X,
    theta = theta,
    alpha = alpha
  )

  tibble(
    method = c("E-process", "HCP"),
    reject = c(e_out$reject, hcp_out[1]),
    n_stop = c(e_out$n, hcp_out[2])
  )
}

############################################################
# Monte Carlo comparison
############################################################

compare_tests <- function(p_c,
                          p_t,
                          N,
                          alpha,
                          m_0,
                          c,
                          theta,
                          B = 5000) {

  map_dfr(
    seq_len(B),
    ~run_one_sim(
      p_c = p_c,
      p_t = p_t,
      N = N,
      alpha = alpha,
      m_0 = m_0,
      c = c,
      theta = theta,
      lambda_star = lambda_star
    ),
    .id = "sim"
  )
}

############################################################
# Type I error
############################################################
#
# Choose parameters satisfying H0
#

null_results <- compare_tests(
  p_c = p_c,
  p_t = p_c,
  N = N,
  alpha = alpha,
  m_0 = m_0,
  c = c,
  theta = theta,
  B = B
)

type1_summary <-
  null_results |>
  group_by(method) |>
  summarise(
    typeI = mean(reject),
    avg_stop = mean(n_stop)
  )

print(type1_summary)

############################################################
# Power curve
############################################################

alt_grid <- tibble(
  p_t = seq(0.35, 0.6, by = 0.05)
)

power_results <-
  map_dfr(
    alt_grid$p_t,
    function(p_t_alt) {

      out <- compare_tests(
        p_c = 0.50,
        p_t = p_t_alt,
        N = N,
        alpha = alpha,
        m_0 = m_0,
        c = c,
        theta = theta,
        B = 2000
      )

      out |>
        group_by(method) |>
        summarise(
          power = mean(reject),
          avg_stop = mean(n_stop),
          .groups = "drop"
        ) |>
        mutate(p_t = p_t_alt)
    }
  )

############################################################
# Plot power
############################################################

ggplot(
  power_results,
  aes(p_t, power, colour = method)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = expression(p[t]),
    y = "Power",
    title = "Power comparison"
  ) +
  theme_bw()

############################################################
# Plot expected stopping time
############################################################

ggplot(
  power_results,
  aes(p_t, avg_stop, colour = method)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = expression(p[t]),
    y = "Average stopping time",
    title = "Expected sample size"
  ) +
  theme_bw()


############################################################
# Example path under alternative
############################################################

plot_sample_paths_overlay <- function(p_c,
                                      p_t,
                                      N,
                                      alpha,
                                      m_0,
                                      c,
                                      theta,
                                      n_rep = 5) {

  lambda_star <-
    (p_t * (1 - p_c) - (1 - p_t) * p_c) /
    (p_t * (1 - p_c) + (1 - p_t) * p_c)

  sim_one <- function(rep_id) {

    X <- sample_patient(N, p_t)

    e_path <- eprocess_run(
      X,
      N,
      alpha,
      lambda_star,
      returnPath = TRUE
    ) |>
      mutate(process = "E-process")

    hcp_path <- tibble(
      n = 1:N,
      value = HCP(
        m_0 = m_0,
        c = c,
        X = X,
        theta = theta,
        alpha = alpha
      ),
      process = "HCP"
    )

    bind_rows(e_path, hcp_path) |>
      mutate(rep = rep_id)
  }

  df <- bind_rows(lapply(1:n_rep, sim_one))

  ggplot(df, aes(n, value, colour = process, group = interaction(process, rep))) +
    geom_line(alpha = 0.5) +
    geom_hline(yintercept = 1 / alpha, linetype = 2) +
    scale_y_log10() +
    labs(
      x = "N",
      y = "Process value (log scale)",
      title = "5 simulated paths (overlay)"
    ) +
    theme_bw()
}


plot_sample_paths_overlay(p_c,
                          0.45,
                          N,
                          alpha,
                          m_0,
                          c,
                          theta,
                          n_rep = 5)


