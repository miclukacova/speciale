library(tidyverse)
source("~/Desktop/Uni/Speciale/speciale/functions/HCP_test.R")

############################################################
# Settings
############################################################

alpha  <- 0.05
N      <- 200

m_0    <- 0.5
c      <- 1/2
theta  <- 3/4

n_sim  <- 5000

############################################################
# E-process
############################################################

eprocess_run <- function(X,
                         N,
                         alpha,
                         lambda,
                         returnPath = FALSE) {

  e_vals <- cumprod(1 + lambda * X)

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

  X <- sample_patient(p_c, p_t, N)

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
    reject = c(e_out$reject, hcp_out$reject),
    n_stop = c(e_out$n, hcp_out$n)
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
                          n_sim = 5000) {

  lambda_star <-
    (p_t * (1 - p_c) - (1 - p_t) * p_c) /
    (p_t * (1 - p_c) + (1 - p_t) * p_c)

  map_dfr(
    seq_len(n_sim),
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

p_c_null <- 0.50
p_t_null <- 0.50

null_results <- compare_tests(
  p_c = p_c_null,
  p_t = p_t_null,
  N = N,
  alpha = alpha,
  m_0 = m_0,
  c = c,
  theta = theta,
  n_sim = n_sim
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
  p_t = seq(0.50, 0.75, by = 0.05)
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
        n_sim = 2000
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
# Illustrative sample path
############################################################

plot_sample_path <- function(p_c,
                             p_t,
                             N,
                             alpha,
                             m_0,
                             c,
                             theta) {

  lambda_star <-
    (p_t * (1 - p_c) - (1 - p_t) * p_c) /
    (p_t * (1 - p_c) + (1 - p_t) * p_c)

  X <- sample_patient(p_c, p_t, N)

  ##########################################################
  # E-process path
  ##########################################################

  e_path <-
    eprocess_run(
      X,
      N,
      alpha,
      lambda_star,
      returnPath = TRUE
    ) |>
    mutate(
      process = "E-process"
    )

  ##########################################################
  # HCP path
  ##########################################################

  hcp_path <-
    HCP(
      m_0 = m_0,
      c = c,
      X = X,
      theta = theta
    ) |>
    mutate(
      process = "HCP"
    )

  bind_rows(
    e_path,
    hcp_path
  ) |>
    ggplot(
      aes(n, value, colour = process)
    ) +
    geom_line(linewidth = 1) +
    geom_hline(
      yintercept = 1 / alpha,
      linetype = 2
    ) +
    scale_y_log10() +
    labs(
      x = "Patient",
      y = "Process value (log scale)",
      title = "Sequential evidence processes"
    ) +
    theme_bw()
}

############################################################
# Example path under alternative
############################################################

plot_sample_path(
  p_c = 0.50,
  p_t = 0.65,
  N = N,
  alpha = alpha,
  m_0 = m_0,
  c = c,
  theta = theta
)

