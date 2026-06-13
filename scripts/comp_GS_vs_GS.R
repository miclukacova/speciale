library(tidyverse)
library(rpact)

source("~/Desktop/Uni/Speciale/speciale/functions/GS_test.R")

# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------

B <- 10^4
N <- 200

p_t_true_grid <- seq(0.30, 0.62, by = 0.04)

p_c <- 0.30
alpha <- 0.05
m_0 <- 0.5
n_looks <- 4

# Soko design
alphas1 <- alphas_soko(
  p_c = p_c,
  n_looks = n_looks,
  Nmax = N,
  B = B,
  alpha = alpha
)

# rpact O'Brien-Fleming design
design_OF <- rpact::getDesignGroupSequential(
  kMax = n_looks,
  alpha = alpha,
  sided = 1,
  typeOfDesign = "OF"
)

alphas2 <- design_OF$criticalValues

# ------------------------------------------------------------------------------
# Data generating mechanism
# ------------------------------------------------------------------------------

sample_patient <- function(N, p_t) {

  pat_c <- rbinom(N, 1, p_c)
  pat_t <- rbinom(N, 1, p_t)

  (pat_t - pat_c + 1) / 2
}

# ------------------------------------------------------------------------------
# Simulation
# ------------------------------------------------------------------------------

compare_tests <- function(N) {

  results <- vector("list", length(p_t_true_grid))

  for (g in seq_along(p_t_true_grid)) {

    cat("Scenario", g, "of", length(p_t_true_grid), "\n")

    p_t_true <- p_t_true_grid[g]

    GS1 <- matrix(NA_real_, nrow = B, ncol = 2)
    GS2 <- matrix(NA_real_, nrow = B, ncol = 2)

    for (i in seq_len(B)) {

      X <- sample_patient(N, p_t_true)

      GS1[i, ] <- gs_run(
        Nmax = N,
        alphas = alphas1,
        n_looks = n_looks,
        X = X,
        m_0 = m_0
      )

      GS2[i, ] <- gs_run(
        Nmax = N,
        alphas = alphas2,
        n_looks = n_looks,
        X = X,
        m_0 = m_0
      )
    }

    colnames(GS1) <- c("Reject", "ESS")
    colnames(GS2) <- c("Reject", "ESS")

    results[[g]] <- tibble(
      p_t_true = p_t_true,

      GS1_power = mean(GS1[, "Reject"]),
      GS1_ESS   = mean(GS1[, "ESS"]),

      GS2_power = mean(GS2[, "Reject"]),
      GS2_ESS   = mean(GS2[, "ESS"])
    )
  }

  bind_rows(results)
}

# ------------------------------------------------------------------------------
# Run simulation
# ------------------------------------------------------------------------------

set.seed(37238493)

res <- compare_tests(N)

# ------------------------------------------------------------------------------
# Power plot
# ------------------------------------------------------------------------------

power_df <- res %>%
  pivot_longer(
    cols = c(GS1_power, GS2_power),
    names_to = "Method",
    values_to = "Power"
  ) %>%
  mutate(
    Method = recode(
      Method,
      GS1_power = "GS (Soko)",
      GS2_power = "GS (rpact OF)"
    )
  )

p_power <- ggplot(
  power_df,
  aes(x = p_t_true, y = Power, colour = Method)
) +
  geom_line(linewidth = 1) +
  theme_minimal() +
  labs(
    x = expression(p[t]),
    y = "Power",
    colour = NULL
  )

p_power

# ------------------------------------------------------------------------------
# ESS plot
# ------------------------------------------------------------------------------

ESS_df <- res %>%
  pivot_longer(
    cols = c(GS1_ESS, GS2_ESS),
    names_to = "Method",
    values_to = "ESS"
  ) %>%
  mutate(
    Method = recode(
      Method,
      GS1_ESS = "GS (Soko)",
      GS2_ESS = "GS (rpact OF)"
    )
  )

p_ESS <- ggplot(
  ESS_df,
  aes(x = p_t_true, y = ESS, colour = Method)
) +
  geom_line(linewidth = 1) +
  theme_minimal() +
  labs(
    x = expression(p[t]),
    y = "Expected Sample Size",
    colour = NULL
  )

p_ESS
