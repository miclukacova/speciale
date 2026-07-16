library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

mu0 <- 0
theta_nu <- 0.2
B <- 100

set.seed(1)

# Simulate under nu
X <- rnorm(B, theta_nu, 1)

# Plug-in estimator
mu_hat <- cumsum(X) / seq_len(B)

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
regret <- (logE - logE_nu) / seq_len(B)

D <- tibble(
  t = 1:B,
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
      "logE_nu" = "black"
    ),
    labels = c(
      "Adaptive",
      "Oracle"
    )
  ) +
  theme_minimal(base_size = 15)

p1

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

p2




mu0 <- 0.5
mu1 <- 0.6
B <- 100

set.seed(1)


