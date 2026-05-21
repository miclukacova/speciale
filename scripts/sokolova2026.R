#-------------------------------------------------------------------------------
# Global Settings
#-------------------------------------------------------------------------------

library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)

source("functions/sokolova_functions.R")

tar_load("skolova_data")

#-------------------------------------------------------------------------------
# LaTeX tables for results
#-------------------------------------------------------------------------------

skolova_data$res_200 %>%
  mutate(
    across(where(is.numeric), round, 3)
  ) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    caption = "Simulation results for $N_{\\max}=200$"
  ) %>%
  kable_styling(
    latex_options = c("hold_position")
  )

skolova_data$res_400 %>%
  mutate(
    across(where(is.numeric), round, 3)
  ) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    caption = "Simulation results for $N_{\\max}=200$"
  ) %>%
  kable_styling(
    latex_options = c("hold_position")
  )


skolova_data$res_noMax %>%
  mutate(
    across(where(is.numeric), round, 3)
  ) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    caption = "Simulation results for SPRT without cap"
  ) %>%
  kable_styling(
    latex_options = c("hold_position")
  )

#-------------------------------------------------------------------------------
# Paths
#-------------------------------------------------------------------------------

n_paths <- 12

set.seed(20840)

e_paths <- map_dfr(1:n_paths, ~eprocess_run(p_c = p_c, p_t = p_t, Nmax = Nmax,
                                            alpha = alpha, returnPath = TRUE) %>%
                     mutate(id = .x))

e_paths$method <- "E-process"

sprt_paths <- map_dfr(1:n_paths, ~sprt_run(p_c = p_c, p_t = p_t, Nmax = Nmax,
                                           alpha = alpha, returnPath = TRUE ) %>%
                        mutate(id = .x))

sprt_paths$method <- "SPRT"

all_paths <- bind_rows(e_paths, sprt_paths, gs_paths)

p1 <- ggplot(bind_rows(e_paths, sprt_paths), aes(n, log(value), group = id)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~method) +
  geom_hline(yintercept = log(1 / alpha), linetype = 2) +
  theme_bw() +
  labs(y = "Statistic (log scale)", title = "E-process and SPRT trajectories under H1")+
  ylim(c(-15,15))

e_paths0 <- map_dfr(1:n_paths, ~eprocess_run(p_c = p_c, p_t = p_t, Nmax = Nmax, null = TRUE,
                                            alpha = alpha, returnPath = TRUE) %>%
                     mutate(id = .x))

e_paths0$method <- "E-process"

sprt_paths0 <- map_dfr(1:n_paths, ~sprt_run(p_c = p_c, p_t = p_t, Nmax = Nmax, null = TRUE,
                                           alpha = alpha, returnPath = TRUE ) %>%
                        mutate(id = .x))

sprt_paths0$method <- "SPRT"

p2 <- ggplot(bind_rows(e_paths0, sprt_paths0), aes(n, log(value), group = id)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~method) +
  geom_hline(yintercept = log(1 / alpha), linetype = 2) +
  theme_bw() +
  labs(y = "Statistic (log scale)", title = "E-process and SPRT trajectories under H0")+
  ylim(c(-15,15))

ggsave("plots/sokolova:trajectories:null.pdf", plot = p2)
ggsave("plots/sokolova:trajectories:alt.pdf", plot = p1)

#-------------------------------------------------------------------------------
# Misspecifications
#-------------------------------------------------------------------------------

# Simulation grid

lambda_star <- (p_t * (1 - p_c) - (1 - p_t) * p_c) / (p_t * (1 - p_c) + (1 - p_t) * p_c)

working_grid <- seq(0.35, 0.60, by = 0.025)
B <- 5000

# Simulations

set.seed(3628303)

e_res <- map_dfr(
  working_grid,
  function(p_work) {

    tmp <- replicate(
      B,
      eprocess_run(
        p_c = 0.30,
        p_t = p_work,
        Nmax = 200,
        alpha = 0.025,
        lambda = lambda_star
      ),
      simplify = FALSE
    )

    tmp <- bind_rows(tmp)

    tibble(
      procedure = "E-process",
      p_working = p_work,
      mean_n = mean(tmp$n),
      power = mean(tmp$reject)
    )
  }
)

gs_res <- map_dfr(
  working_grid,
  function(p_work) {

    tmp <- replicate(
      B,
      gs_run(
        p_c = 0.30,
        p_t = p_work,
        Nmax = 200,
        alpha = alphas,
        n_looks = 4
      ),
      simplify = FALSE
    )

    tmp <- bind_rows(tmp)

    tibble(
      procedure = "GS",
      p_working = p_work,
      mean_n = mean(tmp$n),
      power = mean(tmp$reject)
    )
  }
)


sprt_res <- map_dfr(
  working_grid,
  function(p_work) {

    tmp <- replicate(
      B,
      sprt_run(
        p_c = 0.30,
        p_t = 0.45,
        Nmax = 200,
        alpha = 0.025,
        p_t_misspec = p_work
      ),
      simplify = FALSE
    )

    tmp <- bind_rows(tmp)

    tibble(
      procedure = "SPRT",
      p_working = p_work,
      mean_n = mean(tmp$n),
      power = mean(tmp$reject)
    )
  }
)

plot_data <- bind_rows(
  e_res,
  sprt_res,
  gs_res
)  %>%
  pivot_longer(
    cols = c(mean_n, power),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      mean_n = "Expected sample size",
      power  = "Power"
    )
  )

p <- ggplot(
  plot_data,
  aes(
    x = p_working,
    y = value,
    color = procedure
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(
    ~metric,
    scales = "free_y"
  ) +
  theme_bw() +
  labs(
    x = expression(p[t]^"*"),
    y = NULL,
    color = "Procedure"
  )

ggsave("plots/sokolova:misspec.pdf", plot = p)
