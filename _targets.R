# Load packages required to define the pipeline:
library(targets)

# Set target options:
tar_option_set(
  packages = c("tibble", "ggplot2")
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source("functions")

# Replace the target list below with your own:
list(
  tar_target(
    name = OLS_risk,
    command = get_OLS_risk(tau = 0.5, B = 2*10^4, N = 3000, n = 500)),
  tar_target(
    name = OLS_plot_risk,
    command = get_OLS_plot_risk(OLS_risk, tau = 0.5)),
  tar_target(
    name = OLS_rej_prob,
    command = get_OLS_rej_prob(alpha = 0.05, gammas = seq(0.6,0.9, by = 0.05), tau = 0.5, B = 2*10^4, N = 3000, n = 500)),
  tar_target(
    name = OLS_rej_prob_plot,
    command = get_OLS_rej_prob_plot(data1 = OLS_rej_prob, data2 = OLS_risk)),
  tar_target(
    name = l_q_stability_data,
    command = get_l_q_stability_data()),
  tar_target(
    name = l_q_stability_plot,
    command = get_l_q_stability_plot(l_q_stability_data)),
  tar_target(
    name = cv_test_data,
    command = get_cv_test_data(taus = seq(0.85,1.2, by = 0.005), B = 2*10^4, N = 500, K = 0.5)
  )
)
