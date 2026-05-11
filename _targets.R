# Load packages required to define the pipeline:
library(targets)

# Set target options:
tar_option_set(
  packages = c("tibble", "ggplot2", "tidyr")
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source("functions")

# Replace the target list below with your own:
list(
  tar_target(
    name = OLS_risk,
    command = get_OLS_risk(tau = 0.5, B = 2*10^4, N = 2500, n = 500)),
  tar_target(
    name = OLS_plot_risk,
    command = get_OLS_plot_risk(OLS_risk, tau = 0.5)),
  tar_target(
    name = alg_rej_prob,
    command = get_alg_rej_prob(alpha = 0.05, gammas = seq(0.6, 0.8, by = 0.02),
                               tau = 0.5, B = 3*10^4, N = 2500, n = 500,
                               loss = loss_bin, alg = lm_alg)),
  tar_target(
    name = alg_rej_prob_plot,
    command = get_OLS_rej_prob_plot(data1 = alg_rej_prob, data2 = OLS_risk, tau = 0.5, alpha = 0.05, N = 2500, n = 500)),
  tar_target(
    name = l_q_stability_data,
    command = get_l_q_stability_data(B = 500, q = 2, alg = lm_alg, loss = loss_sq_error)),
  tar_target(
    name = l_q_stability_plot,
    command = get_l_q_stability_plot(l_q_stability_data)),
  tar_target(
    name = cv_test_data,
    command = get_cv_test_data(taus = seq(0.85,1.2, by = 0.005), B = 10^4, N = 500, K = 5,
                               loss = loss_bin, alg = lm_alg)),
  tar_target(
    name = cv_test_plot,
    command = get_cv_test_plot(cv_test_data)),
  tar_target(
    name = rina_bb_test_data,
    command = get_rina_bb_test_data(alpha = 0.05, taus = seq(0.3,1, by = 0.01), B = 10^4, N = 2500,
                                    n = 500, loss = loss_bin, alg = lm_alg)),
  tar_target(
    name = rina_bb_test_plot,
    command = get_rina_bb_test_plot(rina_bb_test_data, B = 10^4, N = 2500, n = 500))
)

tar_target(
  name = e_vs_p_power_data,
  command = get_e_vs_p_power(alphas = c(0.001, 0.01, 0.05),
                                 Bs = c(1000, 5000, 10000),
                                 ns = seq(50, 250, by = 50),
                                 ps = seq(0.1, 0.9, by = 0.1),
                                 n = 100,
                                 B = 5000,
                                 alpha = 0.05,
                                 p = 0.7,
                                 lambda = 0.5))

tar_target(
  name = e_vs_p_power_plot,
  command = get_e_vs_p_power_plot(e_vs_p_power_data)
)



