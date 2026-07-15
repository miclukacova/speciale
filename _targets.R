# Load packages required to define the pipeline:
library(targets)

# Set target options:
tar_option_set(
  packages = c("tibble", "ggplot2", "here", "dplyr", "tidyr", "purrr",
               "kableExtra", "rpact", "patchwork", "future.apply")
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source("functions")

# Replace the target list below with your own:
list(
  # Stability
  tar_target(
    name = OLS_risk,
    command = get_OLS_risk(tau = 0.5,
                           B = 2*10^4,
                           n = 500)),
  tar_target(
    name = OLS_plot_risk,
    command = get_OLS_plot_risk(OLS_risk, tau = 0.5)),
  tar_target(
    name = OLS_rej_prob,
    command = get_OLS_rej_prob(alpha = 0.05,
                               ms = seq(0.6, 0.8, by = 0.02),
                               tau = 0.5,
                               B = 3 * 10^4,
                               n = 500,
                               N = 2500)),
  tar_target(
    name = OLS_rej_prob_plot,
    command = get_OLS_rej_prob_plot(data = OLS_risk,
                                    tau = 0.5,
                                    alpha = 0.05,
                                    N = 2500,
                                    n = 500)),
  tar_target(
    name = cv_test_data,
    command = get_cv_test_data(taus = seq(0.85,1.2, by = 0.005),
                               B = 10^4,
                               N = 500, K = 5,
                               loss = loss_bin, alg = lm_alg)),
  tar_target(
    name = cv_test_plot,
    command = get_cv_test_plot(cv_test_data)),
  tar_target(
    name = rina_bb_test_data,
    command = get_rina_bb_test_data(alpha = 0.05,
                                    taus = seq(0.3,1, by = 0.01),
                                    B = 10^4, N = 2500,
                                    n = 500,
                                    loss = loss_bin,
                                    alg = lm_alg)),
  tar_target(
    name = rina_bb_test_plot,
    command = get_rina_bb_test_plot(rina_bb_test_data,
                                    B = 10^4,
                                    N = 2500,
                                    n = 500)),
  # E-variable
  tar_target(
    name = e_vs_p_power_data,
    command = get_e_vs_p_power(alphas = c(0.001, 0.01, 0.05, 0.08),
                               ns = seq(50, 250, by = 50),
                               n = 100,
                               B = 10^4,
                               alpha = 0.05,
                               sample_fct = sample_binom)),
  tar_target(
    name = e_vs_p_power_plot,
    command = get_e_vs_p_power_plot(e_vs_p_power_data)),
  # Neyman-Pearson test compared to the SPRT and the GS
  tar_target(
    name = NP_vs_SPRT_vs_GS,
    command = get_NP_vs_SPRT_vs_GS()
  ),

  # Comparing sequential test in a Bernoulli or RCT setting
  tar_target(
    name = seq_test_comp_RCT_p_t,
    command = get_seq_test_comp_RCT_p_t(B = 3 * 10^4,
                                        N = 100,
                                        N1 = 200)),
  tar_target(
    name = seq_test_comp_RCT_N,
    command = get_seq_test_comp_RCT_N(B = 3 * 10^4)),
  tar_target(
    name = seq_test_comp_bern,
    command = get_seq_test_comp_bern()
  ),
  # Comparing sequential testing in a normal example
  tar_target(
    name = seq_test_comp_RCT_norm,
    command = get_seq_test_comp_RCT_norm(B = 2 * 10^4,
                                         N = 100,
                                         N1 = 200,
                                         Sigma = matrix(c(1,0,0,1), ncol = 2),
                                         side = 2)
  ),
  # Comparing sequential testing in a normal example with unknown variance
  tar_target(
    name = seq_test_comp_RCT_norm_unknown_sigma,
    command = get_seq_test_comp_RCT_norm(B = 10^2,
                                         N = 100,
                                         N1 = 200,
                                         Sigma = matrix(c(1,0.4,0.4,2), ncol = 2),
                                         side = 2,
                                         sigmaUnknown = TRUE)
  )
  )

