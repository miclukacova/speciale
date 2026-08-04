# Load packages required to define the pipeline:
library(targets)

# Global settings
B <- 5 * 10^4
alpha <- 0.05

# Colors
colss <- c("steelblue", "firebrick", "darkgreen")

# Set target options:
tar_option_set(
  packages = c("tibble", "ggplot2", "here", "dplyr", "tidyr", "purrr",
               "kableExtra", "rpact", "patchwork", "future.apply")
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source("functions")

# Replace the target list below with your own:
list(
  #------------------------------
  # Stability and existence
  #------------------------------
  # MC estimates of the risk of the OLS algorithm in our particular setting
  tar_target(
    name = OLS_risk,
    command = get_OLS_risk(tau = 0.5,
                           B = B,
                           n = 500)),

  # Plots of OLS risk
  tar_target(
    name = OLS_plot_risk,
    command = get_OLS_plot_risk(OLS_risk, tau = 0.5)),

  # Rejection probability of Rina BB test in OLS setting
  tar_target(
    name = OLS_rej_prob,
    command = get_OLS_rej_prob(alpha = alpha,
                               ms = seq(0.6, 0.8, by = 0.02),
                               tau = 0.5,
                               B = B,
                               n = 500,
                               N = 2500)),

  # Plots of rejection probability of Rina BB test
  tar_target(
    name = OLS_rej_prob_plot,
    command = get_OLS_rej_prob_plot(data = OLS_risk,
                                    data2 = OLS_rej_prob,
                                    tau = 0.5,
                                    alpha = alpha,
                                    N = 2500,
                                    n = 500)),

  # CV test
  tar_target(
    name = cv_test_data,
    command = get_cv_test_data(taus = seq(0.85,1.2, by = 0.005),
                               B = B,
                               N = 500,
                               K = 5,
                               loss = loss_sq_error,
                               alg = lm_alg)),

  # CV test plot
  tar_target(
    name = cv_test_plot,
    command = get_cv_test_plot(cv_test_data,
                               B = B,
                               alg = lm_alg,
                               loss = loss_sq_error,
                               N = 500)),

  # BB test data
  tar_target(
    name = rina_bb_test_data,
    command = get_rina_bb_test_data(alpha = alpha,
                                    taus = seq(0.35,1, by = 0.02),
                                    B = B,
                                    N = 500,
                                    n = 400,
                                    loss = loss_bin,
                                    alg = lm_alg)),

  # BB test plot
  tar_target(
    name = rina_bb_test_plot,
    command = get_rina_bb_test_plot(rina_bb_test_data,
                                    taus = seq(0.3,1, by = 0.01),
                                    alpha = alpha,
                                    B = B,
                                    N = 500,
                                    n = 400)),

  #---------------------------------------------
  # Sequential tests
  #---------------------------------------------

  # E-variable and p-variable comparision: used in the motivational example
  tar_target(
    name = e_vs_p_power_data,
    command = get_e_vs_p_power(alphas = c(0.001, 0.01, 0.05, 0.08),
                               ns = seq(50, 250, by = 50),
                               n = 100,
                               B = B,
                               alpha = alpha,
                               sample_fct = sample_binom)),

  # Plot for e-variable and p-variable comparision: used in the motivational example
  tar_target(
    name = e_vs_p_power_plot,
    command = get_e_vs_p_power_plot(e_vs_p_power_data,
                                    alpha = alpha)),

  # Neyman-Pearson test compared to the SPRT and the GS (in the appendix)
  tar_target(
    name = NP_vs_SPRT_vs_GS,
    command = get_NP_vs_SPRT_vs_GS()),

  # Assymptotic log optimal e-process
  tar_target(
    name = ass_log_opt_e_proc,
    command = get_ass_log_opt_e_proc(mu0 = 0,
                                     theta_nu = 0.2,
                                     N = 100)),

  #---------------------------------------------
  # Comparing sequential test in an RCT setting.
  #---------------------------------------------

  # Function outputs plot of power as function of p_T.
  tar_target(
    name = seq_test_comp_RCT_p_t,
    command = get_seq_test_comp_RCT_p_t(B = B,
                                        N = 100,
                                        p_c = 0.3,
                                        m_0 = 1 / 2,
                                        c = 3 / 4,
                                        theta = 1,
                                        alpha = alpha,
                                        gamma = 0.9,
                                        n_looks = 20)),

  # Function outputs plot of power as function of N.
  tar_target(
    name = seq_test_comp_RCT_N,
    command = get_seq_test_comp_RCT_N(B = B,
                                      p_c = 0.3,
                                      m_0 = 1 / 2,
                                      c = 3 / 4,
                                      theta = 1,
                                      alpha = alpha,
                                      gamma = 0.9,
                                      n_looks = 20)),

  # Function outputs plot of power as function of N.
  tar_target(
    name = comp_soko_eproc,
    command = get_comp_soko_eproc(alpha,
                                  B,
                                  N = 200,
                                  p_c = 0.3,
                                  p_t = 0.45,
                                  c = 3 / 4,
                                  theta = 1,
                                  p_t_ms = 0.6)),


  #---------------------------------------------
  # Comparing sequential testing in a gaussian setting
  #---------------------------------------------
  # Unknown variance

  ## Function outputs power and ESS comparisons as function of M_T
  tar_target(
    name = seq_test_comp_RCT_norm,
    command = get_seq_test_comp_RCT_norm(B = 10^4,
                                         N = 100,
                                         N1 = 200,
                                         Sigma = matrix(c(1,0,0,1), ncol = 2),
                                         side = 2,
                                         burnin = 1,
                                         m_init = 0.3,
                                         m = 0.3,
                                         c = 3 / 4,
                                         theta = 1 / 2,
                                         alpha = alpha,
                                         gamma = 0.9,
                                         n_looks = 20)),

  ## Function outputs power and ESS comparisons as function of N
  tar_target(
    name = seq_test_comp_RCT_norm_N,
    command = get_seq_test_comp_RCT_norm_N(B = 10^4,
                                           Sigma = matrix(c(1,0,0,1), ncol = 2),
                                           side = 2,
                                           burnin = 1,
                                           m_init = 0.3,
                                           m = 0.3,
                                           c = 3 / 4,
                                           theta = 1 / 2,
                                           alpha = alpha,
                                           gamma = 0.9,
                                           n_looks = 20)),

  #---------------------------------------------
  # Unknown variance
  ## Function outputs power and ESS comparisons as function of M_T
  tar_target(
    name = seq_test_comp_RCT_norm_unknown,
    command = get_seq_test_comp_RCT_norm_unknown(B = 2000,
                                                 N = 150,
                                                 Sigma = matrix(c(1,0.4,0.4,2), ncol = 2),
                                                 side = 2,
                                                 burnin = 10,
                                                 m_init = 0.3,
                                                 m = 0.3,
                                                 c = 3 / 4,
                                                 theta = 1 / 2,
                                                 alpha = alpha,
                                                 gamma = 0.9,
                                                 n_looks = 20,
                                                 sample_patient = sample_patient_norm)),

  ## Function outputs power and ESS comparisons as function of N
  tar_target(
    name = seq_test_comp_RCT_norm_unknown_N,
    command = get_seq_test_comp_RCT_norm_N_unknown(B = 500,
                                                   Sigma = matrix(c(1,0.4,0.4,2), ncol = 2),
                                                   side = 2,
                                                   burnin = 10,
                                                   m_init = 0.3,
                                                   m = 0.3,
                                                   c = 3 / 4,
                                                   theta = 1 / 2,
                                                   alpha = alpha,
                                                   gamma = 0.9,
                                                   n_looks = 20,
                                                   sample_data = sample_patient_norm)),
  #---------------------------------------------
  # Sampling from contaminated normal
  ## Function outputs power and ESS comparisons as function of M_T
  tar_target(
    name = seq_test_comp_RCT_norm_unknown_mis1,
    command = get_seq_test_comp_RCT_norm_unknown(B = 100,
                                                 N = 150,
                                                 Sigma = matrix(c(1,0.4,0.4,2), ncol = 2),
                                                 side = 2,
                                                 burnin = 10,
                                                 m_init = 0.3,
                                                 m = 0.3,
                                                 c = 3 / 4,
                                                 theta = 1 / 2,
                                                 alpha = alpha,
                                                 gamma = 0.9,
                                                 n_looks = 20,
                                                 sample_patient = sample_patient_contaminated)),

  ## Function outputs power and ESS comparisons as function of N
  tar_target(
    name = seq_test_comp_RCT_norm_unknown_N_mis1,
    command = get_seq_test_comp_RCT_norm_N_unknown(B = 100,
                                                   Sigma = matrix(c(1,0.4,0.4,2), ncol = 2),
                                                   side = 2,
                                                   burnin = 10,
                                                   m_init = 0.3,
                                                   m = 0.3,
                                                   c = 3 / 4,
                                                   theta = 1 / 2,
                                                   alpha = alpha,
                                                   gamma = 0.9,
                                                   n_looks = 20,
                                                   sample_patient = sample_patient_contaminated)),

  #---------------------------------------------
  # Sampling from log normal
  ## Function outputs power and ESS comparisons as function of M_T
  tar_target(
    name = seq_test_comp_RCT_norm_unknown_mis2,
    command = get_seq_test_comp_RCT_norm_unknown(B = 1000,
                                                 N = 150,
                                                 Sigma = matrix(c(1,0.4,0.4,2), ncol = 2),
                                                 side = 2,
                                                 burnin = 10,
                                                 m_init = 0.3,
                                                 m = 0.3,
                                                 c = 3 / 4,
                                                 theta = 1 / 2,
                                                 alpha = alpha,
                                                 gamma = 0.9,
                                                 n_looks = 20,
                                                 sample_patient = sample_patient_lognormal)),

  ## Function outputs power and ESS comparisons as function of N
  tar_target(
    name = seq_test_comp_RCT_norm_unknown_N_mis2,
    command = get_seq_test_comp_RCT_norm_N_unknown(B = 100,
                                                   Sigma = matrix(c(1,0.4,0.4,2), ncol = 2),
                                                   side = 2,
                                                   burnin = 10,
                                                   m_init = 0.3,
                                                   m = 0.3,
                                                   c = 3 / 4,
                                                   theta = 1 / 2,
                                                   alpha = alpha,
                                                   gamma = 0.9,
                                                   n_looks = 20,
                                                   sample_patient = sample_patient_lognormal))



  )

