#-------------------------------------------------------------------------------
# Function to compare power and expected sample size in a scenario where we test the mean of Bernoulli variables

if(FALSE){
  source("~/Desktop/Uni/Speciale/speciale/functions/GS_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/HCP_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/HW_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/UIE_test.R")
}

get_seq_test_comp_RCT_norm_N <- function(B,
                                         Sigma,
                                         side,
                                         sigmaUnknown,
                                         burnin,
                                         m_init,
                                         m,
                                         c,
                                         theta,
                                         alpha,
                                         gamma,
                                         n_looks) {
  # -----------------------------------------------------------
  # Parameters
  # -----------------------------------------------------------
  N_grid <- seq(50, 400, by = 10)
  m_t_true_values <- c(0.6, 0.8)

  # -----------------------------------------------------------
  # Data sampling function
  # -----------------------------------------------------------
  sample_patient <- function(N, m) {
    X <- MASS::mvrnorm(n = N, mu = m, Sigma = Sigma)
    X
  }

  # -----------------------------------------------------------
  # UIE
  # -----------------------------------------------------------
  # Estimator of density in the null
  if(sigmaUnknown){
    log_f0 <- function(X, mu_est, sigma_est) {
      mvtnorm::dmvnorm(X, mean = c(mu_est, mu_est), sigma = sigma_est, log = TRUE)
    }

    # Estimator of density in the alternative
    log_f1 <- function(X, mT_est, mC_est, sigma_est) {
      mvtnorm::dmvnorm(X, mean = c(mT_est, mC_est), sigma = sigma_est, log = TRUE)
    }
  } else {
    # In case Sigma does not change, we can compute the density faster
    R <- chol(Sigma)
    Sigma_inv <- chol2inv(R)
    log_det <- 2 * sum(log(diag(R)))
    d <- nrow(Sigma)

    log_dmvnorm_fast <- function(X, mu){
      Z <- sweep(X, 2, mu)
      quad <- rowSums((Z %*% Sigma_inv) * Z)
      -0.5 * (d * log(2*pi) + log_det + quad)
    }

    # Estimator of density in the alternative
    log_f1 <- function(X, mT_est, mC_est, sigma_est) {
      log_dmvnorm_fast(X, mu = c(mT_est, mC_est))
    }
    log_f0 <- function(X, mu_est, sigma_est) {
      log_dmvnorm_fast(X, mu = c(mu_est, mu_est))
    }
  }

  # -----------------------------------------------------------
  # Precompute HW critical values and Q_n calculating function
  # -----------------------------------------------------------
  cat("Precomputing critical values...\n")

  # If sigma is unknown, we need an alternative approach to calculate Q_n,
  # we follow the approach outlined in Holmes section 11
  if(sigmaUnknown){
    # In this case we use the t-test statistic, and use that it follows a normal distribution
    z_agN <- qt(p = 1 - alpha * gamma / side, df = N - 1)
    z_agN1 <- qt(p = 1 - alpha * gamma / side, df = N1 - 1)

    # We cannot calculate Q_n analytically as we do not have the null distribution, we have to sample from
    # the estimated null
    sample_data_null <- function(N, X) {
      sigma_est <- sd(X)
      X_boot <- rnorm(N, 0, sd = sigma_est)
      T_N_n <- mean(c(X, X_boot)) / sd(c(X, X_boot)) * sqrt(length(c(X, X_boot)))
      if(side == 2)  T_N_n <- abs(T_N_n)
      return(T_N_n)
    }

    calc_q_n <- NULL
  } else {
    # We use the mean of the D's as test statistic
    sigma_D_N <- sqrt(1 / N_grid * (Sigma[1,1] + Sigma[2,2] - 2 * Sigma[1,2]))
    z_agN <- qnorm(p = 1 - alpha * gamma / side, mean = 0, sd = sigma_D_N)
    N_max <- max(N_grid)
    sigmas_lookup <- sqrt(1 / N_max ^ 2 * (N_max - 1:N_max) * (Sigma[1,1] + Sigma[2,2] - 2 * Sigma[1,2]))

    calc_q_n <- function(X, N, z_ag) {

      meanss <- 1 / N * cumsum(X)
      1 - pnorm(z_ag, meanss, sigmas_lookup[1:N]) + pnorm(- z_ag, meanss, sigmas_lookup[1:N])
    }

    sample_data_null <- NULL
  }

  # -----------------------------------------------------------
  # GS critical values and sigma if not unknown
  # -----------------------------------------------------------

  alphas <- rpact::getDesignGroupSequential(kMax = n_looks,
                                            alpha = alpha,
                                            sided = side,
                                            typeOfDesign = "OF")$criticalValues

  sigmaGS <- sqrt(Sigma[1,1] + Sigma[2,2] - 2*Sigma[1,2])

  # -------------------------------------------------
  # Main comparison function
  # -------------------------------------------------

  compare_tests <- function(m_T) {
    mu_true <- c(m_T, m)

    results <- future.apply::future_lapply(
      seq_along(N_grid),
      function(g) {
        print(g)
        N <- N_grid[g]

        HCP_mat <- matrix(0, B, 2)
        HW_mat  <- matrix(0, B, 2)
        UIE_mat <- matrix(0, B, 2)
        GS_mat  <- matrix(0, B, 2)

        for(b in seq_len(B)){
          X <- sample_patient(N, mu_true)

          HCP_mat[b,] <- run_HCP_test(
            m_0 = 1 / 2,
            c = c,
            X = (X[,1] - X[,2] + 5) / 10,
            theta = theta,
            alpha = alpha
            )

          HW_mat[b,] <- run_HW_test(
            N = N,
            X = X[,1] - X[,2],
            calc_q_n = calc_q_n,
            sample_data_null = sample_data_null,
            gamma = gamma,
            quanti = z_agN[N == N_grid],
            B = B
            )

          UIE_mat[b,] <- UIE_test(
            X = X,
            log_f0 = log_f0,
            log_f1 = log_f1,
            N = N,
            Sigma = Sigma,
            sigmaUnknown = sigmaUnknown,
            m_init = m_init,
            burnin = burnin
            )

          GS_mat[b,] <- gs_run(
            Nmax = N,
            alphas = alphas,
            n_looks = n_looks,
            X = X[,1] - X[,2],
            m_0 = 0,
            side = side,
            sigmaUnknown = sigmaUnknown,
            sigma = if(!sigmaUnknown) sigmaGS else NULL
            )
          }

        # ----------------------------
        # AVERAGE OVER B
        # ----------------------------
        tibble::tibble(
          N = N,

          HCP_power = mean(HCP_mat[, 1]),
          HCP_ESS   = mean(HCP_mat[, 2]),

          UIE_power = mean(UIE_mat[, 1]),
          UIE_ESS   = mean(UIE_mat[, 2]),

          HW_power  = mean(HW_mat[, 1]),
          HW_ESS    = mean(HW_mat[, 2]),

          GS_power  = mean(GS_mat[, 1]),
          GS_ESS    = mean(GS_mat[, 2])
        )
      },
      future.seed = TRUE
    )

    dplyr::bind_rows(results)
  }

  # -------------------------------------------------
  # Run simulations
  # -------------------------------------------------

  set.seed(37238493)

  future::plan(
    future::multisession,
    workers = parallel::detectCores() - 1
  )

  res6 <- compare_tests(m_t_true_values[1])
  res8 <- compare_tests(m_t_true_values[2])

  # -------------------------------------------------
  # Clean names
  # -------------------------------------------------

  clean_method_names <- function(x) {

    x <- gsub("HCP_ESS", "HCP", x)
    x <- gsub("HW_ESS", "HW", x)
    x <- gsub("GS_ESS", "GS", x)
    x <- gsub("UIE_ESS", "UIE", x)

    x <- gsub("HCP_power", "HCP", x)
    x <- gsub("HW_power", "HW", x)
    x <- gsub("GS_power", "GS", x)
    x <- gsub("UIE_power", "UIE", x)

    x
  }


  res6$Scenario <- "m_T = 0.6"
  res8$Scenario <- "m_T = 0.8"

  # -------------------------------------------------
  # Power plot
  # -------------------------------------------------

  power_df <- bind_rows(res6, res8) |>
    pivot_longer(
      cols = c(
        HCP_power,
        HW_power,
        GS_power,
        UIE_power
      ),
      names_to = "Method",
      values_to = "Power"
    ) |>
    mutate(
      Method = clean_method_names(Method))

  p_power <- ggplot(power_df, aes(N, Power, colour = Method)) +
    geom_line() +
    facet_wrap(~Scenario, scales = "free") +
    theme_minimal() +
    scale_color_manual(values = c("GS" = "darkgreen",
                                  "HCP" = "firebrick",
                                  "HW" = "steelblue",
                                  "UIE" = "orange"),
                       labels = c("GS" = "GS-test",
                                  "HCP" = "HCP-test",
                                  "HW" = "HW-test",
                                  "UIE" = "UIE-test"))+
    scale_y_continuous(limits = c(0, 1))

  # -------------------------------------------------
  # ESS plot
  # -------------------------------------------------

  ESS_df <- bind_rows(res6, res8) |>
    pivot_longer(
      cols = c(
        HCP_ESS,
        HW_ESS,
        UIE_ESS,
        GS_ESS),
      names_to = "Method",
      values_to = "ESS") |>
    mutate(Method = clean_method_names(Method))

  p_ESS <- ggplot(ESS_df,aes(N, ESS, colour = Method)) +
    geom_line() +
    facet_wrap(~Scenario, scales = "free") +
    theme_minimal() +
    scale_color_manual(values = c("GS" = "darkgreen",
                                  "HCP" = "firebrick",
                                  "HW" = "steelblue",
                                  "UIE" = "orange"),
                       labels = c("GS" = "GS-test",
                                  "HCP" = "HCP-test",
                                  "HW" = "HW-test",
                                  "UIE" = "UIE-test"))

  return(list(
    res6 = res6,
    res8 = res8,
    p_power = p_power,
    p_ESS = p_ESS
  ))
}


