#-------------------------------------------------------------------------------
# Function to compare power and expected sample size in a scenario where we test the mean of Bernoulli variables

if(FALSE){
  source("~/Desktop/Uni/Speciale/speciale/functions/GS_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/HCP_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/HW_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/UIE_test.R")
}

get_seq_test_comp_RCT_norm <- function(B,
                                       N,
                                       N1,
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

  # Parameters
  m_t_true_grid <- seq(0, 1, by = 0.05)

  # GST
  alphas <- rpact::getDesignGroupSequential(kMax = n_looks,
                                            alpha = alpha,
                                            sided = side,
                                            typeOfDesign = "OF")$criticalValues

  sigmaGS <- sqrt(Sigma[1,1] + Sigma[2,2] - 2*Sigma[1,2])

  # Data sampling function
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

  # If sigma is unknown, we need an alternative approach to calculate Q_n,
  # we follow the approach outlined in Holmes section 11
  if(sigmaUnknown){
    # In this case we use the t-test statistic, and use that it follows a normal distribution
    z_agN <- qt(p = 1 - alpha * gamma / side, df = N - 1)
    z_agN1 <- qt(p = 1 - alpha * gamma / side, df = N1 - 1)

    # We cannot calculate Q_n analytically as we do not have the null distribution, we have to sample from
    # the estimated null

    sample_data_null <- function(N, X, B) {

      sigma_est <- sd(X)

      # Generate all bootstrap samples simultaneously
      X_boot <- matrix(
        rnorm(N * B, mean = 0, sd = sigma_est),
        nrow = N
      )

      # Sum of observed data
      sum_X <- sum(X)

      # Bootstrap means
      totals <- colSums(X_boot) + sum_X
      n_total <- length(X) + N

      # Standard deviations
      SDs <- apply(
        rbind(matrix(X, nrow = length(X), ncol = B), X_boot),
        2,
        sd
      )

      if(side == 2) totals <- abs(totals)

      totals / SDs / sqrt(n_total)
    }

    #sample_data_null <- function(N, X) {
    #  sigma_est <- sd(X)
    #  X_boot <- rnorm(N, 0, sd = sigma_est)
    #  T_N_n <- mean(c(X, X_boot)) / sd(c(X, X_boot)) * sqrt(length(c(X, X_boot)))
    #  if(side == 2)  T_N_n <- abs(T_N_n)
    #  return(T_N_n)
    #}

    calc_q_n <- NULL
  } else {
    # We use the mean of the D's as test statistic
    sigma_D_N <- sqrt(1 / N * (Sigma[1,1] + Sigma[2,2] - 2 * Sigma[1,2]))
    sigma_D_N1 <- sqrt(1 / N1 * (Sigma[1,1] + Sigma[2,2] - 2 * Sigma[1,2]))
    z_agN <- qnorm(p = 1 - alpha * gamma / side, mean = 0, sd = sigma_D_N)
    z_agN1 <- qnorm(p = 1 - alpha * gamma / side, mean = 0, sd = sigma_D_N1)


    sigmas <- sqrt(1 / N ^ 2 * (N - 1:N) * (Sigma[1,1] + Sigma[2,2] - 2 * Sigma[1,2]))
    calc_q_n <- function(X, N, z_ag) {
      meanss <- 1 / N * cumsum(X)
      1 - pnorm(z_ag, meanss, sigmas) + pnorm(- z_ag, meanss, sigmas)
    }

    sample_data_null <- NULL
  }

  #-------------------------------------------------------------------------------
  ## Comparisons
  #-------------------------------------------------------------------------------

  compare_tests <- function(N, z_ag) {

    results <- future.apply::future_lapply(
      seq_along(m_t_true_grid),
      function(g) {
        print(g)
        m_t <- m_t_true_grid[g]
        mu_true <- c(m_t, m)

        HCP_rej <- 0
        HW_rej  <- 0
        UIE_rej <- 0
        GS_rej  <- 0
        HCP_ess <- 0
        HW_ess  <- 0
        UIE_ess <- 0
        GS_ess  <- 0

        bigX <- sample_patient(N * B, mu_true)

        for(b in seq_len(B)) {

            # Data
            X <- bigX[(N * (b - 1) + 1):(b * N),]
            D <- X[, 1] - X[, 2]
            D_scaled <- (D + 5) / 10

            out <- run_HCP_test(
              m_0 = 1 / 2,
              c = c,
              X = D_scaled,
              theta = theta,
              alpha = alpha
            )

            HCP_rej <- HCP_rej + out[1]; HCP_ess <- HCP_ess + out[2]

            out <- run_HW_test(
              N = N,
              X = D,
              calc_q_n = calc_q_n,
              sample_data_null = sample_data_null,
              gamma = gamma,
              quanti = z_ag,
              B = 100
            )

            HW_rej <- HW_rej + out[1]; HW_ess <- HW_ess + out[2]

            out <- UIE_test(
              X = X,
              log_f0 = log_f0,
              log_f1 = log_f1,
              N = N,
              Sigma = Sigma,
              sigmaUnknown = sigmaUnknown,
              m_init = m_init,
              burnin = burnin
            )

            UIE_rej <- UIE_rej + out[1]; UIE_ess <- UIE_ess + out[2]

            out <- gs_run(
              Nmax = N,
              alphas = alphas,
              n_looks = n_looks,
              X = D,
              m_0 = 0,
              side = side,
              sigmaUnknown = sigmaUnknown,
              sigma = if(!sigmaUnknown) sigmaGS else NULL
            )

            GS_rej <- GS_rej + out[1]; GS_ess <- GS_ess + out[2]
        }
        # ----------------------------
        # AVERAGE OVER B
        # ----------------------------
        tibble::tibble(
          m_t_true = m_t,

          HCP_power = HCP_rej / B,
          HCP_ESS   = HCP_ess / B,

          UIE_power = UIE_rej / B,
          UIE_ESS   = UIE_ess / B,

          HW_power  = HW_rej / B,
          HW_ESS    = HW_ess / B,

          GS_power  = GS_rej / B,
          GS_ESS    = GS_ess / B
        )
      },
      future.seed = TRUE
    )

    dplyr::bind_rows(results)
  }

  set.seed(37238493)
  # To parallelize
  future::plan(
    future::multisession,
    workers = parallel::detectCores() - 1
  )

  res <- compare_tests(N = N, z_ag = z_agN)
  res1 <- compare_tests(N = N1, z_ag = z_agN1)
  df <- res
  df1 <- res1

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

  # Power plot
  res$Design <- paste0("N = ", N)
  res1$Design <- paste0("N = ", N1)

  power_df <- bind_rows(res, res1) |>
    pivot_longer(
      cols = c(HCP_power,
               HW_power,
               GS_power,
               UIE_power),
      names_to = "Method",
      values_to = "Power"
    ) |>
    mutate(Method = clean_method_names(Method))

  p2 <- ggplot(power_df,
               aes(m_t_true, Power, colour = Method)) +
    geom_line() +
    facet_wrap(~Design, scales = "free_y") +
    theme_minimal()+
    geom_hline(aes(yintercept = alpha), linetype = 2)+
    labs(x = expression(m[T]))+
    scale_color_manual(values = c("GS" = "darkgreen",
                                  "HCP" = "firebrick",
                                  "HW" = "steelblue",
                                  "UIE" = "orange"),
                       labels = c("GS" = "GS-test",
                                  "HCP" = "HCP-test",
                                  "HW" = "HW-test",
                                  "UIE" = "UIE-test"))

  # ESS plot
  ESS_df <- bind_rows(res, res1) |>
    pivot_longer(
      cols = c(HCP_ESS,
               HW_ESS,
               GS_ESS,
               UIE_ESS),
      names_to = "Method",
      values_to = "ESS"
    ) |>
    mutate(Method = clean_method_names(Method))

  p3 <- ggplot(ESS_df,
               aes(m_t_true, ESS, colour = Method)) +
    geom_line() +
    facet_wrap(~Design, scales = "free_y") +
    theme_minimal()+
    labs(x = expression(m[T])) +
    scale_color_manual(values = c("GS" = "darkgreen",
                                  "HCP" = "firebrick",
                                  "HW" = "steelblue",
                                  "UIE" = "orange"),
                       labels = c("GS" = "GS-test",
                                  "HCP" = "HCP-test",
                                  "HW" = "HW-test",
                                  "UIE" = "UIE-test"))

  return(list(df = df, df1 = df1, p2 = p2, p3 = p3))
}
