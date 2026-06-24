#-------------------------------------------------------------------------------
# Function to compare power and expected sample size in a scenario where we test the mean of Bernoulli variables

if(FALSE){
  source("~/Desktop/Uni/Speciale/speciale/functions/GS_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/HCP_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/HW_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/UIE_test.R")
}

get_seq_test_comp_RCT_norm <- function(B = 500,
                                       N = 100,
                                       N1 = 200,
                                       Sigma = matrix(c(1,0,0,1), ncol = 2),
                                       side = 2,
                                       sigmaUnknown = FALSE,
                                       burnin = 1,
                                       m_init = 0.3) {

  # Parameters
  m_t_true_grid <- seq(0.1, 0.7, by = 0.05)
  m <- 0.3
  c = 3 / 4
  theta = 1 / 2
  alpha = 0.05
  gamma = 0.9
  n_looks = 4

  # GST
  alphas <- rpact::getDesignGroupSequential(kMax = n_looks,
                                            alpha = alpha,
                                            sided = side,
                                            typeOfDesign = "OF")$criticalValues

  # Data sampling function
  sample_patient <- function(N, m) {
    X <- MASS::mvrnorm(n = N, mu = m, Sigma = Sigma)
    X
  }

  # Estimator of density in the null
  log_f0 <- function(X, mu_est, sigma_est) {
    mvtnorm::dmvnorm(X, mean = c(mu_est, mu_est), sigma = sigma_est, log = TRUE)
  }

  # Estimator of density in the alternative
  log_f1 <- function(X, mT_est, mC_est, sigma_est) {
    mvtnorm::dmvnorm(X, mean = c(mT_est, mC_est), sigma = sigma_est, log = TRUE)
  }

  # -----------------------------------------------------------
  # Precompute HW critical values and Q_n calculating function
  # -----------------------------------------------------------

  # If sigma is unknown, we need an alternative approach to calculate Q_n, we follow the approach outlined
  # in Holmes section 11
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
    results <- vector("list", length(m_t_true_grid))
    for (g in seq_along(m_t_true_grid)) {
      print(g)
      # True parameter
      mu_true <- c(m_t_true_grid[g], m)
      sample_data <- function(N) sample_patient(N, mu_true)

      # Storage
      HCP_res <- matrix(NA_real_, nrow = B, ncol = 2)
      UIE_res <- HW_res <- GS_res <- HCP_res

      # Simulations
      sim_results <- future_lapply(
        seq_len(B),
        function(b) {

          X <- sample_data(N)

          # HCP
          HCP_res <- run_HCP_test(
            m_0 = 1 / 2,
            c = c,
            X = (X[,1] - X[,2] + 5) / 10,
            theta = theta,
            alpha = alpha
          )

          # HW
          HW <- run_HW_test(
            N = N,
            X = X[,1] - X[,2],
            calc_q_n = calc_q_n,
            sample_data_null,
            gamma = gamma,
            quanti = z_ag,
            B = B
          )

          # UIE
          UIE_res <- UIE_test(
            X = X,
            log_f0 = log_f0,
            log_f1 = log_f1,
            N = N,
            Sigma = Sigma,
            sigmaUnknown = sigmaUnknown,
            m_init = m_init,
            burnin = burnin
            )

          # GS
          if(sigmaUnknown) sigmaGS <- TRUE
          else sigmaGS <- sqrt(1 / N * (Sigma[1,1] + Sigma[2,2] - 2 * Sigma[1,2]))
          GS <- gs_run(
            Nmax = N,
            alphas = alphas,
            n_looks = n_looks,
            X = X[, 1] - X[, 2],
            m_0 = 0,
            side = side,
            sigmaUnknown = sigmaGS
          )

          list(
            HCP = HCP_res,
            HW = HW,
            UIE = UIE_res,
            GS = GS
          )
        },
        future.seed = TRUE
      )

      HCP_res <- matrix(
        unlist(lapply(sim_results, `[[`, "HCP")),
        ncol = 2,
        byrow = TRUE
      )
      colnames(HCP_res) <- c("Reject", "ESS")


      HW_res <- matrix(
        unlist(lapply(sim_results, `[[`, "HW")),
        ncol = 2,
        byrow = TRUE
      )
      colnames(HW_res) <- c("Reject", "ESS")

      UIE_res <- matrix(
        unlist(lapply(sim_results, `[[`, "UIE")),
        ncol = 2,
        byrow = TRUE
      )
      colnames(UIE_res) <- c("Reject", "ESS")

      GS_res <- matrix(
        unlist(lapply(sim_results, `[[`, "GS")),
        ncol = 2,
        byrow = TRUE
      )
      colnames(GS_res) <- c("Reject", "ESS")

      # ----------------------------
      # Summaries
      # ----------------------------

      out <- list(m_t_true = m_t_true_grid[g],
                  HCP_power = mean(HCP_res[, "Reject"]),
                  HCP_ESS   = mean(HCP_res[, "ESS"]),
                  UIE_power = mean(UIE_res[, "Reject"]),
                  UIE_ESS   = mean(UIE_res[, "ESS"]),
                  HW_power = mean(HW_res[, "Reject"]),
                  HW_ESS   = mean(HW_res[, "ESS"]),
                  GS_power = mean(GS_res[, "Reject"]),
                  GS_ESS   = mean(GS_res[, "ESS"]))

      results[[g]] <- tibble::as_tibble(out)
    }

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
    labs(x = expression(m[T]))

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
    labs(x = expression(m[T]))

  return(list(df = df, df1 = df1, p2 = p2, p3 = p3))
}
