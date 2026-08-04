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

  # -----------------------------------------------------------
  # Precompute HW critical values and Q_n calculating function
  # -----------------------------------------------------------
  cat("Precomputing critical values...\n")

  # We use the mean of the D's as test statistic
  sigma_D_N <- sqrt(1 / N_grid * (Sigma[1,1] + Sigma[2,2] - 2 * Sigma[1,2]))
  z_agN <- qnorm(p = 1 - alpha * gamma / side, mean = 0, sd = sigma_D_N)
  N_max <- max(N_grid)
  sigmas_lookup <- sapply(N_grid, function(N) sqrt(1 / N ^ 2 * (N - 1:N) * (Sigma[1,1] + Sigma[2,2] - 2 * Sigma[1,2])))

  calc_q_n <- function(X, N, z_ag) {
    meanss <- 1 / N * cumsum(X)
    sigmass <- sigmas_lookup[[which(N == N_grid)]]
    1 - pnorm(z_ag, meanss, sigmass) + pnorm(- z_ag, meanss, sigmass)
  }

  sample_data_null <- NULL

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

    # Sample all patients
    mu_true <- c(m_T, m)
    maxN <- max(N_grid)
    bigX <- sample_patient(maxN * B, mu_true)

    results <- future.apply::future_lapply(
      seq_along(N_grid),
      function(g) {
        print(g)
        N <- N_grid[g]

        cat("m_T =", m_T,
            " N =", N, "\n")

        HCP_rej <- 0
        HW_rej  <- 0
        UIE_rej <- 0
        GS_rej  <- 0
        HCP_ess <- 0
        HW_ess  <- 0
        UIE_ess <- 0
        GS_ess  <- 0

        for(b in seq_len(B)){

          # Find data
          first <- (b - 1) * maxN + 1
          last  <- first + N - 1
          X <- bigX[first:last,]

          out <- run_HCP_test(
            m_0 = 1 / 2,
            c = c,
            X = (X[,1] - X[,2] + 5) / 10,
            theta = theta,
            alpha = alpha
            )

          HCP_rej <- out[1] + HCP_rej; HCP_ess <- out[2] + HCP_ess

          out <- run_HW_test(
            N = N,
            X = X[,1] - X[,2],
            calc_q_n = calc_q_n,
            sample_data_null = sample_data_null,
            gamma = gamma,
            quanti = z_agN[N == N_grid],
            B = B
            )

          HW_rej <- out[1] + HW_rej; HW_ess <- out[2] + HW_ess

          out <- UIE_test(
            X = X,
            log_f0 = log_f0,
            log_f1 = log_f1,
            N = N,
            Sigma = Sigma,
            sigmaUnknown = FALSE,
            m_init = m_init,
            burnin = burnin
            )

          UIE_rej <- out[1] + UIE_rej; UIE_ess <- out[2] + UIE_ess

          out <- gs_run(
            Nmax = N,
            alphas = alphas,
            n_looks = n_looks,
            X = X[,1] - X[,2],
            m_0 = 0,
            side = side,
            sigmaUnknown = FALSE,
            sigma = sigmaGS
            )

          GS_rej <- out[1] + GS_rej; GS_ess <- out[2] + GS_ess
        }

        # ----------------------------
        # AVERAGE OVER B
        # ----------------------------
        tibble::tibble(
          N = N,

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
                                  "UIE" = "UIE-test"))

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


