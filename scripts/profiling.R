profvis::profvis({
  # The universal inference process
  UIE <- function(X, log_f0, log_f1, N, Sigma, sigmaUnknown, m_init, burnin) {
    # Storage
    logE <- vector(length = N)
    f1s <- vector(length = N)

    # Estimates
    # Mean
    mT_est  <- cumsum(X[,1]) / 1:N
    mC_est  <- cumsum(X[,2]) / 1:N
    m0_est <- (cumsum(X) / (1:(2 * N)))[seq(2, 2 * N, by = 2)]

    # Variance covariance matrix
    if(sigmaUnknown){
      # Moments
      mT2_est <- cumsum(X[,1]^2) / 1:N
      mC2_est <- cumsum(X[,2]^2) / 1:N
      mTC_est <- cumsum(X[,1] * X[,2]) / 1:N

      # Variance covariance
      s2_1 <- mT2_est - (mT_est)^2
      s2_2 <- mC2_est - (mC_est)^2
      s2_12 <- mTC_est - mT_est * mC_est

      # Initialize
      sigma_est <- diag(2)
      sigma_est_1 <-  diag(2)

    } else {
      sigma_est <- Sigma
      sigma_est_1 <- Sigma
    }

    # Since the f1 process should be predictable we add an initial mT and mC value
    mT_est  <- c(m_init, mT_est)
    mC_est  <- c(m_init, mC_est)

    # Calculate the e-process
    for(i in burnin:N){

      if(sigmaUnknown & i > 2){
        # Update variance estimate
        sigma_est_1 <- sigma_est
        sigma_est <- matrix(c(s2_1[i], s2_12[i], s2_12[i], s2_2[i]), ncol = 2)
      }

      # Compute the predictable density estimate evaluated in the newest observation
      f1s[i] <- log_f1(X = X[i,,drop = FALSE],
                       mT_est = mT_est[i],
                       mC_est = mC_est[i],
                       sigma_est = sigma_est_1)

      f0s <- log_f0(X = X[1:i,,drop = FALSE],
                    mu_est = m0_est[i],
                    sigma_est = sigma_est)

      if(!any(is.finite(f0s))) {
        print(i)
        stop("log(f0) is infinite")
      }

      # Compute the log e-values
      logE[i] <- sum(f1s[1:i] - f0s)
    }

    if(!any(is.finite(f1s))) stop("log(f1) is infinite")
    logE
  }

  UIE_test <- function(X, log_f0, log_f1, N, Sigma, sigmaUnknown = FALSE, m_init = 0.3, burnin = 1) {
    logUIE <- UIE(X = X,
                  log_f0 = log_f0,
                  log_f1 = log_f1,
                  N = N,
                  Sigma = Sigma,
                  sigmaUnknown = sigmaUnknown,
                  m_init = m_init,
                  burnin = burnin)

    test_res <- logUIE > log(1 / alpha)
    #test_fut <- UIE_res < alpha / 2

    if(any(test_res)){         #| any(test_fut)
      ESS <- which((test_res) == 1)[1]           # + test_fut
    } else {
      ESS <- N
    }

    return(c(Reject = any(test_res[1:ESS]), ESS = ESS))
  }


  run_HW_test <- function(N,
                          X,
                          sample_data_null = NULL,
                          calc_q_n = NULL,
                          gamma,
                          quanti,
                          B = 1000,
                          return_Q_n = FALSE) {

    if(is.null(calc_q_n)) {

      Q_n <- vector(length = N)
      i <- 2
      while(i <= (N-1) & all(Q_n < gamma)) {
        boot_data <- replicate(B, sample_data_null(N = N - i, X = X[1:i]))
        Q_n[i] <- mean(boot_data >= quanti)
        i <- i + 1
      }
      Q_n[N] <- sample_data_null(N = N, X = X) > quanti

    } else {
      Q_n <- calc_q_n(X, N, quanti)
    }

    Reject = any(Q_n >= gamma)

    if(Reject) ESS <- which(Q_n >= gamma)[1]
    else ESS <- N

    if(return_Q_n) return(Q_n)

    return(c(Reject = Reject, ESS = ESS))
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
        sim_results <- lapply(
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
          })

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
    res <- compare_tests(N = N, z_ag = z_agN)
  }


  get_seq_test_comp_RCT_norm(B = 15,
                             N = 100,
                             N1 = 200,
                             Sigma = matrix(c(1,0.4,0.4,2), ncol = 2),
                             side = 2,
                             sigmaUnknown = TRUE)

})
