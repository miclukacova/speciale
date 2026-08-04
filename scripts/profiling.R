profvis::profvis({
  gs_run <- function(Nmax,
                     alphas,
                     n_looks,
                     X,
                     m_0,
                     side = 1,
                     sigmaUnknown = TRUE,
                     sigma = NULL) {

    # Information times / looks
    look_times <- round(seq(Nmax / n_looks, Nmax, length.out = n_looks))

    cs  <- cumsum(X)
    cs2 <- cumsum(X^2)

    Xbar_k <- cs[look_times] / look_times

    if (sigmaUnknown) {
      var_k <- (cs2[look_times] - cs[look_times]^2 / look_times) /
        (look_times - 1)
      sigma_hat_k <- sqrt(var_k)
    } else {
      sigma_hat_k <- rep(sigma, length(look_times))
    }

    # Z_k_star
    Z_k_star <- (Xbar_k - m_0) / sigma_hat_k * sqrt(look_times)

    # If we do not know sigma, we need to correct with the t distribution
    if(sigmaUnknown) alphas <- qt(p = pnorm(alphas), df = look_times - 1)

    # Boundary crossing
    if(side == 2) test_res <- abs(Z_k_star) >= alphas
    else test_res <- Z_k_star >= alphas

    # Sample size
    n <-  look_times[min(which(test_res), n_looks)]

    c(Reject = any(test_res), ESS = n)

  }

  # Obrian/Pocock boundaries

  alphas_soko <- function(p_c, n_looks, Nmax, B = 10^4, alpha = 0.05) {

    # Define quantities
    look_times <- round(seq(Nmax / n_looks, Nmax, length.out = n_looks))
    info_frac <- look_times/Nmax

    # Sample data
    xT_null_mat <- matrix(stats::rbinom(B * Nmax, 1, p_c),nrow = B)
    xC_null_mat <- matrix(stats::rbinom(B * Nmax, 1, p_c),nrow = B)

    # Sum number of sucesses in each stage for all B replications
    cumT_null <- t(apply(xT_null_mat, 1, cumsum))[, look_times, drop = FALSE]
    cumC_null <- t(apply(xC_null_mat, 1, cumsum))[, look_times, drop = FALSE]

    # Number of observations in each stage for all B replications
    nn_looks <- matrix(look_times, nrow = B, ncol = length(look_times), byrow = TRUE)

    # Estimated effect sizes divided by standard error
    z_null <- evalinger:::.z_matrix(cumT_null, cumC_null, nn_looks)

    # We find the maximal effect size across stages for each of the B trials
    m_null <- apply(z_null * matrix(sqrt(info_frac), nrow = B, ncol = length(info_frac), byrow = TRUE), 1, max)
    # Take the quantile of the maximal effect sizes
    gs_c <- as.numeric(stats::quantile(m_null, probs = 1 - alpha, names = FALSE))
    obf_bounds <- gs_c / sqrt(info_frac)

    return(obf_bounds)
  }

  # Calculating the predictable sequence
  calculate_lambda_tilde <- function(X, alpha){
    # Number of observations
    n <- length(X)
    # mu_t hat
    mu_hat <- (1 /2 + cumsum(X)) / (1:n +1)
    # sigma_t hat
    sigma_hat <- (1 / 4 + cumsum((X - mu_hat) ^ 2)) / (1:n + 1)

    # We shift by 1 to restore predictability
    sigma_hat_t_1 <- c(1 / 4, sigma_hat)[1:n]
    # lambda tilde calculation
    lambda_tilde <- sqrt(2 * log(2 / alpha) / (sigma_hat_t_1 * 1:n * log(1:n + 1)))

    lambda_tilde

  }

  # The hedged capital process
  HCP <- function(m_0, c, X, theta, alpha) {

    # Calculate the lambda's
    lambda_tilde <- calculate_lambda_tilde(X, alpha = alpha)
    lambda_plus <- pmin(abs(lambda_tilde), c / m_0)
    lambda_minus <- pmin(abs(lambda_tilde), c / (1 - m_0))

    # Calculate the capital process
    K_plus <- cumprod(1 + lambda_plus * (X - m_0))
    K_minus <- cumprod(1 - lambda_minus * (X - m_0))

    return(pmax(theta * K_plus, (1 - theta) * K_minus))
  }

  # The HCP test
  run_HCP_test <- function(m_0, c, X, theta, alpha) {
    HCP_res <- HCP(m_0, c, X, theta, alpha)
    test_res <- HCP_res >= 1 / alpha
    #test_fut <- HCP_res < alpha / 2
    if(any(test_res)){         #| any(test_fut)
      ESS <- which((test_res) == 1)[1]           # + test_fut
    } else {
      ESS <- length(HCP_res)
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
        boot_data <- sample_data_null(
          N = N - i,
          X = X[1:i],
          B = B
        )
        Q_n[i] <- mean(boot_data >= quanti)

        #boot_data <- numeric(B)
        #boot_data <- numeric(B)
        #for(b in seq_len(B))
        #  boot_data[b] <- sample_data_null(
        #    N = N - i,
        #    X = X[1:i]
        #  )
        #Q_n[i] <- mean(boot_data >= quanti)
        i <- i + 1
      }

      Q_n[N] <- sample_data_null(N = N, X = X, B) > quanti

    } else {
      Q_n <- calc_q_n(X, N, quanti)
    }

    Reject = any(Q_n >= gamma)

    if(Reject) ESS <- which(Q_n >= gamma)[1]
    else ESS <- N

    if(return_Q_n) return(Q_n)

    return(c(Reject = Reject, ESS = ESS))
  }

  # The universal inference process
  UIE <- function(X, log_f0, log_f1, N, Sigma, sigmaUnknown, m_init, burnin) {
    # Storage
    logE <- vector(length = N)
    f1s <- vector(length = N)

    # Estimates
    # Mean
    i <- 1:N
    mT_est  <- cumsum(X[,1]) / i
    mC_est  <- cumsum(X[,2]) / i
    m0_est <- (cumsum(X) / (1:(2 * N)))[seq(2, 2 * N, by = 2)]

    # Variance covariance matrix
    if(sigmaUnknown){
      # Moments
      mT2_est <- cumsum(X[,1]^2) / i
      mC2_est <- cumsum(X[,2]^2) / i
      mTC_est <- cumsum(X[,1] * X[,2]) / i

      # Variance covariance unbiased
      s2_1 <- (i / (i - 1)) * (mT2_est - mT_est^2)
      s2_2 <- (i / (i - 1)) * (mC2_est - mC_est^2)
      s2_12 <- (i / (i - 1)) * (mTC_est -  mT_est * mC_est)

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

      if(sigmaUnknown && i > 2){
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
      logE[i] <- sum(f1s[burnin:i] - f0s[burnin:i])
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

  get_seq_test_comp_RCT_norm <- function(B,
                                         N,
                                         N1,
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

    #-------------------------------------------------------------------------------
    ## Comparisons
    #-------------------------------------------------------------------------------

    compare_tests <- function(N, z_ag) {

      for(m_t in seq_along(m_t_true_grid)){
          print(m_t)
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
              sigmaUnknown = FALSE,
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
              sigmaUnknown = FALSE,
              sigma = sigmaGS
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
      }
    }

    res <- compare_tests(N = N, z_ag = z_agN)
  }


  get_seq_test_comp_RCT_norm(B = 100,
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
                             n_looks = 20)

})
