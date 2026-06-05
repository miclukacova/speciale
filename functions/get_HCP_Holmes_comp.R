#-------------------------------------------------------------------------------
# Function to compare power and expected sample size in a scenario where we test the mean of Bernoulli variables

get_HCP_Holmes_comp <- function(alpha = 0.05, m_0 = 0.5, m_1 = 0.65, N = 300, c = 1 / 2,
                                theta = 3 / 4, B = 10^4, gamma = 0.9, p_c = 0.3, p_t = 0.45) {

  set.seed(37238493)

  # The optimal N for the Holmes test
  N_opt_Holmes <- 50
  pow <- 0

  while(pow < 1 - alpha) {
    z_ag <- qbinom(p = 1 - (alpha * gamma), size = N_opt_Holmes, prob = m_0)
    pow <- pbinom(z_ag, N_opt_Holmes, m_1, lower.tail = FALSE)
    N_opt_Holmes <- N_opt_Holmes + 1
  }

  #-------------------------------------------------------------------------------
  ## The hedged capital process
  #-------------------------------------------------------------------------------


  # Calculating the predictable sequence
  calculate_lambda_tilde <- function(X, alpha){

    n <- length(X)
    mu_hat <- (1 /2 + cumsum(X)) / (1:n +1)
    sigma_hat <- (1 / 4 + cumsum((X - mu_hat) ^ 2)) / (1:n + 1)
    sigma_hat_t_1 <- c(1 / 4, sigma_hat)[1:n]
    lambda_tilde <- sqrt(2 * log(2 / alpha) / (sigma_hat_t_1 * 1:n * log(1:n + 1)))
    lambda_tilde

  }

  # Data sampling functions
  sample_patient <- function(p_c, p_t, n) {
    pat_c <- rbinom(n,1,p_c)
    pat_t <- rbinom(n,1,p_t)
    (pat_t - pat_c + 1) / 2
  }


  # The hedged capital process
  hedged_cap_proc <- function(m_0, c, sample_data, N, theta) {

    # Sample data
    X <- sample_data(N)

    # Calculate the lambda's
    lambda_tilde <- calculate_lambda_tilde(X, alpha = alpha)
    lambda_plus <- pmin(abs(lambda_tilde), c / m_0)
    lambda_minus <- pmin(abs(lambda_tilde), c / (1 -m_0))

    # Calculate the capital process
    K_plus <- cumprod(1 + lambda_plus * (X - m_0))
    K_minus <- cumprod(1 - lambda_minus * (X - m_0))

    return(pmax(theta * K_plus, (1 - theta) * K_minus))
  }

  HCP_test <- function(HCP, alpha) {
    test_res <- HCP > 1/ alpha
    test_fut <- HCP < alpha
    if(any(test_res)){
      ESS <- which((test_res + test_fut) == 1)[1]
    } else {
      ESS <- length(HCP)
    }
    return(c(Reject = any(test_res[1:ESS]), ESS = ESS))
  }



  # Plot the process
  simulate_paths <- function(p_c, p_t, label) {

    sample_data <- function(N) sample_patient(p_c, p_t, N)

    HCP_mat <- replicate(
      10,
      hedged_cap_proc(
        m_0 = m_0,
        c = c,
        sample_data = sample_data,
        N = N,
        theta = theta
      )
    )

    data.frame(
      time = rep(1:N, 10),
      path = rep(1:10, each = N),
      HCP = as.vector(HCP_mat),
      scenario = label
    )
  }

  df_alt <- simulate_paths(p_c, p_t, "Alternative (m_true ≠ m_0)")
  df_null <- simulate_paths(p_c, p_c, "Null (m_true = m_0)")

  df <- rbind(df_alt, df_null)

  p1 <- ggplot(df, aes(time, HCP, group = path)) +
    geom_line(alpha = 0.4) +
    geom_hline(yintercept = 1 / alpha,
               colour = "red",
               linetype = 2) +
    geom_hline(yintercept = alpha,
               colour = "red",
               linetype = 2) +
    scale_y_log10() +
    theme_minimal() +
    facet_wrap(~ scenario, ncol = 2)

  #-------------------------------------------------------------------------------
  ## SPRT
  #-------------------------------------------------------------------------------

  # Jeg havde to ideer til SPRT'en: 1) lave mixture, men det svarer jo egentlig bare til at vælge et p?
  # 2) lave den adaptive process, men den er lidt mærkelig i det her scenarie, fordi at den ender med at sige
  # at alternativ hypotesen slet ikke er mulig. Derfor gør vi lidt noget andet...

  beta <- alpha
  gamma0 <- beta / (1 - alpha)
  gamma1 <- (1 - alpha) / beta

  sprt_test <- function(N, p_c, p_t, sample_data) {

    X <- sample_data(N)
    f_1 <- p_t * (1 - p_c) * (X == 1) +  (1 - p_t) * p_c * (X == 0) + (X == 1 / 2) * ((1 - p_t) * (1 - p_c) + p_t * p_c)
    f_0 <-  p_c * (1 - p_c) * (X == 1 | X == 0) + (p_c * p_c + (1 - p_c) ^ 2) * (X == 1 / 2)

    L <- cumprod(f_1 / f_0)
    ss <- min(which(L > gamma1 | gamma0 > L), N)

    return(c(L[ss] > gamma1, ss))
  }

  # Adaptive version

  sprt_test_adap <- function(N, sample_data, p_c) {
    X <- sample_data(N)
    # Plug-in estimator
    p_t_est <- c(p_t, 2 * cumsum(X[1:(N-1)]) / 1:(N-1) - 1 + p_c)
    p_t_est[p_t_est < 0] <- 0.1
    # E-process
    f_1 <- p_t_est * (1 - p_c) * (X == 1) +  (1 - p_t_est) * p_c * (X == 0) +
      (X == 1 / 2) * ((1 - p_t_est) * (1 - p_c) + p_t_est * p_c)
    f_0 <-  p_c * (1 - p_c) * (X == 1 | X == 0) + (p_c * p_c + (1 - p_c) ^ 2) * (X == 1 / 2)
    L <- cumprod(f_1 / f_0)
    # Sample size
    ss <- min(which(L > gamma1 | gamma0 > L), N)
    return(c(L[ss] > gamma1, ss))
  }

  sprt_test_adap(N=300, p_c = 0.3, sample_data = sample_data)

  #-------------------------------------------------------------------------------
  ## Holmes method
  #-------------------------------------------------------------------------------

  # Vi bliver nødt til at simulere for at få fat i kvantilerne

  holmes_test <- function(N, m_true, m_0, gamma) {

    z_ag <- qbinom(p = 1 - (alpha * gamma), size = N, prob = m_0)

    X <- sample_bern(N, m_true)
    Q_n <- 1 - pbinom(q = z_ag - cumsum(X), size = seq(N-1, 0), prob = m_0)
    Reject = any(Q_n >= gamma)

    if(Reject) ESS <- which(Q_n >= gamma)[1]
    else ESS <- N

    return(c(Reject = Reject, ESS = ESS))
  }

  #-------------------------------------------------------------------------------
  ## Comparisons
  #-------------------------------------------------------------------------------

  compare_tests <- function(m_true_grid,
                            m1_grid,
                            gamma,
                            N,
                            B) {

    results <- vector("list", length(m_true_grid))

    for (g in seq_along(m_true_grid)) {

      m_true <- m_true_grid[g]

      # ----------------------------
      # Storage
      # ----------------------------

      HCP_res <- matrix(
        NA_real_,
        nrow = B,
        ncol = 2,
        dimnames = list(NULL, c("Reject", "ESS"))
      )

      HOLM_res <- matrix(
        NA_real_,
        nrow = B,
        ncol = 2,
        dimnames = list(NULL, c("Reject", "ESS"))
      )

      n_sprt <- length(m1_grid)

      SPRT_res <- matrix(
        NA_real_,
        nrow = B,
        ncol = 2 * n_sprt + 2
      )

      colnames(SPRT_res) <- c(
        as.vector(rbind(
          paste0("Reject_", m1_grid),
          paste0("ESS_", m1_grid)
        )),
        "Reject_adap",
        "ESS_adap"
      )

      # ----------------------------
      # Simulations
      # ----------------------------

      for (b in seq_len(B)) {

        # HCP
        HCP <- hedged_cap_proc(
          m_0 = m_0,
          m_true = m_true,
          c = c,
          sample_data = sample_bern,
          N = N,
          theta = theta
        )

        HCP_res[b, ] <- HCP_test(HCP, alpha)

        # HOLMES
        HOLM_res[b, ] <- holmes_test(
          N = N,
          m_true = m_true,
          m_0 = m_0,
          gamma = gamma
        )

        # Fixed SPRTs
        for (j in seq_along(m1_grid)) {

          SPRT_res[b, (2 * j - 1):(2 * j)] <-
            sprt_test(
              N = N,
              m_true = m_true,
              m_0 = m_0,
              m_1 = m1_grid[j]
            )
        }

        # Adaptive SPRT
        SPRT_res[b, (2 * n_sprt + 1):(2 * n_sprt + 2)] <-
          sprt_test_adap(
            N = N,
            m_true = m_true,
            m_0 = m_0
          )
      }

      # ----------------------------
      # Summaries
      # ----------------------------

      out <- list(
        m_true = m_true,

        HCP_power = mean(HCP_res[, "Reject"]),
        HCP_ESS   = mean(HCP_res[, "ESS"]),

        HOLM_power = mean(HOLM_res[, "Reject"]),
        HOLM_ESS   = mean(HOLM_res[, "ESS"])
      )

      for (j in seq_along(m1_grid)) {

        out[[paste0("SPRT_power_", m1_grid[j])]] <-
          mean(SPRT_res[, paste0("Reject_", m1_grid[j])])

        out[[paste0("SPRT_ESS_", m1_grid[j])]] <-
          mean(SPRT_res[, paste0("ESS_", m1_grid[j])])
      }

      out$SPRT_power_adap <-
        mean(SPRT_res[, "Reject_adap"])

      out$SPRT_ESS_adap <-
        mean(SPRT_res[, "ESS_adap"])

      results[[g]] <- tibble::as_tibble(out)
    }

    dplyr::bind_rows(results)
  }

  res <- compare_tests(
    m_true_grid = seq(0.5, 0.8, by = 0.03),
    m1_grid = c(0.55, 0.75),
    gamma = 0.95,
    N = N,
    B = B
  )

  res_n_opt_holmes <- compare_tests(
    m_true_grid = seq(0.5, 0.8, by = 0.03),
    m1_grid = c(0.6, 0.75),
    gamma = 0.95,
    N = N_opt_Holmes,
    B = B
  )

  df1 <- res
  df2 <- res_n_opt_holmes

  clean_method_names <- function(x) {

    x <- gsub("SPRT_ESS_", "SPRT(", x)
    x <- gsub("SPRT_power_", "SPRT(", x)
    x <- gsub("$", ")", x)

    x <- gsub("HCP_ESS)", "HCP", x)
    x <- gsub("HOLM_ESS)", "HOLM", x)
    x <- gsub("HCP_power)", "HCP", x)
    x <- gsub("HOLM_power)", "HOLM", x)

    x
  }

  # Power plot
  res$Design <- paste0("N = ", N)
  res_n_opt_holmes$Design <- paste0("N = ", N_opt_Holmes)

  power_df <- bind_rows(res, res_n_opt_holmes) |>
    pivot_longer(
      cols = c(HCP_power,
               HOLM_power,
               starts_with("SPRT_power")),
      names_to = "Method",
      values_to = "Power"
    ) |>
    mutate(Method = clean_method_names(Method))

  p2 <- ggplot(power_df,
         aes(m_true, Power, colour = Method)) +
    geom_line() +
    facet_wrap(~Design, scales = "free_y") +
    theme_minimal()

  # ESS plot
  ESS_df <- bind_rows(res, res_n_opt_holmes) |>
    pivot_longer(
      cols = c(HCP_ESS,
               HOLM_ESS,
               starts_with("SPRT_ESS")),
      names_to = "Method",
      values_to = "ESS"
    ) |>
    mutate(Method = clean_method_names(Method))

  p3 <- ggplot(ESS_df,
         aes(m_true, ESS, colour = Method)) +
    geom_line() +
    facet_wrap(~Design, scales = "free_y") +
    theme_minimal()

  return(list(df1 = df1, df2 = df2, p1 = p1, p2 = p2, p3 = p3))
}
