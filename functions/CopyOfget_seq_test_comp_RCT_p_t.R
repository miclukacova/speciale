#-------------------------------------------------------------------------------
# Function to compare power and expected sample size in a scenario where we test the mean of Bernoulli variables

if(FALSE){
  source("~/Desktop/Uni/Speciale/speciale/functions/GS_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/HCP_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/HW_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/SPRT.R")
}

get_seq_test_comp_RCT_p_t <- function(B,
                                      N = 100,
                                      p_c,
                                      m_0,
                                      c,
                                      theta,
                                      alpha,
                                      gamma,
                                      n_looks) {

  # Parameters
  p_t_true_grid <- seq(0.3, 0.625, by = 0.025)
  sprt_grid <- c(0.45, 0.6)

  # GST
  alphas <- rpact::getDesignGroupSequential(kMax = n_looks,
                                            alpha = alpha,
                                            sided = 1,
                                            typeOfDesign = "OF")$criticalValues

  # The variance of the D_n's under the null (the variance is completely determined by the mean of X_n^T)
  sigma <- sqrt(1 / 2 * (p_c * (1 - p_c)))

  # Data sampling function
  sample_patient <- function(N, p_t) {
    pat_c <- rbinom(N, 1, p_c)
    pat_t <- rbinom(N, 1, p_t)
    (pat_t - pat_c + 1) / 2
  }


  # SPRT
  f0 <- function(X) p_c * (1 - p_c) * (X == 1 | X == 0) + (p_c * p_c + (1 - p_c) ^ 2) * (X == 1 / 2)

  f1_adap <- function(X){
    N <- length(X)
    cum_mean <- c(0.5, cumsum(X[-N]) / seq_len(N-1))
    p_t <- 2 * cum_mean - 1 + p_c
    p_t <- pmin(pmax(p_t, 0.3), 0.999)
    p_t * (1 - p_c) * (X == 1) +  (1 - p_t) * p_c * (X == 0) + (X == (1 / 2)) * ((1 - p_t) * (1 - p_c) + p_t * p_c)
  }

  # List of
  f1_list <- lapply(sprt_grid, function(p_alt){

    function(X){

      ind1 <- (X == 1)
      ind0 <- (X == 0)
      indh <- (X == 1/2)

      ind1 * p_alt * (1 - p_c) +
        ind0 * (1 - p_alt) * p_c +
        indh * ((1 - p_alt) * (1 - p_c) + p_alt * p_c)
    }

  })

  # -----------------------------------------------------------
  # Precompute HW critical values and Q_n calculating function
  # -----------------------------------------------------------

  cat("Precomputing critical values...\n")

  # Function to calculate quantiles under the null
  null_pmf <- function(N){

    probs <- c(p_c * (1 - p_c), p_c^2 + (1 - p_c)^2, p_c * (1 - p_c))
    pmf <- probs

    if(N > 1){
      for(i in 2:N){
        pmf <- convolve(pmf, rev(probs), type="open")
      }
    }
    pmf
  }

  # Precompute the pmf and tail probabilities in order to avoid repeated computations
  pmf_lookup <- vector("list", N)
  tail_lookup <- vector("list", N)
  for(k in 1:N){
    pmf <- null_pmf(k)
    pmf_lookup[[k]] <- pmf

    # We compute the tail probabilities P(X >= x[k])
    tail_lookup[[k]] <- rev(cumsum(rev(pmf)))
  }

  # Function to find quantiles
  z_ag_exact <- function(N){
    cdf <- cumsum(null_pmf(N))
    k <- which(cdf >= 1 - alpha * gamma)[1]
    (k - 1) / 2
  }
  z_ag <- z_ag_exact(N)

  # Function for calculating Q_n
  calc_q_n <- function(X, N, z_ag){

    Q_n <- numeric(N)
    sumX <- cumsum(X)

    for(i in seq_len(N-1)){
      s_obs <- sumX[i]
      # We find the threshold
      threshold <- z_ag - s_obs
      # We find the index in the support of sum_{n+1}^N X_i corresponding to the threshold
      idx <- threshold * 2 + 1
      # If the index exceeds the support, there is 0 probability that the fixed sample size test will reject
      if(idx > (N - i) * 2 + 1) break
      # We calculate the tail probability
      tail_prob <- tail_lookup[[N-i]]
      Q_n[i] <- tail_prob[idx]
      if(Q_n[i] >= gamma) break
    }

    Q_n[N] <- as.numeric(sum(X) >= z_ag)

    Q_n
  }

  #-------------------------------------------------------------------------------
  ## Comparisons
  #-------------------------------------------------------------------------------

  compare_tests <- function(N, z_ag) {

    results <- vector("list", length(p_t_true_grid))

    for (g in seq_along(p_t_true_grid)) {

      print(g)

      p_t_true <- p_t_true_grid[g]

      sample_data <- function(N) sample_patient(N, p_t_true)

      # Simulations
      sim_results <- future_lapply(
        seq_len(B),
        function(b) {

          X <- sample_data(N)

          HCP <- run_HCP_test(
            m_0 = m_0, c = c, X = X,
            theta = theta, alpha = alpha
          )

          HW <- run_HW_test(
            N = N, X = X,
            calc_q_n = calc_q_n,
            gamma = gamma,
            quanti = z_ag
          )

          SPRT_0.45 <- run_sprt_test(
            N, X = X,
            f0 = f0,
            f1 = f1_list[[1]],
            gamma0 = alpha/(1-alpha),
            gamma1 = (1-alpha)/alpha
          )

          SPRT_0.6 <- run_sprt_test(
            N, X = X,
            f0 = f0,
            f1 = f1_list[[2]],
            gamma0 = alpha/(1-alpha),
            gamma1 = (1-alpha)/alpha
          )

          SPRT_adap <- run_sprt_test(
            N, X = X,
            f0 = f0,
            f1 = f1_adap,
            gamma0 = 0,
            gamma1 = (1-alpha)/alpha
          )

          GS <- gs_run(
            Nmax = N,
            alphas = alphas,
            n_looks = n_looks,
            X = X,
            m_0 = m_0,
            sigmaUnknown = FALSE,
            sigma = sigma
          )


          list(
            HCP = HCP,
            HW = HW,
            GS = GS,
            SPRT_0.45 = SPRT_0.45,
            SPRT_0.6 = SPRT_0.6,
            SPRT_adap = SPRT_adap
          )
        },
        future.seed = TRUE
      )


      # Convert simulation output into matrices
      get_results <- function(name) {

        x <- do.call(
          rbind,
          lapply(sim_results, `[[`, name)
        )

        colnames(x) <- c("Reject", "ESS")

        x
      }


      HCP_res <- get_results("HCP")
      HW_res <- get_results("HW")
      GS_res <- get_results("GS")
      SPRT_0.45_res <- get_results("SPRT_0.45")
      SPRT_0.6_res <- get_results("SPRT_0.6")
      SPRT_adap_res <- get_results("SPRT_adap")


      # Helper function for power
      power_summary <- function(x) {

        p <- mean(x[, "Reject"])

        se <- sqrt(
          p * (1-p) / B
        )

        c(
          power = p,
          lower = max(0, p - 1.96*se),
          upper = min(1, p + 1.96*se)
        )
      }


      # Helper function for ESS
      ESS_summary <- function(x) {

        c(
          ESS = mean(x[, "ESS"]),
          lower = quantile(x[, "ESS"], 0.025),
          upper = quantile(x[, "ESS"], 0.975)
        )
      }


      out <- list(
        p_t_true = p_t_true,


        # Power
        HCP_power = power_summary(HCP_res)[1],
        HCP_power_lower = power_summary(HCP_res)[2],
        HCP_power_upper = power_summary(HCP_res)[3],

        HW_power = power_summary(HW_res)[1],
        HW_power_lower = power_summary(HW_res)[2],
        HW_power_upper = power_summary(HW_res)[3],

        GS_power = power_summary(GS_res)[1],
        GS_power_lower = power_summary(GS_res)[2],
        GS_power_upper = power_summary(GS_res)[3],


        SPRT_0.45_power = power_summary(SPRT_0.45_res)[1],
        SPRT_0.45_power_lower = power_summary(SPRT_0.45_res)[2],
        SPRT_0.45_power_upper = power_summary(SPRT_0.45_res)[3],

        SPRT_0.6_power = power_summary(SPRT_0.6_res)[1],
        SPRT_0.6_power_lower = power_summary(SPRT_0.6_res)[2],
        SPRT_0.6_power_upper = power_summary(SPRT_0.6_res)[3],

        SPRT_adap_power = power_summary(SPRT_adap_res)[1],
        SPRT_adap_power_lower = power_summary(SPRT_adap_res)[2],
        SPRT_adap_power_upper = power_summary(SPRT_adap_res)[3],



        # ESS
        HCP_ESS = ESS_summary(HCP_res)[1],
        HCP_ESS_lower = ESS_summary(HCP_res)[2],
        HCP_ESS_upper = ESS_summary(HCP_res)[3],

        HW_ESS = ESS_summary(HW_res)[1],
        HW_ESS_lower = ESS_summary(HW_res)[2],
        HW_ESS_upper = ESS_summary(HW_res)[3],

        GS_ESS = ESS_summary(GS_res)[1],
        GS_ESS_lower = ESS_summary(GS_res)[2],
        GS_ESS_upper = ESS_summary(GS_res)[3],

        SPRT_0.45_ESS = ESS_summary(SPRT_0.45_res)[1],
        SPRT_0.45_ESS_lower = ESS_summary(SPRT_0.45_res)[2],
        SPRT_0.45_ESS_upper = ESS_summary(SPRT_0.45_res)[3],

        SPRT_0.6_ESS = ESS_summary(SPRT_0.6_res)[1],
        SPRT_0.6_ESS_lower = ESS_summary(SPRT_0.6_res)[2],
        SPRT_0.6_ESS_upper = ESS_summary(SPRT_0.6_res)[3],

        SPRT_adap_ESS = ESS_summary(SPRT_adap_res)[1],
        SPRT_adap_ESS_lower = ESS_summary(SPRT_adap_res)[2],
        SPRT_adap_ESS_upper = ESS_summary(SPRT_adap_res)[3]
      )


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

  res <- compare_tests(N = N, z_ag = z_ag)
  df <- res

  # Power plot
  power_bounds <- res |>
    pivot_longer(
      cols = matches("_power_(lower|upper)$"),
      names_to = c("Method", "bound"),
      names_pattern = "(.*)_power_(lower|upper)",
      values_to = "value"
    ) |>
    mutate(Method = paste0(Method, "_power")) |>
    pivot_wider(
      names_from = bound,
      values_from = value
    )

  power_df <- res |>
    pivot_longer(
      cols = matches("_power$"),
      names_to = "Method",
      values_to = "Power"
    ) |>
    left_join(
      power_bounds,
      by = c("p_t_true", "Method")
    )


  p2 <- ggplot(power_df,
              aes(p_t_true, Power,
                  colour = Method,
                  fill = Method)) +

    geom_ribbon(
      aes(ymin = lower, ymax = upper),
      alpha = 0.15,
      colour = NA
    ) +

    geom_line(size = 1) +

    theme_minimal() +

    labs(
      x = expression(p[t]),
      y = "Power"
    ) +

    geom_hline(
      yintercept = alpha,
      linetype = 2
    )

  # ESS plot
  ESS_df <- res |> pivot_longer(cols = c(ends_with("_ESS")),
                                names_to = "Method",values_to = "ESS")

  p3 <- ggplot(ESS_df, aes(p_t_true, ESS, colour = Method)) +
    geom_line() +
    theme_minimal()+
    labs(x = expression(p[t])) +
    scale_color_manual(values = c("GS_ESS" = "darkgreen",
                                  "HCP_ESS" = "firebrick",
                                  "HW_ESS" = "steelblue",
                                  "SPRT_0.45_ESS" = "orange",
                                  "SPRT_0.6_ESS" = "goldenrod",
                                  "SPRT_adap_ESS" = "#009E73"),
                       labels = c("GS_ESS" = "GS-test",
                                  "HCP_ESS" = "HCP-test",
                                  "HW_ESS" = "HW-test",
                                  "SPRT_0.45_ESS" = "SPRT(0.45)",
                                  "SPRT_0.6_ESS" = "SPRT(0.6)",
                                  "SPRT_adap_ESS" = "SPRT(adap)"))

  return(list(df = df, p2 = p2, p3 = p3))
}
