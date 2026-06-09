#-------------------------------------------------------------------------------
# Function to compare power and expected sample size in a scenario where we test the mean of Bernoulli variables

get_seq_test_comp_bern <- function() {

  # Parameters
  m_true_grid <- seq(0.5, 0.8, by = 0.03)
  m1_grid <- c(0.6, 0.75)
  alpha = 0.05
  m_true = 0.6
  m_0 = 0.5
  m_1 = 0.65
  N = 300
  c = 1 / 2
  theta = 3 / 4
  B = 10^4
  gamma = 0.9

  # The optimal N for the Holmes test
  N1 <- 50
  pow <- 0
  while(pow < 1 - alpha) {
    z_ag <- qbinom(p = 1 - (alpha * gamma), size = N1, prob = m_0)
    pow <- pbinom(z_ag, N1, m_1, lower.tail = FALSE)
    N1 <- N1 + 1
  }

  ## SPRT
  f0 <- function(x) dbinom(x, 1,  m_0)
  f1_adap <- function(x) {
    # OBS: estimate has to be predictable
    m_est <- cumsum(c(0.5, x[1:(N-1)])) / seq_len(N)
    dbinom(x, 1,  m_est)
  }


  #-------------------------------------------------------------------------------
  ## Comparisons
  #-------------------------------------------------------------------------------

  compare_tests <- function(N) {

    results <- vector("list", length(m_true_grid))
    for (g in seq_along(m_true_grid)) {

      # True parameter
      m_true <- m_true_grid[g]
      sample_data <- function(N) rbinom(N, 1, m_true)

      # Storage
      HCP_res <- matrix(NA_real_,
                        nrow = B,
                        ncol = 2,
                        dimnames = list(NULL, c("Reject", "ESS")))
      HOLM_res <- HCP_res
      n_sprt <- length(m1_grid)
      SPRT_res <- matrix(NA_real_,
                         nrow = B,
                         ncol = 2 * n_sprt + 2)

      colnames(SPRT_res) <- c(as.vector(rbind(paste0("Reject_", m1_grid), paste0("ESS_", m1_grid))),
                              "Reject_adap", "ESS_adap")

      # Simulations
      for (b in seq_len(B)) {

        # HCP
        HCP_res[b, ] <- run_HCP_test(m_0 = m_0,
                                     c = c,
                                     sample_data = sample_data,
                                     N = N,
                                     theta = theta,
                                     alpha)

        # HOLMES
        z_ag <- qbinom(p = 1 - (alpha * gamma), size = N, prob = m_0)
        calc_q_n <- function(X, N) {
          1 - pbinom(q = z_ag - cumsum(X), size = seq(N-1, 0), prob = m_0)
        }
        HOLM_res[b, ] <- holmes_test(N = N,
                                     sample_data_true = sample_data,
                                     calc_q_n = calc_q_n,
                                     gamma = gamma,
                                     quanti = z_ag,
                                     B = 1000)

        # Fixed SPRTs
        for (j in seq_along(m1_grid)) {
          f1 <- function(x) dbinom(x, 1,  m1_grid[j])

          SPRT_res[b, (2 * j - 1):(2 * j)] <- sprt_test(N = N,
                                                        sample_data = sample_data,
                                                        f0 = f0,
                                                        f1 = f1,
                                                        beta = alpha,
                                                        alpha = alpha)
        }

        # Adaptive SPRT
        SPRT_res[b, (2 * n_sprt + 1):(2 * n_sprt + 2)] <- sprt_test(N,
                                                                    sample_data = sample_data,
                                                                    f0 = f0,
                                                                    f1 = f1_adap,
                                                                    beta = alpha,
                                                                    alpha = alpha)
      }

      # ----------------------------
      # Summaries
      # ----------------------------

      out <- list(m_true = m_true, HCP_power = mean(HCP_res[, "Reject"]),
                  HCP_ESS   = mean(HCP_res[, "ESS"]), HOLM_power = mean(HOLM_res[, "Reject"]),
                  HOLM_ESS   = mean(HOLM_res[, "ESS"]))

      for (j in seq_along(m1_grid)) {
        out[[paste0("SPRT_power_", m1_grid[j])]] <- mean(SPRT_res[, paste0("Reject_", m1_grid[j])])
        out[[paste0("SPRT_ESS_", m1_grid[j])]] <- mean(SPRT_res[, paste0("ESS_", m1_grid[j])])
      }

      out$SPRT_power_adap <- mean(SPRT_res[, "Reject_adap"])
      out$SPRT_ESS_adap <- mean(SPRT_res[, "ESS_adap"])

      results[[g]] <- tibble::as_tibble(out)
    }

    dplyr::bind_rows(results)
  }

  set.seed(37238493)
  res <- compare_tests(N = N)
  res1 <- compare_tests(N = N1)
  df1 <- res
  df2 <- res1

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
  res1$Design <- paste0("N = ", N1)

  power_df <- bind_rows(res, res1) |>
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
  ESS_df <- bind_rows(res, res1) |>
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

  return(list(df1 = df1, df2 = df2, p2 = p2, p3 = p3))
}
