#-------------------------------------------------------------------------------
# Function to compare power and expected sample size in a scenario where we test the mean of Bernoulli variables

if(FALSE){
  source("~/Desktop/Uni/Speciale/speciale/functions/GS_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/HCP_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/Holmes_test.R")
  source("~/Desktop/Uni/Speciale/speciale/functions/SPRT.R")
}
get_seq_test_comp_RCT <- function(B = 500,
                                  B1 = 500,
                                  N = 200,
                                  N1 = 500) {

  # Parameters
  p_t_true_grid <- seq(0.3, 0.7, by = 0.05)
  sprt_grid <- c(0.45, 0.6)
  p_c <- 0.3
  m_0 = 1 / 2
  c = 0.5
  theta = 3 / 2
  alpha = 0.05
  gamma = 0.9
  n_looks = 4
  alphas = alphas_soko(p_c, n_looks, Nmax = N, B = B, alpha = alpha)
  sd0 = sqrt(1 / 2 * p_c * (1 - p_c))

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
    p_t <- c(p_c, 2 * cumsum(X[1:(N - 1)]) / (1:(N - 1)) - 1 + p_c)
    p_t <- pmax(p_t, 0.1)
    p_t * (1 - p_c) * (X == 1) +  (1 - p_t) * p_c * (X == 0) + (X == 1 / 2) * ((1 - p_t) * (1 - p_c) + p_t * p_c)
  }

  # Holmes - quantile estimation for all p_t_true values
  X <- vector(length = B)
  for(i in 1:B) X[i] <- sum(sample_patient(N, p_c))
  z_ag <- quantile(X, 1 - alpha)

  # GST
  #alphas <- rpact::getDesignGroupSequential(kMax = 4,
  #                                          alpha = alpha,
  #                                          sided = 1,
  #                                          typeOfDesign = "OF")$criticalValues

  # Main test comparison function
  compare_tests <- function(N) {
    results <- vector("list", length(p_t_true_grid))
    for (g in seq_along(p_t_true_grid)) {
      print(g)
      # True parameter
      p_t_true <- p_t_true_grid[g]
      sample_data <- function(N) sample_patient(N, p_t_true)

      # Storage
      HCP_res <- matrix(NA_real_, nrow = B, ncol = 2)
      HOLM_res <- matrix(NA_real_, nrow = B, ncol = 2)
      GS_res <- matrix(NA_real_, nrow = B, ncol = 2)
      n_sprt <- length(sprt_grid)
      SPRT_res <- matrix(NA_real_, nrow = B, ncol = 2 * n_sprt + 2)

      # Simulations
      sim_results <- future_lapply(
        seq_len(B),
        function(b) {

          X <- sample_data(N)

          # HCP
          HCP <- run_HCP_test(
            m_0 = 1 / 2,
            c = c,
            X = X,
            theta = theta,
            alpha = alpha
          )

          # HOLMES
          HOLM <- run_holmes_test(
            N = N,
            X = X,
            sample_data_null = function(N) sample_patient(N, p_c),
            gamma = gamma,
            quanti = z_ag,
            B = B1
          )

          # SPRTs
          SPRT <- numeric(2 * n_sprt + 2)

          for (j in seq_along(sprt_grid)) {

            f1 <- function(X) {
              sprt_grid[j] * (1 - p_c) * (X == 1) +
                (1 - sprt_grid[j]) * p_c * (X == 0) +
                (X == 1 / 2) *
                ((1 - sprt_grid[j]) * (1 - p_c) +
                   sprt_grid[j] * p_c)
            }

            SPRT[(2 * j - 1):(2 * j)] <- run_sprt_test(
              N,
              X = X,
              f0 = f0,
              f1 = f1,
              beta = alpha,
              alpha = alpha
            )
          }

          # Adaptive SPRT
          SPRT[(2 * n_sprt + 1):(2 * n_sprt + 2)] <-
            run_sprt_test(
              N,
              X = X,
              f0 = f0,
              f1 = f1_adap,
              beta = alpha,
              alpha = alpha
            )

          # GS
          GS <- gs_run(
            Nmax = N,
            alphas = alphas,
            n_looks = n_looks,
            X = X,
            sd0 = sd0,
            m_0 = m_0
          )

          list(
            HCP = HCP,
            HOLM = HOLM,
            SPRT = SPRT,
            GS = GS
          )
        },
        future.seed = TRUE
      )

      HCP_res <- do.call(
        rbind,
        lapply(sim_results, `[[`, "HCP")
      )
      colnames(HCP_res) <- c("Reject", "ESS")

      HOLM_res <- do.call(
        rbind,
        lapply(sim_results, `[[`, "HOLM")
      )
      colnames(HOLM_res) <- c("Reject", "ESS")

      SPRT_res <- do.call(
        rbind,
        lapply(sim_results, `[[`, "SPRT")
      )
      colnames(SPRT_res) <- c(
        as.vector(rbind(
          paste0("Reject_", sprt_grid),
          paste0("ESS_", sprt_grid)
        )),
        "Reject_adap",
        "ESS_adap"
      )

      GS_res <- do.call(
        rbind,
        lapply(sim_results, `[[`, "GS")
      )
      colnames(GS_res) <- c("Reject", "ESS")

      # ----------------------------
      # Summaries
      # ----------------------------

      out <- list(p_t_true = p_t_true,
                  HCP_power = mean(HCP_res[, "Reject"]),
                  HCP_ESS   = mean(HCP_res[, "ESS"]),
                  HOLM_power = mean(HOLM_res[, "Reject"]),
                  HOLM_ESS   = mean(HOLM_res[, "ESS"]),
                  GS_power = mean(GS_res[, "Reject"]),
                  GS_ESS   = mean(GS_res[, "ESS"]))

      for (j in seq_along(sprt_grid)) {
        out[[paste0("SPRT_power_", sprt_grid[j])]] <- mean(SPRT_res[, paste0("Reject_", sprt_grid[j])])
        out[[paste0("SPRT_ESS_", sprt_grid[j])]] <- mean(SPRT_res[, paste0("ESS_", sprt_grid[j])])
      }

      out$SPRT_power_adap <- mean(SPRT_res[, "Reject_adap"])
      out$SPRT_ESS_adap <- mean(SPRT_res[, "ESS_adap"])

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

  res <- compare_tests(N)
  res1 <- compare_tests(N1)
  df <- res
  df1 <- res1

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
               aes(p_t_true, Power, colour = Method)) +
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
               aes(p_t_true, ESS, colour = Method)) +
    geom_line() +
    facet_wrap(~Design, scales = "free_y") +
    theme_minimal()

  return(list(df = df, df1 = df1, p2 = p2, p3 = p3))
}
