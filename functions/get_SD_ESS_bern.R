get_SD_ESS_bern <- function(){
  if(FALSE){
    source("~/Desktop/Uni/Speciale/speciale/functions/GS_test.R")
    source("~/Desktop/Uni/Speciale/speciale/functions/HCP_test.R")
    source("~/Desktop/Uni/Speciale/speciale/functions/HW_test.R")
    source("~/Desktop/Uni/Speciale/speciale/functions/SPRT.R")
    library(tidyverse)
    library(dplyr)
    library(ggplot2)
    library(ggpattern)
  }

  #==================================================
  # Simulation for one p_t value
  #==================================================

  run_simulation <- function(p_t_true) {

    sample_data <- function(N) sample_patient(N, p_t_true)

    # Storage
    HCP_res <- matrix(NA_real_, nrow = B, ncol = 2)
    HW_res <- matrix(NA_real_, nrow = B, ncol = 2)
    GS_res <- matrix(NA_real_, nrow = B, ncol = 2)
    SPRT_0.45_res <- matrix(NA_real_, nrow = B, ncol = 2)
    SPRT_0.6_res <- matrix(NA_real_, nrow = B, ncol = 2)
    SPRT_adap_res <- matrix(NA_real_, nrow = B, ncol = 2)


    # Simulations
    for(b in seq_len(B)) {

      X <- sample_data(N)

      # HCP
      HCP_res[b,] <- run_HCP_test(
        m_0 = m_0,
        c = c,
        X = X,
        theta = theta,
        alpha = alpha
      )

      # HW
      HW_res[b,] <- run_HW_test(
        N = N,
        X = X,
        calc_q_n = calc_q_n,
        gamma = gamma,
        quanti = z_ag
      )

      # SPRTs
      SPRT_0.45_res[b,] <- run_sprt_test(
        N,
        X = X,
        f0 = f0,
        f1 = f1_list[[1]],
        gamma0 = alpha/(1-alpha),
        gamma1 = (1-alpha)/alpha
      )

      SPRT_0.6_res[b,] <- run_sprt_test(
        N,
        X = X,
        f0 = f0,
        f1 = f1_list[[2]],
        gamma0 = alpha/(1-alpha),
        gamma1 = (1-alpha)/alpha
      )

      SPRT_adap_res[b,] <- run_sprt_test(
        N,
        X = X,
        f0 = f0,
        f1 = f1_adap,
        gamma0 = 0,
        gamma1 = (1-alpha)/alpha
      )

      # GS
      GS_res[b,] <- gs_run(
        Nmax = N,
        alphas = alphas,
        n_looks = n_looks,
        X = X,
        m_0 = m_0,
        sigmaUnknown = FALSE,
        sigma = sigma
      )
    }


    return(list(
      HCP = HCP_res,
      HW = HW_res,
      GS = GS_res,
      SPRT_045 = SPRT_0.45_res,
      SPRT_060 = SPRT_0.6_res,
      SPRT_adap = SPRT_adap_res
    ))

  }


  #==================================================
  # Run simulations
  #==================================================

  res_045 <- run_simulation(p_t_true = 0.45)

  res_060 <- run_simulation(p_t_true = 0.60)



  #==================================================
  # Summarize SD of ESS
  #==================================================

  SD_ESS <- bind_rows(

    tibble(
      p_t = "0.45",
      method = c(
        "GS",
        "HCP",
        "HW",
        "SPRT_0.45",
        "SPRT_0.6",
        "SPRT_adap"
      ),
      SD = c(
        sd(res_045$GS[,2]),
        sd(res_045$HCP[,2]),
        sd(res_045$HW[,2]),
        sd(res_045$SPRT_045[,2]),
        sd(res_045$SPRT_060[,2]),
        sd(res_045$SPRT_adap[,2])
      )
    ),


    tibble(
      p_t = "0.60",
      method = c(
        "GS",
        "HCP",
        "HW",
        "SPRT_0.45",
        "SPRT_0.6",
        "SPRT_adap"
      ),
      SD = c(
        sd(res_060$GS[,2]),
        sd(res_060$HCP[,2]),
        sd(res_060$HW[,2]),
        sd(res_060$SPRT_045[,2]),
        sd(res_060$SPRT_060[,2]),
        sd(res_060$SPRT_adap[,2])
      )
    )

  )



  #==================================================
  # Plot SD of ESS
  #==================================================


  p1 <- ggplot(
      SD_ESS,
      aes(
        x = method,
        y = SD,
        fill = method,
        pattern = p_t
      )
    ) +
    geom_col_pattern(
      position = position_dodge(width = 0.8),
      width = 0.7,
      alpha = 0.85,

      # Pattern settings
      pattern_fill = "white",
      pattern_colour = "white",
      pattern_density = 0.1,
      pattern_spacing = 0.02,
      pattern_angle = 45
    ) +

    scale_fill_manual(values = c(
      "GS" = "darkgreen",
      "HCP" = "firebrick",
      "HW" = "steelblue",
      "SPRT_0.45" = "orange",
      "SPRT_0.6" = "goldenrod",
      "SPRT_adap" = "#009E73"
    ),
    guide = "none") +

    scale_pattern_manual(values = c(
      "0.45" = "none",
      "0.60" = "stripe"
    )) +

    labs(
      x = NULL,
      y = "Standard deviation of ESS",
      fill = "Test",
      pattern = expression(p[t])
    ) +

    theme_bw(base_size = 14) +

    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.minor = element_blank(),

    )

  return()

}





