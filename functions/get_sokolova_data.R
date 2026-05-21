# This function simulates CRT data from different settings similar to the ones from Sokolova2026
# Simulation I seeks to recreate the results from the paper
# Simulation II examines what happens when we let the SPRT run its natural course
# Simulation III examines what happens when Nmax = 300

get_sokolova_data <- function(B = 5*10^4) {
  #-------------------------------------------------------------------------------
  # Global Settings
  #-------------------------------------------------------------------------------

  alpha <- 0.025
  B <- 5*10^4
  set.seed(27940)

  #-------------------------------------------------------------------------------
  # Simulations I
  #-------------------------------------------------------------------------------

  alphas <- alphas_soko(p_c = 0.3,
                        p_t = 0.45,
                        n_looks = 4,
                        Nmax = 200,
                        alpha = alpha)
  #alphas <- getDesignGroupSequential(kMax = 4,
  #                                   alpha = alpha,
  #                                   sided = 1,
  #                                   typeOfDesign = "OF")$criticalValues

  res_200 <- run_scenario(
    p_c = 0.30,
    p_t = 0.45,
    Nmax = 200,
    n_looks = 4,
    alpha = 0.025,
    B = B,
    alphas = alphas,
    lambda = NULL,
    noMax = FALSE
  )

  res_200

  #-------------------------------------------------------------------------------
  # Simulations II
  #-------------------------------------------------------------------------------

  res_noMax <- run_scenario(
    p_c = 0.30,
    p_t = 0.45,
    Nmax = 200,
    n_looks = 4,
    alpha = 0.025,
    B = B,
    alphas = alphas,
    lambda = NULL,
    noMax = TRUE
  )

  res_noMax


  #-------------------------------------------------------------------------------
  # Simulations III
  #-------------------------------------------------------------------------------

  alphas <- alphas_soko(p_c = 0.3,
                        p_t = 0.45,
                        n_looks = 8,
                        Nmax = 400,
                        alpha = alpha)

  res_400 <- run_scenario(
    p_c = 0.30,
    p_t = 0.45,
    Nmax = 400,
    n_looks = 8,
    alpha = 0.025,
    B = B,
    alphas = alphas,
    lambda = NULL,
    noMax = FALSE
  )

  res_400

  return(list(res_200 = res_200, res_noMax = res_noMax, res_400 = res_400))
}
