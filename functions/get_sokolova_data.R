# This function simulates CRT data from different settings similar to the ones from Sokolova2026
# Simulation I seeks to recreate the results from the paper
# Simulation II examines what happens when we let the e-process run longer
# Simulation III examines how the processes react to misspecifications

get_sokolova_data <- function(B = 5*10^4) {
  #-------------------------------------------------------------------------------
  # Global Settings
  #-------------------------------------------------------------------------------

  alpha <- 0.025
  B <- 5*10^4
  set.seed(27940)
  alphas <- alphas_soko(p_c = 0.3,
                        p_t = 0.45,
                        n_looks = 4,
                        Nmax = 200,
                        alpha = alpha)
  alphas <- getDesignGroupSequential(kMax = 4,
                                     alpha = alpha,
                                     sided = 1,
                                     typeOfDesign = "OF")$criticalValues

  #-------------------------------------------------------------------------------
  # Simulations I
  #-------------------------------------------------------------------------------

  res_200 <- run_scenario(
    p_c = 0.30,
    p_t = 0.45,
    Nmax = 200,
    n_looks = 4,
    alpha = 0.025,
    B = B,
    alphas = alphas,
    lambda = NULL
  )

  res_200

  #-------------------------------------------------------------------------------
  # Simulations II
  #-------------------------------------------------------------------------------

  # DO not work yet

  res_500 <- run_scenario(
    p_c = 0.30,
    p_t = 0.45,
    Nmax = 200,
    n_looks = 4,
    alpha = 0.025,
    B = B,
    alphas = alphas,
    lambda = NULL
  )

  res_500


  return(list(res_200 = res_200, res_500 = res_500))
}
