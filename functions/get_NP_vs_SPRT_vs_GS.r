#-------------------------------------------------------------------------------
# Sequential anytime-valid inference using e-processes - Example Ramdas s. 94
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Normal example
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

get_NP_vs_SPRT_vs_GS <- function(B = 5 * 10^4) {
  set.seed(9475)
  alpha <- 10^(-3)
  beta <- alpha

  # Defining parameters
  mu0 <- 0
  mu1 <- 0.3

  #-------------------------------------------------------------------------------
  # Neyman–Pearson
  #-------------------------------------------------------------------------------

  # Definition of the test. It is equivalent to test the mean.
  np_norm_test <- function(n, mu) {
    X <- rnorm(n, mu, 1)
    mean(X) > gamma
  }

  # We find the quantile and needed sample size for the desired power and the desired type I error
  alpha_n <- 1
  beta_n <- 1
  n_np <- 300

  while(beta_n > beta){

    # First we find the needed quantile
    gamma <- mu0 + qnorm(1 - alpha)/sqrt(n_np)
    # Then we calculate what power this results in
    beta_n <- pnorm(gamma, mu1, 1 / sqrt(n_np))

    n_np <- n_np + 1
  }

  n_np <- n_np-1
  print(n_np)

  #-------------------------------------------------------------------------------
  # SPRT
  #-------------------------------------------------------------------------------

  gamma0 <- beta / (1 - alpha)
  gamma1 <- (1 - alpha) / beta

  ratio_norm <- function(mu) {
    ratio_func <- function(){
      X <- rnorm(1, mu, 1)
      prod(dnorm(X,mu1,1)/dnorm(X,mu0,1))
    }
    return(ratio_func)
  }

  sprt_test <- function(ratio) {
    L <- 1
    i <- 0
    while(L <= gamma1 & gamma0 <= L){
      L <- L * ratio()
      i <- i + 1
    }
    return(c(L >= 1, i))
  }
  #-------------------------------------------------------------------------------
  # Group-Sequential trial
  #-------------------------------------------------------------------------------

  gammas <- rpact::getDesignGroupSequential(kMax = 4, alpha = alpha,
                                            sided = 1, typeOfDesign = "OF")$criticalValues

  np_test_gs <- function(mu) {

    gamma <- gammas[1]
    Z_star <- 1
    Z_k <- vector()
    k <- 1

    while((Z_star <= gamma) & (k <= 4)){
      gamma <- gammas[k]
      X <- rnorm(n_gs / 4, mu)
      Z_k[k] <- mean(X) * sqrt(n_gs / 4)
      Z_star <- 1 / sqrt(k) * sum(Z_k)
      k <- k + 1
    }

    return(c(Z_star > gammas[k-1] , n_gs / 4 * (k - 1)))

  }

  n_gs <- 400
  pow_gs <- rep(0,B)
  while(mean(pow_gs) < 1-beta){
    for(j in 1:B) pow_gs[j] <- np_test_gs(mu1)[1]
    n_gs <- n_gs + 4
  }

  print(n_gs)

  #-------------------------------------------------------------------------------
  # Simulation Function
  #-------------------------------------------------------------------------------

  sim_func <- function(mu) {

    test_res <-  matrix(ncol = 4, nrow = B)
    colnames(test_res) <- c("SPRT_power", "SPRT_ESS", "GS_power" , "GS_ESS")
    NP_test_res <- vector()

    for(j in 1:B){
      test_res[cbind(j, 1:2)] <- sprt_test(ratio_norm(mu))
      NP_test_res[j] <- np_norm_test(n_np, mu)
      test_res[cbind(j, 3:4)] <- np_test_gs(mu)
    }

    return(c("NP_power" = mean(NP_test_res), colMeans(test_res)))
  }


  #-------------------------------------------------------------------------------
  # Comparison
  #-------------------------------------------------------------------------------

  mus <- seq(0,0.3,by=0.025)

  # Run simulations over Grid

  results <- matrix(ncol = 5, nrow = length(mus))
  colnames(results) <- c("NP_power","SPRT_power", "SPRT_ESS", "GS_power", "GS_ESS")
  i <- 1

  for(mu in mus){
    results[i,] <- sim_func(mu)
    i <- i + 1
    print(mu)
  }

  # Power curves

  p1 <- ggplot(results, aes(mus)) +
    geom_line(aes(y = NP_power, colour = "Fixed NP"), linewidth = 1.1) +
    geom_line(aes(y = SPRT_power, colour = "SPRT"), linewidth = 1.1) +
    geom_line(aes(y = GS_power, colour = "GS test"), linewidth = 1.1) +
    labs(
      y = "Power",
      colour = "",
      title = "Power curves"
      ) +
    theme_minimal()

  # Expected sample size curves
  p2 <- ggplot(results, aes(mus)) +
    geom_line(aes(y = n_np, colour = "Fixed NP"), linewidth = 1.1) +
    geom_line(aes(y = SPRT_ESS, colour = "SPRT"), linewidth = 1.1) +
    geom_line(aes(y = GS_ESS, colour = "GS test"), linewidth = 1.1) +
    labs(
      y = "Expected sample size",
      colour = "",
      title = "Expected sample size") +
    theme_minimal()

  return(list(results = results, Power = p1, ESS = p2))
}

