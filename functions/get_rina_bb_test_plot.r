get_rina_bb_test_plot <- function(data, B = 10^4, N = 3000, n = 500, alpha = 0.05, tau = 0.5){

  # True Risks
  riskss <- c()
  for(i in 1:10^4){
    X1 <- rnorm(N+1,-10,2)
    X2 <- rnorm(N+1,-1,1)
    X <- cbind(X1[1:N],X2[1:N])

    riskss[i] <- 2* pnorm(-0.5, 0, t(X[n+1,])%*%solve(t(X)%*% X)%*% X[n+1,] +1)
  }

  if(alpha < (1 - tau)^(floor(N / (n + 1)))) {
    print("BB test power can be calculated - bound is exact")
    # calculating bb box power
    g <- function(tau) {
      bound <- alpha * (1 + (tau - risk_gamma)/(1 - tau))^(floor(N / (n+1)))
      return(pmin(bound, 1))
    }
  }
  else{
    # calculating Rinas bound
    g <- function(tau) {
      tau_tilde <- (1 + ((1 / alpha) - 1) / N) * tau
      bound <- alpha * (1 + (tau_tilde - mean(riskss)) / (1 - tau_tilde)) ^ (N / n)
      return(pmin(bound, 1))
    }
    print("Calculating the bound from Rina Thm")
  }

  plot <- ggplot(data)+
    geom_point(aes(x=tau, y = rej_prob), size = 2)+
    geom_function(fun = g, color = "steelblue", size = 2)+
    theme_bw()+
    labs(x = "tau", y = "Rejection Probability")+
    geom_vline(xintercept = mean(riskss), color = "red", linetype = 2)

  return(plot)
}
