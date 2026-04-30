get_rina_bb_test_plot <- function(data, B = 10^4, N = 2500, n = 500, alpha = 0.05,
                                  taus = seq(0.3,1, by = 0.01)){

  # True Risks
  riskss <- c()
  for(i in 1:10^4){
    X1 <- rnorm(N+1,-10,2)
    X2 <- rnorm(N+1,-1,1)
    X <- cbind(X1[1:N],X2[1:N])
    riskss[i] <- 2* pnorm(-0.5, 0, t(X[n+1,])%*%solve(t(X)%*% X)%*% X[n+1,] +1)

  }
  R_p_A <- mean(riskss)
  print(paste("Risk estimate:", R_p_A))

  # The exact power can be calculated for these taus
  exact_idx <- which(alpha < (1 - taus)^(floor(N / (n + 1))))
  bound <- alpha * (1 + (taus[exact_idx] - mean(riskss))/(1 - taus[exact_idx]))^(floor(N / (n+1)))
  exact_pow <- pmin(bound, 1)

  # Calculating Rinas bound
  tau_tilde <- (1 + ((1 / alpha) - 1) / N) * taus
  bound_idx <- which(tau_tilde < 1)
  bound <- alpha * (1 + (tau_tilde[bound_idx] - R_p_A) / (1 - tau_tilde[bound_idx])) ^ (N / n)
  bound <- pmin(bound, 1)

  plot <- ggplot()+
    geom_line(aes(x=taus[bound_idx], y = bound), color = "blue", size = 2)+
    geom_line(aes(x=taus[exact_idx], y = exact_pow), color = "red", size = 3)+
    geom_point(aes(x=data$tau, y = data$rej_prob), size = 2)+
    theme_bw()+
    labs(x = "tau", y = "Rejection Probability")+
    geom_vline(xintercept = R_p_A, linetype = 2)+
    geom_hline(yintercept = alpha, linetype = 2)

  return(plot)
}
