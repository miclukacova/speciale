get_rina_bb_test_plot <- function(data, B = 10^4, N = 3000, n = 500){

  # True Risks
  riskss <- c()
  for(i in 1:10^4){
    X1 <- rnorm(N+1,-10,2)
    X2 <- rnorm(N+1,-1,1)
    X <- cbind(X1[1:N],X2[1:N])

    riskss[i] <- 2* pnorm(-0.5, 0, t(X[n+1,])%*%solve(t(X)%*% X)%*% X[n+1,] +1)
  }

  g <- function(tau) {
    tau_tilde <- (1+((1/0.05)-1)/N)*tau
    bound <- 0.05 * (1 + (tau_tilde - mean(riskss))/(1-tau_tilde))^(N/n)
    return(min(bound, 1))
  }

  vg <- Vectorize(g)

  ggplot(data)+
    geom_point(aes(x=tau, y = rej_prob), size = 2)+
    geom_function(fun = vg, color = "steelblue", size = 2)+
    theme_bw()+
    labs(x = "tau", y = "Rejection Probability")+
    geom_vline(xintercept = mean(riskss), color = "red", linetype = 2)
}
