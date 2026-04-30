get_OLS_rej_prob_plot <- function(data1, data2, tau = 0.5, alpha = 0.05, N = 2500, n = 500) {

  gammas <- data1$gamma
  rej_prob <- data1$rej_prob
  gammas2 <- data2[[2]]$gamma

  # We find the risk of the OLS algorithm for the gammas in range of the rejection probability
  gamma_range <- range(gammas)
  low_idx <- which(abs(gammas2 - gamma_range[1]) < 10^(-5))
  high_idx <- which(abs(gammas2 - gamma_range[2]) < 10^(-5))
  risk_gamma <- data2[[2]]$risk[low_idx: high_idx]

  if(alpha < (1 - tau)^(floor(N / (n + 1)))) {
    print("BB test power can be calculated - bound is exact")
    # calculating bb box power
    bound <- alpha * (1 + (tau - risk_gamma)/(1-tau))^(floor(N/(n+1)))
    bound <- pmin(bound, 1)
  }
  else{
    # calculating Rinas bound
    tau_tilde <- (1+((1/alpha)-1)/N)*tau
    bound3 <- alpha * (1 + (tau_tilde - risk_gamma)/(1-tau_tilde))^(N/n)
    bound3 <- pmin(bound3, 1)
    print("Calculating the bound from Rina Thm")
  }

  # Calculating our bound
  gamma_star <- data2$gamma_star
  bound2 <- alpha + sqrt(1 - sqrt(2*gamma_star*gammas/(gamma_star^2+gammas^2)))*sqrt(2)

  plot <- ggplot() +
    geom_line(aes(x = gammas, y = rej_prob, color = "Rejection probability")) +
    geom_line(aes(x = gammas2[low_idx: high_idx], y = bound, color = "Bound"), linetype = 2) +
    geom_line(aes(x = gammas2[low_idx: high_idx], y = bound3, color = "Bound2"), linetype = 2) +
    geom_line(aes(x = gammas, y = bound2, color = "Bound TV"), linetype = 2) +
    geom_vline(aes(xintercept = gamma_star, color = "gamma*"), linetype = 3) +
    theme_bw() +
    labs(x = expression(gamma),y = "Rejection Probability", color = NULL) +
    scale_color_manual(
      values = c("Rejection probability" = "darkred","Bound" = "blue","gamma*" = "black", "Bound TV" = "darkgreen"),
      labels = c("Bound", expression(gamma^"*"), "Rejection probability", "Bound TV"))

  return(plot)
}
