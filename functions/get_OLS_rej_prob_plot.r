get_OLS_rej_prob_plot <- function(data, tau = 0.5, alpha = 0.05, n = 500, N = 2500) {

  # We find the risk of the OLS algorithm for the ms in an interesting range
  m_star <- data[[1]]
  ms <- data[[2]]$ms
  idx <- which(0.59 <= ms & ms<= 0.81)
  risk_m <- data[[2]]$risk[idx]
  ms <- ms[idx]


  # This is an exact power/rejection probability of the test
  if(alpha < (1-tau)^(floor(N /(n + 1)))) {
    print("BB test power can be calculated - bound is exact")
    # calculating bb box power
    exact_rej_prob <- alpha * (1 + (tau - risk_m) / (1 - tau)) ^ (floor(N / (n + 1)))
    exact_rej_prob <- pmin(exact_rej_prob, 1)
  }
  # calculating Rinas bound
  tau_tilde <- (1+((1/alpha)-1)/N)*tau
  bound <- alpha * (1 + (tau_tilde - risk_m)/(1-tau_tilde))^(N/n)
  bound <- pmin(bound, 1)
  print("Calculating the bound from Rina Thm")

  # Calculating our bound
  #m_star <- data2$m_star
  #bound2 <- alpha + sqrt(1 - sqrt(2*m_star*ms/(m_star^2+ms^2)))*sqrt(2)

  plot <- ggplot() +
    geom_line(aes(x = ms, y = exact_rej_prob, color = "Rejection probability"), size = 0.9) +
    geom_line(aes(x = ms, y = bound, color = "Bound"), linetype = 2, size = 0.9) +
    #geom_line(aes(x = ms, y = bound2, color = "Bound TV"), linetype = 2, size = 0.9) +
    geom_vline(aes(xintercept = m_star), linetype = 3, size = 1, color = "black") +
    annotate(
      "text", label = expression(m^"*"),
      x = 0.75, y = 0.12, size = 5, colour = "black"
    )+
    theme_bw() +
    labs(x = expression(m),y = "Rejection Probability", color = NULL) +
    scale_color_manual(
      values = c("Rejection probability" = "darkred","Bound" = "steelblue", "Bound TV" = "darkgreen"),
      labels = c( "Bound", "Bound TV", "Rejection Rate"))

  #ggsave(here::here("plots/OLS_risk_plot.pdf"), plot = plot)

  return(plot)
}
