get_OLS_plot_risk <- function(OLS_risk, tau = 0.5){

  m_star <- OLS_risk[[1]]
  df <- OLS_risk[[2]]

  plot <- ggplot(df)+
    geom_line(aes(x=ms, y = risk), linewidth = 2, color = "firebrick")+
    xlab(expression(m))+
    labs(x = expression(m), y = "Estimate of Risk")+
    geom_hline(yintercept = tau, color = "black", linetype = 2)+
    geom_vline(xintercept = m_star, color = "black", linetype = 2)+
    annotate(
      "text", label = expression(m^"*"),
      x = 1, y = 0.05, size = 5, colour = "black"
    )+
    annotate(
      "text", label = expression(tau),
      x = 0.05, y = 0.53, size = 5, colour = "black"
    )+
    theme_bw()+
    scale_x_continuous(trans='log10')

  #ggsave(here::here("plots/OLS_risk_plot.pdf"), plot = plot)

  return(plot)
}
