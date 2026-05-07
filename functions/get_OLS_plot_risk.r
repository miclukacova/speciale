get_OLS_plot_risk <- function(OLS_risk, tau = 0.5){

  gamma_star <- OLS_risk[[1]]
  df <- OLS_risk[[2]]

  plot <- ggplot(df)+
    geom_line(aes(x=gamma, y = risk), linewidth = 2, color = "steelblue")+
    xlab(expression(gamma))+
    labs(x = expression(gamma), y = "Estimate of Risk")+
    geom_hline(yintercept = tau, color = "darkred", linetype = 2)+
    geom_vline(xintercept = gamma_star, color = "black", linetype = 2)+
    annotate(
      "text", label = expression(gamma^"*"),
      x = 1, y = 0.05, size = 5, colour = "black"
    )+
    annotate(
      "text", label = expression(tau),
      x = 0.05, y = 0.53, size = 5, colour = "darkred"
    )+
    theme_bw()+
    scale_x_continuous(trans='log10')+
    theme(text = element_text(family = "TT Times New Roman"))

  ggsave(here::here("plots/OLS_risk_plot.pdf"), plot = plot)

  return(plot)
}
