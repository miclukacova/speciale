get_OLS_plot_risk <- function(OLS_risk, tau = 0.5){

  gamma_star <- OLS_risk[[1]]
  df <- OLS_risk[[2]]

  plot <- ggplot(df)+
    geom_line(aes(x=gamma, y = risk), linewidth = 2, color = "steelblue")+
    xlab(expression(gamma))+
    labs(x = expression(gamma), y = "Estimate of Risk")+
    geom_hline(yintercept = tau, color = "darkred", linetype = 2)+
    geom_vline(xintercept = gamma_star, color = "darkred", linetype = 2)+
    theme_bw()+
    scale_x_continuous(trans='log10')

  return(plot)
}
