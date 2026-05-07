get_l_q_stability_plot <- function(l_q_stability_data){
  o_n_2 <- function(n) n^(-1/2)

  plot <- ggplot(l_q_stability_data)+
    geom_point(aes(x=n, y = beta_q), size = 2, color = "darkred")+
    geom_line(aes(x=n, y = beta_q), linetype = 2, color = "darkred")+
    geom_function(fun = Vectorize(o_n_2), color = "steelblue", size = 2)+
    theme_bw()+
    labs(x = "Sample Size", y = "Algorithmic Stability of OLS")

  return(plot)
}

