get_e_vs_p_power_plot <- function(data){

  #Making dataframe suited to plot
  data_long <- pivot_longer(
    data,
    cols = c(p_value_power, e_value_power),
    names_to = "test",
    values_to = "power"
  )

  plot <- ggplot(data_long, aes(x = x, y = power, color = test)) +
    geom_line(linewidth = 1.1) +
    geom_point() +
    facet_wrap(~ variable, scales = "free_x") +
    theme_bw() +
    labs(x = NULL, y = "Empirical power", color = "Test")

  return(plot)

}

