get_e_vs_p_power_plot <- function(data){

  return()

  #Making dataframe suited to plot
  data_long <- pivot_longer(
    data,
    cols = c(p_power, e_power),
    names_to = "test",
    values_to = "power"
  )

  plot <- ggplot(data_long, aes(x = x, y = power, color = test)) +
    geom_point() +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ variable, scales = "free_x") +
    theme_bw() +
    labs(x = NULL, y = "Empirical power", color = "Test") +
    scale_color_manual(values = c("e_power" = "darkred", "p_power" = "steelblue"))+
    theme(
      strip.background = element_rect(fill = "lavender"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      axis.title.y = element_text(size = 8),
    )

    #ggsave(here::here("plots/e_vs_p_power_plot.pdf"), plot = plot, width = 6.2, height = 2.5)

    return(plot)

}
