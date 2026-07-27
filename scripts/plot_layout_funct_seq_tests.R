library(ggplot2)
library(patchwork)

plot_layout_seq_test_results <- function(p1, p2){

common_theme <- theme_bw(base_size = 11) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),

    strip.text = element_text(size = 10),
    strip.background = element_rect(fill = "white"),

    panel.spacing.x = unit(0.7, "lines"),
    panel.background = element_rect(fill = "white"),

    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(1.5, "cm"),

    plot.margin = margin(5, 5, 5, 5)
  )

p1_edit <- p1 +
  common_theme

p2_edit <- p2 +
  common_theme

combined_plot <- (p1_edit | p2_edit) +
  plot_layout(
    widths = c(1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )

return(combined_plot)

}

#HUSK at gemme dem med width = 8.5 og height = 4.71#############################
plot_layout_seq_test_results(seq_test_comp_RCT_p_t$p2, seq_test_comp_RCT_p_t$p3)
plot_layout_seq_test_results(seq_test_comp_RCT_N$p_power, seq_test_comp_RCT_N$p_ESS)
plot_layout_seq_test_results(seq_test_comp_RCT_norm$p2, seq_test_comp_RCT_norm$p3)
plot_layout_seq_test_results(seq_test_comp_RCT_norm_N$p_power, seq_test_comp_RCT_norm_N$p_ESS)
