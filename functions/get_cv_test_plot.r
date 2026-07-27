get_cv_test_plot <- function(cv_test_data){

  plot <- ggplot(cv_test_data) +
    geom_line(aes(x=taus, y = rej_prob), col = "firebrick") +
    geom_point(aes(x=taus, y = rej_prob), col = "firebrick") +
    theme_bw() +
    labs(x = expression(tau), y = "Rejection Probability") +
    geom_vline(xintercept = 1, color ="steelblue")

  return(plot)
}
