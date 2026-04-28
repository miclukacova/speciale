get_cv_test_plot <- function(cv_test_data){

  plot <- ggplot(cv_test_data)+
    geom_point(aes(x=taus, y = rej_prob))+
    theme_bw()+
    labs(x = "tau", y = "Rejection Probability")+
    geom_vline(xintercept = 1)

  return(plot)
}
