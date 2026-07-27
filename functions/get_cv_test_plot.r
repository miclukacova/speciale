get_cv_test_plot <- function(cv_test_data, B, alg, loss, N){

  mse <- c()
  for(b in 1:B){
    #-------------------------------------------------------------------------------
    X1 <- rnorm(N,-10,2)
    X2 <- rnorm(N,-1,1)
    X <- cbind(X1,X2)
    Y <- X %*% beta + rnorm(N)

    X_n_1 <- c(rnorm(1,-10,2), rnorm(1,-1,1))
    Y_n_1 <- t(X_n_1) %*% beta + rnorm(1)
    mse[b] <- loss(alg(X, Y)(X_n_1), Y_n_1)
  }

  true_risk <- mean(mse)
  print(paste("True risk estimate", true_risk))

  plot <- ggplot(cv_test_data) +
    geom_line(aes(x=taus, y = rej_prob), col = "firebrick") +
    geom_point(aes(x=taus, y = rej_prob), col = "firebrick") +
    theme_bw() +
    labs(x = expression(tau), y = "Rejection Probability") +
    geom_vline(xintercept = true_risk, color ="steelblue")

  return(plot)
}
