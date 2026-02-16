# Packages

library(ggplot2)
#install.packages("FNN")                                                        #knn
library(FNN)

#install.packages("MASS")                                                        # ridge
library(MASS)

# Reproducing Liang and Barber 2025: relating m-stability and 1-stability
#-------------------------------------------------------------------------------
# Calculating stability as a function of m

d <- 40
N <- 500
X <- matrix(nrow = N, ncol = d)
Y <- c()
emp.est.knn <- c()
emp.est.ridge <- c()

for(m in 1:25){
  X.m <- matrix(nrow = m+N, ncol = d)
  Y.m <- c()
  diff.knn <- c()
  diff.ridge <- c()

  for(i in 1:1000){

    # Generating the n first data points
    for(j in 1:N){
      X[j,] <- runif(d)
      e <- rbinom(d, 1, 1/3)
      Y[j] <- sum(sin(X[j,])/1:d + e * runif(d, -0.1, 0.1) + (1-e) * runif(d, -1, 1))
    }

    # Generating the next m data points
    X.m[1:N,] <- X
    Y.m[1:N] <- Y

    for(j in 1:m){
      X.m[j + N,] <- runif(d)
      e <- rbinom(d, 1, 1/3)
      Y.m[j + N] <- sum(sin(X.m[j + N,])/1:d + e * runif(d, -0.1, 0.1) + (1-e) * runif(d, -1, 1))
    }

    # Drawing the test point
    X.test <- runif(d)

    # The algorithms
    mu.knn.n <- knn.reg(train = X, test = matrix(X.test, nrow = 1), y = Y, k = 20)$pred
    mu.knn.n.m <- knn.reg(train = X.m, test = matrix(X.test, nrow = 1), y = Y.m, k = 20)$pred

    beta.ridge.n <- solve(t(X) %*% X + 0.01 * diag(ncol(X))) %*% t(X) %*% Y
    beta.ridge.n.m <- solve(t(X.m) %*% X.m + 0.01 * diag(ncol(X.m))) %*% t(X.m) %*% Y.m
    mu.ridge.n <- t(X.test) %*% beta.ridge.n
    mu.ridge.n.m <- t(X.test) %*% beta.ridge.n.m

    # Distance between the algorithm predictions
    diff.knn[i] <- abs(mu.knn.n - mu.knn.n.m)
    diff.ridge[i] <- abs(mu.ridge.n - mu.ridge.n.m)
  }
  # Empirical estimates
  emp.est.knn[m] <- mean(diff.knn)
  emp.est.ridge[m] <- mean(diff.ridge)
}

ggplot()+
  geom_line(aes(x = 1:25, y = emp.est.knn, col = "blue"))+
  geom_line(aes(x = 1:25, y = emp.est.ridge, col = "red3"))

#-------------------------------------------------------------------------------
# Bounds calculated by lemma 5.2
#-------------------------------------------------------------------------------

# To do implementer det her rigtigt

ns <- N:(N + 25)
for(k in 1:25){
  n <- ns[k]

  X.m <- matrix(nrow = 1+n, ncol = d)
  Y.m <- c()
  diff.knn <- c()
  diff.ridge <- c()

  for(i in 1:1000){

    # Generating the n first data points
    for(j in 1:n){
      X[j,] <- runif(d)
      e <- rbinom(d, 1, 1/3)
      Y[j] <- sum(sin(X[j,])/1:d + e * runif(d, -0.1, 0.1) + (1-e) * runif(d, -1, 1))
    }

    # Generating the next m data points
    X.m[1:n,] <- X
    Y.m[1:n] <- Y
    X.m[1 + n,] <- runif(d)
    e <- rbinom(d, 1, 1/3)
    Y.m[1 + n] <- sum(sin(X.m[1 + n,])/1:d + e * runif(d, -0.1, 0.1) + (1-e) * runif(d, -1, 1))

    # Drawing the test point
    X.test <- runif(d)

    # The algorithms
    mu.knn.n <- knn.reg(train = X, test = matrix(X.test, nrow = 1), y = Y, k = 20)$pred
    mu.knn.n.m <- knn.reg(train = X.m, test = matrix(X.test, nrow = 1), y = Y.m, k = 20)$pred

    beta.ridge.n <- solve(t(X) %*% X + 0.01 * diag(ncol(X))) %*% t(X) %*% Y
    beta.ridge.n.m <- solve(t(X.m) %*% X.m + 0.01 * diag(ncol(X.m))) %*% t(X.m) %*% Y.m
    mu.ridge.n <- t(X.test) %*% beta.ridge.n
    mu.ridge.n.m <- t(X.test) %*% beta.ridge.n.m

    # Distance between the algorithm predictions
    diff.knn[i] <- abs(mu.knn.n - mu.knn.n.m)
    diff.ridge[i] <- abs(mu.ridge.n - mu.ridge.n.m)
  }

  # Empirical estimates
  emp.est.knn[k] <- mean(diff.knn)
  emp.est.ridge[k] <- mean(diff.ridge)
}

bound.ridge <- c()
bound.knn <- c()
for(m in 1:25){
  bound.ridge[m] <- sum(emp.est.ridge[N:(N+m-1)])
  bound.knn[m] <- sum(emp.est.knn[N:(N+m-1)])
}

#-------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------

ggplot()+
  geom_line(aes(x = 1:25, y = emp.est.knn, col = "blue"))+
  geom_line(aes(x = 1:25, y = emp.est.ridge, col = "red3"))+
  geom_line(aes(x = 1:25, y = bound.knn, col = "blue"), linetype = 2)+
  geom_line(aes(x = 1:25, y = bound.ridge, col = "red3"), linetype = 2)







