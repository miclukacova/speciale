# Packages

install.packages("caTools")                                                     #knn
library(caTools)

install.packages("MASS")                                                        # ridge
library(MASS)

# Reproducing Liang and Barber 2025: relating m-stability and 1-stability

d <- 40
N <- 500
m <- 1
X <- matrix(nrow = N, ncol = d)
X.m <- matrix(nrow = m+N, ncol = d)
Y <- c()
Y.m <- c()
diff.knn <- matrix(nrow = 1000, ncol = length(m))
diff.ridge <- matrix(nrow = 1000, ncol = length(m))

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
  mu.knn.n <- knn(train = X, test = X.test, cl = Y, k = 20)
  mu.knn.n.m <-  knn(train = X.m, test = X.test, cl = Y.m, k = 20)

  beta.ridge.n <- solve(t(X) %*% X + 0.01 * diag(ncol(X))) %*% t(X) %*% Y
  beta.ridge.n.m <- solve(t(X.m) %*% X.m + 0.01 * diag(ncol(X.m))) %*% t(X.m) %*% Y.m
  mu.ridge.n <- t(X.test) %*% beta.ridge.n
  mu.ridge.n.m <- t(X.test) %*% beta.ridge.n.m

  # Distance between the algorithm predictions
  diff.knn[i] <- abs(mu.knn.n - mu.knn.n.m)
  diff.ridge[i] <- abs(mu.ridge.n - mu.ridge.n.m)
}
