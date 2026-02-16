# Packages

library(ggplot2)
#install.packages("FNN")                                                        #knn
library(FNN)

#install.packages("MASS")                                                        # ridge
library(MASS)

# Reproducing Liang and Barber 2025: Illustrating the bound from THM 3.2
#-------------------------------------------------------------------------------

#Calculating B_m^out

d <- 40
N <- 500
X <- matrix(nrow = N, ncol = d)
Y <- c()
emp.est.knn <- c()
emp.est.ridge <- c()

for(m in 10){
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
	emp.est.knn <- mean(diff.knn)
	emp.est.ridge <- mean(diff.ridge)
}


#Implementing Jackknife as for loop

alpha <- 0.1
B <- 1000
X <- matrix(nrow = N, ncol = d)
Y <- c()
mu.ridge.n_1 <- c()
res <- c()

for(j in 1:N){
	X[j,] <- runif(d)
	e <- rbinom(d, 1, 1/3)
	Y[j] <- sum(sin(X[j,])/1:d + e * runif(d, -0.1, 0.1) + (1-e) * runif(d, -1, 1))
}

indicator <- c()

for (k in 1:B) {

X.test <- runif(d)
e <- rbinom(d, 1, 1/3)
Y.test <- sum(sin(X.test)/1:d + e * runif(d, -0.1, 0.1) + (1-e) * runif(d, -1, 1))

for (i in 1:N) {
	
	beta.ridge.i <- solve(t(X[-i,]) %*% X[-i,] + 0.01 * diag(ncol(X[-i,]))) %*% t(X[-i,]) %*% Y[-i]
	mu.ridge.i <- t(X[i,]) %*% beta.ridge.i
	res[i] <- abs(mu.ridge.i - Y[i])
	mu.ridge.n_1[i] <- t(X.test) %*% beta.ridge.i
		
}

indicator[k] <- (-quantile(-mu.ridge.n_1+res, (1-alpha)*(1+1/N))<= Y.test & Y.test <= quantile(mu.ridge.n_1+res, (1-alpha)*(1+1/N)))

}
indicator
sum(indicator)/B

