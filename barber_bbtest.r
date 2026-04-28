#-------------------------------------------------------------------------------
### Rina Barber paper - Binomial example
#-------------------------------------------------------------------------------
# The idea is to employ the binomial black box algorithm test

set.seed(6389)

taus <- seq(0.3,1, by = 0.01)
B <- 10000
N <- 5000                         # Data size
n <- 500                          # Evaluation size

rej_prob <- c()
i <- 1

# Estimate Risk
loss <- function(y_hat, y) abs(y-y_hat) > 0.5

lm_alg <- function(X, Y) {
  solve(t(X) %*% X) %*% t(X) %*% Y
}

beta <- c(0.3,0.9,-2)
test_res <- c()

for(tau in taus){
  k_star <- 0
  a_star <- -1
  while(a_star < 0 | a_star > 1){
    k_star <- k_star + 1
    a_star <- (0.05 - pbinom(k_star -1, floor(N/(n+1)), prob = tau))/dbinom(k_star, floor(N/(n+1)), prob = tau)
  }

  for(b in 1:B){
    # Data
    #-------------------------------------------------------------------------------
    X1 <- rnorm(N,-10,2)
    X2 <- rnorm(N,-1,1)
    X3 <- rbinom(N, 1, 0.4)
    X <- cbind(X1,X2,X3)
    Y <- X %*% beta + rnorm(N)

    # Black box test
    #-------------------------------------------------------------------------------

    # Fit models
    K <- floor(N/(n+1))
    parti <- sample(1:N, N)

    hat_f_k <- list()
    for(k in 1:K){
      indices <- ((k -1) * (n+1) + 1):(k*(n+1) -1)
      hat_f_k[[k]] <- lm_alg(X[indices,], Y[indices])
    }

    loss_k <- c()

    for(k in 1:K){
      y_hat <- X[k*(n+1),] %*% hat_f_k[[k]]
      loss_k[k] <- loss(y_hat, Y[k*(n+1)])
    }

    S <- sum(loss_k)
    test_res[b] <- (S < k_star) + (S == k_star) * (runif(1) <= a_star)
  }

  rej_prob[i] <- mean(test_res)
  i <- i+1
}


riskss <- c()
for(i in 1:10^4){
  X1 <- rnorm(N+1,-10,2)
  X2 <- rnorm(N+1,-1,1)
  X3 <- rbinom(N+1, 1, 0.4)
  X <- cbind(X1[1:N],X2[1:N],X3[1:N])

  riskss[i] <- 2* pnorm(-0.5, 0, t(X[n+1,])%*%solve(t(X)%*% X)%*% X[n+1,] +1)
}

g <- function(tau) {
  tau_tilde <- (1+((1/0.05)-1)/N)*tau
  bound <- 0.05 * (1 + (tau_tilde - mean(riskss))/(1-tau_tilde))^(N/n)
  return(min(bound, 1))
}


vg <- Vectorize(g)
vg(taus)

ggplot(data.frame(taus, rej_prob))+
  geom_point(aes(x=taus, y = rej_prob), size = 2)+
  geom_function(fun = vg, color = "steelblue", size = 2)+
  theme_bw()+
  labs(x = "tau", y = "Rejection Probability")+
  geom_vline(xintercept = mean(riskss), color = "red", linetype = 2)


#-------------------------------------------------------------------------------
### Rina Barber paper - Binomial example
#-------------------------------------------------------------------------------
# The idea is to employ the binomial black box algorithm test

set.seed(6389)

taus <- seq(0.3,1, by = 0.01)
B <- 10000
N <- 5000                         # Data size
n <- 500                          # Evaluation size

rej_prob <- c()
i <- 1

# Estimate Risk
loss <- function(y_hat, y) abs(y-y_hat) > 0.5

lm_alg <- function(X, Y) {
  solve(t(X) %*% X) %*% t(X) %*% Y
}

beta <- c(0.3,0.9,-2)
test_res <- c()

for(tau in taus){
  k_star <- 0
  a_star <- -1
  while(a_star < 0 | a_star > 1){
    k_star <- k_star + 1
    a_star <- (0.05 - pbinom(k_star -1, floor(N/(n+1)), prob = tau))/dbinom(k_star, floor(N/(n+1)), prob = tau)
  }

  for(b in 1:B){
    # Data
    #-------------------------------------------------------------------------------
    X1 <- rnorm(N,-10,2)
    X2 <- rnorm(N,-1,1)
    X3 <- rbinom(N, 1, 0.4)
    X <- cbind(X1,X2,X3)
    Y <- X %*% beta + rnorm(N)

    # Black box test
    #-------------------------------------------------------------------------------

    # Fit models
    K <- floor(N/(n+1))
    parti <- sample(1:N, N)

    hat_f_k <- list()
    for(k in 1:K){
      indices <- ((k -1) * (n+1) + 1):(k*(n+1) -1)
      hat_f_k[[k]] <- lm_alg(X[indices,], Y[indices])
    }

    loss_k <- c()

    for(k in 1:K){
      y_hat <- X[k*(n+1),] %*% hat_f_k[[k]]
      loss_k[k] <- loss(y_hat, Y[k*(n+1)])
    }

    S <- sum(loss_k)
    test_res[b] <- (S < k_star) + (S == k_star) * (runif(1) <= a_star)
  }

  rej_prob[i] <- mean(test_res)
  i <- i+1
}


riskss <- c()
for(i in 1:10^4){
  X1 <- rnorm(N+1,-10,2)
  X2 <- rnorm(N+1,-1,1)
  X3 <- rbinom(N+1, 1, 0.4)
  X <- cbind(X1[1:N],X2[1:N],X3[1:N])

  riskss[i] <- 2* pnorm(-0.5, 0, t(X[n+1,])%*%solve(t(X)%*% X)%*% X[n+1,] +1)
}

g <- function(tau) {
  tau_tilde <- (1+((1/0.05)-1)/N)*tau
  bound <- 0.05 * (1 + (tau_tilde - mean(riskss))/(1-tau_tilde))^(N/n)
  return(min(bound, 1))
}


vg <- Vectorize(g)
vg(taus)

ggplot(data.frame(taus, rej_prob))+
  geom_point(aes(x=taus, y = rej_prob), size = 2)+
  geom_function(fun = vg, color = "steelblue", size = 2)+
  theme_bw()+
  labs(x = "tau", y = "Rejection Probability")+
  geom_vline(xintercept = mean(riskss), color = "red", linetype = 2)



