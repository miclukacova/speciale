#-------------------------------------------------------------------------------
### Rina Barber paper - d_TV calculation
#-------------------------------------------------------------------------------
# The idea is to find two distributions for which we can calculate the TV distance
# We want to find a distribution in the class P and in the class Q
# We start out by assuming a particular distribution regarding data and considering the OLS algorithm and then find
# two distributions that are very close together

set.seed(6389)
library(ggplot2)

tau <- 0.5
B <- 10^4                         # Bootstrap size
N <- 5000                         # Data size
n <- 500                          # Evaluation size

# The loss function
loss <- function(y_hat, y) abs(y-y_hat) > 0.5

# The algorithm
lm_alg <- function(X, Y) {
  solve(t(X) %*% X) %*% t(X) %*% Y
}

# Our particular choice of betas
beta <- c(0.3,0.9)

# Definition of vectors for the collection of results
riskss <- c()                                                                # This is a vector for every single calculation of the conditional probability
vars <- c(seq(0.05,1, by =0.01), seq(1,10, by = 1), seq(10,100, by = 10))    # This is the sequence standard deviations
risk_est <- c()

for(j in seq_along(vars)){
  for(i in 1:B){
    X1 <- rnorm(N+1,-10,2)
    X2 <- rnorm(N+1,-1,1)
    X <- cbind(X1[1:N],X2[1:N])

    riskss[i] <- 2* pnorm(-0.5, 0, t(X[n+1,])%*%solve(t(X)%*% X)%*% X[n+1,] +vars[j])
  }
  risk_est[j] <- mean(riskss)
}

# We find gamma^*
gamma_star <- vars[min(which(risk_est >= tau))]

ggplot(data.frame(vars, risk_est))+
  geom_line(aes(x=vars, y = risk_est), linewidth = 2, color = "steelblue")+
  xlab(expression(gamma))+
  labs(x = expression(gamma), y = "Estimate of Risk")+
  geom_hline(yintercept = tau, color = "darkred", linetype = 2)+
  geom_vline(xintercept = gamma_star, color = "darkred", linetype = 2)+
  theme_bw()+
  scale_x_continuous(trans='log10')


#-------------------------------------------------------------------------------
### We now want to illustrate that for gamma close to gamma^* the rejection probability is close to alpha
#-------------------------------------------------------------------------------

set.seed(6389)

alpha <- 0.05
gammas <- seq(0.6,0.9, by = 0.05)
risk_gamma <- risk_est[56:86]

rej_prob <- c()
i <- 1
test_res <- c()

for(gamma in gammas){
  k_star <- 0
  a_star <- -1
  while(a_star < 0 | a_star > 1){
    k_star <- k_star + 1
    a_star <- (alpha - pbinom(k_star -1, floor(N/(n+1)), prob = tau))/dbinom(k_star, floor(N/(n+1)), prob = tau)
  }

  for(b in 1:B){
    # Data
    #-------------------------------------------------------------------------------
    X1 <- rnorm(N,-10,2)
    X2 <- rnorm(N,-1,1)
    X3 <- rbinom(N, 1, 0.4)
    X <- cbind(X1,X2,X3)
    Y <- X %*% beta + rnorm(N, sd = gamma)

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

# calculating Rinas bound
tau_tilde <- (1+((1/alpha)-1)/N)*tau
bound <- alpha * (1 + (tau_tilde - risk_gamma)/(1-tau_tilde))^(N/n)
bound <- pmin(bound, 1)

# calculating bb box power

# alpha < (1-tau)^(floor(1000/(n+1)))

tau_tilde <- (1+((1/alpha)-1)/N)*tau
bound <- alpha * (1 + (tau_tilde - risk_gamma)/(1-tau_tilde))^(N/n)
bound <- pmin(bound, 1)

bound2 <- alpha + sqrt(1 - sqrt(2*gamma_star*gammas/(gamma_star^2+gammas^2)))*sqrt(2)


ggplot() +
  geom_line(aes(x = gammas, y = rej_prob, color = "Rejection probability")) +
  geom_line(aes(x = vars[56:86], y = bound, color = "Bound"), linetype = 2) +
  geom_line(aes(x = gammas, y = bound2, color = "Bound TV"), linetype = 2) +
  geom_vline(aes(xintercept = gamma_star, color = "gamma*"), linetype = 2) +
  theme_bw() +
  labs(x = expression(gamma),y = "Rejection Probability", color = NULL) +
  scale_color_manual(
    values = c("Rejection probability" = "darkred","Bound" = "blue","gamma*" = "black", "Bound TV" = "blue"),
    labels = c("Bound", expression(gamma^"*"), "Rejection probability", ))



