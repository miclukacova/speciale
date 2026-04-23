#-------------------------------------------------------------------------------
### Rina Barber paper - d_TV calculation
#-------------------------------------------------------------------------------
# The idea is to find two distributions for which we can calculate the TV distance
# We want to find a distribution in the class P and in the class Q
# We start out by assuming a particular distributin regarding data and considering the OLS algorithm and then find
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
riskss <- c()                                            # This is a vector for every single calculation of the conditional probability
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

ggplot(data.frame(vars, risk_est))+
  geom_line(aes(x=vars, y = risk_est), linewidth = 2, color = "steelblue")+
  xlab(expression(gamma))+
  labs(x = expression(gamma), y = "Estimate of Risk")+
  geom_vline(xintercept = tau, color = "darkred", linetype = 2)+
  theme_bw()+
  scale_x_continuous(trans='log10')



