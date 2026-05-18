#-------------------------------------------------------------------------------
# Sequential anytime-valid inference using e-processes - Example Ramdas s. 94
#-------------------------------------------------------------------------------

# Binomial example

#-------------------------------------------------------------------------------
# Neyman–Pearson
#-------------------------------------------------------------------------------

# Defining parameters
p0 <- 0.5
p1 <- 0.6
alpha <- 0.05
beta <- 0.05
B <- 3*10^4

# Definition of L_n
L_n_binom <- function(n, p) {
  X <- rbinom(n, 1, p)
  prod(dbinom(X,1,p1)/dbinom(X,1,p0))
}

# We empirically find the quantile
Lns <- vector()
n <- 280
for(i in 1:B) Lns[i] <- L_n_binom(n, p0)
hist(Lns)
gamma <- quantile(Lns,0.95)

# We define the Neyman Pearson test
np_test <- function(L_n) as.numeric(L_n >= gamma)

#-------------------------------------------------------------------------------
# SPRT
#-------------------------------------------------------------------------------

gamma0 <- beta
gamma1 <- 1/alpha

ratio <- function(p) {
  ratio_func <- function(){
    X <- rbinom(1, 1, p)
    prod(dbinom(X,1,p1)/dbinom(X,1,p0))
  }
  return(ratio_func)
}

sprt_test <- function(ratio) {
  L <- 1
  i <- 0
  while(L <= gamma1 & gamma0 <= L){
    L <- L * ratio()
    i <- i + 1
  }
  return(c(L >= 1, i))
}

sprt_test(ratio(p1))

#-------------------------------------------------------------------------------
# Comparison
#-------------------------------------------------------------------------------

# Under nul

sprt_test_res <-  matrix(ncol = 2, nrow = B)
NP_test_res <- vector()

for(j in 1:B){
  sprt_test_res[j,] <- sprt_test(ratio(p0))
  L <- L_n_binom(n, p0)
  NP_test_res[j] <- np_test(L)
}

mean(NP_test_res)
colMeans(sprt_test_res)

# Under alternativ

sprt_test_res1 <- matrix(ncol = 2, nrow = B)
NP_test_res1 <- vector()

for(j in 1:B){
  sprt_test_res1[j,] <- sprt_test(ratio(p1))
  L <- L_n_binom(n, p1)
  NP_test_res1[j] <- np_test(L)
}

mean(NP_test_res1)
colMeans(sprt_test_res1)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Normal example
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Neyman–Pearson
#-------------------------------------------------------------------------------

alpha <- 0.05
beta <- alpha
set.seed(9475)

# Defining parameters
mu0 <- 0
mu1 <- 0.3
B <- 3*10^4

# Definition of L_n
L_n_norm <- function(n, mu) {
  X <- rnorm(n, mu, 1)
  prod(dnorm(X,mu1,1)/dnorm(X,mu0,1))
}

# We empirically find the quantile for the desired power and type 1 error
# First we find which n corresponds to the desired significance level and type 2 error
alpha_n <- 1
beta_n <- 1
n <- 120

np_test_res0 <- vector()
np_test_res1 <- vector()
while(alpha_n > alpha | beta_n > beta){

  # First we find the needed quantile
  Lns <- vector()
  for(i in 1:(B)) Lns[i] <- L_n_norm(n, mu0)
  gamma <- quantile(Lns,1-alpha)

  # Then we calculate the type 1 and 2 error
  for(j in 1:(2*B)){
    L <- L_n_norm(n, mu0)
    np_test_res0[j] <- np_test(L)
    L <- L_n_norm(n, mu1)
    np_test_res1[j] <- np_test(L)
  }
  alpha_n <- mean(np_test_res0)
  beta_n <- 1 - mean(np_test_res1)
  n <- n + 1
}

n <- n-1
print(n)

#-------------------------------------------------------------------------------
# SPRT
#-------------------------------------------------------------------------------

gamma0 <- beta
gamma1 <- 1/alpha

ratio_norm <- function(mu) {
  ratio_func <- function(){
    X <- rnorm(1, mu, 1)
    prod(dnorm(X,mu1,1)/dnorm(X,mu0,1))
  }
  return(ratio_func)
}

#-------------------------------------------------------------------------------
# Comparison
#-------------------------------------------------------------------------------

# Under nul

sprt_test_res <-  matrix(ncol = 2, nrow = B)
NP_test_res <- vector()

for(j in 1:B){
  sprt_test_res[j,] <- sprt_test(ratio_norm(mu0))
  L <- L_n_norm(n, mu0)
  NP_test_res[j] <- np_test(L)
}

mean(NP_test_res)
colMeans(sprt_test_res)

# Under alternativ

sprt_test_res1 <- matrix(ncol = 2, nrow = B)
NP_test_res1 <- vector()

for(j in 1:B){
  sprt_test_res1[j,] <- sprt_test(ratio_norm(mu1))
  L <- L_n_norm(n, mu1)
  NP_test_res1[j] <- np_test(L)
}

mean(NP_test_res1)
colMeans(sprt_test_res1)

# Misspecified

sprt_test_res2 <- matrix(ncol = 2, nrow = B)
NP_test_res2 <- vector()

for(j in 1:B){
  sprt_test_res2[j,] <- sprt_test(ratio_norm(0.2))
  L <- L_n_norm(n, 0.2)
  NP_test_res2[j] <- np_test(L)
}

mean(NP_test_res2)
colMeans(sprt_test_res2)

#-------------------------------------------------------------------------------
# Alpha - spending
#-------------------------------------------------------------------------------

set.seed(9372)

n1 <- 80
n2 <- 123

for(i in 1:(2*B)) Lns[i] <- L_n_norm(n1, mu0)
gamma1 <- quantile(Lns,1 - alpha/2)

for(i in 1:(2*B)) Lns[i] <- L_n_norm(n2, mu0)
gamma2 <- quantile(Lns,1 - alpha/2)

np_test_seq <- function(L_n, gamma) as.numeric(L_n >= gamma)
gammas <- c(gamma1, gamma2)

# Under null
ns <- vector()
NP_test_res_seq0 <- vector()
for(j in 1:B){
  L <- L_n_norm(n1, mu0)
  test_res <- np_test_seq(L, gammas[1])
  if(test_res != 1) {
    L <- L_n_norm(n2-n1, mu0) * L
    test_res <- np_test_seq(L, gammas[2])
    ns[j] <- n2
  }
  else{
    ns[j] <- n1
  }
  NP_test_res_seq0[j] <- test_res
}

mean(NP_test_res_seq0)
mean(ns)

# Under alternative
NP_test_res_seq1 <- vector()
for(j in 1:B){
  L <- L_n_norm(n1, mu1)
  test_res <- np_test_seq(L, gammas[1])
  if(test_res != 1) {
    L <- L_n_norm(n2-n1, mu1) * L
    test_res <- np_test_seq(L, gammas[2])
    ns[j] <- n2
  }
  else{
    ns[j] <- n1
  }
  NP_test_res_seq1[j] <- test_res
}

# Under misspecification
NP_test_res_seq2 <- vector()
for(j in 1:B){
  L <- L_n_norm(n1, 0.2)
  test_res <- np_test_seq(L, gammas[1])
  if(test_res != 1) {
    L <- L_n_norm(n2-n1, 0.2) * L
    test_res <- np_test_seq(L, gammas[2])
    ns[j] <- n2
  }
  else{
    ns[j] <- n1
  }
  NP_test_res_seq2[j] <- test_res
}

mean(NP_test_res_seq2)
mean(ns)
