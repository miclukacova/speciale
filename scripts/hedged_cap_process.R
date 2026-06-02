#-------------------------------------------------------------------------------

library(tidyr)
library(ggplot2)
library(dplyr)

# Parameters
alpha = 0.05
m_true = 0.6
m_0 = 0.5
m_1 = 0.6
N = 300
c = 1 / 2
theta = 1 / 2
B = 10^4
gamma = 0.95

# Another way to choose N
N <- 50
pow <- 0
while(pow < 0.95) {
  z_ag <- qbinom(p = 1 - (alpha * gamma), size = N, prob = m_0)
  pow <- pbinom(z_ag, N, 0.6, lower.tail = FALSE)
  N <- N + 1
}

#-------------------------------------------------------------------------------
## The hedged capital process
#-------------------------------------------------------------------------------


# Calculating the predictable sequence
calculate_lambda_tilde <- function(X, alpha){

  n <- length(X)
  mu_hat <- (1 /2 + cumsum(X)) / (1:n +1)
  sigma_hat <- (1 / 4 + cumsum((X - mu_hat) ^ 2)) / (1:n + 1)
  sigma_hat_t_1 <- c(1 / 4, sigma_hat)[1:n]
  lambda_tilde <- sqrt(2 * log(2 / alpha) / (sigma_hat_t_1 * 1:n * log(1:n + 1)))
  lambda_tilde

}

# Data sampling functions
sample_bern <- function(N, m_true) rbinom(N, 1, m_true)

# The hedged capital process
hedged_cap_proc <- function(m_0, m_true, c, sample_data, N, theta) {

  # Sample data
  X <- sample_data(N, m_true)

  # Calculate the lambda's
  lambda_tilde <- calculate_lambda_tilde(X, alpha = alpha)
  lambda_plus <- pmin(abs(lambda_tilde), c / m_0)
  lambda_minus <- pmin(abs(lambda_tilde), c / (1 -m_0))

  # Calculate the capital process
  K_plus <- cumprod(1 + lambda_plus * (X - m_0))
  K_minus <- cumprod(1 - lambda_minus * (X - m_0))

  return(pmax(theta * K_plus, (1 - theta) * K_minus))
}

HCP_test <- function(HCP, alpha) {
  test_res <- HCP > 1/ alpha
  test_fut <- HCP < alpha
  if(any(test_res)){
    ESS <- which((test_res + test_fut) == 1)[1]
  } else {
    ESS <- length(HCP)
  }
  return(c(Reject = any(test_res[1:ESS]), ESS = ESS))
}

# Bootstrap
HCP_test_res <- matrix(nrow = B, ncol = 2)
colnames(HCP_test_res) <- c("Reject", "ESS")

for(i in 1:B) {
  HCP <- hedged_cap_proc(m_0 = m_0, sample_data = sample_bern, N = N, c = c, theta = theta, m_true)
  HCP_test_res[i, ] <- HCP_test(HCP, alpha = alpha)
}

colMeans(HCP_test_res)


# Plot the process

library(ggplot2)

HCP_mat <- replicate(
  10,
  hedged_cap_proc(
    m_0 = m_0,
    c = c,
    sample_data = sample_bern,
    N = N,
    theta = theta,
    m_true = m_true
  )
)

df <- data.frame(
  time = rep(1:N, 10),
  path = rep(1:10, each = N),
  HCP = as.vector(HCP_mat)
)

ggplot(df, aes(time, HCP, group = path)) +
  geom_line(alpha = 0.4) +
  geom_hline(yintercept = 1 / 0.05,
             colour = "red",
             linetype = 2) +
  geom_hline(yintercept = 0.05,
             colour = "red",
             linetype = 2) +
  scale_y_log10() +
  theme_minimal()


#-------------------------------------------------------------------------------
## SPRT
#-------------------------------------------------------------------------------

beta <- alpha
gamma0 <- beta / (1 - alpha)
gamma1 <- (1 - alpha) / beta

ratio_norm <- function(m_true, m_1) {
  ratio_func <- function(){
    X <- rbinom(1, 1, m_true)
    prod(dbinom(X, 1, m_1) / dbinom(X, 1, m_0))
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
  return(c(L >= gamma1, i))
}


# Bootstrap
SPRT_test_res <- matrix(nrow = B, ncol = 2)
colnames(SPRT_test_res) <- c("Reject", "ESS")

for(i in 1:B) {
  SPRT_fun <- ratio_norm(m_true = m_true, m_1 = 0.6)
  SPRT_test_res[i, ] <-sprt_test(SPRT_fun)
}

colMeans(SPRT_test_res)

#-------------------------------------------------------------------------------
## GS t-test
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
## Holmes method
#-------------------------------------------------------------------------------

holmes_test <- function(N, m_true, m_0, z_ag, gamma) {

  X <- sample_bern(N, m_true)
  Q_n <- 1 - pbinom(q = z_ag - cumsum(X), size = seq(N-1, 0), prob = m_0)
  Reject = any(Q_n >= gamma)

  if(Reject) ESS <- which(Q_n >= gamma)[1]
  else ESS <- N

  return(c(Reject = Reject, ESS = ESS))
}


# Bootstrap
holmes_test_res <- matrix(nrow = B, ncol = 2)
colnames(holmes_test_res) <- c("Reject", "ESS")

for(i in 1:B) {
  holmes_test_res[i, ] <- holmes_test(N = N, m_true = m_0, m_0 = m_0,
                                      z_ag = z_ag, gamma = gamma)
}

colMeans(holmes_test_res)


#-------------------------------------------------------------------------------
## Comparisons
#-------------------------------------------------------------------------------

# Overvej hvordan du vælger n??

compare_tests <- function(m_true_grid, m1_grid, z_ag, gamma) {

  results <- vector("list", length(m_true_grid))

  for(g in seq_along(m_true_grid)) {

    m_true <- m_true_grid[g]

    HCP_res   <- matrix(NA, B, 2)
    HOLM_res  <- matrix(NA, B, 2)

    SPRT_res <- matrix(
      NA,
      nrow = B,
      ncol = 2 * length(m1_grid)
    )

    colnames(HCP_res)  <- c("Reject", "ESS")
    colnames(HOLM_res) <- c("Reject", "ESS")

    colnames(SPRT_res) <- as.vector(
      rbind(
        paste0("Reject_", m1_grid),
        paste0("ESS_", m1_grid)
      )
    )

    for(b in seq_len(B)) {

      # ----------------------
      # HCP
      # ----------------------
      HCP <- hedged_cap_proc(
        m_0 = m_0,
        m_true = m_true,
        c = c,
        sample_data = sample_bern,
        N = N,
        theta = theta
      )

      HCP_res[b, ] <- HCP_test(HCP, alpha)

      # ----------------------
      # HOLMES
      # ----------------------
      HOLM_res[b, ] <- holmes_test(
        N = N,
        m_true = m_true,
        m_0 = m_0,
        z_ag = z_ag,
        gamma = gamma
      )

      # ----------------------
      # SPRTs
      # ----------------------
      for(j in seq_along(m1_grid)) {

        SPRT_fun <- ratio_norm(
          m_true = m_true,
          m_1 = m1_grid[j]
        )

        SPRT_res[b, (2*j-1):(2*j)] <- sprt_test(SPRT_fun)
      }
    }

    out <- c(
      m_true = m_true,

      HCP_power  = mean(HCP_res[, "Reject"]),
      HCP_ESS    = mean(HCP_res[, "ESS"]),

      HOLM_power = mean(HOLM_res[, "Reject"]),
      HOLM_ESS   = mean(HOLM_res[, "ESS"])
    )

    for(j in seq_along(m1_grid)) {

      out[paste0("SPRT_power_", m1_grid[j])] <-
        mean(SPRT_res[, paste0("Reject_", m1_grid[j])])

      out[paste0("SPRT_ESS_", m1_grid[j])] <-
        mean(SPRT_res[, paste0("ESS_", m1_grid[j])])
    }

    results[[g]] <- as.data.frame(as.list(out))
  }

  dplyr::bind_rows(results)
}

res <- compare_tests(
  m_true_grid = seq(0.5, 0.8, by = 0.03),
  m1_grid = c(0.55, 0.60, 0.65, 0.70, 0.75),
  z_ag = z_ag,
  gamma = 0.95
)

res


# Power plot
power_df <- res |>
  pivot_longer(
    cols = c(HCP_power, HOLM_power, starts_with("SPRT_power")),
    names_to = "Method",
    values_to = "Power"
  )

ggplot(
  power_df,
  aes(
    x = m_true,
    y = Power,
    colour = Method
  )
) +
  geom_line() +
  geom_point() +
  theme_minimal()

# ESS plot
ess_df <- res |>
  pivot_longer(
    cols = c(HCP_ESS, HOLM_ESS, starts_with("SPRT_ESS")),
    names_to = "Method",
    values_to = "ESS"
  )

ess_df$Method <- gsub("SPRT_ESS_", "SPRT(", ess_df$Method)
ess_df$Method <- gsub("$", ")", ess_df$Method)
ess_df$Method <- gsub("HCP_ESS)", "HCP", ess_df$Method)
ess_df$Method <- gsub("HOLM_ESS)", "HOLM", ess_df$Method)

ggplot(
  ess_df,
  aes(
    x = m_true,
    y = ESS,
    colour = Method
  )
) +
  geom_line() +
  geom_point() +
  theme_minimal()


# Jeg havde to ideer til SPRT'en: 1) lave mixture, men det svarer jo egentlig bare til at vælge et p?
# 2) lave den adaptive process, men den er lidt mærkelig i det her scenarie, fordi at den ender med at sige
# at alternativ hypotesen slet ikke er mulig. Derfor gør vi lidt noget andet...

# Vi skal muligvis have gang i forskellige N'er...


# Simulate under nu
X <- rbinom(N, 1, prob = mu_1)

# Plug-in estimator
mu_hat <- cumsum(c(0.5,X[1:(N-1)])) / seq_len(N)

# Adaptive e-process
E <- cumprod(
  dbinom(X, 1,  mu_hat) / dbinom(X, 1, mu0)
)


