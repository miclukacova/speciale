# Reproducing the results of Sokolova

# Treatment effect under alternative (simpel mod simpel)
# Hvorfor bruger vi ikke bare SPRT her (?)

alpha <- 0.025
p_t <- 0.35
B <- 10^4
p_c <- 0.2
Nmax <- 200

sample_patient <- function(p_c, p_t, n) {
  pat_c <- rbinom(n,1,p_c)
  pat_t <- rbinom(n,1,p_t)
  pat_t - pat_c
}

lambda_star <- (p_t * (1 - p_c) - (1 - p_t) * p_c) / (p_t * (1 - p_c) + (1 - p_t) * p_c)
g <- ((p_t * (1 - p_c) * log(1 + lambda_star)) +  (1 - p_t) * p_c * log(1 - lambda_star))

# Necessary sample size
ns_e_n <- vector()
for(i in 1:B) ns_e_n[i] <- min(which(cumprod(1 + lambda_star * sample_patient(p_c, p_t, n)) >= 1/alpha), n)
mean(ns_e_n)
# Their expected sample size
e_n <- log(1 / alpha) / g

# Type I error and Power
rej1 <- rej0 <- 0
for(i in 1:B) rej1 <- rej1 + any(cumprod(1 + lambda_star * sample_patient(p_c, p_t, n)) >= 1/alpha)
for(i in 1:B) rej0 <- rej0 + any(cumprod(1 + lambda_star * sample_patient(p_c, p_c, n)) >= 1/alpha)
rej1/B; rej0/B

# For fun we compare to the sprt, I have a feeling that it is much better
sprt <- function(p_c, p_t, n, null = FALSE){
  X <- sample_patient(p_c, p_t, n)
  if(null) X <- sample_patient(p_c, p_c, n)
  f_1 <- p_t * (1 - p_c) * (X == 1) +  (1 - p_t) * p_c * (X == -1) + (X == 0) * ((1 - p_t) * (1 - p_c) + p_t * p_c)
  f_0 <-  p_c * (1 - p_c) * (X == 1 | X == -1) + (p_c * p_c + (1 - p_c)^2)* (X == 0)
  cumprod(f_1 / f_0)
}

# It is much better in terms of sample size
ns <- vector()
for(i in 1:B) {
  ratio <- sprt(p_c, p_t, 2*n)
  ns[i] <- min(which(ratio <= alpha | ratio >= 1/alpha))
}
mean(ns[ns != Inf])

# and in terms of Type I and Type II error
rejSPRT1 <- 0; rejSPRT0 <- 0
for(i in 1:B) {
  ratios <- sprt(p_c, p_t, n)
  ratios <- ratios[ratios > alpha]
  rejSPRT1 <- rejSPRT1 + any(ratios >= 1/alpha)
}

for(i in 1:B) {
  ratios <- sprt(p_c, p_t, n, null = TRUE)
  ratios <- ratios[ratios > alpha]
  rejSPRT0 <- rejSPRT0 + any(ratios >= 1/alpha)
}

rejSPRT1/B; rejSPRT0/B

# Obrien/Pocock "conventional" design
n_looks <- 4

look_times <- round(seq(Nmax/n_looks, Nmax, length.out = n_looks))
info_frac <- look_times/Nmax
obf_bounds <- NULL
gs_c <- NULL
xT_null_mat <- matrix(stats::rbinom(B * Nmax, 1, p_c),
                      nrow = B)
xC_null_mat <- matrix(stats::rbinom(B * Nmax, 1, p_c),
                      nrow = B)
xT_alt_mat <- matrix(stats::rbinom(B * Nmax, 1, p_t),
                     nrow = B)
xC_alt_mat <- matrix(stats::rbinom(B * Nmax, 1, p_c),
                     nrow = B)
cumT_null <- t(apply(xT_null_mat, 1, cumsum))[, look_times,
                                              drop = FALSE]
cumC_null <- t(apply(xC_null_mat, 1, cumsum))[, look_times,
                                              drop = FALSE]
cumT_alt <- t(apply(xT_alt_mat, 1, cumsum))[, look_times,
                                            drop = FALSE]
cumC_alt <- t(apply(xC_alt_mat, 1, cumsum))[, look_times,
                                            drop = FALSE]
nn_looks <- matrix(look_times, nrow = B, ncol = length(look_times),
                   byrow = TRUE)

z_null <- evalinger:::.z_matrix(cumT_null, cumC_null, nn_looks)
z_alt <- evalinger:::.z_matrix(cumT_alt, cumC_alt, nn_looks)

m_null <- apply(z_null * matrix(sqrt(info_frac), nrow = B, ncol = length(info_frac), byrow = TRUE), 1, max)
gs_c <- as.numeric(stats::quantile(m_null, probs = 1 - alpha, names = FALSE))
obf_bounds <- gs_c/sqrt(info_frac)

rej_null <- logical(B)
rej_alt <- logical(B)
stop_null <- rep(Nmax, B)
stop_alt <- rep(Nmax, B)
for (scenario in c("null", "alt")) {
  if (scenario == "null") {
    cumT <- cumT_null
    cumC <- cumC_null
  }
  else {
    cumT <- cumT_alt
    cumC <- cumC_alt
  }
  z <- if (scenario == "null") z_null else z_alt
  bounds_mat <- matrix(obf_bounds, nrow = B,
                       ncol = length(look_times), byrow = TRUE)
  rossed <- z >= bounds_mat
  tmp <- evallinger:::.store_results(crossed, look_times, Nmax, scenario,
                                     B, rej_null, rej_alt, stop_null, stop_alt)
  rej_null <- tmp$rej_null
  rej_alt <- tmp$rej_alt
  stop_null <- tmp$stop_null
  stop_alt <- tmp$stop_alt
}



rows <- lapply(methods, function(m) {
  data.frame(method = m, null_rej = mean(raw[[m]]$null_rej),
             alt_rej = mean(raw[[m]]$alt_rej), avg_n_null = mean(raw[[m]]$null_stop),
             avg_n_alt = mean(raw[[m]]$alt_stop), stringsAsFactors = FALSE)
})
results <- do.call(rbind, rows)
structure(list(results = results, design = list(p_C = p_C,
                                                p_T_alt = p_T_alt, Nmax = Nmax, n_looks = n_looks, alpha = alpha,
                                                nrep = B, lambda = lambda, gs_c = gs_c, bayes_thresh = if (exists("bayes_thresh")) bayes_thresh else NULL),
               raw = raw), class = "ecomparison")

#-------------------------------------------------------------------------------
# Their example script
#-------------------------------------------------------------------------------

library(evalinger)

# 1. Design: find the optimal betting fraction
design <- edesign_binary(p_C = 0.30, delta = 0.15, alpha = 0.025, power = 0.80)
print(design)

# 2. Construct an e-process from trial data
x_T <- rbinom(150, 1, 0.45)  # treatment outcomes
x_C <- rbinom(150, 1, 0.30)  # control outcomes
ep <- eprocess_binary(x_T, x_C, lambda = design$lambda, alpha = 0.025)
plot(ep)

# 3. Real-time monitoring (batch-by-batch)
mon <- emonitor(alpha = 0.025, lambda = design$lambda)
mon <- update(mon, x_T = x_T[1:50],  x_C = x_C[1:50])   # interim 1
mon <- update(mon, x_T = x_T[51:100], x_C = x_C[51:100]) # interim 2
print(mon)

# 4. Always-valid confidence sequence
cs <- confseq_binary(x_T, x_C, alpha = 0.05)
plot(cs)

# 5. Compare methods head-to-head
cmp <- simulate_comparison(p_C = 0.30, p_T_alt = 0.45, Nmax = 200,
                           n_looks = 20, alpha = 0.025, nrep = 5000)
print(cmp)
plot(cmp)
