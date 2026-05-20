#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# LaTeX tables
#-------------------------------------------------------------------------------

res_200 %>%
  mutate(
    across(where(is.numeric), round, 3)
  ) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    caption = "Simulation results for $N_{\\max}=200$"
  ) %>%
  kable_styling(
    latex_options = c("hold_position")
  )

res_500 %>%
  mutate(
    across(where(is.numeric), round, 3)
  ) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    caption = "Simulation results for $N_{\\max}=500$"
  ) %>%
  kable_styling(
    latex_options = c("hold_position")
  )

#-------------------------------------------------------------------------------
# Simulations III
#-------------------------------------------------------------------------------

working_grid <- seq(0.35, 0.60, by = 0.025)

misspec_res <- map_dfr(working_grid, function(p_work) {

  tmp <- replicate(
    B,
    eprocess_misspecified(
      p_c = 0.30,
      p_t_true = 0.45,
      p_t_working = p_work,
      Nmax = 200,
      alpha = 0.025
    ),
    simplify = FALSE
  )

  tmp <- bind_rows(tmp)

  tibble(
    p_working = p_work,
    mean_n = mean(tmp$n),
    power = mean(tmp$reject)
  )
})


p1 <- ggplot(
  misspec_res,
  aes(p_working, mean_n)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = expression("Working alternative " ~ p[t]^"*"),
    y = "Mean stopping time",
    title = "Effect of misspecification on stopping time"
  ) +
  theme_bw()


p2 <- ggplot(
  misspec_res,
  aes(p_working, power)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    x = expression("Working alternative " ~ p[t]^"*"),
    y = "Power",
    title = "Effect of misspecification on power"
  ) +
  theme_bw()

p1 / p2


##-------------------------------------------------------------------------------
## The e-process defined in the article
#
## The optimal lambda star and corresponding g
#lambda_star <- (p_t * (1 - p_c) - (1 - p_t) * p_c) / (p_t * (1 - p_c) + (1 - p_t) * p_c)
#g <- ((p_t * (1 - p_c) * log(1 + lambda_star)) +  (1 - p_t) * p_c * log(1 - lambda_star))
#
## Necessary sample size
#ns_e_n <- vector()
#for(i in 1:B) ns_e_n[i] <- min(which(cumprod(1 + lambda_star * sample_patient(p_c, p_t, Nmax)) >= 1/alpha), Nmax)
#mean(ns_e_n)
#
## Their expected sample size
#e_n <- log(1 / alpha) / g
## Type I and type II error
#e_process(p_t, p_c, Nmax)
#
#
##-------------------------------------------------------------------------------
## The SPRT
## It is much better in terms of sample size
#ns <- vector()
#for(i in 1:B) {
#  ratio <- sprt(p_c, p_t, Nmax)
#  ns[i] <- min(which(ratio <= alpha | ratio >= 1/alpha), Nmax)
#}
#mean(ns)
#
## Type I and type II error
#sprt_error(p_t, p_c, Nmax)
#
##-------------------------------------------------------------------------------
## Obrien/Pocock "conventional" design
#
#n_looks <- 4
## Change this so it corresponds to significance level and n_looks
#alphas_book <- getDesignGroupSequential(kMax = 4, alpha = 0.05, sided = 1, typeOfDesign = "OF")$criticalValues
#alphas_soko <- alphas_soko(p_c = p_c, p_t = p_t, n_looks = n_looks, Nmax = Nmax, alpha = alpha)
#
#res1 <- obrien_gs(p_c = p_c, p_t = p_t, n_looks = n_looks, Nmax = Nmax, alphas = alphas_soko)
#res2 <- obrien_gs(p_c = p_c, p_t = p_t, n_looks = n_looks, Nmax = Nmax, alphas = alphas_book)
#
#
#
