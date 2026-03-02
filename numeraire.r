Z <- rnorm(10^6,1,1)
mean(2*dcauchy(Z)/(dnorm(Z,-1,1) + dnorm(Z,1,1)))
mean(dcauchy(Z)/max(dnorm(Z,-1,1), dnorm(Z,1,1)))

Z <- rnorm(10^6,-1,1)
mean(2*dcauchy(Z)/(dnorm(Z,-1,1) + dnorm(Z,1,1)))
mean(dcauchy(Z)/max(dnorm(Z,-1,1), dnorm(Z,1,1)))

Z <- rcauchy(1)
mean(2*dcauchy(Z)/(dnorm(Z,-1,1) + dnorm(Z,1,1)))
mean(dcauchy(Z)/max(dnorm(Z,-1,1), dnorm(Z,1,1)))

hist(Z[Z>1000000])
