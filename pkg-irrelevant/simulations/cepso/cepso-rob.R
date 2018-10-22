setwd("~/Desktop/local")
rm(list = ls())
source("generic.R")
library("parallel")

## 30 samples from f(x) = -1.25x + 2.25x^2 -0.5x^3 + Normal(0, sd = 0.3)
## cubic <- readRDS("./simulations/cubic.RDS")
## plot(y ~ x, cubic)

## 30 samples from f(x) = 1 + tanh(x - 3) + Normal(0, sd = 0.1)
## outliers at 2nd and last 2nd samples (y = 2, 0)
hyper <- readRDS("./simulations/hyper-rob.rds")
## plot(y ~ x, hyper)

## Fit a rational model f(x) = (a + b*x) / (1 + c*x + d*x^2)
## DO NOT REUSE PROPOSAL. THERE IS A COUNTER WITHIN THAT NEEDS TO BE RESET.

indicator <- indicatorAllFac(2L, 2L, min(hyper$x), max(hyper$x), TRUE)
objective <- with(hyper, rationalBiweightFac(y, x, 2L, 2L))
rar.hyper <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), hyper))
names(rar.hyper) <- letters[1:4]
rar.m.hyper <- genMonoRats(rar.hyper, indicator, 1000)

un.hyper <- nls(y ~ (a + b*x) / (1 + c*x + d*I(x^2)), hyper, start = rar.hyper)
## cepso.hyper <- CEPSO(rar.m.hyper, objective, coef(un.hyper), indicator,
##                      200, 2000, 1000)

iter <- rep(2000, len = 100)
cepso.hyper <- mcmapply(CEPSO,
                      iter = iter,
                      MoreArgs = list(starting = rar.m.hyper,
                                      objf = objective,
                                      pilot = coef(un.hyper),
                                      indicator = indicator,
                                      neighbour = 200,
                                      N = 1000),
                      SIMPLIFY = FALSE,
                      mc.cores = 15)
saveRDS(cepso.hyper, "./simulations/cepso/cepso-hyper-rob.rds")

