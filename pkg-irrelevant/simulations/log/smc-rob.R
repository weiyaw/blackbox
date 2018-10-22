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
## smc.hyper <- SMCSA(objective, proposalMonoFac(indicator), rar.m.hyper,
##                    logSchedule, 10000, 200)

prop <- list()
for (i in 1:100) {
    prop[[i]] <- proposalMonoFac(indicator)
}
smc.hyper <- mcmapply(SMCSA,
                      proposal = prop,
                      MoreArgs = list(objf = objective,
                                      schedule = logSchedule,
                                      starting = rar.m.hyper,
                                      N = 1000,
                                      iter = 400),
                      SIMPLIFY = FALSE,
                      mc.cores = 15)
saveRDS(smc.hyper, "./simulations/log/smc-log-hyper-rob.rds")
