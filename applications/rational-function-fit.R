######## We first fit an unconstrained model. ########
rm(list = ls())
## Get data
hyper <- read.csv("../data/hyper.csv")

## Get the algorithms
source("../algorithms/SA-based.R")

## Get helper functions to construct proposal distributions and starting values
## generators
source("../algorithms/hyperparameters.R")

## Get all the functions required for fitting rational functions
source("./rational-function.R")

## Visualise the dataset
plot(y ~ x, data = hyper)

## Fit a rational model f(x) = (a + b*x) / (1 + c*x + d*x^2)
## Construct loss functions and proposal distribution
loss <- rational_rss_fac(hyper$y, hyper$x, 2L, 2L)
rnorm_rw <- rnorm_rw_fac(sd = 2, fraction = 0.99)

## Get starting values using the rearrange and regress method
rar.hyper <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), hyper))
names(rar.hyper) <- letters[1:4]

## Randomly generate 1000 starting values arounf rar.hyper
set.seed(1)
starting <- matrix(rep(rar.hyper, times = 3000), 4, 3000) +
    rcauchy(12000, scale = 1)

## Optimise with SMCSA and 1000 iterations
set.seed(1)
rm1 <- SMCSA(loss, rnorm_rw, starting, log_schedule, 3000, 1000, FALSE, TRUE)
## Optimise with Newton-Gauss
rm2 <- nls(y ~ (a + b*x) / (1 + c*x + d*I(x^2)), hyper, start = rar.hyper)

## Plot the results
plot_rat(rm1$theta[1:2], c(1, rm1$theta[3:4]))
plot_rat(coef(rm2)[1:2], c(1, coef(rm2)[3:4]))


######## Now we fit a constrained model. ########
rm(list = ls())
## Get data
hyper <- read.csv("../data/hyper.csv")

## Get the algorithms
source("../algorithms/SA-based.R")

## Get helper functions to construct proposal distributions and starting values
## generators
source("../algorithms/hyperparameters.R")

## Get all the functions required for fitting rational functions
source("./rational-function.R")

## Fit a rational model f(x) = (a + b*x) / (1 + c*x + d*x^2), subject to constraint
## Construct loss functions, indicator, and proposal distribution
loss <- rational_rss_fac(hyper$y, hyper$x, 2L, 2L)
indicator <- is_monotones_fac(2L, 2L, min(hyper$x), max(hyper$x), TRUE)
rtnorm_rw <- rtnorm_rw_fac(indicator, sd = 2, fraction = 0.98)

## Get starting values using the rearrange and regress method
rar.hyper <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), hyper))
names(rar.hyper) <- letters[1:4]

## Randomly generate 1000 starting values arounf rar.hyper
set.seed(1)
starting <- get_feasibles(rar.hyper, indicator, 3000)

## Optimise with SMCSA and 100 iterations
set.seed(1)
rm1 <- SMCSA(loss, rtnorm_rw, starting, recip_schedule, 3000, 100, FALSE, TRUE)

## Optimise with SMCSA and 1000 iterations. Log schedule is preferable if
## running for an extended period.
set.seed(1)
rm2 <- SMCSA(loss, rtnorm_rw, starting, log_schedule, 3000, 1000, FALSE, TRUE)

## Plot the results
plot_rat(rm1$theta[1:2], c(1, rm1$theta[3:4]))
plot_rat(rm2$theta[1:2], c(1, rm2$theta[3:4]))
