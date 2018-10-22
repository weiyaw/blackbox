## library(microbenchmark)
## microbenchmark(multiSA(objf, proposal_fac(), starting, schedule, 100, 100))

## system.time(multiSA(objf, proposal_fac(), starting, schedule, 100, 100))

## library(profvis)
## profvis( multiSA(objf, proposal_fac(), starting, schedule, 100, 100))

rm(list=ls())
hyper <- readRDS("./simulations/hyper.rds")
plot(y ~ x, hyper)

source("generic.R")

indicator <- indicatorAllFac(2L, 2L, min(hyper$x), max(hyper$x), TRUE)
objective <- with(hyper, rationalResidsFac(y, x, 2L, 2L))
rar.hyper <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), hyper))
names(rar.hyper) <- letters[1:4]
rar.m.hyper <- genMonoRats(rar.hyper, indicator, 1000)

rm0 <- nls(y ~ (a + b*x) / (1 + c*x + d*I(x^2)), hyper, start = rar.hyper)

rm1 <- multiSA(objective, proposalMonoFac(indicator), rar.m.hyper, recipSchedule,
               1000, 400, TRUE, TRUE)

rm2 <- SMCSA(objective, proposalMonoFac(indicator), rar.m.hyper, recipSchedule,
             1000, 400, TRUE, TRUE)

rm3 <- CEPSO(rar.m.hyper[,1:100], objective, coef(rm0), indicator, print = TRUE)

### ROBUST ###
rm(list = ls())
hyper <- readRDS("./simulations/hyper-rob.rds")
plot(y ~ x, hyper.rob)

source("generic.R")

indicator <- indicatorAllFac(2L, 2L, min(hyper$x), max(hyper$x), TRUE)
objective <- with(hyper, rationalBiweightFac(y, x, 2L, 2L))
rar.hyper <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), hyper))
names(rar.hyper) <- letters[1:4]
rar.m.hyper <- genMonoRats(rar.hyper, indicator, 1000)

rm0 <- nls(y ~ (a + b*x) / (1 + c*x + d*I(x^2)), hyper, start = rar.hyper)

rm1 <- multiSA(objective, proposalMonoFac(indicator), rar.m.hyper, recipSchedule,
               1000, 100, TRUE, TRUE)

rm2 <- SMCSA(objective, proposalMonoFac(indicator), rar.m.hyper, recipSchedule,
             1000, 100, TRUE, TRUE)

rm3 <- CEPSO(rar.m.hyper[,1:100], objective, coef(rm0), indicator, N = 1000, print = TRUE)




indicatorAll <- indicatorAllFac(2, 1, xmax = 20)
rm2 <- multiSA(rationalResidsFac(y, x, 2L, 1L), proposalMonoFac(indicator = indicatorAll), starting, schedule, diagnostic = TRUE)

rm3 <- SMCSA(rationalResidsFac(y, x, 2L, 1L), proposalMonoFac(indicator = indicatorAll), starting, schedule, 1000, 100)

starting <- cbind(1:3, c(1,1,1), 4:2)
indicatorAll <- indicatorAllFac(2, 1)
rm4 <- SMCSA(rationalResidFac(y, x, 2L, 1L), proposalMonoFac(indicator = indicatorAll), starting, schedule, 1000, 100)

plotRat(rm4$theta[1:2], c(1, rm4$theta[3]))

## Test residual function
source("generic-subor.R")
rationalResid <- rationalResidFac(y, x, 3L, 2L)
rationalResids <- rationalResidsFac(y, x, 3L, 2L)
theta <- matrix(sample.int(50, 5000, T), 5)
rationalResid(theta) == rationalResids(theta)
microbenchmark(rationalResid(theta), rationalResids(theta))

## Test CEPSO
rm(list = ls())
source("generic.R")
starting1 <- matrix(sample.int(50, 20, T), 5)

delta <- function(x) {colSums(as.matrix(x))}


indicator1 <- function(x) {rep(FALSE, NCOL(x))} # An always false indicator
indicator2 <- function(x) {rep(c(T, F), len = NCOL(x))}   # A random indicator

CEPSO(starting1, delta, c(0,0,0,0,0), indicator2, 3)

sw1 <- initSwarm(starting1, delta, indicator1)
sw2 <- initSwarm(matrix(nrow=5, ncol=0), delta, indicator2)
asw <- exchange(sw1, sw2, 2)

sw1e <- exploit(asw[[1]], matrix(1, 5, 1))
sw2e <- explore(asw[[2]], matrix(1, 5, 1))

sw1a <- assess(sw1e, delta, indicator2)
sw2a <- assess(sw2e, delta, indicator2)

exchange(sw1a, sw2a, 2)

## Test CEPSO with some real data
rm(list=ls())
source("generic.R")
x <- seq(0,20,0.6)
y <- (3*x + 1) / (1.5*x + 4) + rnorm(x, sd = 0.1)
plot(y ~ x)
indicatorAll <- indicatorAllFac(2, 1, xmax = 20)
starting <- cbind(c(1,2,3), c(1,1,1), c(4,3,2))

nl1 <- nls(y ~ (a + b*x) / (1 + c*x), start = list(a=1, b=1, c=2))
crossprod(resid(nl1))
pilot <- coef(nl1)

rm1 <- multiSA(rationalResidsFac(y, x, 2L, 1L), proposalFac(), starting,
               schedule, 2000, 100)

rm2 <- multiSA(rationalResidsFac(y, x, 2L, 1L),
               proposalMonoFac(indicator=indicatorAll,
                               dist.para=list(mean=0, sd=1)),
               starting, schedule, 1000, 100, diagnostic=TRUE)

rm3 <- SMCSA(rationalResidsFac(y, x, 2L, 1L),
             proposalMonoFac(indicator=indicatorAll),
             starting, schedule, 1000, 100)

rm4 <- CEPSO(starting, rationalResidsFac(y, x, 2L, 1L), pilot, indicatorAll, 10, 100, 1000)


## Test with
rm(list=ls())
source("generic.R")
simlidar <- read.csv("data/simlidar.csv")
x <- simlidar$s.range
y <- simlidar$s.y
plot(y ~ x)
indicatorAll <- indicatorAllFac(3, 2, xmax = max(x), xmin = min(x))

nl1 <- nls(y ~ (a1 + a2*x + a3*x^2) / (1 + b1*x + b2*x^2),
           start=list(a1=-0.5, a2=2, a3=1.5, b1=-0.5, b2=0.1))
crossprod(resid(nl1))
pilot <- coef(nl1)
names(pilot) <- NULL
plotRat(pilot[1:3], c(1, pilot[4:5]))

starting <- genMonoRats(cbind(c(1,2,3,4,5), c(3,2,1,2,3), c(5,4,3,2,1),
                              pilot, deparse.level=0),
                        indicatorAll, 20)

## No constraint
rm1 <- multiSA(rationalResidsFac(y, x, 3L, 2L), proposalFac(), starting,
               schedule, 2000, 100)
plotRat(rm1$theta[1:3], c(1, rm1$theta[4:5]), col = 'red')

## Monotone constrained
source("generic.R")
rm2 <- multiSA(rationalResidsFac(y, x, 3L, 2L),
               proposalMonoFac(indicator=indicatorAll,
                               dist.para=list(mean=0, sd=1),
                               fraction=0.97),
               starting, schedule, 2000, 100)
plotRat(rm2$theta[1:3], c(1, rm2$theta[4:5]), col = 'blue')

## Monotone constrained
rm3 <- SMCSA(rationalResidsFac(y, x, 3L, 2L),
             proposalMonoFac(indicator=indicatorAll,
                             dist.para=list(mean=0, sd=2),
                             fraction=0.95),
             starting, schedule, 3000, 100)
plotRat(rm3$theta[1:3], c(1, rm3$theta[4:5]), col = 'green')

## Monotone constrained
rm4 <- CEPSO(starting, rationalResidsFac(y, x, 3L, 2L), pilot,
             indicatorAll, 20, 150, 5000)
plotRat(rm4$theta[1:3], c(1, rm4$theta[4:5]), col = 'brown')

