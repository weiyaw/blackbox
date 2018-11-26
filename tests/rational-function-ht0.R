## Proposal distribution and number of iterations depend on the cooling schedule and problem at hand.

######## We first fit an unconstrained model. ########
## Visualise the dataset
load_all()
scaled_ht0 <- as.data.frame(scale(HT0, FALSE, apply(HT0, 2, max)))
plot(y ~ x, data = scaled_ht0,
     main = "A simulated dataset from a hyperbolic tangent function",
     cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)

## Fit a rational model f(x) = (a + b*x) / (1 + c*x + d*x^2)
## Construct loss functions
loss <- rational_rss_fac(scaled_ht0$y, scaled_ht0$x, 2L, 2L)

## Get starting values using the rearrange and regress method
rar_hyper <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), scaled_ht0))
names(rar_hyper) <- letters[1:4]

## Randomly generate 1000 starting values arounf rar_hyper
set.seed(1)
starting <- matrix(rep(rar_hyper, times = 1000), 4, 1000) +
    rcauchy(4000, scale = 1)

## Optimise with SMCSA and 1000 iterations
set.seed(1)
rnorm_rw <- rnorm_rw_fac(sd = 2, fraction = 0.99)
rm1_rss <- SMCSA(loss, rnorm_rw, starting, log_schedule, 6000, 500, FALSE, TRUE)

## Optimise with Newton-Gauss
rm2_rss <- nls(y ~ (a + b*x) / (1 + c*x + d*I(x^2)), hyper, start = rar_hyper)


## Optimise with SAA and 5000000 iterations
set.seed(1)
rm3_rss <- SAA(loss, NULL, as.matrix(rar_hyper), NULL, 1, 5000000, FALSE, TRUE)


## Plot the results
plot_rat(rm1_rss$state[1:2], c(1, rm1_rss$state[3:4]))
plot_rat(coef(rm2_rss)[1:2], c(1, coef(rm2_rss)[3:4]), xlim = c(0, 6))


######## Now we fit a constrained model. ########
## Fit a rational model f(x) = (a + b*x) / (1 + c*x + d*x^2), subject to constraint
## Construct loss functions and an indicator
plot(y ~ x, data = scaled_ht0,
     main = "A simulated dataset from a hyperbolic tangent function",
     cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)
loss <- rational_rss_fac(hyper$y, hyper$x, 2L, 2L)
is_monotone <- is_monotones_fac(2L, 2L, min(hyper$x), max(hyper$x), TRUE)

## Get 5000 starting values using the rearrange and regress method
rar_hyper <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), hyper))
names(rar_hyper) <- letters[1:4]
set.seed(10)
starting <- get_feasibles(rar_hyper, is_monotone, 5000)

## SMCSA, reciprocal, fraction = 0.97, sd = 1, 500 iter, N = 5000
## Use a more aggresive proposal
load_all()
set.seed(11)
rtnorm_rw <- rtnorm_rw_comp_fac(is_monotone, sd = 1, fraction = 0.97, n_pertub = 2)
rm1_crss <- SMCSA(loss, rtnorm_rw, starting, recip_schedule, 5000, 500, TRUE, TRUE)
plot_rat(rm1_crss$state[1:2], c(1, rm1_crss$state[3:4]))

## SMCSA, log, fraction = 0.97, sd = 1, 1000 iter
## running for an extended period.
set.seed(12)
rtnorm_rw <- rtnorm_rw_fac(is_monotone, sd = 1, fraction = 0.99)
rm2_crss <- SMCSA(loss, rtnorm_rw, starting, log_schedule, 1000, 1000, FALSE, TRUE)
plot_rat(rm2_crss$state[1:2], c(1, rm2_crss$state[3:4]))

## CEPSO, OLS as pilot estimate
set.seed(1)
rar_hyper <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), hyper))
names(rar_hyper) <- letters[1:4]
ht0_pilot <- coef(nls(y ~ (a + b*x) / (1 + c*x + d*I(x^2)), hyper, start = rar_hyper))
plot_rat(ht0_pilot[1:2], c(1, ht0_pilot[3:4]))
set.seed(13)
rm3_crss <- CEPSO(starting, loss, ht0_pilot, is_monotone, 500, 1000, 5000, TRUE)
plot_rat(rm3_crss$theta[1:2], c(1, rm3_crss$theta[3:4]))



## Fit a rational model f(x) = (a + b*x + c*x^2) / (1 + d*x + e*x^2), subject to constraint
## Construct loss functions and an indicator
plot(y ~ x, data = scaled_ht0,
     main = "A simulated dataset from a hyperbolic tangent function",
     cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)
loss <- rational_rss_fac(scaled_ht0$y, scaled_ht0$x, 3L, 2L)
is_monotone <- is_monotones_fac(3L, 2L, min(scaled_ht0$x), max(scaled_ht0$x), TRUE)

## Get starting values using the rearrange and regress method
rar_hyper <- coef(lm(y ~ x + I(x^2) + I(-x*y) + I(-x^2*y), scaled_ht0))
names(rar_hyper) <- letters[1:5]

## Randomly generate 5000 starting values arounf rar_hyper
set.seed(10)
starting <- get_feasibles(rar_hyper, is_monotone, 1000)

## SMCSA, reciprocal, fraction = 0.97, sd = 1, 200 iter, N = 1000, 2 comp
## Use a more aggresive proposal
load_all()
set.seed(11)
rtnorm_rw <- rtnorm_rw_comp_fac(is_monotone, sd = 1, fraction = 0.97, n_pertub = 3)
system.time(rm1_crss <- SMCSA(loss, rtnorm_rw, starting, recip_schedule, 1000, 200, TRUE, TRUE))
plot_rat(rm1_crss$state[1:3], c(1, rm1_crss$state[4:5]))

## SMCSA, log, fraction = 0.97, sd = 1, 1000 iter
## running for an extended period.
set.seed(12)
rtnorm_rw <- rtnorm_rw_comp_fac(is_monotone, sd = 1, fraction = 0.998, n_pertub = 2)
system.time(rm2_crss <- SMCSA(loss, rtnorm_rw, starting, log_schedule, 1000, 3000, TRUE, TRUE))
plot_rat(rm2_crss$state[1:2], c(1, rm2_crss$state[3:4]))

## CEPSO, OLS as pilot estimate
ht0_pilot <- coef(nls(y ~ (a + b*x + c*x^2) / (1 + d*x + e*I(x^2)), hyper,
                      start = rar_hyper))
plot_rat(ht0_pilot[1:3], c(1, ht0_pilot[4:5]))
set.seed(13)
system.time(rm3_crss <- CEPSO(starting, loss, ht0_pilot, is_monotone, 2000, 10000, 5000, TRUE))
plot_rat(rm3_crss$theta[1:3], c(1, rm3_crss$theta[4:5]))

