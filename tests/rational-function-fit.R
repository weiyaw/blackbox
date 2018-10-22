## Proposal distribution and number of iterations depend on the cooling schedule and problem at hand.

######## We first fit an unconstrained model. ########
## Visualise the dataset
load_all()
plot(y ~ x, data = hyper, main = "A simulated dataset from a hyperbolic tangent function",
     cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)

plot(y ~ x, data = hyper_out, main = "A simulated dataset with outliers",
     cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)

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
rm1_rss <- SMCSA(loss, rnorm_rw, starting, log_schedule, 6000, 500, FALSE, TRUE)
## Optimise with Newton-Gauss
rm2_rss <- nls(y ~ (a + b*x) / (1 + c*x + d*I(x^2)), hyper, start = rar.hyper)

## Plot the results
plot_rat(rm1_rss$theta[1:2], c(1, rm1_rss$theta[3:4]))
plot_rat(coef(rm2_rss)[1:2], c(1, coef(rm2_rss)[3:4]), xlim = c(0, 6))


######## Now we fit a constrained model. ########
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
rm1_crss <- SMCSA(loss, rtnorm_rw, starting, recip_schedule, 1000, 100, FALSE, TRUE)

## Optimise with SMCSA and 1000 iterations. Log schedule is preferable if
## running for an extended period.
rm2_crss <- SMCSA(loss, rtnorm_rw, starting, log_schedule, 6000, 500, FALSE, TRUE)

## Plot the results
plot_rat(rm1_crss$theta[1:2], c(1, rm1_crss$theta[3:4]))
plot_rat(rm2_crss$theta[1:2], c(1, rm2_crss$theta[3:4]))


######## Now we fit a constrained model and 'hyper_out'. ########
## Fit a rational model f(x) = (a + b*x) / (1 + c*x + d*x^2), subject to constraint
## Construct loss functions, indicator, and proposal distribution
loss <- rational_rss_fac(hyper_out$y, hyper_out$x, 2L, 2L)
indicator <- is_monotones_fac(2L, 2L, min(hyper_out$x), max(hyper_out$x), TRUE)
rtnorm_rw <- rtnorm_rw_fac(indicator, sd = 1, fraction = 0.97)

## Get starting values using the rearrange and regress method
rar.hyper_out <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), hyper_out))
names(rar.hyper_out) <- letters[1:4]

## Randomly generate 1000 starting values arounf rar.hyper_out
set.seed(1)
starting <- get_feasibles(rar.hyper_out, indicator, 3000)

## Optimise with SMCSA and 100 iterations
rm1_crss_out <- SMCSA(loss, rtnorm_rw, starting, recip_schedule, 1000, 400, FALSE, TRUE)

## Optimise with SMCSA and 1000 iterations. Log schedule is preferable if
## running for an extended period.
rm2_crss_out <- SMCSA(loss, rtnorm_rw, starting, log_schedule, 6000, 500, FALSE, TRUE)

## Plot the results
plot_rat(rm1_crss_out$theta[1:2], c(1, rm1_crss_out$theta[3:4]))
plot_rat(rm2_crss_out$theta[1:2], c(1, rm2_crss_out$theta[3:4]))


######## Now we fit a unconstrained model with Tukey biweight and 'hyper_out'. ########
## Fit a rational model f(x) = (a + b*x) / (1 + c*x + d*x^2), subject to constraint
## Construct loss functions, indicator, and proposal distribution
loss <- rational_biweight_fac(hyper_out$y, hyper_out$x, 2L, 2L)

## Get starting values using the rearrange and regress method
rar_hyper_out <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), hyper_out))
names(rar_hyper_out) <- letters[1:4]

## Randomly generate 1000 starting values arounf rar.hyper_out
set.seed(1)
starting <- matrix(rep(rar_hyper_out, times = 6000), 4, 6000) + rcauchy(4 * 6000, scale = 1)

## Optimise with SMCSA and 100 iterations
rnorm_rw <- rnorm_rw_fac(sd = 2, fraction = 0.97)
rm1_biw_out <- SMCSA(loss, rnorm_rw, starting, recip_schedule, 1000, 400, FALSE, TRUE)

## Optimise with SMCSA and 1000 iterations. Log schedule is preferable if
## running for an extended period.
rnorm_rw <- rnorm_rw_fac(sd = 2, fraction = 0.99)
rm2_biw_out <- SMCSA(loss, rnorm_rw, starting, log_schedule, 6000, 500, FALSE, TRUE)

## Plot the results
plot_rat(rm1_biw_out$theta[1:2], c(1, rm1_biw_out$theta[3:4]))
plot_rat(rm2_biw_out$theta[1:2], c(1, rm2_biw_out$theta[3:4]))


######## Now we fit a constrained model with Tukey biweight and 'hyper_out'. ########
## Fit a rational model f(x) = (a + b*x) / (1 + c*x + d*x^2), subject to constraint
## Construct loss functions, indicator, and proposal distribution
loss <- rational_biweight_fac(hyper_out$y, hyper_out$x, 2L, 2L)
indicator <- is_monotones_fac(2L, 2L, min(hyper_out$x), max(hyper_out$x), TRUE)

## Get starting values using the rearrange and regress method
rar_hyper_out <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), hyper_out))
names(rar_hyper_out) <- letters[1:4]

## Randomly generate 1000 starting values arounf rar.hyper_out
set.seed(1)
starting <- get_feasibles(rar_hyper_out, indicator, 3000)

## Optimise with SMCSA and 100 iterations
rtnorm_rw <- rtnorm_rw_fac(indicator, sd = 2, fraction = 0.97)
rm1_cbiw_out <- SMCSA(loss, rtnorm_rw, starting, recip_schedule, 1000, 400, FALSE, TRUE)

## Optimise with SMCSA and 1000 iterations. Log schedule is preferable if
## running for an extended period.
rtnorm_rw <- rtnorm_rw_fac(indicator, sd = 2, fraction = 0.99)
rm2_cbiw_out <- SMCSA(loss, rtnorm_rw, starting, log_schedule, 3000, 500, FALSE, TRUE)

## Plot the results
plot_rat(rm1_cbiw_out$theta[1:2], c(1, rm1_cbiw_out$theta[3:4]), col = 'red')
plot_rat(rm2_cbiw_out$theta[1:2], c(1, rm2_cbiw_out$theta[3:4]), col = 'blue')


