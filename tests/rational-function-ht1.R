## Proposal distribution and number of iterations depend on the cooling schedule and problem at hand.


load_all()

## fit a rational model f(x) = (a + b*x + c*x^2) / (1 + d*x + e*x^2), subject to constraint
## no scaling on the dataset.

## visualise the dataset
plot(y ~ x, data = ht1, main = "HT1", cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)

## construct the loss functions and an indicator
loss <- rational_biweight_fac(ht1$y, ht1$x, 3L, 2L, cst = 1)
is_monotone <- is_monotones_fac(3L, 2L, min(ht1$x), max(ht1$x), TRUE)

## get a starting value using the rearrange and regress method
rar_ht1 <- coef(lm(y ~ x + I(x^2) + I(-x*y) + I(-x^2*y), ht1))
names(rar_ht1) <- letters[1:5]

library(doParallel)
registerDoParallel(cores = 8)

recip_time <- system.time(recip_ht1 <- foreach(i = 1:40) %dopar% {
    ## randomly generate 1000 starting values arounf rar_ht1
    set.seed(300 + i)
    starting <- get_feasibles(rar_ht1, is_monotone, 1000, dist_para = list(scale = 2))

    ## SMCSA, reciprocal, fraction = 0.97, sd = 1, 500 iter, N = 3000, 2 comp
    ## Use a more aggresive proposal
    rtnorm_rw <- rtnorm_rw_comp_fac(is_monotone, sd = 1, fraction = 0.97, n_pertub = 2)
    runtime <- system.time(rm1 <- SMCSA(loss, rtnorm_rw, starting, recip_schedule, 3000, 1000, FALSE, TRUE))
    rm1$runtime <- runtime
    rm1
})
cat("RECIP DONE.", recip_time[3], "\n")
saveRDS(recip_ht1, "~/Dropbox/honours/extras/ht1/recip_ht1.rds")
## plot_rat(rm1$state[1:3], c(1, rm1$state[4:5]), col = 'red')

recip85_time <- system.time(recip85_ht1 <- foreach(i = 1:40) %dopar% {
    ## randomly generate 1000 starting values arounf rar_ht1
    set.seed(300 + i)
    starting <- get_feasibles(rar_ht1, is_monotone, 1000, dist_para = list(scale = 2))

    ## SMCSA, reciprocal, fraction = 0.97, sd = 1, 500 iter, N = 3000, 2 comp
    ## Use a more aggresive proposal
    rtnorm_rw <- rtnorm_rw_comp_fac(is_monotone, sd = 1, fraction = 0.97, n_pertub = 2)
    runtime <- system.time(rm1 <- SMCSA(loss, rtnorm_rw, starting, recip_schedule_0.85, 3000, 1000, FALSE, TRUE))
    rm1$runtime <- runtime
    rm1
})
cat("RECIP DONE.", recip85_time[3], "\n")
saveRDS(recip85_ht1, "~/Dropbox/honours/extras/ht1/recip85_ht1.rds")
## plot_rat(rm1$state[1:3], c(1, rm1$state[4:5]), col = 'red')



log_time <- system.time(log_ht1 <- foreach(i = 1:40) %dopar% {
    ## randomly generate 1000 starting values arounf rar_ht1
    set.seed(300 + i)
    starting <- get_feasibles(rar_ht1, is_monotone, 1000, dist_para = list(scale = 2))

    ## SMCSA, log, fraction = 0.998, sd = 1, 3000 iter, N = 1000, 2 comp
    rtnorm_rw <- rtnorm_rw_comp_fac(is_monotone, sd = 1, fraction = 0.998, n_pertub = 2)
    runtime <- system.time(rm2 <- SMCSA(loss, rtnorm_rw, starting, log_schedule, 1000, 3000, FALSE, TRUE))
    rm2$runtime <- runtime
    rm2
})
cat("LOG DONE.", log_time[3], "\n")
saveRDS(log_ht1, "~/Dropbox/honours/extras/ht1/log_ht1.rds")
## plot_rat(rm2$state[1:3], c(1, rm2$state[4:5]), col = 'blue')


ht1_pilot <- coef(nls(y ~ (a + b*x + c*x^2) / (1 + d*x + e*I(x^2)),
                        ht1, start = rar_ht1))

cepso_time <- system.time(cepso_ht1 <- foreach(i = 1:40) %dopar% {
    ## randomly generate 1000 starting values arounf rar_ht1
    set.seed(300 + i)
    starting <- get_feasibles(rar_ht1, is_monotone, 1000, dist_para = list(scale = 2))

    ## CEPSO, 2000 iter, N = 3000, neighbour = 3000/5
    ## OLS as pilot estimate
    runtime <- system.time(rm3 <- CEPSO(starting, loss, ht1_pilot, is_monotone, iter = 2000, N = 3000, verbose = TRUE))
    rm3$runtime <- runtime
    rm3
})
cat("CEPSO DONE.", cepso_time[3], "\n")
saveRDS(cepso_ht1, "~/Dropbox/honours/extras/ht1/cepso_ht1.rds")
## plot_rat(rm3$theta[1:3], c(1, rm3$theta[4:5]), col = 'red')


















## ######## Now we fit a constrained model and 'hyper_out'. ########
## ## Fit a rational model f(x) = (a + b*x) / (1 + c*x + d*x^2), subject to constraint
## plot(y ~ x, data = hyper_out, main = "A simulated dataset with outliers",
##      cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)

## ## Construct loss functions and an indicator

## ## Get starting values using the rearrange and regress method
## rar_hyper_out <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), hyper_out))
## names(rar_hyper_out) <- letters[1:4]

## ## Randomly generate 1000 starting values arounf rar_hyper_out
## set.seed(1)
## starting <- get_feasibles(rar_hyper_out, is_monotone, 5000)

## ## Optimise with SMCSA with reciprocal schedule.
## rtnorm_rw <- rtnorm_rw_comp_fac(is_monotone, sd = 1, fraction = 0.97, n_pertub = 2)
## rm1_crss_out <- SMCSA(loss, rtnorm_rw, starting, recip_schedule, 5000, 400, TRUE, TRUE)
## plot_rat(rm1_crss_out$state[1:2], c(1, rm1_crss_out$state[3:4]))


## ## Optimise with SMCSA and 1000 iterations. Log schedule is preferable if
## ## running for an extended period.
## rtnorm_rw <- rtnorm_rw_fac(is_monotone, sd = 1, fraction = 0.97)
## rm2_crss_out <- SMCSA(loss, rtnorm_rw, starting, log_schedule, 1000, 500, TRUE, TRUE)
## plot_rat(rm2_crss_out$state[1:2], c(1, rm2_crss_out$state[3:4]))



## ######## Now we fit a unconstrained model with Tukey biweight and 'ht1'. ########
## ## Fit a rational model f(x) = (a + b*x) / (1 + c*x + d*x^2), subject to constraint
## ## Construct loss functions and an indicator
## loss <- rational_biweight_fac(ht1$y, ht1$x, 3L, 2L)

## ## Get starting values using the rearrange and regress method
## rar_ht1 <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), ht1))
## names(rar_ht1) <- letters[1:4]

## ## Randomly generate 1000 starting values arounf rar_ht1
## set.seed(1)
## starting <- matrix(rep(rar_ht1, times = 6000), 4, 6000) + rcauchy(4 * 6000, scale = 1)

## ## Optimise with SMCSA and 100 iterations
## rnorm_rw <- rnorm_rw_fac(sd = 2, fraction = 0.97)
## rm1_biw_out <- SMCSA(loss, rnorm_rw, starting, recip_schedule, 1000, 400, FALSE, TRUE)

## ## Optimise with SMCSA and 1000 iterations. Log schedule is preferable if
## ## running for an extended period.
## rnorm_rw <- rnorm_rw_fac(sd = 2, fraction = 0.99)
## rm2_biw_out <- SMCSA(loss, rnorm_rw, starting, log_schedule, 6000, 500, FALSE, TRUE)

## ## Plot the results
## plot_rat(rm1_biw_out$state[1:2], c(1, rm1_biw_out$state[3:4]))
## plot_rat(rm2_biw_out$state[1:2], c(1, rm2_biw_out$state[3:4]))




## ht1 <- as.data.frame(scale(HT1, FALSE, apply(HT1, 2, max)))
## plot(y ~ x, data = ht1, main = "HT1",
##      cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)


## ######## Now we fit a constrained model with Tukey biweight and 'ht1'. ########
## ## Fit a rational model f(x) = (a + b*x + c*x^2) / (1 + d*x + e*x^2), subject to constraint
## ## Construct loss functions, indicator, and proposal distribution
## loss <- rational_biweight_fac(ht1$y, ht1$x, 3L, 2L, c = 1)
## is_monotone <- is_monotones_fac(3L, 2L, min(ht1$x), max(ht1$x), TRUE)
## ## is_monotone <- function(x) {rep(TRUE, length = NCOL(x))}

## ## Get starting values using the rearrange and regress method
## rar_ht1 <- coef(lm(y ~ x + I(x^2) + I(-x*y) + I(-x^2*y), ht1))
## names(rar_ht1) <- letters[1:5]
## ## rar_ht1 <- coef(lm(y ~ x + I(-x*y) + I(-x^2*y), ht1))
## ## rar_ht1 <- rep(1, length = 4)
## ## names(rar_ht1) <- letters[1:4]

## ## Randomly generate 1000 starting values arounf rar_ht1
## set.seed(1)
## starting <- get_feasibles(rar_ht1, is_monotone, 1000)

## ## Optimise with SMCSA and 100 iterations
## rtnorm_rw <- rtnorm_rw_fac(is_monotone, sd = 1, fraction = 0.97)
## rm1_cbiw_out <- SMCSA(loss, rtnorm_rw, starting, recip_schedule, 3000, 500, FALSE, TRUE)

## ## Optimise with SMCSA and 1000 iterations. Log schedule is preferable if
## ## running for an extended period.
## rtnorm_rw <- rtnorm_rw_fac(is_monotone, sd = 1, fraction = 0.998)
## rm2_cbiw_out <- SMCSA(loss, rtnorm_rw, starting, log_schedule, 1000, 3000, FALSE, TRUE)

## ht1_pilot <- coef(nls(y ~ (a + b*x) / (1 + c*x + d*I(x^2)),
##                       ht1, start = rar_ht1))
## plot_rat(ht1_pilot[1:2], c(1, ht1_pilot[3:4]), col = 'red')

## ## Plot the results
## plot_rat(rm1_cbiw_out$state[1:3], c(1, rm1_cbiw_out$state[4:5]), col = 'blue')
## plot_rat(rm2_cbiw_out$state[1:3], c(1, rm2_cbiw_out$state[4:5]), col = 'blue')

