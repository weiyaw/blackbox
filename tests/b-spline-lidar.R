## Proposal distribution and number of iterations depend on the cooling schedule and problem at hand.

######## fit a constrained model. ########
## Visualise the dataset
load_all()
rm(list = ls())
## K: number of inner knots
get_bs_design <- function(x, K, deg, EPS = 1e-6) {
    dist <- diff(range(x)) / (K + 1)
    knots <- seq(min(x) - (deg * dist) - EPS, max(x) + (deg * dist) + EPS,
                 len = K + 2 * (deg + 1))
    res <- list()
    res$design <- splines::splineDesign(knots, x, ord = deg + 1)
    res$knots <- knots
    return(res)
}

## RSS for linear model
lm_rss_fac <- function(y, design_mat) {
    if (NROW(y) != NROW(design_mat)) {
        stop("Number of rows for y and design matrix mismatch.")
    }
    function(beta) {
        beta <- as.matrix(beta)
        colSums((y - design_mat %*% beta)^2)
    }
}

#### LIDAR, monotonic decreasing quadratic spline ####
is_beta_decreasing <- function(betas) {
    res <- rep(NA, NCOL(betas))
    for (i in seq_along(res)) {
        res[i] <- all(diff(betas[, i]) <= 0)
    }
    res
}

rtrw_hmc_97 <- function(states, k) {
    if (k == 0) { sd <- 100 } else { sd <- 1 }
    fraction <- 0.97
    sigma <- (sd * fraction^k)^2 * diag(NROW(states))
    F <- diag(NROW(states))
    F[row(F) + 1 == col(F)] <- -1
    F <- F[1:(NROW(states) - 1), ]
    lower <- rep(0, len = NROW(F))
    initial <- seq(NROW(states), 1)

    res <- matrix(NA, NROW(states), NCOL(states))
    for (i in 1:NCOL(states)) {
        res[, i] <- tnorm::rmvtnorm(1, states[, i], cov = sigma, g = lower, F = F,
                                    initial = initial, burn = 10)
    }
    res
}

rtrw_hmc_998 <- function(states, k) {
    if (k == 0) { sd <- 100 } else { sd <- 1 }
    fraction <- 0.97
    sigma <- (sd * fraction^k)^2 * diag(NROW(states))
    F <- diag(NROW(states))
    F[row(F) + 1 == col(F)] <- -1
    F <- F[1:(NROW(states) - 1), ]
    lower <- rep(0, len = NROW(F))
    initial <- seq(NROW(states), 1)

    res <- matrix(NA, NROW(states), NCOL(states))
    for (i in 1:NCOL(states)) {
        res[, i] <- tnorm::rmvtnorm(1, states[, i], cov = sigma, g = lower, F = F,
                                    initial = initial, burn = 10)
    }
    res
}


#### LIDAR, quadratic spline ####
scaled_lidar <- as.data.frame(scale(lidar, FALSE, apply(abs(lidar), 2, max)))
plot(y ~ range, data = scaled_lidar, main = "Light Detection and Ranging (LIDAR)",
     cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)

## Construct loss functions
Bmat <- get_bs_design(scaled_lidar$range, 4, 2)$design
loss <- lm_rss_fac(scaled_lidar$y, Bmat)

## generate an arbitrary starting value
rand_eta <- NCOL(Bmat):1


library(doParallel)
registerDoParallel(cores = 8)

recip_time <- system.time(recip_lidar <- foreach(i = 1:40) %dopar% {
    ## Randomly generate 1000 starting values arounf rar_hyper
    set.seed(500 + i)
    starting <- get_feasibles(rand_eta, is_beta_decreasing, 1000, dist_para = list(scale = 2))

    ## SMCSA, reciprocal, fraction = 0.97, sd = 1, 500 iter, N = 3000, 2 comp
    ## Use a more aggresive proposal
    rtnorm_rw <- rtnorm_rw_comp_fac(is_beta_decreasing, sd = 1, fraction = 0.97, n_pertub = 2)
    runtime <- system.time(bs1 <- SMCSA(loss, rtnorm_rw, starting, recip_schedule, 3000, 1000, FALSE, TRUE))
    bs1$runtime <- runtime
    bs1
})
cat("RECIP DONE.", recip_time[3], "\n")
# plot(y ~ range, data = scaled_lidar, main = "Light Detection and Ranging (LIDAR)",
#      cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)
# lines((Bmat %*% bs1$state) ~ scaled_lidar$range, col = 'red')
saveRDS(recip_lidar, "~/Dropbox/honours/extras/lidar/recip_lidar.rds")

recip85_time <- system.time(recip85_lidar <- foreach(i = 1:40) %dopar% {
  ## Randomly generate 1000 starting values arounf rar_hyper
  set.seed(500 + i)
  starting <- get_feasibles(rand_eta, is_beta_decreasing, 1000, dist_para = list(scale = 2))
  
  ## SMCSA, reciprocal, fraction = 0.97, sd = 1, 500 iter, N = 3000, 2 comp
  ## Use a more aggresive proposal
  rtnorm_rw <- rtnorm_rw_comp_fac(is_beta_decreasing, sd = 1, fraction = 0.97, n_pertub = 2)
  runtime <- system.time(bs1 <- SMCSA(loss, rtnorm_rw, starting, recip_schedule_0.85, 3000, 1000, FALSE, TRUE))
  bs1$runtime <- runtime
  bs1
})
cat("RECIP DONE.", recip85_time[3], "\n")
# plot(y ~ range, data = scaled_lidar, main = "Light Detection and Ranging (LIDAR)",
#      cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)
# lines((Bmat %*% bs1$state) ~ scaled_lidar$range, col = 'red')
saveRDS(recip85_lidar, "~/Dropbox/honours/extras/lidar/recip85_lidar.rds")


log_time <- system.time(log_lidar <- foreach(i = 1:40) %dopar% {
  ## Randomly generate 1000 starting values arounf rar_hyper
  set.seed(500 + i)
  starting <- get_feasibles(rand_eta, is_beta_decreasing, 1000, dist_para = list(scale = 2))

  ## SMCSA, log, fraction = 0.998, sd = 1, 3000 iter, N = 1000, 2 comp
  rtnorm_rw <- rtnorm_rw_comp_fac(is_beta_decreasing, sd = 1, fraction = 0.998, n_pertub = 2)
  runtime <- system.time(bs2 <- SMCSA(loss, rtnorm_rw, starting, log_schedule, 1000, 3000, FALSE, TRUE))
  bs2$runtime <- runtime
  bs2
})
cat("LOG DONE.", log_time[3], "\n")
# plot(y ~ range, data = scaled_lidar, main = "Light Detection and Ranging (LIDAR)",
#      cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)
# lines((Bmat %*% bs1$state) ~ scaled_lidar$range, col = 'red')
saveRDS(log_lidar, "~/Dropbox/honours/extras/lidar/log_lidar.rds")
## plot_rat(bs2$state[1:3], c(1, bs2$state[4:5]), col = 'blue')


lidar_pilot <- solve(crossprod(Bmat)) %*% crossprod(Bmat, scaled_lidar$y)

cepso_time <- system.time(cepso_lidar <- foreach(i = 1:40) %dopar% {
  ## Randomly generate 1000 starting values arounf rar_hyper
  set.seed(500 + i)
  starting <- get_feasibles(rand_eta, is_beta_decreasing, 1000, dist_para = list(scale = 2))

  ## CEPSO, 2000 iter, N = 3000, neighbour = 3000/5
  ## OLS as pilot estimate
  runtime <- system.time(bs3 <- CEPSO(starting, loss, lidar_pilot, is_beta_decreasing, iter = 2000, N = 3000, verbose = TRUE))
  bs3$runtime <- runtime
  bs3
})
cat("CEPSO DONE.", cepso_time[3], "\n")
## plot(y ~ range, data = scaled_lidar, main = "Light Detection and Ranging (LIDAR)",
##      cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)
## lines((Bmat %*% lidar_pilot) ~ scaled_lidar$range, col = 'blue')
saveRDS(cepso_lidar, "~/Dropbox/honours/extras/lidar/cepso_lidar.rds")


















