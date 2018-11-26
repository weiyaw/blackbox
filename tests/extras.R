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


#### LIDAR, quadratic spline ####
scaled_lidar <- as.data.frame(scale(lidar, FALSE, apply(abs(lidar), 2, max)))
plot(y ~ range, data = scaled_lidar, main = "Light Detection and Ranging (LIDAR)",
     cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)

## Construct loss functions
Bmat <- get_bs_design(scaled_lidar$range, 4, 2)$design
loss <- lm_rss_fac(scaled_lidar$y, Bmat)

## generate an arbitrary starting value
rand_eta <- NCOL(Bmat):1

## randomly generate 1000 starting values around an arbitrary starting value
set.seed(1)
starting <- matrix(rep(rand_eta, times = 1000), NCOL(Bmat), 1000) +
    rcauchy(1000 * NCOL(Bmat), scale = 1)

## SMCSA, log schedule
rnorm_rw <- rnorm_rw_fac(sd = 2, fraction = 0.99)
ld1_rss <- SMCSA(loss, rnorm_rw, starting, log_schedule, 6000, 1000, TRUE, TRUE)
lines((Bmat %*% ld1_rss$state) ~ scaled_lidar$range)

## exact solution
ld2_rss <- solve(crossprod(Bmat)) %*% crossprod(Bmat, scaled_lidar$y)
lines((Bmat %*% ld2_rss) ~ scaled_lidar$range, col = 'red')

## SAA and 1000000 iterations
ld3_rss <- SAA(loss, NULL, as.matrix(rand_eta), NULL, 1, 1000000, FALSE, TRUE)
lines((Bmat %*% ld3_rss$state) ~ scaled_lidar$range, col = 'blue')

#### LIDAR, monotonic decreasing quadratic spline ####
is_beta_decreasing <- function(betas) {
    res <- rep(NA, NCOL(betas))
    for (i in seq_along(res)) {
        res[i] <- all(diff(betas[, i]) <= 0)
    }
    res
}

rtrw_hmc_ld <- function(states, k) {
    if (k == 0) { sd <- 100 } else { sd <- 2 }
    fraction <- 0.95
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

## set sd = 100 and k = 0, and use trw_hmc_ld
set.seed(100)
starting <- get_feasibles(7:1, is_beta_decreasing, 1000)
## starting <- rtrw_comp(matrix(rep(7:1, times = 6000), ncol = 6000), 0)

## SMCSA, log schedule, fraction = 0.99, sd = 2
set.seed(101)
ld1_crss <- SMCSA(loss, rtrw_hmc_ld, starting, log_schedule, 6000, 500, TRUE, TRUE)
plot(y ~ range, data = scaled_lidar, main = "Light Detection and Ranging (LIDAR)",
     cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)
lines((Bmat %*% ld1_crss$state) ~ scaled_lidar$range, col = 'blue')

## SMCSA, reciprocal schedule, fraction = 0.95, sd = 2
set.seed(102)
ld2_crss <- SMCSA(loss, rtrw_hmc_ld, starting, recip_schedule, 6000, 300, TRUE, TRUE)
lines((Bmat %*% ld2_crss$state) ~ scaled_lidar$range)

## CEPSO, OLS as pilot estimate
set.seed(103)
ld_pilot <- solve(crossprod(Bmat)) %*% crossprod(Bmat, scaled_lidar$y)
ld3_crss <- CEPSO(starting, loss, ld_pilot, is_beta_decreasing, 500, 1000, 1000, TRUE)
lines(lines((Bmat %*% ld3_crss$theta) ~ scaled_lidar$range))

## SMCSA, log schedule, fraction = 0.99, sd = 2
set.seed(104)
rtrw_comp <- rtnorm_rw_comp_fac(is_beta_decreasing, 2, 0.995, 1)
ld4_crss <- SMCSA(loss, rtrw_comp, starting, log_schedule, 2000, 1000, TRUE, TRUE)
plot(y ~ range, data = scaled_lidar, main = "Light Detection and Ranging (LIDAR)",
     cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)
lines((Bmat %*% ld4_crss$state) ~ scaled_lidar$range, col = 'blue')

## SMCSA, recip schedule, fraction = 0.97, sd = 2
set.seed(105)
load_all()
rtrw_comp <- rtnorm_rw_comp_fac(is_beta_decreasing, 2, 0.97, 1)
ld5_crss <- SMCSA(loss, rtrw_comp, starting, recip_schedule, 2000, 500, TRUE, TRUE)
plot(y ~ range, data = scaled_lidar, main = "Light Detection and Ranging (LIDAR)",
     cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)
lines((Bmat %*% ld5_crss$state) ~ scaled_lidar$range, col = 'blue')





#### Dugong ####
load_all()
scaled_dugong <- as.data.frame(scale(dugong, FALSE, apply(abs(dugong), 2, max)))
plot(Length ~ Age, data = scaled_dugong, main = "Dugong",
     cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)

Bmat <- get_bs_design(scaled_dugong$Age, 2, 2)$design
loss <- lm_rss_fac(scaled_dugong$Length, Bmat)

rtnorm_rw_dg <- function(states, k) {
    sd <- 2
    fraction <- 0.99
    sigma <- (sd * fraction^k)^2 * diag(NROW(states))
    F <- diag(NROW(states)) * -1
    F[row(F) + 1 == col(F)] <- 1
    F <- F[1:(NROW(states) - 1), ]
    lower <- rep(0, len = NROW(F))
    initial <- seq(1, NROW(states))

    res <- matrix(NA, NROW(states), NCOL(states))
    for (i in 1:NCOL(states)) {
        res[, i] <- tnorm::rmvtnorm(1, states[, i], cov = sigma, g = lower, F = F,
                                    initial = initial, burn = 10)
    }
    res
}

## set sd = 100 and k = 0, and use rtnorm_rw_dg
set.seed(100)
starting <- rtnorm_rw_dg(matrix(rep(1:NCOL(Bmat), times = 1000), ncol = 1000), 1)

## SMCSA, log schedule
dg1_crss <- SMCSA(loss, rtnorm_rw_dg, starting, log_schedule, 2000, 1000, TRUE, TRUE)

## SMCSA, log schedule
dg2_crss <- SMCSA(loss, rtnorm_rw_dg, starting, recip_schedule, 2000, 1000, TRUE, TRUE)



dg_pilot <- solve(crossprod(Bmat)) %*% crossprod(Bmat, scaled_dugong$Length)

Bmat_plot <- with(scaled_dugong, get_bs_design(seq(min(Age), max(Age), len = 100), 2, 2)$design)
lines((Bmat_plot %*% dg1_crss$state) ~ with(scaled_dugong, seq(min(Age), max(Age), len = 100)), col = 'black')
lines((Bmat_plot %*% dg_pilot) ~ with(scaled_dugong, seq(min(Age), max(Age), len = 100)), col = 'red')




#### Rabbit ####
load_all()
scaled_rabbit <- as.data.frame(scale(rabbit, FALSE, apply(abs(rabbit), 2, max)))
plot(Lens ~ Days, data = scaled_rabbit, main = "Rabbit",
     cex = 1.3, cex.lab = 1.3, cex.axis = 1.3)
Bmat <- get_bs_design(scaled_rabbit$Days, 2, 2)$design
loss <- lm_rss_fac(scaled_rabbit$Lens, Bmat)

pilot <- solve(crossprod(Bmat)) %*% crossprod(Bmat, scaled_rabbit$Lens)
lines((Bmat %*% pilot) ~ scaled_rabbit$Days, col = 'blue')

## ## get starting values using the rearrange and regress method
## rar_rabbit <- coef(lm(Lens ~ Days + I(-Lens*Days), scaled_rabbit))
## names(rar_rabbit) <- letters[1:3]

## randomly generate 1000 starting values
set.seed(2)
starting <- matrix(rep(rar_rabbit, times = 1000), length(rar_rabbit), 1000) +
    rcauchy(1000 * length(rar_rabbit), scale = 1)

## SMCSA, log schedule
rnorm_rw <- rnorm_rw_fac(sd = 2, fraction = 0.99)
rb1_rss <- SMCSA(loss, rnorm_rw, starting, log_schedule, 1000, 500, FALSE, TRUE)
plot_rat(rb1_rss$state[1:2], c(1, rb1_rss$state[3]))

## Optimise with Newton-Gauss
rb2_rss <- nls(Lens ~ (a + b*Days) / (1 + c*Days), scaled_rabbit, start = rar_rabbit)
plot_rat(coef(rb2_rss)[1:2], c(1, coef(rb2_rss)[3]), col = 'red')

## Optimise with SAA and 5000000 iterations
set.seed(1)
rm3_rss <- SAA(loss, NULL, as.matrix(rar_hyper), NULL, 1, 5000000, FALSE, TRUE)





