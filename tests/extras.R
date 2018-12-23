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





