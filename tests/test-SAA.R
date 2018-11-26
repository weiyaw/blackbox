## a SAA example
ux <- function(x) {
  x <- as.matrix(x)
  x1 <- x[1, ]
  x2 <- x[2, ]
  
  res <- -1 * (x1 * sin(20 * x2) + x2 * sin(20 * x1))^2 * cosh(sin(10 * x1) * x1) -
    (x1 * cos(10 * x2) - x2 * sin(10 * x1))^2
  
  res[abs(x1) > 1.1 | abs(x2) > 1.1] <- Inf
  res
}

load_all()
set.seed(1)
## minimum at
## partitions: u_1 = -8, u_40 = 0.2
## schedule: t_h = 0.5, T_0 = 200, t_* = 0.01
## proposal: random walk N(. , 0.25^2 I )
## sampling dist: zeta = 0.1
## gain: sigma = 1.0, T_0 = 2000
ux_1 <- SAA(ux, NULL, cbind(c(1, 1), c(1, 1)), NULL, 2, 10000, FALSE, F)

rnorm_rw <- rnorm_rw_fac(sd = 2, fraction = 0.99)
starting <- as.matrix(c(1, 1))
ux_2 <- SMCSA(ux, rnorm_rw, starting, log_schedule, 1000, 500, FALSE, F)