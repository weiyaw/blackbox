## get_deriv_nums
rm(list = ls())
source("./rational-function.R")
## one (2, 2) rational function
fs <- matrix(1:3, 3, 1)
gs <- matrix(4:6, 3, 1)
all(get_deriv_nums(fs, gs) == as.matrix(c(3, 12, 3)))

## one horizontal line and two (2, 2) rational functions
fs <- matrix(1:9, 3, 3)
gs <- matrix(c(1:3, 1:6), 3, 3)
all(get_deriv_nums(fs, gs) == cbind(c(0, 0, 0),
                                    c(-3, -12, -3),
                                    c(-3, -12, -3)))

fs <- matrix(0:2, 1, 3)
gs <- matrix(c(1:3, 1:6), 3, 3)
all(get_deriv_nums(fs, gs) == cbind(c(0, 0),
                                    c(-2, -6),
                                    c(-10, -24)))

fs <- matrix(c(0, 1:3, 0, 5:7, 0), 3, 3)
gs <- matrix(c(1, 2, 0, 1, 0, 3, 0, 5, 6), 3, 3)
all(get_deriv_nums(fs, gs) == cbind(c(1, 4, 4),
                                    c(0, -8, 0),
                                    c(-30, -72, -42)))


## eval_pols
rm(list = ls())
source("./rational-function.R")
## one quadratic at 4 points
poly_coefs <- matrix(1:3, 3, 1)
all(eval_pols(poly_coefs, 0:4) == as.matrix(c(1, 6, 17, 34, 57)))

## three quadratic at a single point
poly_coefs <- matrix(1:9, 3, 3)
all(eval_pols(poly_coefs, 4) == matrix(c(57, 120, 183), 1, 3))

## three quadratic at 4 points
poly_coefs <- matrix(c(1, 2, 0, 1, 0, 3, 0, 5, 6), 3, 3)
all(eval_pols(poly_coefs, 0:4) == cbind(c(1, 3, 5, 7, 9),
                                        c(1, 4, 13, 28, 49),
                                        c(0, 11, 34, 69, 116)))


## is_positive
rm(list = ls())
source("./rational-function.R")
## horizontal line
poly_coefs <- matrix(-1:1, 1, 3)
all(is_positives(poly_coefs, 1, 2, TRUE) == c(FALSE, TRUE, TRUE))

## straight line (same slope, different intercepts)
poly_coefs <- matrix(c(-3, 1, -2, 1, -1, 1, 0, 1), 2, 4)
all(is_positives(poly_coefs, 1, 2, TRUE) == c(FALSE, FALSE, TRUE, TRUE))

## quadratic polynomial (roots at x = 1 and x = 2)
poly_coefs <- matrix(c(-2, 3, -1, 2, -3, 1), 3, 2)
all(is_positives(poly_coefs, 1, 2, TRUE) == c(TRUE, FALSE))

poly_coefs <- matrix(c(-2, 3, -1, 2, -3, 1), 3, 2)
all(is_positives(poly_coefs, 1, Inf, FALSE) == c(FALSE, FALSE))

poly_coefs <- matrix(c(-2, 3, -1, 2, -3, 1), 3, 2)
all(is_positives(poly_coefs, 2, Inf, FALSE) == c(TRUE, FALSE))

## quadratic polynomial (roots at x = -1 and x = -2)
poly_coefs <- matrix(c(2, 3, 1, -2, -3, -1), 3, 2)
all(is_positives(poly_coefs, -Inf, Inf, TRUE) == c(FALSE, FALSE))

## quadratic/cubic polynomial (roots at x = -1 and x = 2)
poly_coefs <- matrix(c(1, 2, 1, 0, 4, 0, -3, 1), 4, 2)
all(is_positives(poly_coefs, -1, Inf, TRUE) == c(TRUE, TRUE))

## strange polynomials
poly_coefs <- matrix(c(1, 2, 0, 1, 0, 3, 0, 5, 6), 3, 3)
all(is_positives(poly_coefs, -1, Inf, TRUE) == c(FALSE, TRUE, FALSE))


## is_no_roots_fac
rm(list = ls())
source("./rational-function.R")
## horizontal line
poly_coefs <- matrix(-1:1, 1, 3)
all(is_no_roots(poly_coefs, 1, 2) == c(TRUE, TRUE, TRUE))

## straight line (same slope, different intercepts)
poly_coefs <- matrix(c(-3, 1, -2, 1, -1, 1, 0, 1), 2, 4)
all(is_no_roots(poly_coefs, 1, 2) == c(TRUE, FALSE, FALSE, TRUE))

## quadratic polynomial (roots at x = 1 and x = 2)
poly_coefs <- matrix(c(-2, 3, -1, 2, -3, 1), 3, 2)
all(is_no_roots(poly_coefs, 1, 2) == c(FALSE, FALSE))

## quadratic/quartic polynomial (no roots)
poly_coefs <- matrix(c(2, 2, 2, 0, 0, 2, 3, 1, 3, 10), 5, 2)
all(is_no_roots(poly_coefs, -Inf, Inf) == c(TRUE, TRUE))

## quadratic polynomial (roots at x = -1 and x = -2)
poly_coefs <- matrix(c(2, 3, 1, -2, -3, -1), 3, 2)
all(is_no_roots(poly_coefs, -Inf, Inf) == c(FALSE, FALSE))

## quadratic/cubic polynomial (roots at x = -1 and x = 2)
poly_coefs <- matrix(c(1, 2, 1, 0, 4, 0, -3, 1), 4, 2)
all(is_no_roots(poly_coefs, 0, Inf) == c(TRUE, FALSE))

## strange polynomials
poly_coefs <- matrix(c(1, 2, 0, 1, 0, 3, 0, 5, 6), 3, 3)
all(is_no_roots(poly_coefs, -1, Inf) == c(FALSE, TRUE, FALSE))


## is_monotones_fac
rm(list = ls())
source("./rational-function.R")
## one degree(2, 2) rational function
fs <- matrix(c(1, 2, 1), 3, 1)
gs <- matrix(2:3, 2, 1)
all(is_monotones_fac(3, 2, 0, Inf, FALSE)(rbind(fs, gs)) == TRUE)
all(is_monotones_fac(3, 2, 0, Inf, TRUE)(rbind(fs, gs)) == FALSE)
all(is_monotones_fac(3, 2, -1, 1, TRUE)(rbind(fs, gs)) == FALSE)
all(is_monotones_fac(3, 2, -Inf, -1, FALSE)(rbind(fs, gs)) == TRUE)

## two degree(1, 1), one degree(0, 1) rational functions
fs <- matrix(c(-1, 1, 0, 1, 1, 0), 2, 3)
gs <- matrix(1:3, 1, 3)
all(is_monotones_fac(2, 1, -Inf, -1, TRUE)(rbind(fs, gs)) == c(FALSE, TRUE, FALSE))
all(is_monotones_fac(2, 1, -Inf, -1, FALSE)(rbind(fs, gs)) == c(FALSE, FALSE, TRUE))
all(is_monotones_fac(2, 1, 0, Inf, TRUE)(rbind(fs, gs)) == c(TRUE, TRUE, FALSE))
all(is_monotones_fac(2, 1, -0.5, 1, TRUE)(rbind(fs, gs)) == c(TRUE, FALSE, FALSE))


## strange rational functions
fs <- matrix(c(0, 1:3, 0, 5:7, 0), 3, 3)
gs <- matrix(c(1, 2, 0, 1, 0, 1, 0, 5, 6), 3, 3)
all(is_monotones_fac(3, 3, -Inf, -1, FALSE)(rbind(fs, gs)) == c(TRUE, TRUE, FALSE))
all(is_monotones_fac(3, 3, -Inf, -1, TRUE)(rbind(fs, gs)) == c(FALSE, FALSE, FALSE))
all(is_monotones_fac(3, 3, -0.5, 0, TRUE)(rbind(fs, gs)) == c(FALSE, FALSE, TRUE))
all(is_monotones_fac(3, 3, -0.5, 0, FALSE)(rbind(fs, gs)) == c(FALSE, TRUE, FALSE))
all(is_monotones_fac(3, 3, 0, Inf, TRUE)(rbind(fs, gs)) == c(TRUE, FALSE, FALSE))
all(is_monotones_fac(3, 3, 2, Inf, FALSE)(rbind(fs, gs)) == c(FALSE, TRUE, TRUE))


## simulate datasets for loss functions calculations
rm(list = ls())
source("./rational-function.R")
set.seed(240318)
x1 <- 1:10
y1 <- (1 + 2 * x1) / (1 + 5 * x1) + rnorm(x1, sd = 0.1)
res11 <- y1 - (1 + 2 * x1) / (1 + 5 * x1)
res12 <- y1 - (5 * x1) / (1 + 7 * x1 + 9 * x1^2)
res13 <- y1 - 5 / (1 + 7 * x1 + 9 * x1^2)

x2 <- 1:10
y2 <- (1 + 2 * x2) / (1 + 5 * x2^2) + rnorm(x2, sd = 5)
res21 <- y2 - (1 + 2 * x2) / (1 + 5 * x2)
res22 <- y2 - (5 * x2) / (1 + 7 * x2 + 9 * x2^2)
res23 <- y2 - 5 / (1 + 7 * x2 + 9 * x2^2)

## rational_rss_fac
rss11 <- crossprod(res11)
rss12 <- crossprod(res12)
rss13 <- crossprod(res13)
rss21 <- crossprod(res21)
rss22 <- crossprod(res22)
rss23 <- crossprod(res23)
rational_rss_fac(y1, x1, 2, 1)(matrix(c(1, 2, 5), 3, 1)) - rss11 < 1E-6
all(rational_rss_fac(y1, x1, 2, 2)
(cbind(c(0, 5, 7, 9), c(5, 0, 7, 9))) - c(rss12, rss13) < 1E-6)

## rational_biweight_fac
cst <- 4.685
rho <- function(x, cst = 4.685){
    out <- abs(x) > cst
    out[!out] <- (cst^2)/6 * (1-(1-(x[!out]/cst)^2)^3)
    out
}

biw11 <- sum(rho(res11))
biw12 <- sum(rho(res12))
biw13 <- sum(rho(res13))
biw21 <- sum(rho(res21))
biw22 <- sum(rho(res22))
biw23 <- sum(rho(res23))
rational_biweight_fac(y1, x1, 2, 1)(matrix(c(1, 2, 5), 3, 1)) - biw11 < 1E-6
all(rational_biweight_fac(y1, x1, 2, 2)
(cbind(c(0, 5, 7, 9), c(5, 0, 7, 9))) - c(biw12, biw13) < 1E-6)


## plot_rat
rm(list = ls())
source("./rational-function.R")
## individual rational function
fs <- as.matrix(c(0, 1, 2))
gs <- as.matrix(c(1, 1, 2, 0))
plot_rat(fs, gs, c(-5, 5), FALSE)

fs <- as.matrix(c(3, 0, 5))
gs <- as.matrix(c(1, 1, 0, 1))
plot_rat(fs, gs, c(-5, 5), FALSE, ylim = c(-5, 5))

fs <- as.matrix(c(6, 7, 0))
gs <- as.matrix(c(1, 0, 5, 6))
plot_rat(fs, gs, c(-5, 5), FALSE, ylim = c(-5, 5))

## all in one plot
fs <- matrix(c(0, 1:3, 0, 5:7, 0), 3, 3)
gs <- matrix(c(1, 1, 2, 0, 1, 1, 0, 1, 1, 0, 5, 6), 4, 3)
plot_rat(fs, gs, c(-5, 5), FALSE, ylim = c(-5, 5))
