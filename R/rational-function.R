## ALL CODES FOR FITTING LEAST SQUARES AND TUKEY BIWEIGHT RATIONAL FUNCTIONS

#' Get the numerators of the derivatives of a set of rational functions.
#'
#' Assuming that the rational function is r(x) = f(x) / g(x), this function
#' returns f'(x)g(x) - f(x)g'(x).
#'
#' @param fs a matrix of numerator polynomial in ascending order representing,
#'     each column corresponds to a polynomial. (num mat)
#' @param gs a matrix of denominator coefficients in ascending order, each
#'     column corresponds to a polynomial. (num mat)
#' @param EPS error tolerance. (num)
#'
#' @return a matrix of polynomial coefficients.
get_deriv_nums <- function(fs, gs, EPS = 1e-6) {

    lg <- NROW(gs)
    lf <- NROW(fs)

    if (lg == 1) {
        ## if the denominator is a constant
        stop('This is a polynomial')

    } else if (lf == 1) {
        ## if the numerator is a constant
        ## derivatives of denominator (no reversing)
        gs_deriv <- gs[seq.int(2, lg), , drop = FALSE] * seq_len(lg - 1)
        ## result
        res <- matrix(NA, nrow = lg - 1, ncol = NCOL(fs))
        for (i in seq_len(NROW(gs_deriv))) {
            res[i, ] <- -fs * gs_deriv[i, ]
        }

    } else {
        ## if the numerator is a polynomial
        ## derivatives of denominator (reversed)
        gs_deriv_r <- gs[seq.int(lg, 2), , drop = FALSE] *
            (lg - seq_len(lg - 1))
        ## derivatives of numerator (reversed)
        fs_deriv_r <- fs[seq.int(lf, 2), , drop = FALSE] *
            (lf - seq_len(lf - 1))
        ## result
        res <- matrix(NA, nrow = lg + lf - 2, ncol = NCOL(fs))
        for (i in seq_len(NCOL(fs))) {
            res[, i] <- convolve(gs[, i], fs_deriv_r[, i], type = 'o') -
                convolve(fs[, i], gs_deriv_r[, i], type = 'o')
        }
    }

    ## FAST IMPLEMENTATION, MIGHT NOT WORK IN CERTAIN CASES.
    ## This operation only remove the zero when differentiating rational
    ## function with same degrees of numerator and denominator.
    if (lg == lf) {
        return(res[-NROW(res), , drop = FALSE])
    } else {
        return(res)
    }
}

#' Evaluate polynomials.
#'
#' The polynomials are evaluated at given points using the Horner scheme.
#'
#' @param poly_coefs a matrix of polynomial coefficients in ascending order,
#'     each column corresponds to a polynomial. (num mat)
#' @param x a vector of position where the polynomials are evaluated at. (num vec)
#'
#' @return a matrix, each row and column corresponds to a point in \code{x} and
#'     a polynomial respectively.
eval_pols <- function(poly_coefs, x) {

    ncoefs <- NCOL(poly_coefs)
    nx <- length(x)
    res <- matrix(0, ncoefs, nx)
    X <- matrix(rep(x, each = ncoefs), ncoefs, nx)
    for (bi in seq.int(NROW(poly_coefs), 1)) {
        res <- X * res + poly_coefs[bi, ]
    }
    t(res)
}

#' Are these polynomials positive?
#'
#' Check if a set of polynomial is positive (or negative) inside a given
#' interval.
#'
#' This can be done by checking wherether the roots of the polynomials have even
#' multiplicities. In addition, one of the points in the interval must have the
#' desired sign.
#'
#' @param poly_coefs a matrix of polynomial coefficients in ascending order,
#'     each column corresponds to a polynomial. (num mat)
#' @param a the lower bound. (num)
#' @param b the upper bound. (num)
#' @param positive if it is \code{FALSE}, check if the polynomial is
#'     negative. (bool)
#' @param EPS error tolerance. (num)
#'
#' @return a boolean vector indicating whether the condition in the description
#'     is met.
is_positives <- function(poly_coefs, a, b, positive = TRUE, EPS = 1e-06) {
    if (a >= b) {
        stop("a must be strictly larger than b.")
    }

    ## Choose five points between the limits
    if (is.infinite(a) && is.infinite(b)) {
        chkpoints <- 1:5
    } else if (is.infinite(a)) {
        chkpoints <- b - 1:5
    } else if (is.infinite(b)) {
        chkpoints <- a + 1:5
    } else {
        chkpoints <- seq.int(a, b, length.out = 5)
    }

    ## Check the signs of these points
    if (positive) {
        signs_mat <- eval_pols(poly_coefs, chkpoints) >= -EPS
    } else {
        signs_mat <- eval_pols(poly_coefs, chkpoints) <= EPS
    }

    ## All points should have the same sign
    satisfy <- colSums(signs_mat) == 5

    ## Check the multiplicities of the roots
    for (i in which(satisfy)) {
        roots <- polyroot(poly_coefs[, i])
        re.roots <- Re(roots)[abs(Im(roots)) < EPS]
        re.roots <- re.roots[a + EPS <= re.roots & re.roots <= b - EPS]
        if (length(re.roots) == 0) {
            satisfy[i] <- TRUE
        } else {
            mltplcty <- rowSums(outer(re.roots, re.roots,
                                      function(x, y) abs(x - y) < EPS))
            satisfy[i] <- all(mltplcty%%2 == 0L)
        }
    }
    satisfy
}

#' Is there any roots in these polynomials?
#'
#' For a given set of polynomials, find all roots in a given interval and return
#' \code{TRUE} if there are no roots in it.
#'
#' @param poly_coefs a matrix of polynomial coefficients in ascending order,
#'     each column corresponds to a polynomial. (num mat)
#' @param a the lower bound. (num)
#' @param b the upper bound. (num)
#' @param EPS error tolerance. (num)
#'
#' @return a boolean vector indicating whether the condition in the description
#'     is met.
is_no_roots <- function(poly_coefs, a, b, EPS = 1e-06) {
    satisfy <- rep(NA, NCOL(poly_coefs))
    for (i in seq_len(NCOL(poly_coefs))) {
        roots <- polyroot(poly_coefs[, i])
        re.roots <- Re(roots)[abs(Im(roots)) < EPS]
        re.roots <- re.roots[a - EPS <= re.roots & re.roots <= b + EPS]
        satisfy[i] <- length(re.roots) == 0
    }
    satisfy
}

#' Is there any roots in these polynomials? (factory)
#'
#' Produce a function that checks whether a given set of rational functions is
#' monotone in a given interval.
#'
#' @param dim_f the number of coefficients in the numerator. (num)
#' @param dim_g the number of coefficients in the denominator, excluding the
#'     constant term. (num)
#' @param xmin the lower bound. (num)
#' @param xmax the upper bound. (num)
#' @param increasing set to \code{TRUE} when checking for monotonic increasing,
#'     \code{FALSE} for monotonic decreasing.
#' @param EPS error tolerance. (num)
#'
#' @return a function that takes in a matrix of rational function coefficients
#'     (numerator precede denominator, in ascending order), each column
#'     corresponds to a rational function, and returns a boolean vector
#'     indicating whether the condition in the description is met.
is_monotones_fac <- function(dim_f, dim_g, xmin = 0, xmax = Inf, increasing = TRUE,
                            EPS = 1e-06) {
    force(xmin)
    force(xmax)
    force(increasing)
    force(EPS)
    idx_f <- seq.int(1, length.out = dim_f)
    idx_g <- seq.int(dim_f + 1, length.out = dim_g)
    n_dims <- dim_f + dim_g

    function(rat_coefs) {
        if (NROW(rat_coefs) != n_dims) {
            stop("Input's dimension mismatch with specified dimensions.")
        }
        ## Get the numerators and denominators of the rational functions
        rat_f <- rat_coefs[idx_f, , drop = FALSE]
        rat_g <- rbind(rep(1, NCOL(rat_coefs)),
                       rat_coefs[idx_g, , drop = FALSE])

        ## Continuity condition
        satisfy <- is_no_roots(rat_g, xmin, xmax, EPS)

        if (any(satisfy)) {
            ## Monotone condition
            derivs <- get_deriv_nums(rat_f[, satisfy, drop = FALSE],
                                     rat_g[, satisfy, drop = FALSE],
                                     EPS)
            satisfy[satisfy] <- is_positives(derivs, xmin, xmax, increasing,
                                             EPS)
        }
        satisfy
    }
}

#' Residual sum of squares of rational function models
#'
#' Produce a function that calculates the residual sum of squares of a set of
#' rational function models.
#'
#' The produced function can be minimised using SMC-SA, multiSA or CEPSO.
#'
#' @param y a response vector. (num vec)
#' @param x a predictor vector. (num vec)
#' @param dim_f the number of coefficients in the numerator. (num)
#' @param dim_g the number of coefficients in the denominator, excluding the
#'     constant term. (num)
#'
#' @return a function that takes in a matrix of rational function coefficients
#'     (numerator precede denominator, in ascending order), each column
#'     corresponds to a rational function, and returns a vector of redisual sum
#'     of squares.
rational_rss_fac <- function(y, x, dim_f, dim_g){

    ## Evaluate function's arguments eagerly
    force(x)
    force(y)
    force(dim_f)
    force(dim_g)

    ## Type-checking
    if (!(is.vector(x) & is.vector(y))) {
        stop("Non-vector dataset")
    }

    ## Sequence of indices for both numerator and denominator
    idx_f <- seq.int(1, dim_f)
    idx_g <- seq.int(dim_f + 1, dim_f + dim_g)
    total_dim <- dim_f + dim_g

    ## Check the dimension of the response
    if (is.matrix(y) && (NCOL(y) == 1 || NROW(y) == 1)) {
        y <- as.vector(y)
    } else if (!is.vector(y)) {
        stop("Response data must be a vector / 1-dim matrix.")
    }

    ## A function that returns the residuals for each column of "theta" that
    ## must be supplied as a matrix
    function(theta){
        n_row <- NROW(theta)
        n_col <- NCOL(theta)

        if (is.vector(theta)) {
            warning("Argument of residual function supplied as vector.")
            theta <- as.matrix(theta)
        }

        if (n_row != total_dim) {
            warning("Argument of residual funciton mismatch with dimension supplied.")
        }

        fit <- eval_pols(theta[idx_f, , drop=FALSE], x) /
            eval_pols(rbind(rep(1, len=n_col), theta[idx_g, , drop=FALSE]), x)
        return(colSums((y - fit)^2))
       }
}

#' Tukey's biweight of rational function models
#'
#' Produce a function that calculates the Tukey's biweight of a set of rational
#' function models.
#'
#' The explicit formula used is rho(x) = (c^2)/6 * (1-(1-(x/c)^2)^3) if x < |c|,
#' or rho(x) = 1 if x > |c|, where c is an arbitrary constant. The sum of
#' rho(y_i - r(x_i)), where r(.) is the rational function, is then
#' minimised. The produced function can be minimised using SMC-SA, multiSA or
#' CEPSO.
#'
#' @param y a response vector. (num vec)
#' @param x a predictor vector. (num vec)
#' @param dim_f the number of coefficients in the numerator. (num)
#' @param dim_g the number of coefficients in the denominator, excluding the
#'     constant term. (num)
#' @param cst the constant in Tukey's biweight. (num)
#'
#' @return a function that takes in a matrix of rational function coefficients
#'     (numerator precede denominator, in ascending order), each column
#'     corresponds to a rational function, and returns a vector Tukey's
#'     biweight.
rational_biweight_fac <- function(y, x, dim_f, dim_g, cst = 4.685){
    force(x)
    force(y)
    force(dim_f)
    force(dim_g)
    force(cst)

    ## Sequence of indices for both numerator and denominator
    idx_f <- seq.int(1, dim_f)
    idx_g <- seq.int(dim_f + 1, dim_f + dim_g)
    total_dim <- dim_f + dim_g

    ## Check the dimension of the response
    if (is.matrix(y) && (NCOL(y) == 1 || NROW(y) == 1)) {
        y <- as.vector(y)
    } else if (!is.vector(y)) {
        stop("Response data must be a vector / 1-dim matrix.")
    }

    ## The biweight function (integration of phi)
    ## x is a matrix
    rho <- function(x, cst = 4.685){
        res <- matrix((cst^2)/6, NROW(x), NCOL(x))
        outside <- abs(x) > cst
        res[!outside] <- (cst^2)/6 * (1-(1-(x[!outside]/cst)^2)^3)
        return(res)
    }

    ## A function that return biweight for each column
    ## "theta" must be supplied as matrix
    function(theta) {
        n_row <- NROW(theta)
        n_col <- NCOL(theta)

        if (is.vector(theta)) {
            warning("Argument of residual function supplied as vector.")
            theta <- as.matrix(theta)
        }

        if (n_row != total_dim) {
            stop("Argument of objective mismatch with the dimension supplied.")
        }

        fit <- eval_pols(theta[idx_f, , drop=FALSE], x) /
            eval_pols(rbind(rep(1, len=n_col), theta[idx_g, , drop=FALSE]), x)
        return(colSums(rho(y - fit)))
    }

}



#' Line plots of rational functions
#'
#' Plot one or multiple rational functions, given a set of coefficients.

#' @param fs a matrix of numerator polynomial in ascending order representing,
#'     each column corresponds to a polynomial. (num mat)
#' @param gs a matrix of denominator coefficients in ascending order, each
#'     column corresponds to a polynomial. (num mat)
#' @param xlim the x-axis range to plot
#' @param sketch if set to \code{TRUE}, this function add a rational function
#'     onto an existing plot, rather than creating a new graphical device.
#' @param fine how many points between xlims?
#' @param ... optional graphical parameters to be passed to plot function.
plot_rat <- function(fs, gs, xlim = c(0, 10), sketch = TRUE, fine = 1000, ...) {

    ## fs and gs must be matrices
    fs <- as.matrix(fs)
    gs <- as.matrix(gs)

    x <- seq(xlim[1], xlim[2], length.out = fine)
    y <- eval_pols(fs, x) / eval_pols(gs, x)

    if (sketch) {
        lines(y[, 1] ~ x, ...)
    } else {
        plot(y[, 1] ~ x, type = 'l', ...)
    }

    ## If there are multiple rational functions, plot the rest
    if (NCOL(fs) > 1) {
        for (i in seq_len(NCOL(fs))[-1]) {
            lines(y[, i] ~ x, col = i+1, ...)
        }
    }
}
