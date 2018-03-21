numDeriv <- function(f, g, EPS = 1e-6) {

  # Calculate numerator of derivative
  # f: coefficients vector of numerator (ascending)
  # g: coefficients vector of denominator (ascending)
  #
  # Return: 1 x n array

    lg <- length(g)
    lf <- length(f)
    if (lg == 1) {
        stop('This is a polynomial')
    } else {
        ## Derivatives of denominator
        g.deriv <- g[seq.int(2, lg)] * seq_len(lg - 1)
    }

    if (lf == 1) {
        ## Derivatives of numerator (if it's a constant)
        temp <- -f * g.deriv
    } else {
        ## If it's not a constant
        f.deriv <- f[seq.int(2, lf)] * seq_len(lf - 1)
        temp <- convolve(f.deriv, rev(g), type = 'o') -
            convolve(f, rev(g.deriv), type = 'o')
    }

    ## Get rid of trailing zeros
    ## FAST IMPLEMENTATION, MIGHT NOT WORK IN CERTAIN CASES
    if (lg == lf) {
        return(temp[-length(temp)])
    } else {
        return(temp)
    }

    ## while (abs(temp[length(temp)]) < EPS) {
    ##     temp <- temp[-length(temp)]
    ## }
}

numDerivs <- function(fs, gs, EPS = 1e-6) {

  # Calculate numerator of derivative
  # f: coefficients matrix of numerator (ascending)
  # g: coefficients matrix of denominator (ascending)
  #
  # Return: 1 x n array

    lg <- NROW(gs)
    lf <- NROW(fs)

    if (lg == 1) {
        ## if the denominator is a constant
        stop('This is a polynomial')

    } else if (lf == 1) {
        ## if the numerator is a constant
        ## derivatives of denominator (no reversing)
        gs_deriv <- gs[seq.int(2, lg), , drop = FALSE]
        for (i in seq_len(lg - 1)[-1]) {
            gs_deriv[i, ] <- gs_deriv[i, ] * i
        }
        ## result
        temp <- for (i in seq_len(NROW(gs_deriv))) {-fs * gs_deriv[i, ]}

    } else {
        ## if the numerator is a polynomial
        ## derivatives of denominator (reversed)
        gs_deriv_r <- gs[rev(seq.int(2, lg)), , drop = FALSE]
        for (i in seq_len(lg - 2)) {
            gs_deriv_r[i, ] <- gs_deriv_r[i, ] * (lg - i)
        }
        ## derivatives of numerator (reversed)
        fs_deriv_r <- fs[rev(seq.int(2, lf)), , drop = FALSE]
        for (i in seq_len(lf - 2)) {
            fs_deriv_r[i, ] <- fs_deriv_r[i, ] * (lf - i)
        }
        ## result
        temp <- matrix(NA, nrow = lg + lf - 2, ncol = NCOL(fs))
        for (i in seq_len(NCOL(fs))) {
            temp[, i] <- convolve(gs[, i], fs_deriv_r[, i], type = 'o') -
                convolve(fs[, i], gs_deriv_r[, i], type = 'o')
        }
    }

    ## Get rid of trailing zeros
    ## FAST IMPLEMENTATION, MIGHT NOT WORK IN CERTAIN CASES
    if (lg == lf) {
        return(temp[-NROW(temp), , drop = FALSE])
    } else {
        return(temp)
    }
}

evalPol <- function (x, beta) {

  # Routine in 'MonoPoly'
  # Evaluate a polynomial at 'x' using Horner scheme

  res <- 0
  for (bi in rev(beta)) res <- res * x + bi
  res
}

evalPols <- function (x, beta) {

    ## Modified routine in 'MonoPoly'
    ## Evaluate multiple polynomials at 'x' using Horner scheme

    ## res <- 0
    ## for (bi in NROW(beta):1) res <- res * x + beta[bi, , drop = F]
    ## res
    nbeta <- NCOL(beta)
    nx <- length(x)
    res <- matrix(0, nx, nbeta)
    X <- matrix(rep(x, times = nbeta), nx, nbeta)
    for (bi in seq.int(NROW(beta), 1)) {
        res <- X*res +
            matrix(rep(beta[bi, ], times = nx), nx, nbeta, byrow = TRUE)
    }
    res
}

isPositive <- function(poly.coef, a = 0, b = Inf, positive = TRUE, EPS = 1e-06) {

  # Modified Kev & Berwin's method in 'MonoPoly'
  # Check whether polynomial 'poly.coef' is positive over the region (a, b)
  # Assumption: 'poly.coef' is an ascending polynomial coefficient vector
  # Return: boolean

  if (is.infinite(b) & (poly.coef[length(poly.coef)] > -EPS) != positive) {
    # cat('infinite \n')
    return(FALSE)
  }
  if (is.finite(b) & (evalPol(b, poly.coef) > -EPS) != positive) {
    # cat('finite \n')
    return(FALSE)
  }

  roots <- polyroot(poly.coef)
  re.roots <- Re(roots)[abs(Im(roots)) < EPS]
  re.roots <- re.roots[a < re.roots & re.roots < b]
  if (length(re.roots) == 0) {
    return(TRUE)
  } else {
    mltplcty <- rowSums(outer(re.roots, re.roots,
        function(x, y) abs(x - y) < EPS))
    return(all(mltplcty%%2 == 0))
  }
}

isPositives <- function(poly_coefs, a, b, positive = TRUE, EPS = 1e-06) {
    if (is.infinite(b)) {
        k_going <- (poly_coefs[NROW(poly_coefs), ] >= -EPS) == positive
    } else {
        k_going <- (evalPols(b, poly_coefs) >= -EPS) == positive
    }
    if (any(k_going)) {
        singlePositive <- function(one_poly) {
              roots <- polyroot(one_poly)
              re.roots <- Re(roots)[abs(Im(roots)) < EPS]
              re.roots <- re.roots[a < re.roots & re.roots < b]
              if (length(re.roots) == 0) {
                  return(TRUE)
              } else {
                  mltplcty <- rowSums(outer(re.roots, re.roots,
                                            function(x, y) abs(x - y) < EPS))
                  return(all(mltplcty%%2 == 0))
              }
        }
        k_going[k_going]<- apply(poly_coefs[, k_going, drop = FALSE], 2, singlePositive)
    }
    return(k_going)
}

## microbenchmark(apply(matrix(rnorm(5000, sd = 100), 5),2,isPositive), isPositives(matrix(rnorm(5000, sd = 100), 5),0,Inf))

noRoots <- function(poly_coef, a, b, EPS = 1e-06) {

    ## Check whether polynomial 'poly_coef' has roots over the region
    ## [a, b].
    ## Assumption: 'poly_coef' is an ascending polynomial coefficient vector
    ## Return: boolean
    roots <- polyroot(poly_coef)
    re.roots <- Re(roots)[abs(Im(roots)) < EPS]
    re.roots <- re.roots[a <= re.roots & re.roots <= b]
    return(length(re.roots) == 0)
}


noRootsFac <- function(a, b, EPS = 1e-06) {

  # Check whether polynomial 'poly_coef' has roots over the region
  # [a, b].
  # Assumption: 'poly_coef' is an ascending polynomial coefficient vector
                                        # Return: boolean
    force(a)
    force(b)
    force(EPS)
    function(poly_coef) {
        roots <- polyroot(poly_coef)
        re.roots <- Re(roots)[abs(Im(roots)) < EPS]
        re.roots <- re.roots[a <= re.roots & re.roots <= b]
        return(length(re.roots) == 0)
    }
}


oracle <- function(f, g, a = 0, b = Inf, increasing = TRUE,
                   EPS = 1e-06) {

    return(noRoots(g, a, b, EPS) &
           isPositive(numDeriv(f, g), a, b, positive = increasing, EPS))
}

indicatorFac <- function(dim_f, dim_g, xmin = 0, xmax = Inf, increasing = TRUE,
                          EPS = 1e-06) {
    force(xmin)
    force(xmax)
    force(increasing)
    force(EPS)
    idx_f <- seq.int(1, length.out = dim_f)
    idx_g <- seq.int(dim_f + 1, length.out = dim_g)
    noRoots <- noRootsFac(xmin, xmax, EPS)

    function(column)
        return(noRoots(column[idx_g]) &
               isPositive(numDeriv(column[idx_f], column[idx_g]),
                          xmin, xmax, increasing, EPS))
}

##
indicatorAllFac <- function(dim_f, dim_g, xmin = 0, xmax = Inf, increasing = TRUE,
                            EPS = 1e-06) {
    force(xmin)
    force(xmax)
    force(increasing)
    force(EPS)
    idx_f <- seq.int(1, length.out = dim_f)
    idx_g <- seq.int(dim_f + 1, length.out = dim_g)
    n_dims <- dim_f + dim_g
    noRoots <- noRootsFac(xmin, xmax, EPS)

    function(mat) {
        if (NROW(mat) != n_dims) {
            stop("Input's dimension mismatch with specified dimensions")
        }
        mat_g <- rbind(rep(1, NCOL(mat)), mat[idx_g, , drop = FALSE])
        k_going <- apply(mat_g, 2, noRoots)
        if (any(k_going)) {
            dervs <- numDerivs(mat[idx_f, k_going, drop = FALSE],
                               mat_g[, k_going, drop = FALSE],
                               EPS)
            k_going[k_going] <- isPositives(dervs, xmin, xmax, increasing, EPS)
            return(k_going)
        }
        else {
            return(k_going)
        }
    }
}

revmat <- function(mat) {
    mat[,NCOL(mat):1L]
}

## First test
## indicator <- indicatorFac(3,4)
## indicatorAll <- indicatorAllFac(3,4)
## xx <- matrix(rnorm(7000, sd = 100), 7)
## microbenchmark(apply(xx, 2, indicator), indicatorAll(xx))
## all(apply(xx, 2, indicator) == indicatorAll(xx))

## err1 <- xx[, ncol(xx)]
## err2 <- xx[, ncol(xx)-1]
## err <- cbind(err1, err2)

## apply(err, 2, indicator)
## indicatorAll(err)
