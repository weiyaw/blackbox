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


evalPol <- function (x, beta) {

  # Routine in 'MonoPoly'
  # Evaluate a polynomial at 'x' using Horner scheme

  res <- 0
  for (bi in rev(beta)) res <- res * x + bi
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
