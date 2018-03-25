source("oracle.R")


#' Generate feasible states
#'
#' Generate feasible states by constantly adding random noises to supplied
#' (non-feasible) states.
#'
#' This is a rejection sampler. It samples a state from \code{dist} centred on
#' \code{theta} repeatedly until a feasible state is attained (i.e. when
#' \code{indicator} returns 1.
#'
#' @param theta a matrix with each column corresponding to a current state. (num
#'     mat)
#' @param indicator an indicator function that takes a matrix and return a
#'     boolean vector, each element of which corresponds to a column of the
#'     input matrix, indicating its feasibility . (num mat -> bool vec)
#' @param N the number of feasible states. If this number is greater than the
#'     column size of \code{theta}, \code{theta} will be recycled. (num)
#' @param dist a sampler built-in in R (e.g. rnorm, rcauchy, rt). (function)
#' @param dist.para a list of parameters used in \code{dist}. (list)
#'
#' @return a list consisting of a matrix of transitioned states, a vector of
#'     objective function values corresponding to the transitioned states, and
#'     the acceptance proportion (useful for diagnostic purpose).
genMonoRats <- function(theta, indicator, N=NCOL(theta), dist=rcauchy,
                           dist.para=list(scale=2)) {

  ## Generate a neigbouring monotone rational function given a non-monotone
  ## rational function.

    force(indicator)
    force(dist)

    n_row <- NROW(theta)
    n_col <- NCOL(theta)

    if (N < n_col) {
        stop("Supplied estimates more than number of desired estimated")
    }

    theta <- as.matrix(theta)
    idx <- sample(seq(1, n_col), N - n_col, TRUE)
    prem <- cbind(theta, theta[, idx, drop=FALSE])
    final <- prem
    bad <- !indicator(final)

    for (k in 1:10000) {
        dist.para$n <- n_row * sum(bad)
        pertub <- matrix(do.call(dist, dist.para), nrow = n_row)

        ## add noise into the matrix
        final[, bad] <- prem[, bad, drop = FALSE] + pertub
        bad[bad] <- !(indicator(final[, bad, drop = FALSE]))

        if (!any(bad)) {
            return(final)
        }
    }

    stop('No monotone function found. Please give a new initial value.')
}


squaredNorm <- function(mat) {

    ## Compute the squared l2 norm for each column (only works for matrix!)

    if (is.matrix(mat)) {
        colSums(mat^2)
    } else {
        stop("Non-matrix in squaredNorm")
    }
}

vecToMat <- function(vec, n, colwise = TRUE) {

    ## Create a matrix with identical columns (or rows)

    if (colwise) {
        matrix(rep(vec, times = n), ncol = n, nrow = length(vec))
    } else {
        matrix(rep(vec, times = n), nrow = n, ncol = length(vec), byrow = TRUE)
    }
}

proposalFac <- function(dist=rnorm, dist.para=list(mean=0, sd=1),
                        fraction=0.97) {

    ## A random walk proposal distributionn, with pertubation distribution specific
    ## in "dist" and its parameters list "dist.para"

    ## dist : a R function that generates random numbers
    ## dist.para : the parameters of the distribution
    ## RETURN : a random walk proposal

    ## dist : R function
    ## dist.para : list
    ## RETURN : d x N numeric matrix -> d x N numeric matrix

    ## evaluate function's arguments eagerly
    force(dist)
    force(dist.para)

    ## initiate counter for decreasing standard deviation (for normal dist.)
    k <- 0

    function(theta) {
        ## vector must be supplied as a column matrix
        if (is.vector(theta)) {
            theta <- as.matrix(theta)
            cat("proposal input is a vector. \n")
        }
        dist.para$n <- prod(dim(theta))

        ## decreasing standard deviations (only valid for normal dist.)
        dist.para$sd <- dist.para$sd * (fraction)^k

        pertub <- matrix(do.call(dist, dist.para), nrow(theta), ncol(theta))
        k <<- k + 1
        return(theta + pertub)
    }
}
proposalMonoFac <- function(indicator, dist=rnorm, dist.para=list(mean=0, sd=1),
                            fraction=0.97) {

    ## A constrained random walk proposal distribution, with pertubation
    ## distribution specific in "dist" and its parameters list "dist.para"

    ## dist : a R function that generates random numbers
    ## dist.para : the parameters of the distribution
    ## indicator : an indicator function that implements constraints
    ## RETURN : a random walk proposal

    ## dist : R function
    ## dist.para : list
    ## indicator : d x N numeric matrix -> d x N numeric matrix
    ## RETURN : d x N numeric matrix -> d x N numeric matrix

    ## evaluate function's arguments eagerly
    force(dist)
    force(dist.para)
    force(indicator)

    ## initiate counter for decreasing standard deviation (for normal dist.)
    k <- 0

    function(theta) {
        ## vector must be supplied as a column matrix
        if (is.vector(theta)) {
            theta <- as.matrix(theta)
            cat("proposal input is a vector. \n")
        }
        nr <- NROW(theta)

        ## initialise final output
        final <- theta

        ## every column is bad in the beginning
        bad <- rep(TRUE, NCOL(theta))

        ## decreasing standard deviations (only valid for normal dist.)
        dist.para$sd <- dist.para$sd * fraction^k
        while (any(bad)) {
            nbad <- sum(bad)
            dist.para$n <- nr * nbad
            pertub <- matrix(do.call(dist, dist.para), nr, nbad)

            ## add
            final[, bad] <- theta[, bad, drop = FALSE] + pertub
            bad[bad] <- !(indicator(final[, bad, drop = FALSE]))
        }

        k <<- k + 1
        return(final)
    }
}

## BENCHMARKING
## xx <- matrix(rnorm(7000), 7)
## proposalMono <- proposalMonoFac(indicator = indicatorAll)
## yy <- proposalMono(xx)
## microbenchmark(proposalMono(xx))

recipSchedule <- function(k, T_init) {
    ## 0.9^k    # 0.8 < alpha < 0.9
    ## 1 / (1 + 1.5 * log(1 + k))    # alpha > 1
    ## 1 / (1 + 0.85 * k)    # alpha > 0
    T_init / (1 + 0.85 * (k-1)^2)    # alpha > 0
}

logSchedule <- function(k, T_init) {
    T_init / log(k + 1)
}

boltzmannFac <- function(f, T) {

    ## Create a Boltzmann kernel

    force(f)
    force(T)
    function(theta) { exp(- f(theta) / T) }
}

## CEPSO ##







