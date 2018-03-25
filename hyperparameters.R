#' Generate feasible states
#'
#' Generate feasible states by constantly adding random noises to supplied
#' (non-feasible) states.
#'
#' This is a rejection sampler that repeatedly adds noises to \code{theta}
#' feasible states are attained (i.e. when \code{indicator} returns 1). If
#' \code{N < ncol(theta)}, N states will be randomly sampled from \code{theta},
#' whereas if \code{N > ncol(theta)}, additional states will be sampled from
#' \code{theta} to make up the number.
#'
#' This function will terminate after 10000 of failed trials.
#'
#' @param theta a matrix with each column corresponding to an (unfeasible)
#'     state. (num mat)
#' @param indicator an indicator function that takes a matrix and return a
#'     boolean vector, each element of which corresponds to a column of the
#'     input matrix, indicating its feasibility . (num mat -> bool vec)
#' @param N the number of feasible states.
#' @param dist the distribution of the random noise. This must be a sampler such
#'     as rnorm, rcauchy and rt. (function)
#' @param dist.para a list of parameters used in \code{dist}. (list)
#'
#' @return a list consisting of a matrix of transitioned states, a vector of
#'     objective function values corresponding to the transitioned states, and
#'     the acceptance proportion (useful for diagnostic purpose).
get_feasibles <- function(theta, indicator, N=NCOL(theta), dist=rcauchy,
                          dist_para=list(scale=2)) {

    force(indicator)
    force(dist)

    theta <- as.matrix(theta)
    n_dims <- nrow(theta)
    n_theta <- ncol(theta)

    if (N < n_theta) {
        warning("Supplied estimates more than number of desired estimated")
        idx <- sample(1:n_theta, N, TRUE)
    } else {
        idx <- c(1:n_theta, sample(seq(1, n_theta), N - n_theta, TRUE))
    }
    theta <- theta[, idx]
    n_theta <- ncol(theta)

    fea_states <- theta
    bad <- !indicator(fea_states)

    for (k in 1:10000) {
        ## Add noise into the matrix
        dist_para$n <- n_dims * sum(bad)
        pertub <- matrix(do.call(dist, dist_para), nrow = n_dims)
        fea_states[, bad] <- theta[, bad, drop = FALSE] + pertub

        ## Are they all feasible now?
        bad[bad] <- !(indicator(fea_states[, bad, drop = FALSE]))

        if (!any(bad)) {
            return(fea_states)
        }
    }
    stop('No monotone function found. Please give a new initial value.')
}


#' Random walk normal distribution (factory)
#'
#' This is a factory that creates a random walk normal distribution that adds
#' some noise to the provided values.
#'
#' This proposal is designed to decrease after each iteration. The default
#' option is to decrese 3 percentage from the standard deviation at the previous
#' iteration.
#'
#' @param sd the standard deviation of the normal distribution at the first
#'     iteration, that is when \code{k = 1}. (num)
#' @param fraction the percentage of standard deviation (relative to the
#'     previous iteration) left after each iteration. (num)
#'
#' @return a function that takes a matrix of current states (each column
#'     corresponds to a state) (num mat) and the current iteration count (num),
#'     and produces a matrix of proposed states (num mat).
rnorm_rw_fac <- function(sd = 1, fraction = 0.97) {

    ## Evaluate the function's argument eagerly
    force(fraction)

    function(theta, k) {
        ## theta must be a matrix
        theta <- as.matrix(theta)

        ## Decrease standard deviations and generate noise
        pertub <- rnorm(length(theta), sd = sd * (fraction^(k - 1)))

        theta + pertub
    }
}



#' Random walk truncated normal distribution (factory)
#'
#' This is a factory that creates a random walk normal distribution that adds
#' some noise to the provided values, subject to an indicator constraint.
#'
#' This factory essentially creates a rejection sampler that samples from a
#' truncated normal distribution.
#'
#' This proposal is designed to decrease after each iteration. The default
#' option is to decrese 3 percentage from the "standard deviation" at the
#' previous iteration. The "standard deviation" mentioned in this section refers
#' to the standard deviation of the untruncated normal random walk, not the
#' actual standard deviation of the proposal distribution.
#'
#' @param indicator an indicator function that takes a matrix and return a
#'     boolean vector, each element of which corresponds to a column of the
#'     input matrix, indicating its feasibility . (num mat -> bool vec)
#' @param sd the "standard deviation" of the normal distribution at the first
#'     iteration, that is when \code{k = 1}. (num)
#' @param fraction the percentage of "standard deviation" (relative to the
#'     previous iteration) left after each iteration. (num)
#'
#' @return a function that takes a matrix of current states (each column
#'     corresponds to a state) (num mat) and the current iteration count (num),
#'     and produces a matrix of proposed states (num mat).
rtnorm_rw_fac <- function(indicator, sd = 1, fraction = 0.97) {

    ## Evaluate the function's arguments eagerly
    force(indicator)
    force(fraction)

    function(theta, k) {
        ## theta must be a matrix
        theta <- as.matrix(theta)
        n_dims <- NROW(theta)

        ## Initialise final output
        fea_states <- theta

        ## Every column is bad in the beginning, to make sure that noise is
        ## added at least once.
        bad <- rep(TRUE, NCOL(theta))

        while (any(bad)) {
            n_bad <- sum(bad)
            ## Decrease standard deviations and generate noise
            pertub <- rnorm(n_dims * n_bad, sd = sd * (fraction^(k - 1)))

            ## Add pertubation to the current states
            fea_states[, bad] <- theta[, bad, drop = FALSE] + pertub
            bad[bad] <- !(indicator(fea_states[, bad, drop = FALSE]))
        }
        fea_states
    }
}

recip_schedule <- function(k, T_init) {
    ## alpha = 0.85
    T_init / (1 + 0.85 * (k-1)^2)
}

log_schedule <- function(k, T_init) {
    T_init / log(k + 1)
}








