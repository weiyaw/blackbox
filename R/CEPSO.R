vecToMat <- function(vec, n, colwise = TRUE) {

    ## Create a matrix with identical columns (or rows)

    if (colwise) {
        matrix(rep(vec, times = n), ncol = n, nrow = length(vec))
    } else {
        matrix(rep(vec, times = n), nrow = n, ncol = length(vec), byrow = TRUE)
    }
}

squaredNorm <- function(mat) {

    ## Compute the squared l2 norm for each column (only works for matrix!)

    if (is.matrix(mat)) {
        colSums(mat^2)
    } else {
        stop("Non-matrix in squaredNorm")
    }
}


#' @title Update the inertia "w" for v1.
#'
#' @param in_boundary a boolean vector specifying whether the particles are
#'     within the boundary. (bool vec)
updateInertia <- function(in_boundary) {

    rho <- mean(in_boundary)
    if (rho < 0.9) 0.5
    else if (rho < 1) 0.73
    else 1
}

#' Factory of closeness function.
#'
#' The "closeness" function described in the paper. It is defined as the
#' difference of the objective function values between pilot and current
#' states.
#'
#' @param pilot an pilot estimate. (num vec)
#' @param objf the objective function. See CEPSO for more details. (num mat ->
#'     num vec)
#'
#' @return a 'closeness' function that takes in a particle matrix and returns a
#'     vector of delta values.
deltaFac <- function(pilot, objf) {
    force(objf)
    pilot_objv <- objf(pilot)
    function(s) {
        return(objf(s) - pilot_objv)
    }
}

#' Initialise particle swarm.
#'
#' Initiate a swarm of particles with the given starting values.
#'
#' @param starting a matrix of starting values. (num mat)
#' @param delta the closeness function that takes a matrix and return a
#'     vector. (num mat -> num vec)
#' @param indicator an indicator function that takes a matrix, and returns a
#'     boolean vector that indicates the feasibility of each column of the
#'     matrix. (num mat -> bool vec)
#'
#' @return a list with various properties of the swarm.
initSwarm <- function(starting, delta, indicator) {

    ## Return a list with the following properties.
    ## particles   : the particles. (num mat)
    ## vp          : the initial velocities for each particles, which are 0.
    ##               (num mat)
    ## objv        : the delta values associated to each particles. (num vec)
    ## is_feasible : specifiying the feasibility of each particles. (bool vec)
    ## pbest       : the best states visited by each particles so far. (num mat)
    ## pbest_objv  : the delta values associated to the particles' best state.
    ##               (num vec)
    ## worst_pbest : the worst state among the particles' bests, in the current
    ##               swarm. (num vec)
    ## wpbest_objv : the delta value of the worst particles' bests. (num)
    ## gbest       : the best state seen so far from the start of the algorithm.
    ##               (num vec)
    ## gbest_objv  : the delta value of the current global minimum. (num)

    force(delta)
    force(indicator)

    sw <- list()
    n_dims <- NROW(starting)
    n_parts <- NCOL(starting)

    ## Starting particles
    if (!(is.matrix(starting))) {
        if (is.vector(starting)) {
            starting <- as.matrix(starting)
        } else {
            stop("Non-matrix starting particles")
        }
    }
    sw$particles <- starting

    ## Initial velocities and objective values
    sw$vp <- matrix(0, n_dims, n_parts)
    sw$objv <- delta(sw$particles)

    ## Determime the feasibility of a particle
    sw$is_feasible <- indicator(sw$particles)

    ## Initialise presonal bests, if possible
    sw$pbest <- matrix(NA, n_dims, n_parts)
    sw$pbest_objv <- rep(Inf, n_parts)
    sw$pbest[, sw$is_feasible] <- sw$particles[, sw$is_feasible]
    sw$pbest_objv[sw$is_feasible] <- sw$objv[sw$is_feasible]

    ## Find the worst and best personal best (a.k.a. global best)
    if (all(is.na(sw$pbest[1, ]))) {
        ## If feasible states haven't been visited
        sw$worst_pbest <- matrix(NA, n_dims, 1)
        sw$wpbest_objv <- Inf
        sw$gbest <- matrix(NA, n_dims, 1)
        sw$gbest_objv <- Inf
    } else {
        temp_objv <- sw$pbest_objv
        temp_objv[is.infinite(sw$pbest_objv)] <- NA
        max_index <- which.max(temp_objv)
        min_index <- which.min(temp_objv)
        sw$worst_pbest <- sw$pbest[, max_index, drop=FALSE]
        sw$gbest <- sw$pbest[, min_index, drop=FALSE]
        sw$wpbest_objv <- sw$objv[max_index]
        sw$gbest_objv <- sw$objv[min_index]
    }
    return(sw)
}


#' Update a particle swarm.
#'
#' Update the properties of an particle swarm.
#'
#' @param sw the swarm produced by 'initSwarm'. (list)
#' @param delta the closeness function that takes a matrix and return a
#'     vector. (num mat -> num vec)
#' @param indicator an indicator function that takes a matrix, and returns a
#'     boolean vector that indicates the feasibility of each column of the
#'     matrix. (num mat -> bool vec)
#'
#' @return a list with various properties of the swarm.
assess <- function(sw, delta, indicator) {

    force(delta)
    force(indicator)

    if (!all(c("particles", "vp", "nature", "pbest", "pbest_objv", "gbest",
               "gbest_objv") %in%
             names(sw))) {
        stop("Assess required information not available")
    }
    n_dims <- NROW(sw$particles)
    n_parts <- NCOL(sw$particles)

    ## Update objective values for each particle
    sw$objv <- delta(sw$particles)

    ## Determime the feasibility of a particle
    sw$is_feasible <- indicator(sw$particles)

    ## Update the personal best of each particle
    update_pbest <- (sw$objv < sw$pbest_objv) & sw$is_feasible
    sw$pbest[, update_pbest] <- sw$particles[, update_pbest]
    sw$pbest_objv[update_pbest] <- sw$objv[update_pbest]

    ## Find the worst and best personal best (a.k.a. global best)
    ## FIX THIS
    if (all(is.na(sw$pbest[1, ]))) {
        sw$worst_pbest <- matrix(NA, n_dims, 1)
        sw$wpbest_objv <- Inf
        sw$gbest <- matrix(NA, n_dims, 1)
        sw$gbest_objv <- Inf
    } else {
        temp_objv <- sw$pbest_objv
        temp_objv[is.infinite(sw$pbest_objv)] <- NA
        max_index <- which.max(temp_objv)
        min_index <- which.min(temp_objv)
        sw$worst_pbest <- sw$pbest[, max_index, drop=FALSE]
        sw$wpbest_objv <- sw$pbest_objv[max_index]
        if (sw$pbest_objv[min_index] < sw$gbest_objv) {
            sw$gbest <- sw$pbest[, min_index, drop=FALSE]
            sw$gbest_objv <- sw$pbest_objv[min_index]
        }
    }

    return(sw)
}

#' Do an exploration step.
#'
#' Do an exploration step for a given particle swarm.
#'
#' @param sw the swarm produced by 'initSwarm'. (list)
#' @param delta the closeness function that takes a matrix and return a
#'     vector. (num mat -> num vec)
#' @param indicator an indicator function that takes a matrix, and returns a
#'     boolean vector that indicates the feasibility of each column of the
#'     matrix. (num mat -> bool vec)
#'
#' @return a list with various properties of the swarm.
explore <- function(sw, pilot) {

    if (all(c("particles", "vp", "nature", "in_explore", "gbest", "gbest_objv",
               "pbest", "pbest_objv") %in% names(sw))) {
        if (sw$nature != "explore") {
            warning("Non-exploit swarm running on exploit") }
    } else {
        stop("Not enough information for explore swarm")
    }

    n_dims <- NROW(sw$particles)
    n_parts <- NCOL(sw$particles)

    if (n_parts == 0L) {
        return(sw)
    }

    ## Evaluate v1
    w <- updateInertia(sw$in_explore)
    v1 <- sw$vp
    v1[, !sw$in_explore] <- 0

    ## Evaluate v2
    pilot_col <- vecToMat(pilot, n_parts)
    pilot_to_parts <- pilot_col - sw$particles
    v2 <- matrix(runif(n_dims * n_parts), n_dims, n_parts) * pilot_to_parts

    ## Create a random global best for each particles if it does not exist
    if (is.na(sw$gbest[1])) {
        direction <- matrix(runif(n_dims * n_parts), n_dims, n_parts)
        dist <- sqrt(colSums((pilot_to_parts)^2) / colSums(direction^2))
        gbest_col <- pilot_col + (direction * vecToMat(dist, n_dims, FALSE))
    } else {
        gbest_col <- vecToMat(sw$gbest, n_parts)
    }

    ## Evaluate v3
    v3 <- matrix(runif(n_dims * n_parts), n_dims, n_parts) *
        (gbest_col - sw$particles)

    ## Update position
    sw$vp <- w*v1 + 1.5*v2 + 1.5*v3
    sw$particles <- sw$particles + sw$vp

    sw$in_explore <- NULL
    sw$k_pbest <- NULL
    return(sw)
}

#' Do an exploitation step.
#'
#' Do an exploitation step for a given particle swarm.
#'
#' @param sw the swarm produced by 'initSwarm'. (list)
#' @param delta the closeness function that takes a matrix and return a
#'     vector. (num mat -> num vec)
#' @param indicator an indicator function that takes a matrix, and returns a
#'     boolean vector that indicates the feasibility of each column of the
#'     matrix. (num mat -> bool vec)
#'
#' @return a list with various properties of the swarm.
exploit <- function(sw, pilot) {

    if (all(c("particles", "vp", "nature", "in_exploit", "pbest", "pbest_objv",
              "gbest", "gbest_objv", "k_pbest") %in% names(sw))) {
        if (sw$nature != "exploit") {
            warning("Non-exploit swarm running on exploit") }
    } else {
        stop("Not enough information for exploit swarm")
    }

    n_dims <- NROW(sw$particles)
    n_parts <- NCOL(sw$particles)

    if (n_parts == 0L) {
        return(sw)
    }

    ## Evaluate v1
    w <- updateInertia(sw$in_exploit)
    v1 <- sw$vp
    v1[, !sw$in_exploit] <- 0

    ## Evaluate v2
    v2 <- matrix(runif(n_dims * n_parts), n_dims, n_parts) *
        (sw$pbest - sw$particles)

    ## Evaluate v3
    r2 <- matrix(runif(n_dims * n_parts), n_dims, n_parts)
    v3 <- r2 * (sw$k_pbest - sw$particles)

    ## If neither pbest nor k_pbest exist
    neither <- is.na(sw$pbest) & is.na(sw$k_pbest)
    v3[neither] <- r2[neither] * (vecToMat(sw$gbest, n_parts) - sw$pbest)[neither]

    v2[is.na(v2)] <- 0
    v3[is.na(v3)] <- 0

    ## Update position
    sw$vp <- w*v1 + 1.5*v2 + 1.5*v3
    sw$particles <- sw$particles + sw$vp

    sw$in_exploit <- NULL
    sw$k_pbest <- NULL
    return(sw)
}



#' Exchange between two swarm.
#'
#' Exchange the particle between the exploitation and exploration swarms.
#'
#' @param sw the swarm produced by 'initSwarm'. (list)
#' @param delta the closeness function that takes a matrix and return a
#'     vector. (num mat -> num vec)
#' @param indicator an indicator function that takes a matrix, and returns a
#'     boolean vector that indicates the feasibility of each column of the
#'     matrix. (num mat -> bool vec)
#'
#' @return a list with various properties of the swarm.
exchange <- function(sw1, sw2, k) {

    if (!all(c("particles", "vp", "objv", "is_feasible", "pbest", "pbest_objv",
              "worst_pbest", "wpbest_objv", "gbest", "gbest_objv")
            %in% names(sw1))) {
        stop("Not enough information for the first swarm")
    }

    if (!all(c("particles", "vp", "objv", "is_feasible", "pbest", "pbest_objv",
              "worst_pbest", "wpbest_objv", "gbest", "gbest_objv")
            %in% names(sw2))) {
        stop("Not enough information for the second swarm")
    }

    exploit_sw <- list()
    explore_sw <- list()

    ## Update global best and worst personal best
    if (is.finite(sw1$wpbest_objv) && is.finite(sw2$wpbest_objv)) {
        if (sw1$wpbest_objv > sw2$wpbest_objv) {
            worst_pbest <- sw1$worst_pbest
            wpbest_objv <- sw1$wpbest_objv
        } else {
            worst_pbest <- sw2$worst_pbest
            wpbest_objv <- sw2$wpbest_objv
        }
    } else if (is.finite(sw1$wpbest_objv)) {
        ## If sw2 doesn't have a feasible state
        worst_pbest <- sw1$worst_pbest
        wpbest_objv <- sw1$wpbest_objv
    } else if (is.finite(sw2$wpbest_objv)) {
        ## If sw1 doesn't have a feasible state
        worst_pbest <- sw2$worst_pbest
        wpbest_objv <- sw2$wpbest_objv
    } else {
        ## If both swarm don't have a feasible state
        worst_pbest <- NA
        wpbest_objv <- Inf
    }

    ## If both swarm don't have a feasible state, gbest will be unavailable.
    if (sw1$gbest_objv < sw2$gbest_objv) {
        gbest  <- sw1$gbest
        gbest_objv <- sw1$gbest_objv
    } else {
        gbest <- sw2$gbest
        gbest_objv <- sw2$gbest_objv
    }


    ## Determine the promising particles
    sw1_n_parts <- NCOL(sw1$particles)
    sw2_n_parts <- NCOL(sw2$particles)

    if (is.na(gbest[1])) {
        ## If a feasible state is yet to be found
        sw1_promising <- rep(FALSE, sw1_n_parts)
        sw2_promising <- rep(FALSE, sw2_n_parts)
    } else {
        sw1_in_explore <- sw1$objv <= gbest_objv
        sw2_in_explore <- sw2$objv <= gbest_objv
        sw1_parts_to_gbest <- sw1$particles - vecToMat(gbest, sw1_n_parts)
        sw2_parts_to_gbest <- sw2$particles - vecToMat(gbest, sw2_n_parts)
        sw1_wpbest_to_gbest <- vecToMat(worst_pbest - gbest, sw1_n_parts)
        sw2_wpbest_to_gbest <- vecToMat(worst_pbest - gbest, sw2_n_parts)
        sw1_promising <- (squaredNorm(sw1_parts_to_gbest) <=
                          squaredNorm(sw1_wpbest_to_gbest)) & sw1_in_explore
        sw2_promising <- (squaredNorm(sw2_parts_to_gbest) <=
                          squaredNorm(sw2_wpbest_to_gbest)) & sw2_in_explore
    }

    ## Determine the particles to be exchanged
    exploit_sw1 <- sw1$is_feasible | sw1_promising
    exploit_sw2 <- sw2$is_feasible | sw2_promising

    ## Exchange particles
    exploit_sw$particles <- cbind(sw1$particles[, exploit_sw1, drop=FALSE],
                                  sw2$particles[, exploit_sw2, drop=FALSE])
    explore_sw$particles <- cbind(sw1$particles[, !exploit_sw1, drop=FALSE],
                                  sw2$particles[, !exploit_sw2, drop=FALSE])

    ## Exchange vp

    exploit_sw$vp <- cbind(sw1$vp[, exploit_sw1, drop=FALSE],
                           sw2$vp[, exploit_sw2, drop=FALSE])
    explore_sw$vp <- cbind(sw1$vp[, !exploit_sw1, drop=FALSE],
                           sw2$vp[, !exploit_sw2, drop=FALSE])

    exploit_sw$nature <- "exploit"
    explore_sw$nature <- "explore"

    ## Update in_exploit and in_explore
    exploit_objv <- c(sw1$objv[exploit_sw1], sw2$objv[exploit_sw2])
    explore_objv <- c(sw1$objv[!exploit_sw1], sw2$objv[!exploit_sw2])

    exploit_sw$in_exploit <- exploit_objv <= wpbest_objv
    explore_sw$in_explore <- explore_objv <= gbest_objv

    ## Insert current global minimum
    exploit_sw$gbest <- explore_sw$gbest <- gbest
    exploit_sw$gbest_objv <- explore_sw$gbest_objv <- gbest_objv

    ## Always keep track of personal best for each particles
    exploit_sw$pbest <- cbind(sw1$pbest[, exploit_sw1, drop=FALSE],
                              sw2$pbest[, exploit_sw2, drop=FALSE])
    explore_sw$pbest <- cbind(sw1$pbest[, !exploit_sw1, drop=FALSE],
                              sw2$pbest[, !exploit_sw2, drop=FALSE])

    exploit_sw$pbest_objv <- c(sw1$pbest_objv[exploit_sw1],
                               sw2$pbest_objv[exploit_sw2])
    explore_sw$pbest_objv <- c(sw1$pbest_objv[!exploit_sw1],
                               sw2$pbest_objv[!exploit_sw2])

    ## Calculate the kth nearest best particle. Only calculated for the exploit
    ## group.
    n_dims <- NROW(exploit_sw$particles)
    n_parts <- NCOL(exploit_sw$particles)

    if (k > n_parts) {
        ## warning("Number of nearest neighbour greater than number of particles")
        k <- n_parts
    }
    exploit_sw$k_pbest <- matrix(NA, n_dims, n_parts)

    if (n_parts != 0L && k != 0) {
        ## Calculate the indices of the kth nearest particles
        indices <- (vecToMat(seq_len(n_parts), k, FALSE) +
                    vecToMat(seq(0, k-1), n_parts)) %% n_parts
        indices[indices == 0] <- n_parts

        for (i in seq_len(n_parts)) {
            k_objv <- exploit_sw$pbest_objv[indices[, i]]
            if (any(is.finite(k_objv))) {
                relative <- which.min(k_objv)
                absolute <- indices[relative, i]
                exploit_sw$k_pbest[, i] <- exploit_sw$particles[, absolute]
            }
        }
    }

    return(list(exploit_sw, explore_sw))
}



#' Constrained minimisation with CEPSO
#'
#' Minimise an objective, subject to the constraint specified by an indicator
#' function.
#'
#' @param starting a matrix of starting values with each column corresponds to a
#'     particle. (num mat)
#' @param objf the objective function that takes a matrix with each column
#'     corresponds to a particle, and return a vector with respective objective
#'     function values. (num mat -> num mat)
#' @param pilot a pilot estimate vector. (num vec)
#' @param indicator an indicator function that takes a particle matrix and
#'     return a boolean vector, each columnm of which corresponds to a column of
#'     the input matrix. (num mat -> bool vec)
#' @param neighbour the number of neighbour (num)
#' @param iter iteration number. (num)
#' @param N number of total particles. (num)
#' @param verbose show the algorithm progression? (bool)
#'
#' @return a list that contains the best particle and its objective function
#'     value.
CEPSO <- function(starting, objf, pilot, indicator, neighbour = N/5, iter = 2000,
                  N = 100, verbose = FALSE) {

    force(objf)
    force(pilot)
    force(indicator)
    force(neighbour)

    names(starting) <- NULL
    names(pilot) <- NULL

    if (N < 2) {
        stop("The number of particles must be greater than 2")
    }

    ## if the supplied starting values is less than N
    if (NCOL(starting) < N) {
        starting <- starting[, rep_len(1:NCOL(starting), N)]
    }

    ## construct the 'closeness' function
    delta <- deltaFac(as.matrix(pilot), objf)

    ## calculate the objective function
    pilot_objv <- objf(as.matrix(pilot))

    ## initialise exploit and explore swarms
    sw1 <- initSwarm(starting[, 1, drop=FALSE], delta, indicator)
    sw2 <- initSwarm(starting[, -1, drop=FALSE], delta, indicator)
    sw <- exchange(sw1, sw2, neighbour)


    ## records best objf at 50, 100, 200, 400, 600, 1000, 2000, 3000
    records <- rep(NA, length = 8)
    names(records) <- c(50, 100, 200, 400, 600, 1000, 2000, 3000)

    ## iterates
    for (i in seq(1, iter)) {
        sw1 <- exploit(sw[[1]], pilot)
        sw1 <- assess(sw1, delta, indicator)
        sw2 <- explore(sw[[2]], pilot)
        sw2 <- assess(sw2, delta, indicator)
        sw <- exchange(sw1, sw2, neighbour)
        if (verbose && (i %% 200 == 0)) {
            cat(i, 'iterations done, current best:',
                pilot_objv + sw[[1]]$gbest_objv, '\n')
        }

        ## records best objf at 50, 100, 200, 400, 600, 1000, 2000, 3000
        if (i %in% c(50, 100, 200, 400, 600, 1000, 2000, 3000)) {
            records[paste(i)] <- pilot_objv + sw[[1]]$gbest_objv
    }

    }
    final_objv <- pilot_objv + sw[[1]]$gbest_objv
    cat('Final objective value:', final_objv, '\n')
    return(list(theta = sw[[1]]$gbest, objv = final_objv, records = records))
}
