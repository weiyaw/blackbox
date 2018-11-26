## ALL CODE FOR SA-BASED OPTIMISATION

#' SA-Move (Markov kernel)
#'
#' Performs a Markov transition step for each state. This function is intended
#' for minimisation.
#'
#' @param state a matrix with each column corresponding to a current state. (num
#'     mat)
#' @param objf an objective function that takes a matrix with each column
#'     corresponds to a state, and return a vector with respective objective
#'     function values. (num mat -> num vec)
#' @param proposal a random walk proposal that takes a matrix, and produces a
#'     matrix of proposed states, centred on the corresponding columns of the
#'     input matrix. It should also take in the current iteration count. ((num
#'     mat, num) -> num mat)
#' @param temp temperature of the current iteration. (num)
#' @param k current iteration count. (num)
#' @param objv an optional vector specifying the objective function values for
#'     each column in \code{state}. (num vec)
#'
#' @return a list consisting of a matrix of transitioned states, a vector of
#'     objective function values corresponding to the transitioned states, and
#'     the acceptance proportion (useful for diagnostic purpose).
SAmove <- function(state, objf, proposal, temp, k, objv = NULL){

    ## WARNING: THIS IS A MINIMISATION ALGORITHM
    ## WARNING: THIS IS A MINIMISATION ALGORITHM
    ## WARNING: THIS IS A MINIMISATION ALGORITHM


    if (is.vector(state)) {
        state <- t(as.matrix(state))
    }

    if (is.null(objv)) {
        objv <- objf(state)
    }

    state_t <- proposal(state, k)
    objv_t <- objf(state_t)

    ## Acceptance step
    u <- runif(NCOL(state), 0, 1)
    evolve <- log(u) < ((objv - objv_t) / temp)

    ## Update accepted states and thier objective values
    state[, evolve] <- state_t[, evolve]
    objv[evolve] <- objv_t[evolve]

    list(state = state, objv = objv, acceptance = mean(evolve))
}

#' SMC-SA
#'
#' Minimise an objective function using SMC-SA with supplied starting states, a
#' proposal distribution and a cooling schedule.
#'
#' The constraints are enforced by using an appropriate proposal
#' distribution. It is advisable that the proposal variance is decreased after
#' each iteration (hence the reason the proposal distribution should take the
#' current iteration count in addition to the current state). The iterator in
#' this code starts from 1 to N.
#'
#' Diagnostic plots include a boxplot of the objective function
#' values and a scatterplot of the acceptance probabilities, both of which are
#' plotted for each iteration. The supplied starting values are recycled to
#' construct \code{N} starting values.
#'
#' @param objf an objective function that takes a matrix with each column
#'     corresponds to a state, and return a vector with respective objective
#'     function values. (num mat -> num vec)
#' @param proposal a random walk proposal that takes a matrix, and produces a
#'     matrix of proposed states, centred on the corresponding columns of the
#'     input matrix. It should also take in the current iteration count. ((num
#'     mat, num) -> num mat)
#' @param starting a matrix with each column corresponding to a starting
#'     state. This code assumes that the starting values are feasible. (num mat)
#' @param schedule a cooling schedule, taking the current iteration count and
#'     the absolute value of the current lowest objective function value as its
#'     first and second arguments respectively, and return a temperature. (num,
#'     num -> num)
#' @param N the number of initial states. If this number is greater than the
#'     column size of \code{starting}, \code{starting} will be recycled. (num)
#' @param iter total number of iterations. (num)
#' @param diagnostic print out diagnostic plots? (bool)
#' @param verbose show the algorithm progression? (bool)
#'
#' @return a list that contains the best state, its objective function
#'     value, and the acceptance rate at each iteration.
#' @export
SMCSA <- function(objf, proposal, starting, schedule, N = 1000, iter = 100,
                  diagnostic = FALSE, verbose = FALSE) {

    ## This is a SMC-SA starting from different states specified in 'starting'.

    ## objf : objective function
    ## proposal : a random walk proposal
    ## starting : starting values matrix, will be padded s.t. ncol(starting) = N
    ## schedule : temperature schedule, current interation count as its argument
    ## RETURN : the global minimum and its objective value

    ## objf : d x N numeric matrix -> N numeric vector
    ## proposal : d x N numeric matrix -> d x N numeric matrix
    ## starting : d x _ numeric matrix
    ## schedule : int -> numeric
    ## RETURN : d numeric vector, numeric

    force(proposal)
    force(objf)
    force(schedule)

    ## Starting values should be feasible, otherwise the proposal should return
    ## infinity for infeasible starting values.
    d <- NROW(starting)
    state_ls <- list(objv = objf(starting))
    minimum <- list(state = starting[which.min(state_ls$objv)],
                    objv = min(state_ls$objv))


    ## Recycle the starting values if it is less than N
    if (NCOL(starting) < N) {
        idx <- rep_len(1:NCOL(starting), N)
        state_ls$state <- starting[, idx, drop = FALSE]
        state_ls$objv <- rep_len(state_ls$objv, N)
    } else {
        if (NCOL(starting) > N) {
            warning("Starting values supplied exceed desired number of MC chain")
        }
        state_ls$state <- starting
        N <- NCOL(state_ls$state)
    }

    ## The previous temperature is set to infinity to make sure that the weight
    ## calculation is correct.
    temp_prev <- Inf

    ## Diagnostic
    acceptance_vec <- rep(0.5, iter)
    track_log_rss <- vector("list", iter)

    for (k in seq.int(1, iter)) {

        ## Get the current temperature
        temp <- schedule(k, abs(minimum$objv))

        ## if (k > 200) {cat(k)}

        ## Importance sampling
        log_w <- -state_ls$objv * (1/temp - 1/temp_prev)
        ## w <- exp(-state_ls$objv * (1/temp - 1/temp_prev))

        index <- sample.int(length(log_w), size = N, replace = TRUE,
                            prob = exp(log_w + min(abs(log_w))))

        ## SA move
        state_ls <- SAmove(state_ls$state[,index, drop = FALSE], objf, proposal,
                         temp, k, state_ls$objv[index])

        ## Look up the current best minimum
        min_col <- which.min(state_ls$objv)
        if (state_ls$objv[min_col] < minimum$objv) {
            minimum$state <- state_ls$state[, min_col]
            minimum$objv <- state_ls$objv[min_col]
            ## if (verbose) {
            ##     cat("At ", k, " iteration, ", minimum$objv, "\n")
            ## }
        }

        ## Current temperature becomes the previous temperature
        temp_prev <- temp

        ## Diagnostic
        acceptance_vec[k] <- state_ls$acceptance
        track_log_rss[[k]] <- log(state_ls$objv)

        if (verbose && (k %% 100 == 0)) {
            cat(k, 'iterations done, current best:', minimum$objv, '\n')
        }
    }

    ## Diagnostic
    if (diagnostic) {
        plot(acceptance_vec ~ seq.int(1, iter), ylim = c(0, 1),
             pch = 20, cex = 0.5)
        ## boxplot(track_log_rss, pch = 20, cex = 0.5)
        ## minimum$track_log_rss <- track_log_rss
    }
    minimum$acc <- acceptance_vec

    if (verbose) {
        cat('Final objective value:', minimum$objv, '\n')
    }
    minimum

}

#' multi SA
#'
#' Minimise an objective function using simulated annealing starting at various
#' locations, with a supplied proposal distribution and cooling schedule.
#'
#' The constraints are enforced by using an appropriate proposal
#' distribution. It is advisable that the proposal variance is decreased after
#' each iteration (hence the reason the proposal distribution should take the
#' current iteration count in addition to the current state). The iterator in
#' this code starts from 1 to N.
#'
#' Diagnostic plots include a boxplot of the objective function values and a
#' scatterplot of the acceptance probabilities, both of which are plotted for
#' each iteration. The supplied starting values are recycled to construct
#' \code{N} starting values.
#'
#' @param objf an objective function that takes a matrix with each column
#'     corresponds to a state, and return a vector with respective objective
#'     function values. (num mat -> num vec)
#' @param proposal a random walk proposal that takes a matrix, and produces a
#'     matrix of proposed states, centred on the corresponding columns of the
#'     input matrix. It should also take in the current iteration count. ((num
#'     mat, num) -> num mat)
#' @param starting a matrix with each column corresponding to a starting
#'     state. This code assumes that the starting values are feasible. (num mat)
#' @param schedule a cooling schedule, taking the current iteration count and
#'     the absolute value of the current lowest objective function value as its
#'     first and second arguments respectively, and return a temperature. (num,
#'     num -> num)
#' @param N the number of initial states. If this number is greater than the
#'     column size of \code{starting}, \code{starting} will be recycled. (num)
#' @param iter total number of iterations. (num)
#' @param diagnostic print out diagnostic plots? (bool)
#' @param verbose show the algorithm progression? (bool)
#'
#' @return a list that contains the best state, its objective function
#'     value, and the acceptance rate at each iteration.
#' @export
multiSA <- function(objf, proposal, starting, schedule, N = 1000, iter = 100,
                    diagnostic = FALSE, verbose = FALSE){
    ## this is a multi-start SA starting from different states specified in
    ## 'starting'.

    ## objf : objective function
    ## proposal : a random walk proposal
    ## starting : starting values matrix, will be padded s.t. ncol(starting) = N
    ## schedule : temperature schedule, current interation count as its argument
    ## RETURN : the global minimum and its objective value

    force(proposal)
    force(objf)
    force(schedule)

    ## Starting values should be feasible, otherwise the proposal should return
    ## infinity for infeasible starting values.
    d <- NROW(starting)
    state_ls <- list(objv = objf(starting))
    minimum <- list(state = starting[which.min(state_ls$objv)],
                    objv = min(state_ls$objv))

    ## Recycle the starting values if it is less than N
    if (NCOL(starting) < N) {
        idx <- rep_len(1:NCOL(starting), N)
        state_ls$state <- starting[, idx, drop = FALSE]
        state_ls$objv <- rep_len(state_ls$objv, N)
    } else if (NCOL(starting) > N) {
        warning("Starting values supplied exceed N")
        state_ls$state <- starting
    }

    ## Diagnostic
    acceptance_vec <- rep(NA, iter)
    track_log_rss <- vector("list", iter)

    for (k in seq.int(1, iter)) {

        ## Get the current temperature
        temp <- schedule(k, abs(minimum$objv))

        ## SA move
        state_ls <- SAmove(state_ls$state, objf, proposal, temp, k, state_ls$objv)

        ## Look up the current best minimum
        min_col <- which.min(state_ls$objv)
        if (state_ls$objv[min_col] < minimum$objv) {
            minimum$state <- state_ls$state[, min_col]
            minimum$objv <- state_ls$objv[min_col]
        }

        ## Diagnostic
        acceptance_vec[k] <- state_ls$acceptance
        track_log_rss[[k]] <- log(state_ls$objv)

        if (verbose && (k %% 100 == 0)) {
            cat(k, 'iterations done, current best:', minimum$objv, '\n')
        }
    }

    ## Diagnostic
    if (diagnostic) {
        plot(acceptance_vec ~ seq.int(1, iter), ylim = c(0, 1),
             pch = 20, cex = 0.5)
        ## boxplot(track_log_rss, pch = 20, cex = 0.5)
        ## minimum$track_log_rss <- track_log_rss
    }
    minimum$acc <- acceptance_vec

    if (verbose) {
        cat('Final objective value:', minimum$objv, '\n')
    }
    minimum
}


