get_partition_fac <- function(wall) {

    idx <- 1:length(wall)
    get_partition_single <- function(x) {
        for (i in idx) {
            if (x <= wall[i]) {
                return(i)
            }
        }
        stop("No partition found.")
    }

    get_partition <- function(objvs) {
        vapply(objvs, get_partition_single, 0L)
    }
    get_partition
}


## This is a Metropolis-Hastings kernel specifically for SAA algorithms.
SAA_MH <- function(state, partition, get_partition, objf, theta, proposal, temp,
                   k, objv = NULL){

    ## WARNING: THIS IS A MINIMISATION ALGORITHM
    ## WARNING: THIS IS A MINIMISATION ALGORITHM
    ## WARNING: THIS IS A MINIMISATION ALGORITHM

    if (is.vector(state)) {
        state <- t(as.matrix(state))
    }

    if (is.null(objv)) {
        objv <- objf(state)
    }

    ## Proposal move
    state_t <- proposal(state, 1)
    ## state_t <- matrix(NA, NROW(state), NCOL(state))
    ## for (i in 1:NCOL(state)) {
    ##     state_t[, i] <- proposal(state[, i], 1000 - partition[i])
    ## }
    objv_t <- objf(state_t)
    partition_t <- get_partition(objv_t)

    ## Acceptance step
    u <- runif(NCOL(state), 0, 1)
    theta_idx <- cbind(partition, 1:NCOL(state))
    theta_t_idx <- cbind(partition_t, 1:NCOL(state))
    evolve <- log(u) < ((objv - objv_t) / temp +
                        (theta[theta_idx] - theta[theta_t_idx]))

    ## Update accepted states and thier objective values
    state[, evolve] <- state_t[, evolve]
    objv[evolve] <- objv_t[evolve]
    partition[evolve] <- partition_t[evolve]

    list(state = state, objv = objv, partition = partition, acceptance = mean(evolve))
}



## Stochastic approximation annealing (SAA) The partitions, schedule, proposal,
## sampling distributions and gain functions are all hard-coded.

SAA <- function(objf, proposal, starting, schedule, N = 1000, iter = 100,
                diagnostic = FALSE, verbose = FALSE) {

    ## This is a SMC-SA starting from different states specified in 'starting'.

    ## objf : objective function
    ## proposal : a random walk proposal
    ## starting : starting values matrix, will be padded s.t. ncol(starting) = N
    ## schedule : temperature schedule, current interation count as its argument
    ## RETURN : the global minimum and its objective value

    ## objf : d x N numeric matrix -> N numeric vector
    ## proposal : d x N numeric matrix -> d x N numeric matrix
    ## starting : d x N numeric matrix
    ## schedule : int -> numeric
    ## RETURN : d numeric vector, numeric

    force(proposal)
    force(objf)
    force(schedule)

    ## suggested temperature schedule for SAA
    schedule <- function(k, defunct) {
        t_h <- 0.5
        T_0 <- 200
        t_star <- 0.01
        t_h * sqrt(T_0 / max(k, T_0)) + t_star
    }

    ## special 1-component proposal
    proposal <- function(state, partition) {
        n_dims <- 2
        dims <- sample(1:NROW(state), n_dims)

        stdnorm <- matrix(rnorm(n_dims * NCOL(state)), n_dims)
        pertub <- stdnorm * 2
        state[dims, ] <- state[dims, ] + pertub
        state
    }

    ## set partitions and thetas
    wall <- c(seq(0, min(objf(starting)), len = 1000), Inf) # upper bounds of
                                                          # partitions
    get_partition <- get_partition_fac(wall)

    ## Starting values should be feasible, otherwise the proposal should return
    ## infinity for infeasible starting values.
    d <- NROW(starting)
    state_ls <- list(state = starting,
                     objv = objf(starting))
    state_ls$partition <- get_partition(state_ls$objv)
    minimum <- list(state = starting[which.min(state_ls$objv)],
                    objv = min(state_ls$objv))

    ## Recycle the starting values if it is less than N
    if (NCOL(starting) < N) {
        idx <- rep_len(1:NCOL(starting), N)
        state_ls$state <- starting[, idx, drop = FALSE]
        state_ls$objv <- rep_len(state_ls$objv, N)
        state_ls$partition <- rep_len(state_ls$partition, N)
    } else if (NCOL(starting) > N) {
        warning("Starting values exceed desired number of MC chains")
        N <- NCOL(state_ls$state)
    }

    ## initialise theta
    theta <- matrix(0, length(wall), N)

    ## sampling distribution, pi
    tune <- 0.1                         # the tuning parameter for gain, zeta
    samp_dist <- exp(-tune * 0:(length(wall) - 1)) /
        sum(exp(-tune * 0:(length(wall) - 1)))
    samp_dist <- matrix(samp_dist, length(samp_dist), N)

    ## gain sequence
    gain <- function(k) {
        sig <- 1.0
        T_0 <- 2000
        (T_0 / max(k, T_0))^sig
    }

    ## Diagnostic
    acceptance_vec <- rep(NA, iter)
    track_rss <- vector("list", iter)

    for (k in seq.int(1, iter)) {

        ## Get the current temperature
        temp <- schedule(k, abs(minimum$objv))

        ## SAA MH update
        state_ls <- SAA_MH(state_ls$state, state_ls$partition, get_partition,
                           objf, theta, proposal, temp, k, state_ls$objv)

        if (k %% 10000 == 0) {
            ## browser()
        }

        ## theta-updating (assuming no truncation)
        e_idx <- cbind(state_ls$partition, 1:N)
        H <- -1 * samp_dist
        H[e_idx] <- H[e_idx] + 1
        theta <- theta + gain(k) * H


        ## Look up the current best minimum
        min_col <- which.min(state_ls$objv)
        if (state_ls$objv[min_col] < minimum$objv) {
            minimum$state <- state_ls$state[, min_col]
            minimum$objv <- state_ls$objv[min_col]
        }

        ## Diagnostic
        acceptance_vec[k] <- state_ls$acceptance
        track_rss[[k]] <- state_ls$objv

        if (verbose && (k %% 100000 == 0)) {
            cat(k, 'iterations done, current best:', minimum$objv, ' ', state_ls$partition, '\n')
        }
    }

    ## Diagnostic
    if (diagnostic) {
        plot(acceptance_vec ~ seq.int(1, iter), ylim = c(0, 1),
             pch = 20, cex = 0.5)
        ## boxplot(track_rss, pch = 20, cex = 0.5)
        minimum$track_rss <- track_rss
    }
    minimum$acc <- acceptance_vec

    if (verbose) {
        cat('Final objective value:', minimum$objv, '\n')
    }
    minimum

}
