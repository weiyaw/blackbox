source('oracle.R')
source('subor.R')

monoRat <- function(y, x, f, g, f.chg = NA, g.chg = NA, limits = NA, 
                    increasing = TRUE, dist = rcauchy, dist.p = list(scale = 2), 
                    SA.step = 1, iter = 100, N = 1000, diagnt = FALSE,
                    EPS = 1e-06, schedule = 'quad', ...) {
  
  # Fit a least square rational function that is pole-less and monotone across
  # region [a,b].
  # 
  # y: response data
  # x: predictor data
  # f: numerator coefficients on the same scale as data (ascending)
  # g: denominator coefficients on the same scale as data (ascending)
  # f.chg, g.chg: specify which coefficients to be optimised
  # N: number of starting values, artifically generated around 'f' and 'g'
  # limit: lower and upper bounds of region of interest
  # inreasing: can be monotone increasing or decreasing
  # SA.only: skip SMC step and perform N times SA simutaneously
  # dist: a symmetrical distribution for generating starting values
  # dist.p: parameters of 'dist' distribution
  # SA.step: step size of SA move. (i.e. 'sd' of normal distribution)
  # iter: number of iteration
  # 
  # 
  # 
  # Only 'SA.only' if input rational function is monotone.
  
  if (length(f) == length(g)) {
    if (all(f == g)) {
      stop('Do not input identical numerator and denominator')
    }
  }
  
  if (is.data.frame(y) | is.data.frame(x)) {
    stop('Transform x and y into vectors or matrices.')
  }  

  ### Default settings below ###
  if (is.na(limits[1])) {
    limits <- range(x)
  }
  a <- limits[1]
  b <- limits[2]
  
  if (is.na(f.chg[1])) {
    f.chg <- rep(TRUE, length(f))
  }
  if (is.na(g.chg[1])) {
    g.chg <- rep(TRUE, length(g))
    g.chg[1] <- FALSE
  }
  ### Default settings above ###
  
  para <- optimal <- list()
  para$R <- rep(NA, length = N)
  para$f <- matrix(NA, nrow = length(f), ncol = N)
  para$g <- matrix(NA, nrow = length(g), ncol = N)
  
  for (i in seq_len(N)) {
    temp <- genMonoRat(y, x, f, g, f.chg, g.chg, a, b, increasing, EPS, 
                       dist = dist, dist.p = dist.p)
    para$R[i] <- temp$R
    para$f[, i] <- temp$f
    para$g[, i] <- temp$g
  }
  
  if (oracle(f, g, a, b, increasing, EPS)) {
    para$R[1] <- ratRSS(y, x, f, g)
    para$f[, 1] <- f
    para$g[, 1] <- g
  }
  
  min.col <- which.min(para$R)
  optimal$R <- para$R[min.col]
  optimal$f <- para$f[, min.col]
  optimal$g <- para$g[, min.col]
  
  # cat(optimal$R, '\n')
  # cat('tag \n')
  celsius_0 <- abs(optimal$R) / log(2)    # Cooling schedule T_0
  celsius <- c(celsius_0, NA)
  
  acc.rate <- rep(NA, iter) # Diagnostic
  track.rss <- vector("list", iter)
  
  for (k in seq_len(iter)) {
    if (k == 1) {
      index <- weight(para$R, k, N, celsius[1])
    } else {
      index <- weight(para$R, k, N, rev(celsius))
    }
    
    para <- SA(y, x, para$f[, index], para$g[, index], f.chg, g.chg, 
               para$R[index], celsius[1], SA.step*(0.95^k), a, b, increasing,
               EPS)
    
    min.col <- which.min(para$R)
    acc.rate[k] <- para$acc.rate.indv  # Diagnostic
    track.rss[[k]] <- para$R
    
    if (para$R[min.col] < optimal$R) {
      optimal$R <- para$R[min.col]
      optimal$f <- para$f[, min.col]
      optimal$g <- para$g[, min.col]
      # para.deriv <- ratRSSDeriv(y, x, optimal$f, optimal$g)
      # cat(optimal$R, '\n')
      # cat('da = ', para.deriv$da, '\n')
      # cat('db = ', para.deriv$db, '\n \n')
    }
    
    # Cooling schedule
    if (schedule == 'quad') {
      celsius <- cooling(k, celsius_0, celsius[1])
    } else if (schedule == 'log') {
      celsius <- c(abs(optimal$R) / log(k + 2), celsius[1])
    }
    else {
      stop('Please choose a cooling schedule \'log\' or \'quad\'.')
    }
    
    if (k %% 10 == 0) {
      N <- N + 50
      cat(k, 'iterations done, RSS:', optimal$R, '\n')
    }
  }
  if (diagnt) {
    plot(acc.rate, ylim = c(0, 1))  # Diagnostic
    boxplot(track.rss, ...)
  }
  return(optimal)
}

SAonly <- function(y, x, f, g, f.chg = NA, g.chg = NA, limits = NA, 
                   increasing = TRUE, dist = rcauchy, dist.p = list(scale = 2), 
                   SA.step = 1, iter = 50, N = 1000, diagnt = FALSE,
                   EPS = 1e-06, ...) {
  # This is a multi-start SA starting from same location.
  
  
  if (length(f) == length(g)) {
    if (all(f == g)) {
      stop('Do not input identical numerator and denominator')
    }
  }
  
  if (is.data.frame(y) | is.data.frame(x)) {
    stop('Transform x and y into vectors or matrices.')
  }  
  
  ### Default settings below ###
  if (is.na(limits[1])) {
    limits <- range(x)
  }
  a <- limits[1]
  b <- limits[2]
  
  if (is.na(f.chg[1])) {
    f.chg <- rep(TRUE, length(f))
  }
  if (is.na(g.chg[1])) {
    g.chg <- rep(TRUE, length(g))
    g.chg[1] <- FALSE
  }
  ### Default settings above ###
  
  if (!oracle(f, g, a, b, increasing, EPS)) {
    stop('SA.only only takes monotone rational function')
  }

  para <- list()
  para$f <- matrix(rep(f, times = N), ncol = N)
  para$g <- matrix(rep(g, times = N), ncol = N)
  para$R <- rep(ratRSS(y, x, f, g), length = N)
  optimal <- list(f = f, g = g, R = para$R[1])

  cat('tag \n')
  celsius_0 <- abs(optimal$R) / log(2)    # Cooling schedule T_0
  celsius <- c(celsius_0, NA)
  
  acc.rate <- rep(NA, iter) # Diagnostic
  track.rss <- vector("list", iter)
  
  for (k in seq_len(iter)) {
    para <- SA(y, x, para$f, para$g, f.chg, g.chg, para$R, celsius[1], 
               SA.step*(0.95^k), a, b, increasing, EPS)
    min.col <- which.min(para$R)
    acc.rate[k] <- para$acc.rate.indv  # Diagnostic
    track.rss[[k]] <- para$R
    
    if (para$R[min.col] < optimal$R) {
      optimal$R <- para$R[min.col]
      optimal$f <- para$f[, min.col]
      optimal$g <- para$g[, min.col]
      para.deriv <- ratRSSDeriv(y, x, optimal$f, optimal$g)
      cat(optimal$R, '\n')
      cat('da = ', para.deriv$da, '\n')
      cat('db = ', para.deriv$db, '\n \n')
    }
    # celsius <- cooling(k, celsius_0, celsius[1])
    celsius <- c(abs(optimal$R) / log(k + 2), celsius[1])    # Cooling schedule
    if (k %% 10 == 0) {
      cat(k, 'iterations done \n')
    }
  }
  if (diagnt) {
    plot(acc.rate, ylim = c(0, 1))  # Diagnostic
    boxplot(track.rss, ...)
  }
  return(optimal)
}

multiSA <- function(y, x, f, g, f.chg = NA, g.chg = NA, limits = NA, 
                   increasing = TRUE, dist = rcauchy, dist.p = list(scale = 2), 
                   SA.step = 1, iter = 100, N = 1000, diagnt = FALSE,
                   EPS = 1e-06, schedule = 'quad', ...) {
  # This is a multi-start SA starting from different states.
  
  
  if (length(f) == length(g)) {
    if (all(f == g)) {
      stop('Do not input identical numerator and denominator')
    }
  }
  
  if (is.data.frame(y) | is.data.frame(x)) {
    stop('Transform x and y into vectors or matrices.')
  }  
  
  ### Default settings below ###
  if (is.na(limits[1])) {
    limits <- range(x)
  }
  a <- limits[1]
  b <- limits[2]
  
  if (is.na(f.chg[1])) {
    f.chg <- rep(TRUE, length(f))
  }
  if (is.na(g.chg[1])) {
    g.chg <- rep(TRUE, length(g))
    g.chg[1] <- FALSE
  }
  ### Default settings above ###
  
  para <- optimal <- list()
  para$f <- matrix(NA, nrow = length(f), ncol = N)
  para$g <- matrix(NA, nrow = length(g), ncol = N)
  
  for (i in seq_len(N)) {
    temp <- genMonoRat(y, x, f, g, f.chg, g.chg, a, b, increasing, EPS, 
                       dist = dist, dist.p = dist.p)
    para$R[i] <- temp$R
    para$f[, i] <- temp$f
    para$g[, i] <- temp$g
  }
  
  if (oracle(f, g, a, b, increasing, EPS)) {
    para$R[1] <- ratRSS(y, x, f, g)
    para$f[, 1] <- f
    para$g[, 1] <- g
  }
  
  min.col <- which.min(para$R)
  optimal$R <- para$R[min.col]
  optimal$f <- para$f[, min.col]
  optimal$g <- para$g[, min.col]
  
  # cat('tag \n')
  celsius_0 <- abs(optimal$R) / log(2)    # Cooling schedule T_0
  celsius <- c(celsius_0, NA)
  
  acc.rate <- rep(NA, iter) # Diagnostic
  track.rss <- vector("list", iter)
  
  for (k in seq_len(iter)) {
    para <- SA(y, x, para$f, para$g, f.chg, g.chg, para$R, celsius[1], 
               SA.step*(0.95^k), a, b, increasing, EPS)
    min.col <- which.min(para$R)
    acc.rate[k] <- para$acc.rate.indv  # Diagnostic
    track.rss[[k]] <- para$R
    
    if (para$R[min.col] < optimal$R) {
      optimal$R <- para$R[min.col]
      optimal$f <- para$f[, min.col]
      optimal$g <- para$g[, min.col]
      para.deriv <- ratRSSDeriv(y, x, optimal$f, optimal$g)
      # cat(optimal$R, '\n')
      # cat('da = ', para.deriv$da, '\n')
      # cat('db = ', para.deriv$db, '\n \n')
    }
    
    # Cooling schedule
    if (schedule == 'quad') {
      celsius <- cooling(k, celsius_0, celsius[1])
    } else if (schedule == 'log') {
      celsius <- c(abs(optimal$R) / log(k + 2), celsius[1])
    }
    else {
      stop('Please choose a cooling schedule \'log\' or \'quad\'.')
    }
    
    if (k %% 10 == 0) {
      cat(k, 'iterations done, RSS:', optimal$R, '\n')
    }
  }
  if (diagnt) {
    plot(acc.rate, ylim = c(0, 1))  # Diagnostic
    boxplot(track.rss, ...)
  }
  return(optimal)
}

SA <- function(y, x, f, g, f.chg, g.chg, R, celsius, sd = 1, a = 0, b = Inf, 
               increasing = TRUE, EPS = 1e-06) {
  
  # WARNING: THIS IS A MINIMISATION ALGORITHM
  # WARNING: THIS IS A MINIMISATION ALGORITHM
  # WARNING: THIS IS A MINIMISATION ALGORITHM
  #
  # Minimise R = RSS to get least square curve. This is the SA move.
  # 
  # x, y: dataset
  # f, g: parameters matries, ncol = N
  # R   : RSS of corresponding column
  # k   : iteration count for cooling schedule
  if (is.vector(f)) {
    f <- t(as.matrix(f))
  }
  if (is.vector(g)) {
    g <- t(as.matrix(g))
  }
  
  c <- 0
  f.temp <- f
  g.temp <- g
  f.size <- sum(f.chg)
  g.size <- sum(g.chg)
  chg.size <- f.size + g.size
  bad <- rep(TRUE, length(R))
  acc.count <- 0  # Diagnostic
  
  while(any(bad)) {
    bad.size <- sum(bad)
    norm.rv <- matrix(rnorm(chg.size * bad.size, sd = sd), nrow = chg.size)
    f.temp[f.chg, bad] <- f[f.chg, bad] + head(norm.rv, f.size, addrownums = F)
    g.temp[g.chg, bad] <- g[g.chg, bad] + tail(norm.rv, g.size, addrownums = F)
    for (i in which(bad)) {
      if (!oracle(f.temp[, i], g.temp[, i], a, b, increasing, EPS)) {
        next
      }
      bad[i] <- FALSE
      R.temp <- ratRSS(y, x, f.temp[, i], g.temp[, i])
      if (acceptReject(R[i], R.temp, celsius)){
        f[, i] <- f.temp[, i]
        g[, i] <- g.temp[, i]
        R[i] <- R.temp
        acc.count <- acc.count + 1  # Diagnostic
      }
    }
  }
  
  # for (i in 1:N) {
  #   f.temp <- f[, i]
  #   g.temp <- g[, i]
  #   repeat {    # Generate a valid move
  #     norm.rv <- rnorm(f.size + g.size, 0, sd)
  #     f.temp[f.chg] <- f[f.chg, i] + norm.rv[1:f.size]  # NOISE DISTRIBUTION
  #     g.temp[g.chg] <- g[g.chg, i] + norm.rv[-(1:f.size)]  # NOISE DISTRIBUTION
  #     # c <- c + 1
  #     if (oracle(f.temp, g.temp, a, b, increasing, EPS)) {   # ORACLE HERE
  #       if (c > 20){
  #         cat(c, ' times \n')
  #         c <- 0
  #       }
  #       break
  #     }
  #   }
  #   
  #   H.temp <- -ratRSS(y, x, f.temp, g.temp)
  #   # cat('RSS ', -H[i], ' RSS.temp', -H.temp, '\n')
  #   
  #   if (acceptReject(H[i], H.temp, celsius)){
  #     f[, i] <- f.temp
  #     g[, i] <- g.temp
  #     H[i] <- H.temp
  #     acc.count <- acc.count + 1  # Diagnostic
  #   }
  # }
  return(list(R = R, f = f, g = g, acc.rate.indv = (acc.count / length(R))))
}