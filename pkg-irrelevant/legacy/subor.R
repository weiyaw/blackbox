# Dependence: 'oracle.R'
source('oracle.R')

ratRSS <- function(y, x, f, g, value = TRUE) {

  # Calculate residuals sum of square (RSS) of fitted rational function.
  #
  # y: response data
  # x: predictor data
  # f: numerator coefficients (ascending)
  # g: denominator coefficients (ascending)
  # value: if TRUE, return RSS (scalar). Return resid vector otherwise.
  #
  # return:
  # RSS, if 'value' is TRUE
  # residuals vector, if 'value' is FALSE

  fit <- evalPol(x, f) / evalPol(x, g)
  if (value) {
    return(crossprod(y - fit))
  } else {
    return(y - fit)
  }
}

ratRSSDeriv <- function(y, x, f, g) {

  # Return derivative w.r.t. parameters.
  #
  # y: response data
  # x: predictor data
  # f: numerator coefficients (ascending)
  # g: denominator coefficients (ascending)
  #
  # return:
  # da: derivative w.r.t. (alpha_0, alpha_1 ..., alpha_p)
  # db: derivative w.r.t. (beta_0, beta_1 ..., beta_q)

  num <- evalPol(x, f)
  denom <- evalPol(x, g)
  RS <- (num / denom) - y
  alpha.term <- RS / denom
  beta.term <- (RS * num) / (denom^2)
  da <- rep(NA, length(f))
  db <- rep(NA, length(g))

  # alpha part
  for (k in 0:(length(f) - 1)) {
    da[k + 1] <- 2 * sum(alpha.term * (x^k))
  }

  # beta part
  for (k in 0:(length(g) - 1)) {
    db[k + 1] <- -2 * sum(beta.term * (x^k))
  }

  return(list(da = da, db = db))
}

weight <- function(R, k, N, celsius) {

  # Calculate the weights of every elements in R according to their objective
  # value and current temperature, and sample N elements from H with
  # replacement.
  #
  # R: objective values vector
  # k: current iteration number
  # N: total number of next input
  # celsius: (T_previous, T)
  #
  # return:
  # a vector of size N with index number of H vector.


  if (k == 1) {
    w <- exp(-R / celsius)
  } else {
    w <- exp(-R * (1 / celsius[2] - 1 / celsius[1]))
  }
  return(sample.int(length(R), N, replace = TRUE, prob = w))
}

acceptReject <- function(R, R.next, celsius) {

  # TRUE: accept H.next
  # FALSE: reject H.next

  rho <- min(exp((R - R.next) / celsius), 1)
  return(runif(1) <= rho)
}

cooling <- function(k, celsius_0, celsius) {
  # celsius.next <- celsius_0 * 0.9^k    # 0.8 < alpha < 0.9
  # celsius.next <- celsius_0 / (1 + 1.5 * log(1 + k))    # alpha > 1
  # celsius.next <- celsius_0 / (1 + 0.85 * k)    # alpha > 0
  celsius.next <- celsius_0 / (1 + 0.85 * k^2)    # alpha > 0
  return(c(celsius.next, celsius))
}

genMonoRat <- function(y, x, f, g, f.chg, g.chg, a, b,increasing = T, EPS = 1e-6, dist = rcauchy, dist.p = list(scale = 2)) {

  # Generate a neigbouring monotone rational function given a non-monotone
  # rational function.
  # Only coefficients that are marked TRUE will be changed.

  f.temp <- f
  g.temp <- g
  dist.f <- dist.p
  dist.g <- dist.p
  dist.f$n <- sum(f.chg)
  dist.g$n <- sum(g.chg)

  for (k in 1:10000) {
    f.temp[f.chg] <- f[f.chg] + do.call(dist, dist.f)
    g.temp[g.chg] <- g[g.chg] + do.call(dist, dist.g)
    # k <- k + 1
    if (oracle(f.temp, g.temp, a, b, increasing, EPS)) {
      R <- ratRSS(y, x, f.temp, g.temp)
      return(list(R = R, f = f.temp, g = g.temp))
    }
  }
  stop('No monotone function found. Check initial value.')
}

genMonoRatOnly <- function(f, g, f.chg, g.chg, a, b,increasing = T, EPS = 1e-6, dist = rcauchy, dist.p = list(scale = 2)) {

  # Generate a neigbouring monotone rational function given a non-monotone
  # rational function.

  # k <- 0
  f.temp <- f
  g.temp <- g
  dist.f <- dist.p
  dist.g <- dist.p
  dist.f$n <- sum(f.chg)
  dist.g$n <- sum(g.chg)

  for (k in 1:10000) {
    f.temp[f.chg] <- f[f.chg] + do.call(dist, dist.f)
    g.temp[g.chg] <- g[g.chg] + do.call(dist, dist.g)
    # k <- k + 1
    if (oracle(f.temp, g.temp, a, b, increasing, EPS)) {
      return(list(f = f.temp, g = g.temp))
    }
  }
  stop('No monotone function found. Please give a new initial value.')
}


plotRat <- function(f, g, xlim = c(0, 10), sketch = T, ...) {

  # Plot rational function

  for (i in seq_len(NCOL(f))) {
    x <- seq(xlim[1], xlim[2], length.out = 1000)
    y <- evalPol(x, f) / evalPol(x, g)
    if (NCOL(f) > 1) {
      lines(y ~ x, col = i+1, ...)
      next
    } else {
      if(sketch) {
        lines(y ~ x, ...)
      } else {
        plot(y ~ x, type = 'l', ...)
      }
    }
  }
}

fiveMonoRat <- function(y, x, f, g, iter, N, schedule) {
  result <- list()
  rss <- vector()
  for (i in 1:5) {
    cat('Run ', i, '\n')
    # result[[i]] <- monoRat(y, x, f, g, limits = c(0, 60), iter = iter, N = N, schedule = schedule)
    result[[i]] <- monoRat(y, x, f, g, c(F,T,T), limits = c(0, 2), iter = iter, N = N, schedule = schedule)
    rss[i] <- result[[i]]$R
  }
  result[[6]] <- list(average = mean(rss), diff = diff(range(rss)), min = min(rss), which = which.min(rss))
  cat('Average = ', result[[6]]$average, '\n')
  cat('Difference = ', result[[6]]$diff, '\n')
  cat('Minimum = ', result[[6]]$min, '\n')
  return(result)
}

fiveMultiSA <- function(y, x, f, g, iter, N, schedule) {
  result <- list()
  rss <- vector()
  for (i in 1:5) {
    cat('Run ', i, '\n')
    # result[[i]] <- multiSA(y, x, f, g, limits = c(0, 60), iter = iter, N = N, schedule = schedule)
    result[[i]] <- multiSA(y, x, f, g, c(F,T,T), limits = c(0, 2), iter = iter, N = N, schedule = schedule)
    rss[i] <- result[[i]]$R
  }
  result[[6]] <- list(average = mean(rss), diff = diff(range(rss)), min = min(rss), which = which.min(rss))
  cat('Average = ', result[[6]]$average, '\n')
  cat('Difference = ', result[[6]]$diff, '\n')
  cat('Minimum = ', result[[6]]$min, '\n')
  return(result)
}
