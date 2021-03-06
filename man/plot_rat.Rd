% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rational-function.R
\name{plot_rat}
\alias{plot_rat}
\title{Line plots of rational functions}
\usage{
plot_rat(fs, gs, xlim = c(0, 10), sketch = TRUE, fine = 1000, ...)
}
\arguments{
\item{fs}{a matrix of numerator polynomial in ascending order representing,
each column corresponds to a polynomial. (num mat)}

\item{gs}{a matrix of denominator coefficients in ascending order, each
column corresponds to a polynomial. (num mat)}

\item{xlim}{the x-axis range to plot}

\item{sketch}{if set to \code{TRUE}, this function add a rational function
onto an existing plot, rather than creating a new graphical device.}

\item{fine}{how many points between xlims?}

\item{...}{optional graphical parameters to be passed to plot function.}
}
\description{
Plot one or multiple rational functions, given a set of coefficients.
}
