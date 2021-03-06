% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rational-function.R
\name{rational_biweight_fac}
\alias{rational_biweight_fac}
\title{Tukey's biweight of rational function models}
\usage{
rational_biweight_fac(y, x, dim_f, dim_g, cst = 4.685)
}
\arguments{
\item{y}{a response vector. (num vec)}

\item{x}{a predictor vector. (num vec)}

\item{dim_f}{the number of coefficients in the numerator. (num)}

\item{dim_g}{the number of coefficients in the denominator, excluding the
constant term. (num)}

\item{cst}{the constant in Tukey's biweight. (num)}
}
\value{
a function that takes in a matrix of rational function coefficients
    (numerator precede denominator, in ascending order), each column
    corresponds to a rational function, and returns a vector Tukey's
    biweight.
}
\description{
Produce a function that calculates the Tukey's biweight of a set of rational
function models.
}
\details{
The explicit formula used is rho(x) = (c^2)/6 * (1-(1-(x/c)^2)^3) if x < |c|,
or rho(x) = 1 if x > |c|, where c is an arbitrary constant. The sum of
rho(y_i - r(x_i)), where r(.) is the rational function, is then
minimised. The produced function can be minimised using SMC-SA, multiSA or
CEPSO.
}
