% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CEPSO.R
\name{exchange}
\alias{exchange}
\title{Exchange between two swarm.}
\usage{
exchange(sw1, sw2, k)
}
\arguments{
\item{sw}{the swarm produced by 'initSwarm'. (list)}

\item{delta}{the closeness function that takes a matrix and return a
vector. (num mat -> num vec)}

\item{indicator}{an indicator function that takes a matrix, and returns a
boolean vector that indicates the feasibility of each column of the
matrix. (num mat -> bool vec)}
}
\value{
a list with various properties of the swarm.
}
\description{
Exchange the particle between the exploitation and exploration swarms.
}
