% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/var_sim.R
\name{var_sim}
\alias{var_sim}
\title{Simulate a stationary Vector Autoregressive (VAR) time series.}
\usage{
var_sim(T, A, Sigma)
}
\arguments{
\item{T}{An integer giving the number of timepoints.}

\item{A}{A d x d transition matrix.}

\item{Sigma}{A d x d innovation covariance matrix.}
}
\description{
Simulate a stationary Vector Autoregressive (VAR) time series.
}
\examples{

theta    <- diag(c(.7,.8,.9,.6,.7,.9))
data     <- t(var_sim(100, theta, diag(.1,6)))
datalag  <- embed(data, 2)
b        <- datalag[,1:6]
A        <- datalag[,7:12]
A_est    <- fista_sparse(A, b, 1, theta, niter = 10, backtrack = TRUE)$out.x
var_forecast(t(b), 2, A_est)

}
\keyword{simulate}
\keyword{var}
