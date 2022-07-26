% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variation_functions.R
\name{logVarMC}
\alias{logVarMC}
\title{Monte Carlo Logx Variance-Covariance}
\usage{
logVarMC(mu, Sigma, K = 1e+05)
}
\arguments{
\item{mu}{The mean vector of a dataset following the multinomial logit-Normal model}

\item{Sigma}{The variance matrix of a dataset following the multinomial logit-Normal model}

\item{K}{Number of samples}
}
\value{
The estimated variance-covariance matrix, \code{logx}
}
\description{
Estimates the variance-covariance of the log of the data, using Monte
Carlo integration.
}
\examples{
mu <- mle$est.min$mu
Sigma <- mle$est.min$Sigma
logVarMC(mu, Sigma)

}
