% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mle_multinom_logit.R
\name{ebicPlot}
\alias{ebicPlot}
\title{Extended Bayesian Information Criterion Plot}
\usage{
ebicPlot(fit, xlog = FALSE)
}
\arguments{
\item{fit}{The model fit}

\item{xlog}{TRUE or FALSE. Renders plot with the x-axis in the log-scale if \code{TRUE}}
}
\value{
Plot of the EBIC (y-axis) against each lambda (x-axis).
}
\description{
Plots the extended Bayesian information criterion (EBIC) of the model fit for
various penalization parameters \code{lambda}.
}
\examples{
mle$ebic #y-axis of ebic plot
log(mle$lambda.gl) #x-axis of ebic plot

ebicPlot(mle, xlog = TRUE) 

}
