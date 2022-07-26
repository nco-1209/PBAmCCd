% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variation_functions.R
\name{g}
\alias{g}
\title{Log of the Inverse Additive Log-ratio}
\usage{
g(x)
}
\arguments{
\item{x}{Compositional data vector which has already been transformed by the 
additive log-ratio}
}
\value{
A vector which is the log of the inverse of the data which has been
transformed by the additive logratio transformation.
}
\description{
Calculates the log of a vector after applying the
inverse additive log-ratio transformation to it. The transformation is 
\code{g(x)=log(alrInx(x))}, where \code{x} is a vector transformed by the 
additive log-ratio.
}
