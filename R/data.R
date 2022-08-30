#' dat.ss
#'
#' The single-cell mRNA dataset of 96 samples (rows), and 100 genes (columns). 
#'
#' @docType data
#'
#' @usage dat.ss
#'
#' @format single-cell mRNA table. 
#'
#' @keywords datasets
#'
#' @examples
#' dat.ss
#' dim(dat.ss)
#' 
#' 
"dat.ss"

#' MLE
#'
#' Maximum likelihood estimates from dat.ss, such as the \code{mean} and 
#' \code{Sigma}. MLE based on the multinomial logit-Normal model.
#'
#' @docType data
#'
#' @usage mle
#'
#' @format Gives mle estimates of \code{mu} and \code{Sigma} of \code{dat.ss}.
#'
#' @keywords datasets
#'
#' @examples
#' mu <- mle$est.min$mu
#' Sigma <- mle$est.min$Sigma
#' 
#' 
"mle"