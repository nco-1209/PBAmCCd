#' MLE Data
#'
#' Dataset of maximum likelihood estimated data used for simulating single-cell mRNA data.
#'
#' @docType data
#'
#' @usage data(mle_results_p100)
#'
#' @format Gives mle estimates such as \code{mu} and \code{Sigma} (\code{mle}), as well as a single-cell
#' dataset of 96 samples and 100 genes (\code{dat.ss}).  
#'
#' @keywords datasets
#'
#' @examples
#' data(mle_results_p100)
#' mu <- mle$est.min$mu
#' Sigma <- mle$est.min$Sigma
#' dat.ss
#' 
#' 
"mle_results_p100"

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