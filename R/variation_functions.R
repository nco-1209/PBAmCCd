#' Estimated Variation Matrix
#' 
#' Estimates the variation matrix of compositional count data.
#'
#' @param counts count compositional dataset
#' @param p.model Probability model for the counts: (\code{logitNormal}, \code{dirichlet},
#' or \code{plugin}) 
#' @param type Type of variation metric to be calculated: \code{standard}, 
#' \code{phi}, or \code{rho}
#' @param refColumn The reference column to be used for the counts matrix
#'
#' @return The estimated variation matrix for the counts. May be a \code{standard},
#' \code{phi}, or \code{rho} variation matrix depending on \code{type} specified. 
#' 
#' @examples 
#' 
#' @export
#'
varEst <- function(counts, p.model=c("logitNormal", "dirichlet", "plugin"), type=c("standard","phi","rho"), refColumn=NULL) {
  if (any(counts<0) | any(counts!=floor(counts)) | !is.matrix(counts)) stop("'counts' must be a matrix containing non-negative integers")
  p.model <- match.arg(p.model)
  type <- match.arg(type)
  
  if (is.null(refColumn)&p.model=="logitNormal") { 
    refColumn <- J
  }
  
  if (p.model=="logitNormal") {
    result <- logitNormalVariation(mu, Sigma, type=type)
  } else if (p.model=="dirichlet") {
    result <- dirichletVariation(counts)
  } else {
    result <- pluginVariation(counts)
  }
  
  return(result)
}


#' Logit Normal Variation
#' 
#' Estimates the empirical metrics of association of count compositional data 
#' which follows a multinomial logit-Normal distribution. Estimate is performed using 
#' only the parameters of the distribution. 
#'
#' @param mu The mle estimate of the mu matrix 
#' @param Sigma The mle estimate of the sigma matrix
#' @param type Type of variation metric to be calculated: \code{standard}, \code{phi},
#' \code{phis} (a symmetrical version of \code{phi}), or \code{rho}
#' @param order The order of the Taylor-series approximation to be used in the 
#' estimation
#'
#' @return An estimation of the variation matrix, \code{V}.
#'
#' @examples
#' 
#' @export
logitNormalVariation <- function(mu, Sigma, type=c("standard","phi", "phis","rho"), 
                                 order=c("first", "second")) {
  type <- match.arg(type)
  J <- length(mu)
  
  ones <- rep(1, J)
  d.S <- diag(Sigma)
  V <- tcrossprod(d.S, ones) + tcrossprod(ones, d.S) - 2*Sigma
  
  if (type=="phi") {
    lv <- logVarTaylor(mu, Sigma, order=order)
    lv.row <- tcrossprod(diag(lv)[-(J+1)], ones)
    
    V <- V/lv.row
  } else if (type=="phis") {
    lv <- logVarTaylor(mu, Sigma, order=order)
    lv.d <- diag(lv)
    den <- outer(lv.d, lv.d, "+") + 2*lv
    V <- V/den
  } else if (type=="rho") {
    lv <- logVarTaylor(mu, Sigma, order=order)
    lv.d <- diag(lv)
    den <- outer(lv.d, lv.d, "+")
    V <- 2*lv/den
  }
  
  return(V)
}

#' Naive Variation
#' 
#' Estimates the variation matrix of counte compositional data using strictly the 
#' dataset.
#'
#' @param counts Dataset of compositional data
#' @param pseudo.count Scaler value added to the data matrix to prevent infinite 
#' values caused by taking the log of the counts
#' @param type Type of variation metric to be calculated: \code{standard}, \code{phi},
#'  \code{phis} (a symmetrical version of \code{phi}), \code{rho}, or \code{logx}
#' @param use If equal to \code{"everything"} and there are no infinite values after
#' taking the log the calculation will use all data. If equal to \code{"everything"} 
#' and there are infinite values after taking the log it is recommended to  run 
#' the function again, instead setting the \code{use} parameter to \code{"pairwise.complete.obs"}
#' @param set.inf.na If \code{TRUE}, sets any infinite values in \code{counts} to \code{NA}
#' @param already.log If \code{FALSE}, the counts have not been transformed by 
#' by the log. This transformation is of the form is \eqn{log(frac{X_{ij}}{s_{i}})}, where 
#' \eqn{s_{i} = \sum{n=1}^{j} X_{in}}, where \eqn{X_{ij}} is element \eqn{counts[i,j]}
#'
#' @return The naive variation matrix, \code{v}.
#' 
#' @examples
#' 
#' @export
#'
naiveVariation <- function(counts, pseudo.count=0, type=c("standard","phi", "phis","rho", "logx"), 
                           use="everything",
                           set.inf.na=TRUE, already.log=FALSE) {
  type <- match.arg(type)
  J <- NCOL(counts)
  l <- counts
  if (!already.log) {
    l <- l + pseudo.count
    l <- l/rowSums(l)
    l <- log(l)
    get.inf <- is.infinite(l)
    l[get.inf] <- NA
    if (any(get.inf) & use=="everything") {
      warning("There are infinities after taking log.  Consider setting paramter use='pairwise.complete.obs'")
    }
  }
  
  v <- matrix(0,J,J)
  for (i in 1:J) {
    for (j in 1:J){
      if (type=="standard") {
        v[i,j] <- var(l[,i]-l[,j], use=use)
      } else if (type=="phi") {
        v[i,j] <- var(l[,i]-l[,j], use=use)/var(l[,i], use=use)
      } else if (type=="phis") {
        v[i,j] <- var(l[,i]-l[,j], use=use)/(var(l[,i]+l[,j], use=use))
      } else if (type=="rho") {
        v[i,j] <- 2*cov(l[,i],l[,j], use=use)/(var(l[,i], use=use)+var(l[,j], use=use))
      }
    }
  }
  
  if (type=="logx") v <- cov(l, use=use)
  
  return(v)
}

#' Monte Carlo Sample
#' 
#' Generates a Monte Carlo sample based on the multinomial logit-normal model. 
#'
#' @param mu Mean matrix of the underlying distribution
#' @param Sigma Variance matrix of the underlying distribution
#' @param K Number of samples to generate
#'
#' @return K samples from the multinomial logit-normal model. The number of features
#' in the sample is of length(\code{mu})+1.
#'
#' @examples
#' 
#' @export
#'
MCSample <- function(mu, Sigma, K=1) {
  x.norm <- mvtnorm::rmvnorm(K, mu, Sigma)
  ex <- cbind(exp(x.norm), 1)
  return(ex/rowSums(ex))
}

#' Monte Carlo Variation
#' 
#' Estimates the "true" values of the metrics of association for count compositional 
#' data, using large-sample Monte-Carlo methods.  
#'
#' @param mu The mean matrix of the underlying distribution 
#' @param Sigma The variance matrix of the underlying distribution 
#' @param x A sample from the multinomial logit-Normal model. If \code{mu=NULL} 
#' and \code{Sigma=NULL}, then x must be provided
#' @param K Number of Monte Carlo samples to generate
#' @param type Type of variation metric to be calculated: \code{standard}, \code{phi},
#'  \code{phis} (a symmetrical version of \code{phi}), \code{rho}, or \code{logx}
#'
#' @return The variance matrix, \code{v}.
#' 
#' @examples
#' 
#' @export
#'
MCVariation <- function(mu=NULL, Sigma=NULL, x=NULL, K=1e6, 
                        type=c("standard","phi", "phis","rho","logx")) {
  
  if (is.null(mu) & is.null(Sigma) & is.null(x)) {
    stop("If x is missing then mu and Sigma must both not be missing.")
  }
  
  type <- match.arg(type)
  
  if (is.null(x)) {
    x <- MCSample(mu, Sigma, K=K)
  }
  
  l <- log(x)
  J <- NCOL(x)
  
  v <- matrix(0,J,J)
  for (i in 1:J) {
    for (j in 1:J){
      if (type=="standard") {
        v[i,j] <- var(l[,i]-l[,j])
      } else if (type=="phi") {
        v[i,j] <- var(l[,i]-l[,j])/var(l[,i])
      } else if (type=="phis") {
        v[i,j] <- var(l[,i]-l[,j])/(var(l[,i])+var(l[,j]))
      } else if (type=="rho") {
        v[i,j] <- 2*cov(l[,i],l[,j])/(var(l[,i]+l[,j]))
      }
    }
  }
  
  if (type=="logx") v <- cov(l)
  
  return(v)
}

#' Log of the Inverse Additive Log-ratio 
#' 
#' Function which takes the logarithm of a vector, after the inverse additive 
#' log-ratio transformation has been applied to it. 
#'
#' @param x Compositional data vector which has already been transformed by the 
#' additive log-ratio
#'
#' @return A vector which is the log of the inverse of the data which has been
#' transformed by the additive logratio transformation. 
#' 
#' @examples
#' vec <- sample(1:250, 9)
#' ref <- vec[length(vec)]
#' alr.vec <- log(vec[-length(vec)]/ref)
#' 
#' g(alr.vec)
#' 
#'
g <- function(x) {
  ls <- log(1+sum(exp(x)))
  p1 <- x - ls
  c(p1, -ls)
}

#' Logx Variance-Covariance
#' 
#' Function which estimates the variance-covariance of the log of a count 
#' compositional dataset, using a Taylor-series approximation. 
#'
#' @param mu The mean vector of the underlying distribution
#' @param Sigma The sigma matrix of the underlying distribution
#' @param transf The desired transformation. If \code{transf="alr"} the inverse 
#' additive log-ratio transformation is applied. If \code{transf="clr"} the
#' inverse centered log-ratio transformation is applied. 
#' @param order The desired order of the Taylor Series approximation
#'
#' @return The estimated variance-covariance matrix, \code{logx}. 
#' 
#' @examples
#' 
#' @export
#' 
logVarTaylor <- function(mu, Sigma, transf=c("alr", "clr"), order=c("first","second")) {
  transf <- match.arg(transf)
  order <- match.arg(order)
  
  D <- length(mu)
  ones <- rep(1, D)
  emu <- exp(mu)
  if (transf=="alr") {
    ainv <- emu/(1+sum(emu))
  } else {
    ainv <- emu/sum(emu)
  }
  M <- diag(D)-tcrossprod(ones, ainv)
  t2 <- 0
  if (order=="second") {
    mat <- Sigma%*%(tcrossprod(ainv)-diag(ainv))
    t2 <- sum(diag(mat%*%mat))
    #print(t2)
  }
  M%*%tcrossprod(Sigma, M) + 0.5*t2
}

#' Monte Carlo Logx Variance-Covariance
#' 
#' Function which estimates the variance0covariance of the log of the data, using Monte
#' Carlo integration. 
#'
#' @param mu The mean vector of the underlying distribution
#' @param Sigma The variance matrix of the underlying distribution
#' @param K Number of samples 
#'
#' @return The estimated variance-covariance matrix, \code{logx}
#'
#' @examples
#' 
#' @export
#' 
logVarMC <- function(mu, Sigma, K=100000) {
  x.norm <- mvtnorm::rmvnorm(K, mu, Sigma)
  ex <- cbind(exp(x.norm), 1)
  x <- ex/rowSums(ex)
  lx <- log(x)
  cov(lx)
}

#' Full Logx Variance-Covariance
#' 
#' Function which estimates the variance-covariance of the log of the data, using a 
#' Taylor-series approximation. This function differs from \code{Logx Variance} in
#' that the resultant matrix includes a reference category. 
#'
#' @param mu The mean vector of the underlying distribution
#' @param Sigma The sigma matrix of the underlying distribution
#' @param transf The desired transformation. If \code{transf="alr"} the inverse 
#' additive logratio transformation is applied. If \code{transf="clr"} the
#' inverse centered logratio transformation is applied. 
#' @param order The desired order of the Taylor Series approximation
#'
#' @return The estimated variance-covariance matrix, \code{logx}.
#' 
#' @examples
#'
#' @export
#'
logVarTaylorFull <- function(mu, Sigma, transf=c("alr", "clr"), order=c("first", "second")) {
  transf <- match.arg(transf)
  order <- match.arg(order)
  
  D <- length(mu)
  ones <- rep(1, D+1)
  emu <- exp(mu)
  if (transf=="alr") {
    ainv <- emu/(1+sum(emu))
  } else {
    ainv <- emu/sum(emu)
  }
  #cat("note: function has been changed here!\n")
  M <- rbind(diag(D),0)-tcrossprod(ones, ainv)
  #M <- matrix(1, J+1, J) - tcrossprod(ones, ainv)
  t2 <- 0
  if (order=="second") {
    mat <- Sigma%*%(tcrossprod(ainv)-diag(ainv))
    t2 <- sum(diag(mat%*%mat))
    #print(t2)
  }
  M%*%tcrossprod(Sigma, M) + 0.5*t2
}

