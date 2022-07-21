#' Function to calculate variation matrix under various models for lattice compositional data.
#'
#' @param counts 
#'
#' @return
#' @export
#'
#' @examples
varEst <- function(counts, p.model=c("logitNormal", "dirichlet", "plugin"), type=c("standard","phi","rho"), refColumn=NULL) {
  if (any(counts<0) | any(counts!=floor(counts)) | !is.matrix(counts)) stop("'counts' must be a matrix containing non-negative integers")
  p.model <- match.arg(p.model)
  type <- match.arg(type)
  
  if (is.null(refColumn)&p.model=="logitNormal") { 
    refColumn <- J
  }
  
  if (p.model=="logitNormal") {
    result <- logitNormalVariation(mu, Sigma, lmu, lsigma, type=type)
  } else if (p.model=="dirichlet") {
    result <- dirichletVariation(counts)
  } else {
    result <- pluginVariation(counts)
  }
  
  return(result)
}

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

logitNormalVariationOLD <- function(mu, Sigma, lmu, lsigma, type=c("standard","phi","rho")) {
  type <- match.arg(type)
  J <- length(mu)
  
  d.S <- diag(Sigma)
  ones <- rep(1, J)
  
  tcp.dS <- tcrossprod(d.S, ones)
  Elr <- tcp.dS + t(tcp.dS) - 2*Sigma
  En <- exp(-lmu+lsigma^2/2)
  
  tcp.e.mu.dS <- tcrossprod(exp(-mu+d.S/2), ones)
  sum.rows <- rowSums(tcp.e.mu.dS+tcp.dS/2-Sigma)
  
  Ep <- tcp.e.mu.dS*(1+tcrossprod(sum.rows, ones))
  
  V <- Elr + En*(2+Ep+t(Ep))
  diag(V) <- 0
  return(V)
}

dirichletVariation <- function(counts) {
  return(NULL)
}

pluginVariation <- function(counts) {
  return(NULL)
}

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

# Function to get MCMC samples from multinom logit-normal model
# K is the number of samples to return. Number of features returned
# is length(mu)+1
MCSample <- function(mu, Sigma, K=1) {
  x.norm <- mvtnorm::rmvnorm(K, mu, Sigma)
  ex <- cbind(exp(x.norm), 1)
  return(ex/rowSums(ex))
}

# Calculating variation using Monte Carlo integration
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

meanLN <- function(mu, Sigma, ind) {
  sd <- diag(Sigma)[-ind]
  sr <- Sigma[ind,-ind]
  si <- Sigma[ind, ind]
  mu.sub <- mu[-ind]
  e <- exp(mu.sub - mu[ind] + 0.5*(sd+si-2*sr))
  1 + exp(-mu[ind]+si/2) + sum(e)
}

meanLOGN <- function(lmu, lsigma) {
  exp(-lmu+lsigma^2/2)
}

logVarMC <- function(mu, Sigma, K=100000) {
  x.norm <- mvtnorm::rmvnorm(K, mu, Sigma)
  ex <- cbind(exp(x.norm), 1)
  x <- ex/rowSums(ex)
  lx <- log(x)
  cov(lx)
}

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

# Includes reference category
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

# Estimating variance using unscented transformation
# https://users.isy.liu.se/en/rt/fredrik/reports/07SSPut.pdf
logVarUnscented <- function(mu, Sigma, transf=c("alr", "clr"), alpha=1e-3, beta=2, kappa=0) {
  transf <- match.arg(transf)
  
  D <- length(mu)
  lambda <- alpha^2*(D+kappa)-D
  
  ones <- rep(1, 2*D+1)
  ones.short <- rep(1, D)
  
  # svd.Sigma <- svd(Sigma)
  # v <- svd.Sigma$v # This is the transpose of u from Hendeby/Gustafsson paper
  # s <- svd.Sigma$d
  #eig <- eigen(Sigma+D+lambda)
  #eig <- eigen(Sigma+D)
  #sqrt.Sigma <- eig$vectors%*%tcrossprod(diag(sqrt(eig$values)), eig$vectors)
  L <- t(chol(D*Sigma))
  
  x.pm <- tcrossprod(ones, mu)
  # pm.term <- sqrt(D+lambda)*tcrossprod(s, ones.short)*v
  # pm.term <- tcrossprod(s, ones.short)*v
  # x.pm[1:D,] <- x.pm[1:D,] + pm.term
  # x.pm[(D+2):(2*D+1),] <- x.pm[(D+2):(2*D+1),] - pm.term
  x.pm[1:D,] <- x.pm[1:D,] + t(L) #sqrt.Sigma
  x.pm[(D+2):(2*D+1),] <- x.pm[(D+2):(2*D+1),] - L #sqrt.Sigma
  
  #omega.pm <- rep(1/(2*(D+lambda)), 2*D+1)
  #omega.pm[D+1] <- lambda/(D+lambda)
  omega.pm <- rep(1, D+1)
  
  z <- t(apply(x.pm, 1, g))
  mu.z <- colSums(z*omega.pm)
  zmm <- z - mu.z
  
  P <- matrix(0, D+1, D+1)
  for (i in 1:(2*D+1)) {
    #a <- ifelse(i==D+1, (1-alpha^2+beta), 0)
    #P <- P + (omega.pm[i]+a)*tcrossprod(zmm[i,])
  }
  
  P <- cov(z)
  
  return(P)
}

logVarUnscentedW <- function(mu, Sigma, transf=c("alr", "clr")) {
  transf <- match.arg(transf)
  
  D <- length(mu)
  
  ones <- rep(1, 2*D+1)
  ones.short <- rep(1, D)
  
  # svd.Sigma <- svd(Sigma)
  # v <- svd.Sigma$v # This is the transpose of u from Hendeby/Gustafsson paper
  # s <- svd.Sigma$d
  #eig <- eigen(Sigma+D+lambda)
  #eig <- eigen(Sigma+D)
  #sqrt.Sigma <- eig$vectors%*%tcrossprod(diag(sqrt(eig$values)), eig$vectors)
  w0 <- 0.05
  L <- t(chol(D/(1-w0)*Sigma))
  
  x.pm <- tcrossprod(ones, mu)
  # pm.term <- sqrt(D+lambda)*tcrossprod(s, ones.short)*v
  # pm.term <- tcrossprod(s, ones.short)*v
  # x.pm[1:D,] <- x.pm[1:D,] + pm.term
  # x.pm[(D+2):(2*D+1),] <- x.pm[(D+2):(2*D+1),] - pm.term
  x.pm[1:D,] <- x.pm[1:D,] + t(L) #sqrt.Sigma
  x.pm[(D+2):(2*D+1),] <- x.pm[(D+2):(2*D+1),] - L #sqrt.Sigma
  
  omega.pm <- rep((1-w0)/(2*(D)), 2*D+1)
  omega.pm[D+1] <- w0
  
  z <- t(apply(x.pm, 1, g))
  mu.z <- colSums(z*omega.pm)
  zmm <- z - mu.z
  
  P <- matrix(0, D+1, D+1)
  for (i in 1:(2*D+1)) {
    P <- P + omega.pm[i]*tcrossprod(zmm[i,])
  }
  
  return(P)
}

# Transformation of interest: g(x)=log(alr.inv(x))
g <- function(x) {
  ls <- log(1+sum(exp(x)))
  p1 <- x - ls
  c(p1, -ls)
}
