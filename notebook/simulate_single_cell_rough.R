# Simulate from Multinom LN model using parameters from single-cell dataset

library(zCompositions)
source("~/Dropbox/Documents/research/code/variation_estimation/variation_functions.R")
source("~/Dropbox/Documents/research/code/variation_estimation/mle_multinom_logit/mle_multinom_logit_normal.R")

# Load in MLE data
load("~/Dropbox/Documents/research/code/variation_estimation/single_cell/data/mle_results_p100.RData")

n.g <- length(mle$est.min$mu)+1
mu <- mle$est.min$mu
Sigma <- mle$est.min$Sigma

# Need to calculate "true" (MC-based) values before subsetting
# True value of v using exact formula
# "True" values of phi, phis, and rho using large MC sample
K <- 1e5 # Number of MC samples
mc <- MCSample(mu, Sigma, K)
v.true <- logitNormalVariation(mu, Sigma)
phi.true <- MCVariation(x=mc, type="phi")#[-n.g,-n.g]
phis.true <- MCVariation(x=mc, type="phis")#[-n.g,-n.g]
rho.true <- MCVariation(x=mc, type="rho")#[-n.g,-n.g]
logx.true <- MCVariation(x=mc, type="logx")#[-n.g,-n.g]

# Simulating data
#set.seed(66573486)
n <- 100 #n.g+1
x.logit <- mvtnorm::rmvnorm(n, mu, Sigma)
xl.exp <- cbind(exp(x.logit), 1)
x <- xl.exp/rowSums(xl.exp)

rd <- rowSums(dat.ss)
ln.mean <- mean(log(rd))
ln.var <- var(log(rd))
ni <- rlnorm(n, ln.mean, ln.var)

y <- mc2d::rmultinomial(n, ni, x)

# Filtering out genes with too many zeros
thresh <- 0.2
pzero <- apply(y, 2, function(x){mean(x==0)})
incl <- pzero<thresh
y <- y[,incl]

# Subsetting true values
v.true <- v.true[incl[-n.g], incl[-n.g]]
phi.true <- phi.true[incl, incl]
phis.true <- phis.true[incl, incl]
rho.true <- rho.true[incl, incl]
logx.true <- logx.true[incl, incl]

# Looking at mu and Sigma using subset of genes
#mu <- mu[incl[-n.g]]
#Sigma <- Sigma[incl[-n.g],incl[-n.g]]

# Updating number of genes after filtering
n.g <- NCOL(y)

# Bayesian imputation of zeros
y.no0 <- as.matrix(cmultRepl(y, output = "p-count"))

n.lam <- 8
lmr <- 0.01
#lam <- c(0.2, 0.5)
mle.sim <- mlePath(y, tol=1e-4, tol.nr=1e-4, n.lambda = n.lam, 
                   lambda.min.ratio = lmr,
                   gamma = 0.1)
# mle.sim <- mlePath(y, tol=1e-4, tol.nr=1e-4, lambda.gl = 0,
#                    gamma = 0.1)
ebicPlot(mle.sim, xlog = TRUE)

mu.hat <- mle.sim$est.min$mu
Sigma.hat <- mle.sim$est.min$Sigma

# Naive versions
v.naive <- naiveVariation(y.no0)[-n.g,-n.g]
phi.naive <- naiveVariation(y.no0, type="phi")[-n.g,-n.g]
phis.naive <- naiveVariation(y.no0, type="phis")[-n.g,-n.g]
rho.naive <- naiveVariation(y.no0, type="rho")[-n.g,-n.g]
logx.naive <- naiveVariation(y.no0, type="logx")[-n.g,-n.g]

# Estimated versions
v.est <- logitNormalVariation(mu.hat, Sigma.hat)
phi.est <- logitNormalVariation(mu.hat, Sigma.hat, type="phi", order="second")
phis.est <- logitNormalVariation(mu.hat, Sigma.hat, type="phis", order="second")
rho.est <- logitNormalVariation(mu.hat, Sigma.hat, type="rho", order="second")
logx.est <- logVarTaylor(mu.hat, Sigma.hat)

# Plots
# Plot for v
lt <- lower.tri(v.est)
rge <- range(v.naive, v.est)
plot(v.true[lt], v.naive[lt], pch=19, cex=0.5, col="red",ylim=rge)
abline(0,1)
points(v.true[lt], v.est[lt], pch=19, cex=0.5, col="blue")

# Plot for phi
lt <- lower.tri(phi.true[-n.g,-n.g])
rge <- range(phi.naive[lt], phi.est[lt])
plot(phi.true[-n.g,-n.g][lt], phi.naive[lt], pch=19, cex=0.5, col="red",ylim=rge)
abline(0,1)
points(phi.true[-n.g,-n.g][lt], phi.est[lt], pch=19, cex=0.5, col="blue")

# Plot for phis
lt <- lower.tri(phis.true[-n.g,-n.g])
rge <- range(phis.naive[lt], phis.est[lt])
plot(phis.true[-n.g,-n.g][lt], phis.naive[lt], pch=19, cex=0.5, col="red",ylim=rge)
abline(0,1)
points(phis.true[-n.g,-n.g][lt], phis.est[lt], pch=19, cex=0.5, col="blue")

# Plot for rho
lt <- lower.tri(rho.true[-n.g,-n.g])
rge <- range(rho.naive[lt], rho.est[lt])
plot(rho.true[-n.g,-n.g][lt], rho.naive[lt], pch=19, cex=0.5, col="red",ylim=rge)
abline(0,1)
points(rho.true[-n.g,-n.g][lt], rho.est[lt], pch=19, cex=0.5, col="blue")

# Plot for logx
lt <- lower.tri(logx.true[-n.g,-n.g])
rge <- range(logx.naive[lt], logx.est[lt])
plot(logx.true[-n.g,-n.g][lt], logx.naive[lt], pch=19, cex=0.5, col="red",ylim=rge)
abline(0,1)
points(logx.true[-n.g,-n.g][lt], logx.est[lt], pch=19, cex=0.5, col="blue")

# Plot for diagonal of logx
di.true <- diag(logx.true[-n.g,-n.g])
di.naive <- diag(logx.naive)
di.est <- diag(logx.est)
rge <- range(di.naive, di.est)
plot(di.true, di.naive, pch=19, cex=0.5, col="red",ylim=rge)
abline(0,1)
points(di.true, di.est, pch=19, cex=0.5, col="blue")

# Sums of pairs of diagonals of logx
spair.true <- outer(di.true, di.true, "+")
spair.naive <- outer(di.naive, di.naive, "+")
spair.est <- outer(di.est, di.est, "+")
rge <- range(spair.naive, spair.est)
plot(spair.true, spair.naive, pch=19, cex=0.5, col="red",ylim=rge)
abline(0,1)
points(spair.true, spair.est, pch=19, cex=0.5, col="blue")

# RMSE
# v
rmse(v.true[lt], v.naive[lt])
rmse(v.true[lt], v.est[lt])

# phi
lt <- lower.tri(phi.est)
rmse(phi.true[-n.g,-n.g][lt], phi.naive[lt])
rmse(phi.true[-n.g,-n.g][lt], phi.est[lt])

# phis
rmse(phis.true[-n.g,-n.g][lt], phis.naive[lt])
rmse(phis.true[-n.g,-n.g][lt], phis.est[lt])

# rho
rmse(rho.true[-n.g,-n.g][lt], rho.naive[lt])
rmse(rho.true[-n.g,-n.g][lt], rho.est[lt])

# logx
rmse(logx.true[-n.g,-n.g][lt], logx.naive[lt])
rmse(logx.true[-n.g,-n.g][lt], logx.est[lt])

# Diagonal logx
rmse(di.true, di.naive)
rmse(di.true, di.est)


