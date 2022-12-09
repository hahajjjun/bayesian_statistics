##########################
### JHPark HW 03 Recap ###
##########################

install.packages("mvtnorm")
install.packages("nimble")
install.packages("adaptMCMC")

library(mvtnorm)
library(coda)
library(nimble)
library(adaptMCMC)

diagnosis <- function(vec, name){
  ts.plot(vec, main=name)
  plot(density(vec), main=name)
  acf(vec, main=name)
  #x <- seq(-100, 100, by=0.1)
  #lines(x, dcauchy(x, eta, theta), col="red")
  cat("[",name,"]\n")
  cat("Acceptance rate:", length(unique(vec))/length(vec), "\n")
  cat("HPD interval:", HPDinterval(as.mcmc(vec)), "\n")
  cat("Posterior mean:", mean(vec), "\n")
  cat("Effective samle size:", effectiveSize(vec), "\n")
}

# 1. Simulation

X <- rmvnorm(1000, mean=rep(0, 4), sigma = diag(4))
beta <- c(0.5, -0.5, 0, 1)
mu <- exp(X%*%beta)
y <- rep(0, 1000)
for(i in 1:length(mu)){
  y[i] <- rpois(n=1, lambda=mu[i])
}
glm(y~X, family=poisson) # Check

# 2. MCMC from scratch
prior <- function(beta){
  ll <- 0
  for(i in 1:length(beta)){
    ll <- ll + dnorm(beta[i], mean=0, sd=sqrt(10), log=TRUE)
  }
  return(ll)
}

likelihood <- function(X, y, beta){
  ll <- 0
  lambda <- exp(X%*%beta)
  for(i in 1:length(y)){
    ll <- ll + dpois(y[i], lambda[i], log=TRUE)
  }
  return(ll)
}

posterior <- function(X,y,beta){
  ll <- likelihood(X,y,beta)+prior(beta)
  return(ll)
}

beta_init <- c(0,0,0,0)
chain <- beta_init
indicator <- 0
for(idx in 1:1500){
  beta_prop <- t(rmvnorm(1, mean=beta_init, sigma=0.001*diag(4)))
  thr <- posterior(X,y,beta_prop)-posterior(X,y,beta_init)
  if(!is.nan(thr)){
    if(log(runif(1))<thr){
    beta_init <- beta_prop
    }
  }
  if(indicator%%5==0){
    chain <- cbind(chain, beta_init)
  }
  indicator <- indicator + 1
}
chain <-chain[,-1]
diagnosis(chain[1,],"beta1")
diagnosis(chain[2,],"beta2")
diagnosis(chain[3,],"beta3")
diagnosis(chain[4,],"beta4")

diagnosis(t(chain), "beta")


# 3. Nimble

model_glm <- nimbleCode({
  # Data Model
  for(i in 1:n){
    lambda[i] <- exp(XB[i])
    Z[i] ~ dpois(lambda=lambda[i])
  }
  XB[1:n] <- X[1:n, 1:p]%*%beta[1:p]
  
  # Parameter Model
  beta[1:p] ~ dmnorm(mean=M[1:p], cov=Cov[1:p,1:p])
})

niter <- 2e3
consts <- list(n=dim(X)[1], p=dim(X)[2], X=X, M=rep(0,dim(X)[2]), Cov=10*diag(dim(X)[2]))
dat <- list(Z=y)
inits <- list(beta=rep(0,dim(X)[2]))

samples_glm <- nimbleMCMC(model_glm, data=dat, inits=inits, constants=consts, monitors=c("beta"), samplesAsCodaMCMC=TRUE, WAIC=FALSE, summary=TRUE, niter=niter, nburnin=1000, thin=1, nchains=1)

diagnosis(samples_glm$samples[,1], "beta1")
diagnosis(samples_glm$samples[,2], "beta2")
diagnosis(samples_glm$samples[,3], "beta3")
diagnosis(samples_glm$samples[,4], "beta4")
# 4. adaptMCMC

init.pars <- c(beta_1=0, beta_2=0, beta_3=0, beta_4=0)

logPost <- function(pars){
  with(as.list(pars), {
    beta <- pars
    lambda <- exp(X%*%beta)
    sum(dpois(y,lambda,log=TRUE))+sum(dnorm(beta,mean=0,sd=sqrt(10),log=TRUE))
  })
}

adaptMCMC_modify <- function(burn_in, thinning, niter){
  out.mcmc <- MCMC(p=logPost, n=niter*thinning+burn_in, init=init.pars, scale=c(0.1,0.1,0.1,0.1),adapt=TRUE,acc.rate=.5)
  samples <- NULL
  for(i in 1:niter){
    samples <- cbind(samples, out.mcmc$samples[i*thinning+burn_in,])
  }
  return(samples)
}

adaptMCMC_samples <- adaptMCMC_modify(burn_in=500, thinning=2, niter=1e3)
diagnosis(adaptMCMC_samples[1,], "beta1")
