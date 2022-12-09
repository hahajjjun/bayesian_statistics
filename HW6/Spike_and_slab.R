#################
### Problem 1.###
#################
install.packages("nimble")
library(nimble)
library(mvtnorm)
library(invgamma)

X <- rmvnorm(10000, rep(0,10), diag(10))
beta <- c(0.5, -0.5, 1, -1, 0.7, 0,0,0,0,0)
mu <- exp(X%*%beta)/(exp(X%*%beta)+1)
y <- rep(0, 10000)
for(i in 1:length(mu)){
  y[i] <- rbinom(1, 1, mu[i])
}

##################
### Problem 2. ###
##################

log_likelihood <- function(X, y, beta){
  ll <- 0
  mu <- exp(X%*%beta)/(1+exp(X%*%beta))
  for(i in 1:length(y)){
    ll <- ll + dbinom(y[i], 1, mu[i], log=TRUE)
  }
  return(ll)
}
log_beta_prior <- function(beta, lambda, ss1, ss2){
  lp <- 0
  for(i in 1:length(beta)){
    lp <- lp + dnorm(beta[i], 0, sqrt(ss1*lambda[i]^2+ss2*(1-lambda[i])^2), log=TRUE)
  }
  return(lp)
}
log_ss1_prior <- function(ss1){
  return(dinvgamma(ss1, 1, scale=20, log=TRUE)*10)
}
log_ss2_prior <- function(ss2){
  return(dgamma(ss2, 1, 20, log=TRUE)*10)
}
log_beta_posterior <- function(X, y, beta, lambda, ss1, ss2){
  return(log_likelihood(X, y, beta)+log_beta_prior(beta, lambda, ss1, ss2))
}
log_lambda_posterior <- function(beta, lambda, ss1, ss2){
  return(log_beta_prior(beta, lambda, ss1, ss2))
}
log_ss1_posterior <- function(beta, lambda, ss1, ss2){
  return(log_beta_prior(beta, lambda, ss1, ss2)+log_ss1_prior(ss1))
}
log_ss2_posterior <- function(beta, lambda, ss1, ss2){
  return(log_beta_prior(beta, lambda, ss1, ss2)+log_ss2_prior(ss2))
}


beta_chain = matrix(nrow=1, ncol=10)
lambda_chain = matrix(nrow=1, ncol=10)
ss1_chain = c(1)
ss2_chain = c(0.1)
beta_chain[1,] <- rep(1, 10)
lambda_chain[1,] <- rep(1, 10)
step_size = 0.1

for(i in 1:2000){
  if(i%%100==0){
    cat(i, "th iteration \n")
  }
  beta_init <- beta_chain[i,]
  lambda_init <- lambda_chain[i,]
  ss1_init <- ss1_chain[i]
  ss2_init <- ss2_chain[i]
  for(j in 1:10){
    beta_prop <- beta_init
    beta_prop[j] <- beta_init[j]+rnorm(1, 0, sd=step_size)
    a <- log_beta_posterior(X, y, beta_prop, lambda_init, ss1_init, ss2_init)
    b <- log_beta_posterior(X, y, beta_init, lambda_init, ss1_init, ss2_init)
    thr <- a-b
    if(!is.nan(thr)){
      if(log(runif(1))<thr){
        beta_init <- beta_prop
      }
    }
  }
  for(j in 1:10){
    lambda_prop <- lambda_init
    lambda_prop[j] <- rbinom(1, 1, 0.5)
    a <- log_lambda_posterior(beta_init, lambda_prop, ss1_init, ss2_init)
    b <- log_lambda_posterior(beta_init, lambda_init, ss1_init, ss2_init)
    thr <- a-b
    if(!is.nan(thr)){
      if(log(runif(1))<thr){
        lambda_init <- lambda_prop
      }
    }
  }
  ss1_prop <- ss1_init+rnorm(1, 0, sd=step_size*50)
  a <- log_ss1_posterior(beta_init, lambda_init, ss1_prop, ss2_init)
  b <- log_ss1_posterior(beta_init, lambda_init, ss1_init, ss2_init)
  thr <- a-b
  if(!is.nan(thr)){
    if(log(runif(1))<thr){
      ss1_init <- ss1_prop
    }
  }
  ss2_prop <- ss2_init+rnorm(1, 0, sd=step_size*0.05)
  a <- log_ss2_posterior(beta_init, lambda_init, ss1_init, ss2_prop)
  b <- log_ss2_posterior(beta_init, lambda_init, ss1_init, ss2_init)
  thr <- a-b
  if(!is.nan(thr)){
    if(log(runif(1))<thr){
      ss2_init <- ss2_prop
    }
  }
  beta_chain <- rbind(beta_chain, beta_init)
  lambda_chain <- rbind(lambda_chain, lambda_init)
  ss1_chain <- c(ss1_chain, ss1_init)
  ss2_chain <- c(ss2_chain, ss2_init)
}

diagnosis <- function(chain, name){
  ts.plot(chain, main=name)
  plot(density(chain), main=name)
  acf(chain, main=name)
  cat("[",name,"]\n")
  cat("Posterior mean:",mean(chain),"\n")
  cat("Acceptance prob:",length(unique(chain))/length(chain), "\n")
}


for(i in 1:10){
  diagnosis(beta_chain[,i], paste("beta",as.character(i)))
}
for(i in 1:10){
  diagnosis(lambda_chain[,i], paste("lambda",as.character(i)))
}
diagnosis(ss1_chain, "ss1")
diagnosis(ss2_chain, "ss2")



for(i in 1:10){
  diagnosis(samples_nimble[,i], paste("beta",as.character(i)))
}
for(i in 1:10){
  diagnosis(samples_nimble[,10+i], paste("lambda",as.character(i)))
}
diagnosis(samples_nimble[,21], "ss1")
diagnosis(samples_nimble[,22], "ss2")

###################
### NIMBLE CODE ###
###################
model <- nimbleCode({
 XB[1:N] <- X[1:N, 1:p]%*%beta[1:p]
 for(i in 1:N){
   mu[i] <- exp(XB[i])/(1+exp(XB[i]))
   y[i] ~ dbern(mu[i])
 }
 
 ss1 ~ dinvgamma(shape=1, scale=20)
 ss2 ~ dgamma(shape=1, rate=20)
 for(i in 1:p){
   lambda[i] ~ dbern(0.5)
   beta[i] ~ dnorm(mean=0, sd=sqrt(ss1*(lambda[i]^2)+ss2*(1-lambda[i])^2))
 }
})

consts <- list(N=dim(X)[1], p=dim(X)[2])
dat <- list(X=X, y=y)
inits <- list(beta=rep(0, dim(X)[2]), lambda=rep(0, dim(X)[2]), ss1=0.5, ss2=0.5)
samples_nimble <- nimbleMCMC(model, data=dat, inits=inits, constants=consts, monitors=c("beta","lambda","ss1","ss2"), samplesAsCodaMCMC = TRUE, WAIC=FALSE, summary=FALSE, niter=1e3, nburnin=0, thin=4, nchains=1)
##################################
### POSTERIOR PREDICTIVE CHECK ###
##################################
y_pred <- c()
for(i in 1:1000){
  if(i%%100==0){
    cat(i,"th iteration \n")
  }
  y_pred_t <- c()
  for(j in 1:dim(X)[1]){
    y_pred_t <- c(y_pred_t, rbinom(1, 1, exp(X[j,]%*%beta_chain[i,])/(exp(X[j,]%*%beta_chain[i,])+1))) 
  }
  y_pred <- c(y_pred, mean(y_pred_t))
}
cat("Bayesian p-value:", sum(y_pred>mean(y))/length(y_pred), "\n")

