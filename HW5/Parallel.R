####################
### HW 05 JHPARK ###
####################

###### PROBLEM1 ######
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(coda)
library(parallel)

X <- rmvnorm(300000, rep(0,4), diag(4))
beta <- c(0.5, -0.5, 0, 1)
mu <- exp(X%*%beta)/(1+exp(X%*%beta))
y <- rep(0, length(mu))
for(i in 1:length(y)){
  y[i] <- rbinom(1, 1, mu[i])
}
# VERIFY:
# glm(y~X, family=binomial)

###### PROBLEM 2 ######
sourceCpp("utils.cpp")

############################################
#R_log_likelihood <- function(beta, Y, X){
#  mu <- exp(X%*%beta)/(1+exp(X%*%beta))
#  return(sum(dbinom(Y, 1, mu, log=TRUE)))
#}
############################################
Rcpp_samples <- Rcpp_MCMC(2000, c(1,1,1,1), y, X, 0.005)

diagnosis <- function(chain, name){
  ts.plot(chain, main=name)
  plot(density(chain), main=name)
  acf(chain, main=name)
  cat("[",name,"]\n")
  cat("Acceptance rate:", length(unique(chain))/length(chain) ,"\n")
  cat("Effective sample size:", effectiveSize(as.mcmc(chain)), "\n")
  cat("HPD interval:", HPDinterval(as.mcmc(chain)), "\n")
  cat("Posterior mean:", mean(chain), "\n")
}
for(i in 1:4){
  diagnosis(Rcpp_samples[1001:2000,i], paste("beta",as.character(i)))
}

n_cores = detectCores()

ptm <- proc.time()
Rcpp_parallel_samples <- Rcpp_MCMC_parallel(2000, n_cores, y, X, 0.005)
rcpp_time <- proc.time()-ptm

weights <- array(0, dim=c(4, 4, n_cores))
weighted_sum <- array(0, dim=c(2000, 4, n_cores))
for(i in 1:n_cores){
  weights[,,i] <- solve(cov(Rcpp_parallel_samples[i,,]))
  weighted_sum[,,i] <- Rcpp_parallel_samples[i,,]%*%weights[,,i]
}
sum_weights <- rowSums(weights, dims=2)
weighted_samples <- rowSums(weighted_sum, dims=2)%*%solve(sum_weights)


#install.packages("SBmedian")
library(SBmedian)
c = list()
for(i in 1:4){ # 16개는 무리...
  c[[i]] = Rcpp_parallel_samples[i,1001:1500,]
}

run <- mpost.euc(c)
ts.plot(run$med.atoms)
ts.plot(run$weiszfeld.history)

for(j in 1:4){
  for(i in 1:n_cores){
    cat("Core",i,"\n")
    cat(length(unique(Rcpp_parallel_samples[i,,j]))/2000,"\n")
    cat(effectiveSize(Rcpp_parallel_samples[i,,j]),"\n")
    cat(mean(Rcpp_parallel_samples[i,,j]),"\n")
    cat(HPDinterval(as.mcmc(Rcpp_parallel_samples[i,,j])),"\n")
  }
  plot(x=seq(1:2000), y=Rcpp_parallel_samples[1,,j], col=rgb(0, 1,1), type='l',  main=paste("TS plot of beta",j, sep=""), xlab="iteration",ylab="value")
  for(i in 2:n_cores){
    lines(x=seq(1:2000), y=Rcpp_parallel_samples[i,,j], col=rgb(0, (n_cores-i)/n_cores, (n_cores-i)/n_cores))
  }
  lines(x=seq(1:2000), y=weighted_samples[,j], col=rgb(1,0,0))
  
  plot(density(Rcpp_parallel_samples[1,,j]), col=rgb(0,1,1), main=paste("Density plot of beta",j, sep=""), xlab="value")
  for(i in 2:n_cores){
    lines(density(Rcpp_parallel_samples[i,,j]), col=rgb(0, (n_cores-i)/n_cores, (n_cores-i)/n_cores))
  }
  lines(density(weighted_samples[,j]), col=rgb(1,0,0))
}

for(j in 1:4){
  c = list()
  for(i in 1:n_cores){
    c[[i]] = as.mcmc(Rcpp_parallel_samples[i,1001:2000,j])
  }
  print(gelman.diag(c))
}