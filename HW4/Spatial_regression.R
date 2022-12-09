rm(list=ls())
library(classInt)
library(fields)
library(maps)
library(sp)
library(gstat)
library(geoR)
library(mvtnorm)
library(MCMCpack)
library(coda)

#############################
## California temperatures ##
#############################

load("CAtemps.RData")

ploteqc <- function(spobj, z, breaks, ...){
  pal <- tim.colors(length(breaks)-1)
  fb <- classIntervals(z, n = length(pal), 
                       style = "fixed", fixedBreaks = breaks)
  col <- findColours(fb, pal)
  plot(spobj, col = col, ...)
  image.plot(legend.only = TRUE, zlim = range(breaks), col = pal)
}

## Plotting

range(CAtemp$avgtemp)
breaks <- seq(40, 75, by = 5)
ploteqc(CAtemp, CAtemp$avgtemp, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Average Annual Temperatures, 1961-1990, Degrees F")

range(CAgrid$elevation)
breaks <- seq(-100, 3600, by = 100)
ploteqc(CAgrid, CAgrid$elevation, breaks, pch = 19)
map("county", region = "california", add = TRUE)
title(main = "Elevations at prediction locations, m")

## Prior parameters
linmod <- lm(avgtemp~lon+lat+elevation, data = CAtemp)
summary(linmod)

m.beta <- rep(0, 4); V.beta <- 1000 * diag(4)
a.s2 <- 0.001; b.s2 <- 0.001
a.t2 <- 0.001; b.t2 <- 0.001

rhoseq <- seq(0.01, 300, length = 100)
a.rho <- 2; b.rho <- 50
plot(rhoseq, dgamma(rhoseq, shape = a.rho, scale = b.rho), type = "l") # prior for rho

## Setup, storage, and starting values

y <- CAtemp$avgtemp
n <- nrow(CAtemp); m <- nrow(CAgrid)
d <- rdist.earth(coordinates(CAtemp))
X <- cbind(rep(1, n), CAtemp$lon, CAtemp$lat, CAtemp$elevation)
Xpred <- cbind(rep(1, m), CAgrid$lon, CAgrid$lat, CAgrid$elevation)

B <- 2000

beta.samps <- matrix(NA, nrow = 4, ncol = B)
beta.samps[,1] <- coef(linmod) # INITIALIZE VALUES

s2.samps <- t2.samps <- rho.samps <- rep(NA, B)
s2.samps[1] <- 5
rho.samps[1] <- 10
t2.samps[1] <- 2

eta.obs.samps <- matrix(NA, nrow = n, ncol = B)

v.prop <- 0.1^2 # step size

## MCMC sampler

Gamma <- exp(-d/rho.samps[1]) # initalize Gamma matrix
Ginv <- solve(Gamma)

for(i in 2:B){
  
  if(i%%100==0) print(i)
  
  ## eta_obs | Rest
  V <- solve((1/t2.samps[i-1])*diag(dim(X)[1])+(1/s2.samps[i-1])*Ginv)
  m <- V%*%((1/t2.samps[i-1])*y+(1/s2.samps[i-1])*Ginv%*%X%*%beta.samps[,i-1])
  eta.obs.samps[,i] <- rmvnorm(1, mean = m, sigma = V, method = "svd")
  
  ## beta | Rest
  V <- solve((t(X)%*%Ginv%*%X)/s2.samps[i-1]+(diag(4)/1000))
  m <- V%*%(t(X)%*%Ginv%*%eta.obs.samps[,i]/s2.samps[i-1])
  beta.samps[,i] <- rmvnorm(1, mean = m, sigma = V, method = "svd")

  ## s2 | Rest
  a <- a.s2+length(y)/2
  b <- b.s2+t(y-X%*%beta.samps[,i])%*%Ginv%*%(y-X%*%beta.samps[,i])/2
  s2.samps[i] <- rinvgamma(1, a, b)
  
  ## t2 | Rest
  a <- a.t2+length(y)/2
  b <- b.t2+t(y-eta.obs.samps[,i])%*%(y-eta.obs.samps[,i])/2
  t2.samps[i] <- rinvgamma(1, a, b)
  
  ## rho | Rest   
  # MH update
  rho.prop <- rnorm(1, mean=rho.samps[i-1], sd=sqrt(v.prop))
  if(rho.prop<0){
    rho.samps[i] <- rho.samps[i-1]
  }else{
    gamma.prop <- exp(-d/rho.prop)
    ll.post.prop <- log(rho.prop)-0.5*log(det(gamma.prop))-50*rho.prop-t(eta.obs.samps[,i]-X%*%beta.samps[,i])%*%solve(gamma.prop)%*%(eta.obs.samps[,i]-X%*%beta.samps[,i])/(2*s2.samps[i])
    ll.post.cur <- log(rho.samps[i-1])-0.5*log(det(Gamma))-50*rho.samps[i-1]-t(eta.obs.samps[,i]-X%*%beta.samps[,i])%*%Ginv%*%(eta.obs.samps[,i]-X%*%beta.samps[,i])/(2*s2.samps[i])
    thr <- ll.post.prop-ll.post.cur
    if(!is.nan(thr)){
      if(log(runif(1))<thr){
        rho.samps[i] <- rho.prop
        Gamma <- gamma.prop
        Ginv <- solve(gamma.prop)
      }
      else{
        rho.samps[i] <- rho.samps[i-1]
      }
    }
  }
}

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

for(i in 1:4){
  diagnosis(beta.samps[i,1001:2000], paste("beta", as.character(i)))
}
diagnosis(eta.obs.samps[1, 1001:2000], "eta")
diagnosis(s2.samps[1001:2000], "s2")
diagnosis(t2.samps[1001:2000], "t2")
diagnosis(rho.samps[1001:2000], "rho")

## Prediction

dcross <- rdist.earth(coordinates(CAtemp), coordinates(CAgrid)) # gamma_cross
dpred <- rdist.earth(coordinates(CAgrid)) # Gamma_pred


eta.pred <- matrix(NA, nrow = nrow(CAgrid), ncol = B)

for(j in 1:B){
  print(j)
  Gcross <- exp(-dcross/rho.samps[j])
  Gpred <- exp(-dpred/rho.samps[j])
  Gobs <- exp(-d/rho.samps[j])
  m <- Xpred%*%beta.samps[,j]+t(Gcross)%*%solve(Gobs)%*%(y-X%*%beta.samps[,j])
  V <- s2.samps[j]*(Gpred-t(Gcross)%*%solve(Gobs)%*%Gcross)
  eta.pred[,j] <- rmvnorm(1, m, V, method = "svd")
}

range(rowMeans(eta.pred[,1001:2000]))
breaks <- seq(30, 75, by = 5)
ploteqc(CAgrid, rowMeans(eta.pred[,1001:2000]) , breaks, pch = 19)
map("county", region = "california", add = TRUE)