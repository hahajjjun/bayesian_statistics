###########################################
### JHPARK, Cauchy distribution example ###
###########################################
install.packages("dgof")
install.packages("coda")
library(dgof)
library(coda)

custom_dcauchy <- function(x, eta, theta){
  f <- 1/(theta*pi*(1+((x-eta)/theta)^2))
  return(log(f))
}

theta <- 2
eta <- 0
init <- 0
chain <- c(init)
for(i in 1:1e4){
  prop <- init+rnorm(n=1, mean=0, sd=10)
  thr <- custom_dcauchy(prop, eta, theta)-custom_dcauchy(init, eta, theta)
  if(log(runif(n=1))<thr){
    init <- prop
  }
  chain <- c(chain, init)
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

diagnosis(chain, "sample")

sort(chain)[2500];sort(chain)[5000];sort(chain)[7500]
qcauchy(0.25, eta, theta);qcauchy(0.5, eta, theta);qcauchy(0.75, eta, theta)

cauchy_dist <- rcauchy(n=10000, location=0, scale=2)

t.test(cauchy_dist, chain[2:10001])

###########################
### theta, eta sampling ###
###########################

prior <- function(eta, theta){
  f <- dnorm(eta, mean=0, sd=10, log=TRUE)
  f <- f+dunif(theta, min=0, max=10, log=TRUE)
  return(f)
}
likelihood <- function(eta, theta, x){
  ll <- 0
  for(i in x){
    ll <- ll+custom_dcauchy(i, eta, theta)
  }
  return(ll)
}
posterior <- function(eta, theta, x){
  return(likelihood(eta, theta, x)+prior(eta, theta))
}

mean(x)

eta_init <- 3
theta_init <- 1
chain <- list(eta=c(eta_init), theta=c(theta_init))
idx <- 0
while(length(chain$eta)<=1e3){
  idx <- idx + 1
  if(idx%%100 == 0){
    cat("Iteration",idx,"\n")
  }
  eta_prop <- eta_init + rnorm(1, 0, sd=0.2)
  thr <- posterior(eta_prop, theta_init, x)-posterior(eta_init, theta_init, x)
  if(!is.nan(thr)){
    if(log(runif(1))<thr){
      eta_init <- eta_prop
    }
  }
  
  theta_prop <- theta_init + rnorm(1, 0, sd=0.2)
  thr <- posterior(eta_init, theta_prop, x)-posterior(eta_init, theta_init, x)
  
  if(!is.nan(thr)){
    if(log(runif(1))<thr){
      theta_init <- theta_prop 
    }
    
  }
  if(idx>1e3 & idx%%2 == 0){ # Burn-in
    chain$eta <- c(chain$eta, eta_init)
    chain$theta <-c(chain$theta, theta_init)
  }
}

diagnosis(chain$theta, "theta")
diagnosis(chain$eta, "eta")

#######
