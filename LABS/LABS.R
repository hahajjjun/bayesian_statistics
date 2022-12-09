###########################################
## LAB1: EXAMPLE of parameter estimation ##
###########################################

mu <- 5
ss <- 10
X <- rnorm(1000, mean=mu, sd=sqrt(ss))
# PRIOR:
# mu ~ N(0,100)
# ss ~ IG(0.1, 0.1)
# TRUE:
# mu = 5, ss = 10

a <- 0
b <- 100
c <- 1
d <- 10


sampler_mu <- function(X, ss){
  mu_mean <- (b/(b+ss/length(X)))*mean(X)+(1-(b/(b+ss/length(X))))*a
  mu_var <- 1/(length(X)/ss + 1/b)
  return(rnorm(1, mu_mean, sqrt(mu_var)))
}

sampler_ss <- function(X, mu){
  ss_shape <- c+length(X)/2
  ss_scale <- d+sum((X-mu)^2)/2
  return(rinvgamma(1, shape=ss_shape, scale=ss_scale))
}

mu_chain <- c(1)
ss_chain <- c(1)
for(i in 1:1000){
  mu_prop <- sampler_mu(X, ss_chain[i])
  mu_chain <- c(mu_chain, mu_prop)
  ss_prop <- sampler_ss(X, mu_chain[i+1])
  ss_chain <- c(ss_chain, ss_prop)
}

############################################
### LAB2: EXAMPLE OF BAYESIAN REGRESSION ###
############################################

### GIBBS SAMPLING ###
mu <- 0
tau <- 100
a <- 0.01
b <- 0.01

x <- rnorm(1000, 0, 1)
beta0 = 0.5
beta1 = 1
ss = 1
Y <- rep(0, 1000)
for(i in 1:1000){
  Y[i] <- rnorm(1, beta0+beta1*x[i], sqrt(ss))
}

sample_beta0 <- function(beta1, ss, x, Y){
  M <- (sum(Y-beta1*x)/ss + mu/tau)/(1/tau + length(x)/ss)
  V <- 1/((1/tau + length(x)/ss))
  return(rnorm(1, M, V))
}

sample_beta1 <- function(beta0, ss, x, Y){
  M <- (sum((Y-beta0)*x)/ss + mu/tau)/(sum(x^2)/ss + 1/tau)
  V <- 1/(sum(x^2)/ss + 1/tau)
  return(rnorm(1, M, V))
}

sample_ss <- function(beta0, beta1, x, Y){
  shape <- a+length(x)/2
  scale <- b+sum((Y-beta0-beta1*x)^2)/2
  return(rinvgamma(1, shape, scale))
}

beta0_chain <- c(1)
beta1_chain <- c(1)
ss_chain <- c(1)
for(i in 1:1000){
  beta0_prop <- sample_beta0(beta1_chain[i], ss_chain[i], x, Y)
  beta0_chain <- c(beta0_chain, beta0_prop)
  beta1_prop <- sample_beta1(beta0_chain[i+1], ss_chain[i], x, Y)
  beta1_chain <- c(beta1_chain, beta1_prop)
  ss_prop <- sample_ss(beta0_chain[i+1], beta1_chain[i+1], x, Y)
  ss_chain <- c(ss_chain, ss_prop)
}

##############################################
### Nested model with model comparison: BF ###
##############################################

X <- rnorm(1000, 0, 1)
beta <- c(0.5, 1)
mean <- beta[1]+beta[2]*X
y <- rep(0, 1000)
for(i in 1:1000){
  y[i] <- rnorm(1, mean[i], sd=0.1)
}

beta_init <- c(1, 1, 1)
gamma_init <- 1

# Prior of gamma : Bernoulli(0.5)
# Prior of beta : beta_i ~ N(0, 100)
# Prior of ss : ss ~ G(1, 20)

log_ss_prior <- function(ss){
  return(dgamma(ss, 1, 20, log=TRUE))
}
log_beta_prior <- function(beta){
  lp <- 0
  for(i in 1:length(beta)){
    lp <- lp + dnorm(beta[i], mean=0, sd=10, log=TRUE)
  }
  return(lp)
}
log_likelihood <- function(beta, gamma, ss, X, y){
  ll <- 0
  for(i in 1:length(X)){
    mean <- beta[1]+beta[2]*X[i]+beta[3]*(X[i]^2)*gamma
    ll <- ll + dnorm(y[i], mean, sqrt(ss), log=TRUE)
  }
  return(ll)
}
log_beta_posterior <- function(beta, gamma, ss, X, y){
  return(log_beta_prior(beta)+log_likelihood(beta, gamma, ss, X, y))
}
log_ss_posterior <- function(beta, gamma, ss, X, y){
  return(log_ss_prior(ss)+log_likelihood(beta, gamma, ss, X, y))
}


beta_chain <- beta_init <- c(1, 1, 1)
gamma_chain <- gamma_init <- 1
ss_chain <- ss_init <- 1

for(i in 1:1000){
  for(j in 1:3){
    beta_prop <- beta_init
    beta_prop[j] <- beta_init[j]+rnorm(1, 0, sd=0.1)
    thr <- log_beta_posterior(beta_prop, gamma_init, ss_init, X, y)-log_beta_posterior(beta_init, gamma_init, ss_init, X, y)
    if(!is.nan(thr)){
      if(log(runif(1))<thr){
        beta_init <- beta_prop
      }
    }
  }
  gamma_prop <- rbinom(1, 1, 0.5)
  thr <- log_likelihood(beta_init, gamma_prop, ss_init, X,y)-log_likelihood(beta_init, gamma_init, ss_init, X, y)
  if(!is.nan(thr)){
    if(log(runif(1))<thr){
      gamma_init <- gamma_prop
    }
  }
  ss_prop <- ss_init+rnorm(1, 0, sd=0.1)
  thr <- log_ss_posterior(beta_init, gamma_init, ss_prop, X,y)-log_ss_posterior(beta_init, gamma_init, ss_init, X, y)
  if(!is.nan(thr)){
    if(log(runif(1))<thr){
      ss_init <- ss_prop
    }
  }
  beta_chain = cbind(beta_chain, beta_init)
  gamma_chain = c(gamma_chain, gamma_init)
  ss_chain = c(ss_chain, ss_init)
}
