# prior theta ~ G(a,b)
# likelihood X|theta ~ Poi(theta)
# posterior theta|X ~ G(sum(X)+a, n+b) whre a and b are shape & rate parameters

# G(5,scale = 6) (prior) --> Implemented as G(5, rate = 1/6)
# X = 52 (likelihood) --> Implemented as G(53, rate = 1)
# G(57, rate = 7/6) (posterior) 

prior_theta <- rgamma(n=1000, shape=5, rate=1/6)
posterior_theta <- rgamma(n=1000, shape= 57, rate=7/6)
likelihood_theta <- rgamma(n=1000, shape= 53, rate=1)
plot(density(prior_theta), col="black", ylim=c(0, 0.1))
lines(density(likelihood_theta), col="blue")
lines(density(posterior_theta), col="red")