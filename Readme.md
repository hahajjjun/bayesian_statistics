# Practical implementation details

### Notation
- G(a,b) in lecture note indicates distribution with **shape=a, rate=b**.
- IG(a,b) in lecture note indicates distribution with **shape=a, scale=b**.

- Kernel desnity: $G(a,b): x^{a-1}\exp(-bx), IG(a,b): x^{-a-1}\exp(-b/x)$
- Realization of such distribution in *R*
    - `G(a, rate=b)`: R base function
    - `IG(a, scale=b)`: R nimble package function
    - `IG(a, rate=b)`: R invgamma package function
- Conjugacy relationship
    - Always based on setting **b=rate**
    
 ### Conjugacy relationships
 |Prior|Likelihood|Posterior|
 |:-:|:-:|:-:|
 |$\theta \sim N(a,b)$|$X\mid\theta \sim N(\theta, \sigma^2)$| $\theta\mid X \sim N(\frac{b}{b+\frac{\sigma^2}{n}}\bar{X}+(1-\frac{b}{b+\frac{\sigma^2}{n}})a, (\frac{n}{\sigma^2}+\frac{1}{b})^{-1})$|
 |$\theta \sim IG(a, rate=b)$|$X\mid\theta \sim N(\mu, \theta)$|$\theta\mid X \sim IG(a+\frac{1}{2}, rate=b+\frac{(X-\mu)^2}{2})$|
 |$\theta \sim G(a, rate=b)$|$X\mid\theta \sim Poisson(\theta)$|$\theta\mid X \sim G(a+\sum{X_i},b+n)$|
 |$\theta \sim Beta(a, b)$|$X\mid\theta \sim Binomial(n, \theta)$|$\theta\mid X \sim Beta(a+X, b+n-X)$|
 
 ### Kernel of distributions and statistics
 |Implementation|Kernel|Mean|
 |:-:|:-:|:-:|
 |`dbinom(x, n, prob)`|$nCx\cdot p^x(1-p)^{n-x}$|
 |`dbernoulli(x, prob)`|$p^x(1-p)^{n-x}$|
 |`dpois(x, lambda)`|$\lambda^x\cdot e^{-\lambda} \over x!
 |`dbeta(x, a, b)`|$x^{a-1}(1-x)^{b-1}$|$\frac{a}{a+b}$|
 |IG(a, b) > `dinvgamma(x, shape=a, rate=b)` - rinvgamma package|$x^{-a-1}e^{-b/x}$|$\frac{b}{a-1}$|
 |IG(a, b) > `dinvgamma(x, shape=a, scale=b)` - nimble package|$x^{-a-1}e^{-b/x}$|$\frac{b}{a-1}$|
 |G(a, b) > `dgamma(x, shape=a, rate=b)`|$x^{a-1}e^{-bx}$|$\frac{a}{b}$ ( $ab$ if b is scale parameter)|
