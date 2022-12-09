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
