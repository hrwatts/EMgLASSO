# EMgLASSO
## Simulate a precision matrix using graphical LASSO and EM

Consider a d-dimensional mixture distribution with k mixture components associated with
the parameter vector Θ = (τ, µ, Σ),
representing respectively the mixture probability, the mean vector, and the covariance matrix
for each component in the MMN distribution. For a random sample X
0 = (X
0
1
, ..., X
0
n
) of size
n drawn independently from the population, the unobserved data are the random variables
Z
0 = (Z1, ..., Zn) that determine with probability τj from which mixture the observation
originates:
Xi
|(Zi = j) ∼ Nd(µj
, Σj), P(Zi = j) = τj
, j = 1, ..., k,
subject to X
k
j=1
τj = 1.
These variables constitute a random sample on the multinomial random variable Z. Therefore Z1, ..., Zn are independently and identically distributed multinomial random variables
with probabilities τ. Let x
0 = (x
0
1
, ..., x
0
n
) be the observed values of X given by the sample.
And z
0 = (z1, ...,zn) are the latent values associated with the realization of the random sample. The observed likelihood function is L(Θ|x) =
Qn
i=1
Pk
j=1
τj
f(xi
; µjΣj), where where f
Date: September 2018.
1
2 HARRISON WATTS
is the pdf of the jth multivariate normal distribution in the mixture and is given by
f(x; µ, Σ) =
exp(−
1
2
(x − µ)
0Σ
−1
(x − µ))
p
(2π)
d
|Σ|
,
where |Σ| is the determinant of the covariance matrix Σ. The complete likelihood function
is the product
L
c
(Θ; x, z) =
Yn
i=1
Y
k
j=1
τj
f(xn; µj
, Σj)
1(z
i=j)
,
where 1 is the indicator function. The product is used with an indicator function in the
exponent because the latent values are given in the complete likelihood function, which
is unknown in practice. 
