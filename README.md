# WeibullParetoDist

[![Build Status](https://github.com/toshiakiasakura/WeibullParetoDist.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/toshiakiasakura/WeibullParetoDist.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package enables you to implement the Weibull distribution parameterised to include the Pareto-approximated parameter.


The following parameterisation of the Weibull distribution is implemented in Distribution.jl package, and see the source file [here](https://github.com/JuliaStats/Distributions.jl/blob/master/src/univariate/continuous/weibull.jl).  
$f(x;\alpha, \theta) = \frac{\alpha}{\theta} \left( \frac{x}{\theta} \right)^{\alpha-1} e^{-(x/\theta)^\alpha}$

The Weibull distribution has a property to be power-law if $\alpha$ takes less than 1 value. Following the definition [Ivan et al., 2018](https://par.nsf.gov/servlets/purl/10202111) , the power-law distribution is power-law if the complementary cumulative function (that is a survival function), $\bar{F(k)}$, follows $\bar{F}(k)=l(k)k^{-\alpha}$.  Since the form of the complementary cumulative function for the Weibull distribution is $\bar{F(k)} = e^{-(x/\theta)^\alpha}$, which satisfies the definition of power-law. 

However, when the Weibull distribution is fitted to the data which follows the power-law, the order of estimates of $\theta$ becomes $10{^-3}$ and $10^{-8}$. Also, the strong correlation between $\alpha$ and $\theta$ is observed, which hampers the stable estimation. 

As [Akira et al., 2022](https://pubmed.ncbi.nlm.nih.gov/36137054/) proved in the supplementary material, let $\kappa$ be $\alpha/\theta^\alpha$, and replace $\kappa$ with $\alpha$. At the region of 0 in log-log plot of the distribution, $ \log \bar{F}(k)=-\kappa/\alpha e^{\alpha \log x} $, the slope becomes $-\kappa$, which does not vary so much, and empirically the correlation between $\alpha$ and $\kappa$ becomes very week. 

In this package, we implemented the above-parameterised version of the Weibull distribution. 


### Usage
```
α = 0.5
κ = 1.1
wb = WeibullPareto(α, κ)
x = 1:100
pdf.(wb, x)
cdf.(wb, x)
ccdf.(wb, x)
mean(wb)
var(wb)
sd(wb)
```

Also, truncated version is surpported.
```
α = 0.5
κ = 1.1
wb = WeibullPareto(α, κ)
l = 1
u = 10
wb_doubly = truncated(wb, l, u)
wb_left = truncated(wb, l, nothing)
wb_left = truncated(wb; lower=l)
wb_right = truncated(wb, nothing, u)
wb_right = truncated(wb; upper=u)
```

## Mathematical expression 
The relationship between $\kappa$ and $\theta$.   
$ \theta = (\frac{\alpha}{\kappa})^{1/\alpha} $ and $ \kappa=\frac{\alpha}{\theta^\alpha}$.

### Weibull distribution

$ pdf(x) = \kappa x^{\alpha-1}e^{-x^\alpha \kappa/\alpha}$     
$ cdf(x) = 1-e^{-x^\alpha\kappa/\alpha}$  

### Left truncated Weibull distribution
Let $l$ be the left truncation point.  

$ pdf(x) = \kappa x^{\alpha-1}e^{-(x^\alpha - l^\alpha)\kappa/\alpha} $  
$ cdf(x) = 1 - e^{-(x^\alpha - l^\alpha)\kappa/\alpha} $

### Doubly truncated Weibull distribution
Let $r$ be the right truncation point. 

$ pdf(x) = (\kappa x^{\alpha-1}e^{-x^\alpha\kappa/\alpha})/(
{e^{-l^\alpha \kappa/\alpha} - e^{-r^\alpha \kappa/\alpha}}) $  
$ cdf(x) = (1 - e^{-x^\alpha\kappa/\alpha})/({e^{-l^\alpha \kappa/\alpha} - e^{-r^\alpha \kappa/\alpha}}) $


