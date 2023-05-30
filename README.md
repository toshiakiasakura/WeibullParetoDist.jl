# WeibullParetoDist

[![Build Status](https://github.com/toshiakiasakura/WeibullParetoDist.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/toshiakiasakura/WeibullParetoDist.jl/actions/workflows/CI.yml?query=branch%3Amain)


This package enables you to implement the Pareto approximated version of the Weibull distribution. 
The scale parameter, θ, is replaced with κ (Pareto approximated) using the shape parameter, α, as θ = (α/κ)^(1/α).

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
wb_doubly = trunc(wb, l, u)
wb_left = trunc(wb, l, nothing)
wb_left = trunc(wb; lower=l)
wb_right = trunc(wb, nothing, u)
wb_right = trunc(wb; upper=u)
```