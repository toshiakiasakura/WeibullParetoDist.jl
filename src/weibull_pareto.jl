using Distributions
using Random
using SpecialFunctions

"""Weibull distribution with kappa approximation.
"""
struct WeibullPareto{T<:Real} <: ContinuousUnivariateDistribution 
    α::T # shape
    κ::T # Pareto parameter
    function WeibullPareto{T}(α::T, κ::T) where {T <: Real}
        new{T}(α, κ)
    end
end

WeibullPareto(α::T, κ::T) where {T <: Real} = WeibullPareto{T}(α, κ)

# get parameters
params(d::WeibullPareto) = (d.α, d.κ)

Distributions.minimum(d::WeibullPareto) = 0.
Distributions.maximum(d::WeibullPareto) = Inf

Distributions.logpdf(d::WeibullPareto, x::Real) = log(d.κ) + (d.α-1)*log(x) - x^d.α *d.κ/d.α
Distributions.pdf(d::WeibullPareto, x::Real) = exp(logpdf(d, x))
Distributions.logccdf(d::WeibullPareto, x::Real) = -x^d.α*d.κ/d.α
Distributions.ccdf(d::WeibullPareto, x::Real) = exp(logccdf(d, x))
Distributions.cdf(d::WeibullPareto, x::Real) = 1 - ccdf(d, x)

function Distributions.mean(d::WeibullPareto)
    α, κ = params(d)
    θ = (α/κ)^(1/α)
    return θ * gamma(1+1/α)
end

function Distributions.var(d::WeibullPareto)
    α, κ = params(d)
    θ = (α/κ)^(1/α)
    return θ^2 * (gamma(1+2/α) - gamma(1+1/α)^2)
end
Distributions.std(d::WeibullPareto) = sqrt(var(d))

function Distributions.rand(rng::AbstractRNG, d::WeibullPareto)
    α, κ = params(d)
    z = randexp(rng)
    return (α/κ)^(1/α)  * z^(1/α)
end
