using Distributions
import Distributions.value_support
import Distributions.truncated

"""
    truncated(d0::WeibullPareto; [lower::Real], [upper:real])
    truncated(d0::WeibullPareto, lower::Real, upper::Real)

```julia
truncated(d0; lower=l)     # left-truncated to the interval [l, Inf)
truncated(d0; upper=u)     # right-truncated to the interval [0, u)
truncated(d0; lower=l, upper=u)     # doubly-truncated to the interval [l, u]
```


"""
function truncated end
function truncated(d::WeibullPareto; 
    lower::Union{Real, Nothing}=nothing, 
    upper::Union{Real, Nothing}=nothing,
    )
    return truncated(d, lower, upper)
end

truncated(d::WeibullPareto, ::Nothing, ::Nothing) = d
truncated(d::WeibullPareto, l::Real, ::Nothing) = LeftTruncatedWeibull(d, l)
function truncated(d::WeibullPareto, ::Nothing, u::Real)
    u_ccdf = ccdf(d, u)
    tp = 1 - u_ccdf
    logtp = log(tp)
    RightTruncatedWeibull(d, u, u_ccdf, tp, logtp)
end
function truncated(d::WeibullPareto, l::Real, u::Real) 
    l_ccdf = ccdf(d, l)
    u_ccdf = ccdf(d, u)
    tp = l_ccdf - u_ccdf
    logtp = log(tp)
    DoublyTruncatedWeibull(d, l, u, l_ccdf, u_ccdf, tp, logtp)
end

"""
    Truncated

Generic wrapper for a truncated distribution
"""
struct LeftTruncatedWeibull{T<:Real} <:ContinuousUnivariateDistribution 
    untruncated::WeibullPareto
    lower::T
    function LeftTruncatedWeibull(d::D, l::T) where {D<:WeibullPareto, T<:Real}
        new{T}(d, l)
    end
end
struct RightTruncatedWeibull{T<:Real} <: ContinuousUnivariateDistribution 
    untruncated::WeibullPareto
    upper::T
    u_ccdf::T
    tp::T
    logtp::T
    function RightTruncatedWeibull(d::D, u::T, u_ccdf::T, tp::T, logtp::T
        ) where {D<:WeibullPareto, T<:Real}
        new{T}(d, u, u_ccdf, tp, logtp)
    end
end
struct DoublyTruncatedWeibull{T<:Real} <: ContinuousUnivariateDistribution
    untruncated::WeibullPareto
    lower::T
    upper::T
    l_ccdf::T # complementary cdf of lower bound
    u_ccdf::T # complementary cdf of upper bound
    tp::T # The probability of the truncated part, i.e. l_ccdf - u_ccdf = u_cdf - l_cdf
    logtp::T # log(tp), i.e. log(l_ccdf-u_ccdf)
    function DoublyTruncatedWeibull(d::D, l::T, u::T, l_ccdf::T, u_ccdf::T, tp::T, logtp::T
        ) where {D<:WeibullPareto, T<:Real}
        new{T}(d, l, u, l_ccdf, u_ccdf, tp, logtp)
    end
end
convert_to_θ(α::Real, κ::Real) = (α/κ)^(1/α)


"""
    Left Truncated Weibull 

Functions  for a left-truncated weibull distribution
"""
params(d::LeftTruncatedWeibull) = (d.untruncated.α, d.untruncated.κ, d.lower)

function Distributions.logpdf(d::LeftTruncatedWeibull, x::Real)
    α,κ,l = params(d)
    log(κ) + (α-1)*log(x) - (x^α - l^α)*κ/α
end
Distributions.pdf(d::LeftTruncatedWeibull, x::Real) = exp(logpdf(d, x))

function Distributions.logccdf(d::LeftTruncatedWeibull, x::Real) 
    α,κ,l = params(d)
    -(x^α - l^α)*κ/α
end
Distributions.ccdf(d::LeftTruncatedWeibull, x::Real) = exp(logccdf(d, x))
Distributions.cdf(d::LeftTruncatedWeibull, x::Real) = 1 - ccdf(d, x) 

function Distributions.mean(d::LeftTruncatedWeibull)
    α, κ, l = params(d)
    θ =  convert_to_θ(α, κ)
    t1 = θ / exp(-l^α*κ/α)
    log_t2_1 = loggamma(1/α+1)
    log_t2_2 = loggamma(1/α+1) + log(gamma_inc(1/α+1, l^α*κ/α, 0)[1])
    return t1*(exp(log_t2_1) - exp(log_t2_2))
end

function Distributions.var(d::LeftTruncatedWeibull)
    α, κ, l = params(d)
    θ =  convert_to_θ(α, κ)
    t1 = θ^2 / exp(-l^α*κ/α)
    log_t2_1 = loggamma(2/α + 1)
    log_t2_2 = loggamma(2/α + 1) + log(gamma_inc(2/α + 1, l^α*κ/α,0)[1])
    t3 = θ^2 / exp(-2*l^α*κ/α)
    log_t4_1 = loggamma(1/α + 1)
    log_t4_2 = loggamma(1/α + 1) + log(gamma_inc(1/α + 1, l^α*κ/α,0)[1]) 
    return t1*(exp(log_t2_1) - exp(log_t2_2)) - t3*(exp(log_t4_1) - exp(log_t4_2))^2
end
Distributions.std(d::LeftTruncatedWeibull) = sqrt(var(d))

function Distributions.rand(rng::AbstractRNG, d::LeftTruncatedWeibull)
    α, κ, l = params(d)
    z = randexp(rng)
    return (z*α/κ + l^α)^(1/α)
end

"""
    Right Truncated Weibull 

Functions  for a left-truncated weibull distribution
"""
params(d::RightTruncatedWeibull) = (d.untruncated.α, d.untruncated.κ, d.upper)

Distributions.logpdf(d::RightTruncatedWeibull, x::Real) = logpdf(d.untruncated, x) - d.logtp
Distributions.pdf(d::RightTruncatedWeibull, x::Real) = exp(logpdf(d, x))
Distributions.cdf(d::RightTruncatedWeibull, x::Real) = cdf(d.untruncated, x)/d.tp
Distributions.ccdf(d::RightTruncatedWeibull, x::Real) = (ccdf(d.untruncated, x) - d.u_ccdf)/d.tp
Distributions.logccdf(d::RightTruncatedWeibull, x::Real) = log(ccdf(d, x)) 

function Distributions.mean(d::RightTruncatedWeibull)
    α, κ, u = params(d)
    θ =  convert_to_θ(α, κ)
    t1 = θ / d.tp
    log_t2 = loggamma(1/α+1) + log(gamma_inc(1/α+1, u^α*κ/α, 0)[1])
    return t1*exp(log_t2)
end

function Distributions.var(d::RightTruncatedWeibull)
    α, κ, u = params(d)
    θ =  convert_to_θ(α, κ)
    t1 = θ^2 / d.tp
    log_t2 = loggamma(2/α + 1) + log(gamma_inc(2/α + 1, u^α*κ/α,0)[1])
    t3 = θ^2 / d.tp^2
    log_t4 = loggamma(1/α + 1) + log(gamma_inc(1/α + 1, u^α*κ/α,0)[1]) 
    return t1*exp(log_t2) - t3*exp(2*log_t4)
end
Distributions.std(d::RightTruncatedWeibull) = sqrt(var(d))

"""
    Doubly Truncated Weibull 

Functions  for a left-truncated weibull distribution
"""
params(d::DoublyTruncatedWeibull) = (d.untruncated.α, d.untruncated.κ, d.lower, d.upper)

Distributions.logpdf(d::DoublyTruncatedWeibull, x::Real) = logpdf(d.untruncated, x) - d.logtp
Distributions.pdf(d::DoublyTruncatedWeibull, x::Real) = exp(logpdf(d, x))
Distributions.cdf(d::DoublyTruncatedWeibull, x::Real) = (cdf(d.untruncated, x) - (1 - d.l_ccdf))/d.tp
Distributions.ccdf(d::DoublyTruncatedWeibull, x::Real) = (ccdf(d.untruncated, x) - d.u_ccdf)/d.tp
Distributions.logccdf(d::DoublyTruncatedWeibull, x::Real) = log(ccdf(d, x))

function Distributions.mean(d::DoublyTruncatedWeibull)
    α, κ, l, u = params(d)
    θ =  (α/κ)^(1/α)
    t1 = θ/(d.l_ccdf - d.u_ccdf)
    log_t2_1 = loggamma(1/α+1) + log(gamma_inc(1/α+1, u^α*κ/α, 0)[1])
    log_t2_2 = loggamma(1/α+1) + log(gamma_inc(1/α+1, l^α*κ/α, 0)[1])
    return t1*(exp(log_t2_1) - exp(log_t2_2))
end

function Distributions.var(d::DoublyTruncatedWeibull)
    α, κ, l, u = params(d)
    θ =  convert_to_θ(α, κ)
    t1 = θ^2 / (d.l_ccdf - d.u_ccdf)
    log_t2_1 = loggamma(2/α + 1) + log(gamma_inc(2/α + 1, u^α*κ/α,0)[1])
    log_t2_2 = loggamma(2/α + 1) + log(gamma_inc(2/α + 1, l^α*κ/α,0)[1])
    t3 = θ^2 / (d.l_ccdf - d.u_ccdf)^2
    log_t4_1 = loggamma(1/α + 1) + log(gamma_inc(1/α + 1, u^α*κ/α,0)[1]) 
    log_t4_2 = loggamma(1/α + 1) + log(gamma_inc(1/α + 1, l^α*κ/α,0)[1]) 
    return t1*(exp(log_t2_1) - exp(log_t2_2)) - t3*(exp(log_t4_1) - exp(log_t4_2))^2
end
Distributions.std(d::DoublyTruncatedWeibull) = sqrt(var(d))
