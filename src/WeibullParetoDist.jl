module WeibullParetoDist

using Distributions
using Random
using SpecialFunctions
export 
    WeibullPareto,
    params,
    minimum,
    maximum,
    logpdf,
    pdf,
    logccdf,
    ccdf,
    cdf,
    mean,
    var,
    std,
    truncated,
    LeftTruncatedWeibull,
    RightTruncatedWeibull,
    DoublyTruncatedWeibull
        
include("weibull_pareto.jl")
include("trunc_weibull.jl")

end
