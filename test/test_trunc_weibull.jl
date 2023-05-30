module TestTruncatedWeibull
using Distributions
using Test
using Random
using WeibullParetoDist
Random.seed!(10)
#include("../src/trunc_weibull.jl")

ϵ1 = 1e-7

α = 0.5
κ = 1.1
wb = WeibullPareto(α, κ)
#Test the declaration of Truncated Weibull distributions.
@test isa(truncated(wb, nothing, nothing), WeibullPareto)
@test isa(truncated(wb, nothing, 1.0), RightTruncatedWeibull)
@test isa(truncated(wb, 1.0, nothing), LeftTruncatedWeibull)
@test isa(truncated(wb, 1.0, 1.0), DoublyTruncatedWeibull)
@test isa(truncated(wb; lower=1.0, upper=1.0), DoublyTruncatedWeibull)

#Test truncated Weibull distribution.
α = 0.5
κ = 1.1
wb = WeibullPareto(α, κ)

# Test to match with Ordinal Weibull
wb_left = truncated(wb, 0.0, nothing)
wb_right = truncated(wb, nothing, Inf)
wb_doubly = truncated(wb, 0.0, Inf)

# Test mean, var
@test abs(mean(wb) - mean(wb_left)) < ϵ1
@test abs(mean(wb) - mean(wb_right)) < ϵ1
@test abs(mean(wb) - mean(wb_doubly)) < ϵ1

@test abs(var(wb) - var(wb_left)) < ϵ1
@test abs(var(wb) - var(wb_right)) < ϵ1
@test abs(var(wb) - var(wb_doubly)) < ϵ1

function test_two_scalers(y1, y2; ϵ2=0.01)
    dif = (y1 .- y2) .|> abs
    @test (dif .< ϵ2) |> all
end

# Test left-truncated weibull
for l in [1.0, 5.0, 10.0]
    x_ = l:100
    wb = WeibullPareto(α, κ)
    wb_trunc = truncated(wb; lower=l)

    # Test logpdf
    y = logpdf.(wb, x_)
    tp= ccdf(wb, l)
    y_trunc = logpdf.(wb_trunc, x_)
    test_two_scalers(y.-log(tp), y_trunc)

    # Test ccdf 
    y = ccdf.(wb, x_)
    tp = ccdf(wb, l)
    y_trunc = ccdf.(wb_trunc, x_)
    test_two_scalers(y./tp, y_trunc)

    # Test random variables
    rand_vars = rand(wb, 100_000_000)
    rand_vars = rand_vars[rand_vars .> l] 
    rand_trunc = rand(wb_trunc, 10_000_000)
    @test abs(mean(rand_vars)  - mean(rand_trunc))  < 0.1

    # Test mean, var, std
    @test abs(mean(rand_vars)  - mean(wb_trunc))  < 0.1
    @test abs(var(rand_vars) - var(wb_trunc)) < 0.1
    @test abs(std(rand_vars) - std(wb_trunc)) < 0.1
end 

# Test Right-truncated weibull
for u in [100.0, 50, 10]
    x_ = 1:u
    wb = WeibullPareto(α, κ)
    wb_trunc = truncated(wb; upper=u)

    # Test logpdf
    y = logpdf.(wb, x_)
    tp = cdf(wb, u)
    y_trunc = logpdf.(wb_trunc, x_)
    test_two_scalers(y.-log(tp), y_trunc)

    # Test cdf 
    y = cdf.(wb, x_)
    tp = cdf(wb, u)
    y_trunc = cdf.(wb_trunc, x_)
    test_two_scalers(y./tp, y_trunc)

    # Test ccdf 
    y = ccdf.(wb, x_) .- ccdf(wb, u)
    tp = cdf(wb, u)
    y_trunc = ccdf.(wb_trunc, x_)
    test_two_scalers(y./tp, y_trunc)

    # Test random variables
    rand_vars = rand(wb, 100_000_000)
    rand_vars = rand_vars[rand_vars .< u] 

    # Test mean, var, std
    @test abs(mean(rand_vars)  - mean(wb_trunc))  < 0.1
    @test abs(var(rand_vars) - var(wb_trunc)) < 0.1
    @test abs(std(rand_vars) - std(wb_trunc)) < 0.1
end 

# Test doubly-truncated weibull
for l in [1.0,3.0], u in [100.0, 10.0]
    x_ = l:u
    wb = WeibullPareto(α, κ)
    wb_trunc = truncated(wb, l, u)

    # Test logpdf
    y = logpdf.(wb, x_)
    tp = ccdf(wb, l) - ccdf(wb, u)
    y_trunc = logpdf.(wb_trunc, x_)
    test_two_scalers(y.-log(tp), y_trunc)

    # Test cdf 
    y = cdf.(wb, x_) .- cdf(wb, l)
    tp = cdf(wb, u) - cdf(wb, l)
    y_trunc = cdf.(wb_trunc, x_)
    test_two_scalers(y./tp, y_trunc)

    # Test ccdf 
    y = ccdf.(wb, x_) .- ccdf(wb, u)
    tp = cdf(wb, u) - cdf(wb, l)
    y_trunc = ccdf.(wb_trunc, x_)
    test_two_scalers(y./tp, y_trunc)

    # Test random variables
    rand_vars = rand(wb, 100_000_000)
    cond = (rand_vars .> l) .& (rand_vars .< u)
    rand_vars = rand_vars[cond] 

    # Test mean, var, std
    @test abs(mean(rand_vars)  - mean(wb_trunc))  < 0.1
    @test abs(var(rand_vars) - var(wb_trunc)) < 0.1
    @test abs(std(rand_vars) - std(wb_trunc)) < 0.1
end 

end