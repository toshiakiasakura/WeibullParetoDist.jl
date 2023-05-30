module TestWeibullPareto
using Distributions
using Test
using Random
using WeibullParetoDist
Random.seed!(10)
#include("../src/weibull_pareto.jl")

ϵ = 1e-7

"""Test Weibull Pareto distribution.
"""
α = 0.5
κ = 1.1
θ = (α/κ)^(1/α)

wb = Weibull(α, θ)
wb_pareto = WeibullPareto(α, κ)



# Test log pdf.
x = 1:100
y = logpdf.(wb, x)
y_pareto = logpdf.(wb_pareto, x)
dif = (y .- y_pareto) .|> abs
@test (dif .< ϵ) |> all

# Test pdf
y = pdf.(wb, x)
y_pareto = pdf.(wb_pareto, x)
dif = (y .- y_pareto) .|> abs
@test (dif .< ϵ) |> all

# Test logccdf, ccdf, cdf
y = cdf.(wb, x)
y_pareto = cdf.(wb_pareto, x)
dif = (y .- y_pareto) .|> abs
@test (dif .< ϵ) |> all

# Test mean, var
@test mean(wb) == mean(wb_pareto) 
@test abs(var(wb) - var(wb_pareto)) < ϵ

# Test mean, var, sd by random variables
rand_vars = rand(wb_pareto, 10_000_000)
@test abs(mean(rand_vars)  - mean(wb_pareto))  < 0.01
@test abs(var(rand_vars) - var(wb_pareto)) < 0.01
@test abs(std(rand_vars) - std(wb_pareto)) < 0.01

end