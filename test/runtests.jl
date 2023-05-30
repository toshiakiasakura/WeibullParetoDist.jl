using WeibullParetoDist
using Test

@testset "WeibullParetoDist.jl" begin
    include("test_weibull_pareto.jl")
    include("test_trunc_weibull.jl")
end
