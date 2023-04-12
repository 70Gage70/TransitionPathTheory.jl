using Test
using TransitionPathTheory, HDF5

include("test-functions.jl")
include("test-cases.jl")

@testset "tpt" begin
    for t in test_cases
        res, ftest = tpt_infinite(t[1]), t[2]
        @test tpt_test(ftest, t[1], res)
    end
end
