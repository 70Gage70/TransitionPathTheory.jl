using Test
using TransitionPathTheory, HDF5

include("test-functions.jl")
include("test-cases.jl")

test_cases_stationary, test_cases_nonstationary = test_cases()

@testset "tpt_homog_stationary" begin
    for t in test_cases_stationary
        res, ftest = tpt_stationary_statistics(t[1]), t[2]
        @test tpt_test(ftest, res)
    end
end

@testset "tpt_homog_nonstationary" begin
    for t in test_cases_nonstationary
        res, ftest = tpt_nonstationary_statistics(t[1]), t[2]
        @test tpt_test(ftest, res)
    end
end
