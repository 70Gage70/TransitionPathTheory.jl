using Test
using TransitionPathTheory, HDF5

include("test-functions.jl")
include("test-cases.jl")

test_cases_stationary, test_cases_nonstationary = test_cases()

@testset "tpt_homog_stationary" begin
    for t in test_cases_stationary
        tpt_stat, ftest = tpt_stationary_statistics(t[1]), t[2]
        parts = minimal_partitions(t[1], tpt_stat)
        @test tpt_test(ftest, tpt_stat)
        @test parts_test(ftest, parts)
    end
end

@testset "tpt_homog_nonstationary" begin
    for t in test_cases_nonstationary
        tpt_nonstat, ftest = tpt_nonstationary_statistics(t[1]), t[2]
        parts = minimal_partitions(t[1], tpt_nonstat)
        @test tpt_test(ftest, tpt_nonstat)
        @test parts_test(ftest, parts)
    end
end
