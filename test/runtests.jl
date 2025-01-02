using Test, TransitionPathTheory
import Random
Random.seed!(1234)

@testset "Homogeneous" begin
    A = [1, 2, 3]
    B = [3, 4, 5]
    P = TransitionMatrix(10)
    
    tpt = HomogeneousTPTProblem(P, A, B)
    st_stats = stationary_statistics(tpt)
    ns_stats = nonstationary_statistics(tpt, 10, initial_dist = :uniform)
    ns_stats = nonstationary_statistics(tpt, 10, initial_dist = :stat)

    f = st_stats.forward_current
    current2arrows(f, rand(size(f, 1) - 1), rand(size(f, 1) - 1))
end

