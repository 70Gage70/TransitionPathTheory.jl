using Test, TransitionPathTheory
import Random
Random.seed!(1234)

@testset "Homogeneous" begin
    A = [1, 2, 3]
    B = [3, 4, 5]
    P = TransitionMatrix(10)
    
    tpt = HomogeneousTPTProblem(P, A, B)
    st_stats = stationary_statistics(tpt)
    ns_stats = nonstationary_statistics(tpt, 10, B_to_S = :interior)
    ns_stats = nonstationary_statistics(tpt, 10, B_to_S = :uniform)
    ns_stats = nonstationary_statistics(tpt, 10, B_to_S = :balanced)
end

