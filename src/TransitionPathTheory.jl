module TransitionPathTheory

using Graphs: SimpleDiGraph, is_strongly_connected, is_weakly_connected
using ArgCheck
using LinearAlgebra: I, normalize, eigvecs
using StatsBase: sample
using Random: seed!
using PrecompileTools: @compile_workload

include("main.jl")
export TPTProblem
export Stochasticity, Stochastic, SuperStochastic, NonStochastic
export Connectivity, StronglyConnected, WeaklyConnected, Disconnected
export TransitionMatrix

include("homogeneous.jl")
export HomogeneousTPTProblem
export 𝒫, 𝒜, ℬ, 𝒮, Ω, 𝒜_true, ℬ_true, 𝒞
export stationary_distribution, 𝒫_backwards, forward_committor, backward_committor
export 𝒮_plus, 𝒫_plus, remaining_time, hitting_location_distribution
export stationary_statistics, nonstationary_statistics

@compile_workload begin
    A = [1, 2, 3]
    B = [3, 4, 5]
    P = TransitionMatrix(10)

    tpt = HomogeneousTPTProblem(P, A, B)
    st_stats = stationary_statistics(tpt)
    ns_stats = nonstationary_statistics(tpt, 10)
end

end # module
