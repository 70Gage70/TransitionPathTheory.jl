import LinearAlgebra
import Graphs

include("../statistics/infinite/committors.jl")

"""
    TPTHomog{T, U}

A container for all the fundamental TPT quantities. 
    
All advanced statistics can be calculated from this.
"""
struct TPTHomog{T<:Real, U<:Integer} 
    sets::TPTSets{U}
    P::Matrix{T}
    pi_stat::Vector{T}
    q_plus::Vector{T}
    q_minus::Vector{T}
end

"""
    TPTHomog(P, A, B)

Initialize a homogenous TPT calculation with the transition matrix `P` and state subsets `A` and `B`.

Fundamental TPT statistics are computed automatically such as `q_plus` and `q_minus`.

### Constraints
The following constraints are imposed on `P`:
- It must be a square matrix of size at least 2.
- Its entries must be positive.
- It must be stochastic or substochastic.
- It must be strongly connected.

The following constraints are imposed on `A` and `B`:
- They must be non-empty vectors of positive elements.
- They can not be subsets of each other.
- They CAN intersect, but TPT avoids these indices (this is a feature.)
"""
function TPTHomog(
    P::Matrix{T},
    A::Vector{U},
    B::Vector{U}) where {T<:Real, U<:Integer}

    @assert size(P, 1) >= 2 "P is too small." 
    @assert size(P, 1) == size(P, 2) "P must be square."
    @assert all(P .>= 0.0) "Entries of P can not be negative."

    # ensure that P is at least sub-stochastic
    substoc = false
    for i = 1:1:size(P, 1)
        rowsum = sum(P[i, :])

        if !(rowsum ≈ 1.0)
            @assert rowsum < 1.0 "row $i of P has a sum greater than 1."
            substoc = true
        end
    end

    if substoc @info "P is sub-stochastic." end

    # ensure P represents an ergodic Markov chain
    Padj::Matrix{Int64} = [P[i,j] != 0.0 ? 1 : 0 for i in 1:size(P, 1), j in 1:size(P, 1)]
    @assert Graphs.is_strongly_connected(Graphs.SimpleDiGraph(Padj)) "P must be strongly connected, i.e. there must be a path between every two states that P is the matrix of."

    # create sets
    if !isempty(intersect(A, B)) @info "A and B intersect; intersection states are avoided by TPT." end

    S = collect(1:size(P, 1))
    sets = TPTSets(S, A, B)

    # calculate fundamental TPT quantities
    pi_stat = abs.(LinearAlgebra.normalize(LinearAlgebra.eigvecs(transpose(P))[:,size(P)[1]], 1))
    P_minus = [(pi_stat[j]/pi_stat[i])*P[j,i] for i in sets.S, j in sets.S]
    qp = q_plus(P, sets.A_true, sets.B_true, sets.C) 
    qm = q_minus(P_minus, sets.A_true, sets.B_true, sets.C)

    # C_plus is the set of indices i outside of B such that there exists a j in S
    # such that i can reach j in one step and j is reactive-connected to B
    C_plus = [i for i in setdiff(sets.S, sets.B) if sum(P[i, j]*qp[j] for j in sets.S) > 0.0] 
    @assert !isempty(C_plus) "C_plus can not be empty. This happens if there is no reactive path from A to B."
    sets = TPTSets(sets.S, sets.A, sets.B, C_plus = C_plus)

    return TPTHomog(sets, P, pi_stat, qp, qm)
end