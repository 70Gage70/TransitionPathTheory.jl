import LinearAlgebra
using Graphs:is_strongly_connected, SimpleDiGraph

include("../statistics/homogeneous/committors.jl")

"""
    TPTHomog{T, U}

A container for all the fundamental TPT quantities. 
    
All advanced statistics can be calculated from this.
"""
struct TPTHomog{T<:Real, U<:Integer} 
    sets::TPTSets{U}
    P::Matrix{T}
    P_plus::Matrix{T}
    pi_stat::Vector{T}
    q_plus::Vector{T}
    q_minus::Vector{T}
end

"""
    TPTHomog(P, A, B)

Initialize a homogenous TPT calculation with the transition matrix `P`, source set `A` and target set `B`.

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

    # ensure P has the correct size
    @assert size(P, 1) >= 2 "P is too small." 
    @assert size(P, 1) == size(P, 2) "P must be square."

    # ensure entries of P are positive numbers
    bad_inds = findfirst(x->!(x>=0), P)
    if bad_inds !== nothing
        bad_entry = P[bad_inds]
        bad_inds = bad_inds.I
        error("Entries of P must be positive numbers. Found that P$(bad_inds) = $(bad_entry).")
    end

    # ensure that P is stochastic
    for i = 1:size(P, 1)
        @assert sum(P[i, :]) ≈ 1.0 "Row $i of P is not stochastic to within float tolerance."
    end

    # ensure P represents an ergodic Markov chain
    Padj::Matrix{Int64} = [P[i,j] != 0.0 ? 1 : 0 for i in 1:size(P, 1), j in 1:size(P, 1)]
    @assert is_strongly_connected(SimpleDiGraph(Padj)) "P must be strongly connected, i.e. there must be a path between every two states that P is the matrix of."

    # create sets
    if !isempty(intersect(A, B)) @info "A and B intersect; intersection states are avoided by TPT." end

    S = collect(1:size(P, 1))
    sets = TPTSets(S, A, B)

    # calculate fundamental TPT quantities
    pi_stat = abs.(LinearAlgebra.normalize(LinearAlgebra.eigvecs(transpose(P))[:,size(P)[1]], 1))
    P_minus = [(pi_stat[j]/pi_stat[i])*P[j,i] for i in sets.S, j in sets.S]
    qp = q_plus(P, sets.A_true, sets.B_true, sets.C) 
    qm = q_minus(P_minus, sets.A_true, sets.B_true, sets.C)

    # S_plus is the set of indices i such that
    # if i is in B_true, then i is in S_plus
    # if i is outside B_true, then i is in S_plus if both q_minus[i] > 0.0 and sum(P[i, j]*qp[j] for j in sets.S) > 0.0
    # intuitively: i is either in B_true, or could have come from A_true and is connected to a state that is reactive_connected to B
    S_plus = [i for i in sets.S if (i in sets.B_true) || (qm[i] > 0.0 && sum(P[i, j]*qp[j] for j in sets.S) > 0.0)]
    sets = TPTSets(sets.S, sets.A, sets.B, S_plus = S_plus)

    # P_plus is the reactive alalogue of P.
    P_plus = zeros(size(P))
    for i in sets.S_plus
        if i in sets.B_true
            P_plus[i, sets.B_true] = P[i, sets.B_true]/sum(P[i, k] for k in sets.B_true)
        else
            P_plus[i, sets.S_plus] =[P[i, j]*qp[j]/sum(P[i, k]*qp[k] for k in sets.S) for j in sets.S_plus]
        end
    end

    return TPTHomog(sets, P, P_plus, pi_stat, qp, qm)
end

function Base.show(io::IO, x::TPTHomog)
    print(io, "TPTHomog[S(")
    show(io, length(x.sets.S))
    print(io, "), A(")
    show(io, length(x.sets.A))
    print(io, "), B(")
    show(io, length(x.sets.B))
    print("), AB_int(")
    show(io, length(intersect(x.sets.A, x.sets.B)))
    print(")]")
end