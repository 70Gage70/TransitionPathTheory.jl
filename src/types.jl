import LinearAlgebra

include("statistics/infinite/committors.jl")

# abstract type AbstractTPTStatistic end

struct TPTSets{U<:Integer}
    S::Vector{U}
    A::Vector{U}
    B::Vector{U}
    C::Vector{U}
    A_true::Vector{U}
    B_true::Vector{U}
    C_plus::Union{Vector{U}, Nothing}
end

function TPTSets(
    S::Vector{U}, 
    A::Vector{U}, 
    B::Vector{U};
    C_plus::Union{Vector{T}, Nothing} = nothing) where {T<:Real, U<:Integer}

    # Scrub any accidental repeated indices
    S = unique(S)
    A = unique(A)
    B = unique(B)    

    # ensure sets are well-defined
    @assert length(A) > 0 "A can not be empty."
    @assert length(B) > 0 "B can not be empty."
    @assert minimum(A) > 0 "Elements of A must be positive."
    @assert minimum(B) > 0 "Elements of B must be positive."
    @assert minimum(S) > 0 "Elements of S must be positive."
    @assert issubset(A, S) "A must be a subset of S."
    @assert issubset(B, S) "B must be a subset of S."
    @assert !issubset(A, B) "A can not be a subset of B."
    @assert !issubset(B, A) "B can not be a subset of A."    

    C = setdiff(S, union(A, B))

    # the "true" A and B are indices unique to A and B
    # where they intersect is avoided by TPT
    A_true = setdiff(A, intersect(A, B))
    B_true = setdiff(B, intersect(A, B))

    return TPTSets(S, A, B, C, A_true, B_true, C_plus)
end

########################################################################
########################################################################
########################################################################

struct TPTHomog{T<:Real, U<:Integer} 
    sets::TPTSets{U}
    P::Matrix{T}
    pi_stat::Vector{T}
    q_plus::Vector{T}
    q_minus::Vector{T}
end

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

    # create sets
    S = collect(1:size(P, 1))
    sets = TPTSets(S, A, B)

    # calculate fundamental TPT quantities
    pi_stat = abs.(LinearAlgebra.normalize(LinearAlgebra.eigvecs(transpose(P))[:,size(P)[1]], 1))
    P_minus = [(pi_stat[j]/pi_stat[i])*P[j,i] for i in sets.S, j in sets.S]
    qp = q_plus(P, sets.A_true, sets.B_true, sets.C) 
    qm = q_minus(P_minus, sets.A_true, sets.B_true, sets.C)

    # C_plus is the set of indices i outside of B such that there exists a j in S
    # such that i can reach j in one step and j is reactive-connected to B
    C_plus = [i for i in setdiff(S, B) if sum(P[i, j]*qp[j] for j in S) > 0.0] 
    sets = TPTSets(S, A, B, C_plus = C_plus)

    return TPTHomog(sets, P, pi_stat, qp, qm)
end

########################################################################
########################################################################
########################################################################

struct TPTStats{T<:Real, U<:Integer}
    reactive_density::Vector{T}
    normalized_reactive_density::Vector{T}
    reactive_current::Matrix{T}
    forward_current::Matrix{T}
    reactive_rate_vanE::T
    transition_time_vanE::T
    remaining_time::Vector{T}
    hitting_distribution::Matrix{T}
    hitting_locations::Vector{U}
end


########################################################################
########################################################################
########################################################################

