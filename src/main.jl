import Pkg
Pkg.activate(".")

using Graphs: SimpleDiGraph, is_strongly_connected, is_weakly_connected
using ArgCheck
using LinearAlgebra: I, normalize, eigvecs
using StatsBase: sample

abstract type Stochasticity end

struct SubStochastic <: Stochasticity end
struct Stochastic <: Stochasticity end
struct SuperStochastic <: Stochasticity end 
struct NonStochastic <: Stochasticity end 

abstract type Connectivity end 

struct StronglyConnected <: Connectivity end 
struct WeaklyConnected <: Connectivity end 
struct Disconnected <: Connectivity end 

struct TransitionMatrix{S<:Stochasticity, C<:Connectivity}
    P::Matrix{Float64}
    stochasticity::S
    connectivity::C
end

function TransitionMatrix(P::Matrix{<:Real})
    # ensure P has the correct size
    @argcheck size(P, 1) == size(P, 2) "P must be square."
    @argcheck size(P, 1) >= 2 "P is too small."
    n_states = size(P, 1) 
    
    # ensure entries of P are positive numbers
    bad_inds = findfirst(x -> isnan(x) || x < 0 , P)
    if bad_inds !== nothing
        bad_entry = P[bad_inds]
        bad_inds = bad_inds.I
        error("Entries of P must be positive numbers. Found that P$(bad_inds) = $(bad_entry).")
    end

    # decide stochasticity
    row_sums = [sum(P[i, :]) for i = 1:size(P, 1)]
    if all(row_sums .≈ 1.0)
        stochasticity = Stochastic()
    elseif all(row_sums .<= 1.0)
        stochasticity = SubStochastic()
    elseif all(row_sums .>= 1.0)
        stochasticity = SuperStochastic()
    else
        stochasticity = NonStochastic()
    end

    # decide connectivity
    P_con = [iszero(P[i,j]) ? 0 : 1 for i in 1:n_states, j in 1:n_states] |> SimpleDiGraph
    if is_strongly_connected(P_con)
        connectivity = StronglyConnected()
    elseif is_weakly_connected(P_con)
        connectivity = WeaklyConnected()
    else
        connectivity = Disconnected()
    end

    return TransitionMatrix(P, stochasticity, connectivity)
end

function TransitionMatrix(P_size::Integer; n_zeros::Integer = 0, normalize::Bool = true)
    @argcheck P_size >= 2
    @argcheck n_zeros >= 0 

    P = rand(P_size, P_size) 
    
    if n_zeros > 0
        P[sample(eachindex(P), n_zeros, replace = false)] .= 0
    end

    row_sums = [sum(P[i, :]) for i = 1:P_size]
    for i = 1:P_size
        if iszero(row_sums[i])
            P[i, rand(1:P_size)] = 1
        else
            if normalize
                P[i, :] .= P[i, :]/row_sums[i]
            end
        end
    end

    return TransitionMatrix(P)
end

abstract type TPTProblem end

struct HomogeneousTPTProblem{TM<:TransitionMatrix} <: TPTProblem
    P::TM
    source::Vector{Int64}
    target::Vector{Int64}
end

function HomogeneousTPTProblem(
    P::TransitionMatrix, 
    source::Vector{<:Integer}, 
    target::Vector{<:Integer}; 
    avoid::Vector{<:Integer} = Int64[])

    # ensure no repeated indices
    @argcheck allunique(source)
    @argcheck allunique(target)
    @argcheck allunique(avoid)

    # ensure indices take appropriate values
    n_states = size(P.P, 1)
    @argcheck all(1 .<= source .<= n_states)
    @argcheck all(1 .<= target .<= n_states)
    @argcheck all(1 .<= avoid .<= n_states)

    source = [source ; avoid] |> sort
    target = [target ; avoid] |> sort

    return HomogeneousTPTProblem(P, source, target)
end

𝒫(tpt::HomogeneousTPTProblem) = tpt.P.P                          
𝒜(tpt::HomogeneousTPTProblem) = tpt.source                        
ℬ(tpt::HomogeneousTPTProblem) = tpt.target                        
𝒮(tpt::HomogeneousTPTProblem) = collect(1:size(𝒫(tpt),1))    
Ω(tpt::HomogeneousTPTProblem) = intersect(𝒜(tpt), ℬ(tpt))          
𝒜_true(tpt::HomogeneousTPTProblem) = setdiff(𝒜(tpt), Ω(tpt))      
ℬ_true(tpt::HomogeneousTPTProblem) = setdiff(ℬ(tpt), Ω(tpt))      
𝒞(tpt::HomogeneousTPTProblem) = setdiff(𝒮(tpt), union(𝒜(tpt), ℬ(tpt)))

function stationary_distribution(tpt::HomogeneousTPTProblem)
    P = 𝒫(tpt)

    if tpt.P.connectivity === StronglyConnected()
        return normalize(abs.(eigvecs(P')[:,end]), 1)
    else
        error("Transition matrix is not strongly connected, so the stationary distribution is not unique.")
    end
end

function 𝒫_backwards(tpt::HomogeneousTPTProblem)
    P, S = 𝒫(tpt), 𝒮(tpt)
    pi_stat = stationary_distribution(tpt)

    return [(pi_stat[j]/pi_stat[i])*P[j,i] for i in S, j in S]
end

function forward_committor(tpt::HomogeneousTPTProblem)
    P, A, B, C = 𝒫(tpt), 𝒜(tpt), ℬ(tpt), 𝒞(tpt)

    # solve the linear algebra problem q = P q + b on the C subspace
    M = I - P[C, C] 
    b = [sum(P[i, k] for k in B) for i in C]
    sol = M\b
    
    # assign the values to a vector
    q = zeros(size(P, 1))
    q[B] .= 1.0
    q[A] .= 0.0 # handle A after B since intersection states need to be 0
    q[C] = sol
    
    # trim very small (possibly negative due to floats) values
    q[abs.(q) .< 1e-16] .= 0.0
 
    return q
end

function backward_committor(tpt::HomogeneousTPTProblem)
    P_back, A, B, C = 𝒫_backwards(tpt), 𝒜(tpt), ℬ(tpt), 𝒞(tpt)

    # The calculation is identical to q_plus, except we use P_back instead of P 
    # and switch the roles of A and B
    M = I - P_back[C, C] 
    b = [sum(P_back[i, k] for k in A) for i in C]
    sol = M\b
    
    # assign the values to a vector; note that q[B] = 0.0 is already handled
    q = zeros(size(P_back, 1))
    q[A] .= 1.0
    q[B] .= 0.0 # handle B after A since intersection states need to be 0
    q[C] = sol
    
    # trim very small (possibly negative due to floats) values
    q[abs.(q) .< 1e-16] .= 0.0

    return q      
end

"""
    𝒮_plus(tpt)

The set of indices `i` such that
- if `i` is in `B_true`, then `i` is in `S_plus`
- if `i` is outside `B_true`, then `i` is in `S_plus` if both 
 - `q_minus[i] > 0.0 `
 - `sum(P[i, j]*qp[j] for j in sets.S) > 0.0`

Intuitively: `i` is either in `B_true`, or could have come from `A_true` and is connected to a state that is reactively connected to `B`.
"""
function 𝒮_plus(tpt::HomogeneousTPTProblem)
    S, B_true = 𝒮(tpt), ℬ_true(tpt)
    P = 𝒫(tpt)
    qp, qm = forward_committor(tpt), backward_committor(tpt)

    return [i for i in S if (i in B_true) || (qm[i] > 0.0 && sum(P[i, j]*qp[j] for j in S) > 0.0)]
end

"""
    𝒫_plus(tpt)

The forward-reactive analogue of `𝒫`.
"""
function 𝒫_plus(tpt::HomogeneousTPTProblem)
    S, S_plus, B_true = 𝒮(tpt), 𝒮_plus(tpt), ℬ_true(tpt)
    P = 𝒫(tpt)
    qp = forward_committor(tpt)

    P_plus = zeros(size(P))
    for i in S_plus
        if i in B_true
            P_plus[i, B_true] = P[i, B_true]/sum(P[i, k] for k in B_true)
        else
            P_plus[i, S_plus] = [P[i, j]*qp[j]/sum(P[i, k]*qp[k] for k in S) for j in S_plus]
        end
    end

    return P_plus
end

function remaining_time(tpt::HomogeneousTPTProblem)
    S, S_plus, B_true, P_plus = 𝒮(tpt), 𝒮_plus(tpt), ℬ_true(tpt), 𝒫_plus(tpt)
    outside_B = setdiff(S_plus, B_true)

    # solve the linear algebra problem t = P_plus t + b restricted to outside_B
    M = I - P_plus[outside_B, outside_B]
    b = fill(1.0, length(outside_B))
    sol = M\b

    # assign the values to a vector; note that t_rem[B] = 0.0 is already handled
    t_rem = zeros(length(S))
    t_rem[outside_B] = sol
    
    # trim very small (possibly negative due to floats) values
    t_rem[abs.(t_rem) .< 1e-16] .= 0.0

    return t_rem       
end

function stationary_statistics(tpt::HomogeneousTPTProblem)
    P = 𝒫(tpt)
    pi_stat, qp, qm = stationary_distribution(tpt), forward_committor(tpt), backward_committor(tpt)
    S = 𝒮(tpt)

    # reactive density
    muAB = [qm[i]*pi_stat[i]*qp[i] for i in S]

    # reactive current
    fij = [qm[i]*pi_stat[i]*P[i, j]*qp[j] for i in S, j in S]

    # forward current
    fplusij = [max(fij[i, j] - fij[j, i], 0) for i in S, j in S]

    # vanE rate
    kAB = sum(fij[i, j] for i in A, j in S)

    # vanE time
    tAB = sum(muAB)/kAB

    # remaining time
    t_rem = remaining_time(tpt)

    return (
        stationary_distribution = pi_stat, 
        forward_committor = qp, 
        backward_committor = qm, 
        reactive_density = muAB, 
        reactive_current = fij, 
        forward_current = fplusij, 
        reactive_rate = kAB, 
        reactive_time = tAB, 
        remaining_time = t_rem)
end

############################################################

A = [1, 2, 3]
B = [3, 4, 5]
P = TransitionMatrix(10)

tpt = HomogeneousTPTProblem(P, A, B)
stats = stationary_statistics(tpt)