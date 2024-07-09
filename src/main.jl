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

    ### PARSE SOURCE/TARGET/AVOID

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

function forward_committor(tpt::HomogeneousTPTProblem)
    P, B, C = 𝒫(tpt), ℬ(tpt), 𝒞(tpt)

    # solve the linear algebra problem q = P q + b on the C subspace
    M = I - P[C, C] 
    b = [sum(P[i, k] for k in B) for i in C]
    sol = M\b
    
    # assign the values to a vector; note that q[A] = 0.0 is already handled
    q = zeros(size(P, 1))
    q[B] .= 1.0
    q[C] = sol
    
    # trim very small (possibly negative due to floats) values
    q[abs.(q) .< 1e-16] .= 0.0
 
    return q
end

############################################################

A = [1, 2, 3]
B = [3, 4, 5]
P = TransitionMatrix(10)

tpt = HomogeneousTPTProblem(P, A, B)
pi_stat = stationary_distribution(tpt)
q_plus = forward_committor(tpt)


