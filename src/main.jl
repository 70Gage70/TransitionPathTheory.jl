"""
    abstract type TPTProblem

Supertype for all TPT problems.
"""
abstract type TPTProblem end

"""
    abstract type Stochasticity

Supertype for all transition matrix stochasticity types.
"""
abstract type Stochasticity end

"""
    struct SubStochastic

A `Stochasticity` such that the row sums of the transition matrix are all strictly less than 1.
"""
struct SubStochastic <: Stochasticity end

"""
    struct Stochastic

A `Stochasticity` such that the row sums of the transition matrix are all equal to 1.
"""
struct Stochastic <: Stochasticity end

"""
    struct SuperStochastic

A `Stochasticity` such that the row sums of the transition matrix are all strictly greater than 1.
"""
struct SuperStochastic <: Stochasticity end 

"""
    struct NonStochastic

A `Stochasticity` such that the row sums of the transition matrix are incomparable to 1.
"""
struct NonStochastic <: Stochasticity end 

"""
    abstract type Connectivity

Supertype for all transition matrix connectivity types.
"""
abstract type Connectivity end 

"""
    struct StronglyConnected

A `Connectivity` such that the directed graph implied by the transition matrix is [strongly connected](https://en.wikipedia.org/wiki/Connectivity_(graph_theory)#Connected_vertices_and_graphs).
"""
struct StronglyConnected <: Connectivity end 

"""
    struct WeaklyConnected

A `Connectivity` such that the directed graph implied by the transition matrix is [weakly connected](https://en.wikipedia.org/wiki/Connectivity_(graph_theory)#Connected_vertices_and_graphs).
"""
struct WeaklyConnected <: Connectivity end 

"""
    struct Disconnected

A `Connectivity` such that the directed graph implied by the transition matrix is [disconnected](https://en.wikipedia.org/wiki/Connectivity_(graph_theory)#Connected_vertices_and_graphs).
"""
struct Disconnected <: Connectivity end 

"""
    struct TransitionMatrix{S, C}

A transition probability matrix with [`Stochasticity`](@ref) `S` and [`Connectivity`](@ref) `C`.

### Fields 

- `P`: The transition probability `Matrix` itself.
- `stochasticity`: `S`
- `connectivity`: `C`

### Constructors

    TransitionMatrix(P_size; n_zeros = 0, normalize = true, seed = 1234)

Construct a randomly generated transition matrix.

Arguments 
- `P_size`: The number of rows (== columns) of the matrix.

Optional Arguments 
- `n_zeros`: This many zeros will be placed in the matrix at random. If this results in a row of the matrix having all zeros, a `1` will be placed randomly in that row.
- `normalize`: Whether to normalize the row sums of the resulting matrix.
- `seed`: A seed for reproducible randomness.

#

    TransitionMatrix(P)

Construct a `TransitionMatrix` from a matrix `P`. The properties of `P` are inferred automatically.

Arguments 
- `P`: The transition matrix.
"""
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

function TransitionMatrix(P_size::Integer; n_zeros::Integer = 0, normalize::Bool = true, seed::Integer = 1234)
    @argcheck P_size >= 2
    @argcheck n_zeros >= 0 

    seed!(seed)

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
