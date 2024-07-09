"""
    struct HomogeneousTPTProblem

A [`TPTProblem`](@ref) where the transition matrix is homogeneous, i.e. time-independent.

### Fields

- `transition_matrix`: A [`TransitionMatrix`](@ref).
- `source`: A vector of indices defining the TPT source, aka `A`.
- `target`: A vector of indices defining the TPT target, aka `B`.

### Constructor

    HomogeneousTPTProblem(P, source, target; avoid = Int64[])

Build a `HomogeneousTPTProblem` where each argument of the constructor maps to the corresponding field.

Optional Arguments 

- `avoid`: Indices provided here will be appended to both `A` and `B` and will hence be avoided by transition path theory.
"""
struct HomogeneousTPTProblem{TM<:TransitionMatrix} <: TPTProblem
    transition_matrix::TM
    source::Vector{Int64}
    target::Vector{Int64}
end

function HomogeneousTPTProblem(
    transition_matrix::TransitionMatrix, 
    source::Vector{<:Integer}, 
    target::Vector{<:Integer}; 
    avoid::Vector{<:Integer} = Int64[])

    # ensure no repeated indices
    @argcheck allunique(source)
    @argcheck allunique(target)
    @argcheck allunique(avoid)

    # ensure indices take appropriate values
    n_states = size(transition_matrix.P, 1)
    @argcheck all(1 .<= source .<= n_states)
    @argcheck all(1 .<= target .<= n_states)
    @argcheck all(1 .<= avoid .<= n_states)

    source = [source ; avoid] |> sort
    target = [target ; avoid] |> sort

    return HomogeneousTPTProblem(transition_matrix, source, target)
end

"""
    𝒫(tpt)

Return `tpt.transition_matrix.P`.
"""
𝒫(tpt::HomogeneousTPTProblem) = tpt.transition_matrix.P

"""
    𝒜(tpt)

Return `tpt.source`.
"""
𝒜(tpt::HomogeneousTPTProblem) = tpt.source    

"""
    ℬ(tpt)

Return `tpt.target`.
"""
ℬ(tpt::HomogeneousTPTProblem) = tpt.target

"""
    𝒮(tpt)

Return all possible indices implied by the transition matrix, i.e. `1:size(P, 1)`.
"""
𝒮(tpt::HomogeneousTPTProblem) = collect(1:size(𝒫(tpt),1))

"""
    Ω(tpt)

Return the set of avoided states.
"""
Ω(tpt::HomogeneousTPTProblem) = intersect(𝒜(tpt), ℬ(tpt))

"""
    𝒜_true(tpt)

Return the set of states that are members of the `source`, but are not avoided.
"""
𝒜_true(tpt::HomogeneousTPTProblem) = setdiff(𝒜(tpt), Ω(tpt))      

"""
    ℬ_true(tpt)

Return the set of states that are members of the `target`, but are not avoided.
"""
ℬ_true(tpt::HomogeneousTPTProblem) = setdiff(ℬ(tpt), Ω(tpt))     

"""
    𝒞(tpt)

Return the set of states that are not in the source or target.
"""
𝒞(tpt::HomogeneousTPTProblem) = setdiff(𝒮(tpt), union(𝒜(tpt), ℬ(tpt)))

"""
    stationary_distribution(tpt)

Compute the stationary distribution of `tpt.transition_matrix.P`.

Errors when `P` is not strongly connected.
"""
function stationary_distribution(tpt::HomogeneousTPTProblem)
    P = 𝒫(tpt)

    if tpt.transition_matrix.connectivity === StronglyConnected()
        return normalize(abs.(eigvecs(P')[:,end]), 1)
    else
        error("Transition matrix is not strongly connected, so the stationary distribution is not unique.")
    end
end

"""
    𝒫_backwards(tpt)

Compute the "backwards" transition matrix; return a `Matrix`.
"""
function 𝒫_backwards(tpt::HomogeneousTPTProblem)
    P, S = 𝒫(tpt), 𝒮(tpt)
    pi_stat = stationary_distribution(tpt)

    return [(pi_stat[j]/pi_stat[i])*P[j,i] for i in S, j in S]
end

"""
    forward_committor(tpt)

Compute the forward committor.
"""
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

"""
    backward_committor(tpt)

Compute the backward committor.
"""
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

"""
    remaining_time(tpt)

Compute the remaining time, aka the lead time.
"""
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

"""
    hitting_location_distribution(tpt)

Compute the distribution of target-hitting locations. 

Return a matrix `R` such that `R[i, j]` is the probability that - starting in state `i` - the first visit to `B` occurs in state `j`. Note that 
`R[i, j] == 0` for all `j ∉ B`.
"""
function hitting_location_distribution(tpt::HomogeneousTPTProblem)
    S, S_plus, B_true, P_plus = 𝒮(tpt), 𝒮_plus(tpt), ℬ_true(tpt), 𝒫_plus(tpt)

    # solve the linear algebra problem r = P_plus r + b restricted to outside_B
    outside_B = setdiff(S_plus, B_true)
    M = I - P_plus[outside_B, outside_B]

    # the probability that state i hits B at index j
    rij = zeros(length(S), length(S))
    for Bind in B_true
        b = P_plus[outside_B, Bind]
        rij[outside_B, Bind] = M\b
    end

    # starting at any given B_true index, guaranteed to hit that 
    rij[B_true, B_true] = [i == j ? 1.0 : 0.0 for i = 1:length(B_true), j = 1:length(B_true)]
    
    # trim very small (possibly negative due to floats) values
    rij[abs.(rij) .< 1e-16] .= 0.0

    return rij
end

"""
    stationary_statistics(tpt)

Compute and return the following statistics in a `NamedTuple`:

- `stationary_distribution`
- `forward_committor`
- `backward_committor`
- `reactive_density`
- `reactive_current`
- `forward_current`
- `reactive_rate`
- `reactive_time`
- `remaining_time`
- `hitting_location_distribution`
"""
function stationary_statistics(tpt::HomogeneousTPTProblem)
    P = 𝒫(tpt)
    pi_stat, qp, qm = stationary_distribution(tpt), forward_committor(tpt), backward_committor(tpt)
    S, A = 𝒮(tpt), 𝒜(tpt)

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

    # hitting locations
    rij = hitting_location_distribution(tpt)

    return (
        stationary_distribution = pi_stat, 
        forward_committor = qp, 
        backward_committor = qm, 
        reactive_density = muAB, 
        reactive_current = fij, 
        forward_current = fplusij, 
        reactive_rate = kAB, 
        reactive_time = tAB, 
        remaining_time = t_rem,
        hitting_location_distribution = rij)
end


"""
    nonstationary_statistics(tpt, horizon)

### Arguments

- `tpt`: The [`TPTProblem`](@ref).
- `horizon`: The time step at which to cut off the calculation. Note that the `horizon` value does NOT enforce that trajectories leaving A hit B by that time.

Compute and return the following statistics in a `NamedTuple`:

- `density`
- `reactive_density`
- `t_cdf`
- `t_cdf_AB`
"""
function nonstationary_statistics(tpt::HomogeneousTPTProblem, horizon::Integer)
    @argcheck horizon >= 1

    P, P_plus = 𝒫(tpt), 𝒫_plus(tpt)
    A_true, B_true, S, S_plus = 𝒜_true(tpt), ℬ_true(tpt), 𝒮(tpt), 𝒮_plus(tpt)

    i0 = [i in A_true ? 1.0/length(A_true) : 0.0 for i in S] # uniform distribution supported on A_true

    density = Vector{Float64}[]
    push!(density, i0)

    reactive_density = Vector{Float64}[]
    push!(reactive_density, i0)

    outside_B = setdiff(S_plus, B_true)
    t_cdf = Vector{Float64}[]
    push!(t_cdf, [i in outside_B ? sum(P_plus[i, j] for j in B_true) : i in B_true ? 1.0 : 0.0 for i in S])

    for _ = 2:horizon
        # density
        push!(density, transpose(P) * density[end])

        # normalized reactive density
        push!(reactive_density, transpose(P_plus) * reactive_density[end])

        # hitting time distribution
        t_cdf_next = [i in B_true ? 1.0 : 0.0 for i in S]
        t_cdf_next[outside_B] .= t_cdf[1][outside_B] + P_plus[outside_B, outside_B] * t_cdf[end][outside_B]
        push!(t_cdf, t_cdf_next)
    end

    # create time cdf for A to B
    t_cdf_AB = [sum(i0[i]*t_cdf[n][i] for i in A_true) for n = 1:horizon]

    return (
        density = density, 
        reactive_density = reactive_density, 
        t_cdf = t_cdf, 
        t_cdf_AB = t_cdf_AB)
end
