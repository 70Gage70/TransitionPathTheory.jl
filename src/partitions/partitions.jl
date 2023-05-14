using Clustering:kmeans, kmedoids
using Graphs:strongly_connected_components, SimpleDiGraph
using HDF5
using StatsBase:countmap, mode
using LinearAlgebra
using Random

"""
    partitionable_submatrix(P, tpt_sets)

Return a tuple `(P_part, inds)` where 

- `P_part` is a submatrix of `P` such that `P_part` is the transition \
probability matrix of an ergodic Markov chain defined on the reactive subset of `tpt_sets`
- `inds` is a list of integers such that `P[inds, inds] = P_part`
"""
function partitionable_submatrix(P::Matrix{<:Real}, tpt_sets::TPTSets)
    # C_part avoids non-reactive states
    AB = union(tpt_sets.A_true, tpt_sets.B_true)
    C_part = setdiff(tpt_sets.S_plus, AB)
    P_part = P[C_part, C_part]

    # find the largest strongly connected component 
    Padj::Matrix{Int64} = [P_part[i,j] != 0.0 ? 1 : 0 for i = 1:size(P_part,1), j = 1:size(P_part,1)]
    scc = sort(sort(strongly_connected_components(SimpleDiGraph(Padj)), by = length)[end])
    P_part = P_part[scc, scc]

    # stochasticize
    P_part = P_part ./ sum(P_part,  dims = 2)
    inds = C_part[scc]

    return (P_part, inds)
end


"""
    partition_spectral(tpt_homog)

Partition the space minimally based on Gary Froyland's spectral method applied to the matrix `P`.

This algorithm uses K-means rather than fuzzy C-means.
"""
function partition_spectral(tpt_homog::TPTHomog; rseed = 123)
    # the matrix and reversed matrix
    P = tpt_homog.P
    pi_stat = tpt_homog.pi_stat
    P_minus = [(pi_stat[j]/pi_stat[i])*P[j,i] for i in 1:length(pi_stat), j in 1:length(pi_stat)]

    # the "reversibilized" matrrix
    R = (P + P_minus)/2

    # feature extraction of the subdominant eigenvector with k means
    Reigvecs = real(eigvecs(R)[:, end - 1]')
    Random.seed!(rseed)
    parts = kmeans(Reigvecs, 2).assignments
    
    return parts
end


"""
    partition_current(tpt_res)

Partition the space minimally based on [`partition_spectral`](@ref) applied to `tpt_res.reactive_current`.
"""
function partition_current(tpt_res::TPTHomogStatResult; rseed = 123)
    fij, inds = partitionable_submatrix(tpt_res.reactive_current, tpt_res.sets)
    tpth_parts = TPTHomog(fij, [1], [2]) # A and B are arbitrary
    parts_spectral = partition_spectral(tpth_parts, rseed = rseed)

    # add zeros to A/B indices
    parts = zeros(Int64, length(tpt_res.sets.S))
    parts[inds] = parts_spectral

    return parts
end


"""
    partition_hitting_location(tpt_res)

Partition the space minimally based on the B-hitting location of the states.

States are in the same partition if they tend to hit B with similar distributions.
"""
function partition_hitting_location(tpt_res::TPTHomogStatResult; rseed = 123)
    # recall that rij[i, j] is the probability that state i hits B at state j
    # Of course, rij[i, X] = 0 for X not in B
    rij = tpt_res.hitting_location_distribution

    # the distance matrix between states, where "distance" is the L2 norm between their distributions
    # C_part avoids non-reactive states
    AB = union(tpt_res.sets.A_true, tpt_res.sets.B_true)
    C_part = setdiff(tpt_res.sets.S_plus, AB)
    dist = [norm(rij[i,:] - rij[j, :]) for i in C_part, j in C_part]

    # construct partitions using kmedoids, since kmeans requires actual coordinates, not just distances
    Random.seed!(rseed)
    km = kmedoids(dist, 2).assignments

    # add 0's in A/B locations
    parts = zeros(Int64, length(tpt_res.sets.S))
    parts[C_part] = km

    return parts
end


"""
    partition_P_plus(tpt_homog::TPTHomog; rseed = 123)

    Partition the space minimally based on [`partition_spectral`](@ref) applied to `tpt_res.reactive_current`.
"""
function partition_P_plus(tpt_homog::TPTHomog; rseed = 123)
    P_plus, inds = partitionable_submatrix(tpt_homog.P_plus, tpt_homog.sets)
    tpth_parts = TPTHomog(P_plus, [1], [2]) # A and B are arbitrary
    parts_spectral = partition_spectral(tpth_parts, rseed = rseed)

    # add zeros to A/B indices
    parts = zeros(Int64, length(tpt_homog.sets.S))
    parts[inds] = parts_spectral

    return parts
end


"""
    minimal_partitions(tpt_homog, tpt_res::TPTHomogStatResult)

Construct all the TPT minimal partitions and return a [`PartitionsStatResult`](@ref).
"""
function minimal_partitions(tpt_homog::TPTHomog, tpt_res::TPTHomogStatResult; rseed = 123)
    spectral_P = partition_spectral(tpt_homog, rseed = rseed)
    spectral_f = partition_current(tpt_res, rseed = rseed)
    hitting_location = partition_hitting_location(tpt_res, rseed = rseed)
    
    return PartitionsStatResult(tpt_homog, spectral_P, spectral_f, hitting_location)
end


"""
    minimal_partitions(tpt_homog, tpt_res::TPTHomogNonStatResult)

Construct all the TPT minimal partitions and return a [`PartitionsNonStatResult`](@ref).
"""
function minimal_partitions(tpt_homog::TPTHomog, tpt_res::TPTHomogNonStatResult; rseed = 123)
    P_plus = partition_P_plus(tpt_homog, rseed = rseed)
    
    return PartitionsNonStatResult(tpt_homog, P_plus)
end