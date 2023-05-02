using Clustering:kmeans, kmedoids
using HDF5
using StatsBase:countmap, mode
using LinearAlgebra


"""
    partition_spectral(tpt_homog)

Partition the space minimally based on Gary Froyland's spectral method applied to the matrix `P`.

This algorithm uses K-means rather than fuzzy C-means.
"""
function partition_spectral(tpt_homog::TPTHomog)
    # the matrix and reversed matrix
    P = tpt_homog.P
    pi_stat = tpt_homog.pi_stat
    P_minus = [(pi_stat[j]/pi_stat[i])*P[j,i] for i in 1:length(pi_stat), j in 1:length(pi_stat)]

    # the "reversibilized" matrrix
    R = (P + P_minus)/2

    # feature extraction of the subdominant eigenvector with k means
    Reigvecs = real(eigvecs(R)[:, end - 1]')
    parts = kmeans(Reigvecs, 2).assignments
    
    return parts
end


"""
    partition_current(tpt_res)

Partition the space minimally based on [`partition_spectral`](@ref) applied to `tpt_res.reactive_current`.
"""
function partition_current(tpt_res::TPTHomogStatResult)
    # the spectral partition using currents
    # we use the subset fij[C, C] since the currents aren't defined at A and B and then stochasticize
    C = tpt_res.sets.C
    fij = tpt_res.reactive_current[C, C]
    fij = fij ./ sum(fij,  dims = 2)

    tpt_homog_C = TPTHomog(fij, [1], [2]) # A and B are arbitrary
    parts_spectral = partition_spectral(tpt_homog_C)

    # add zeros to A/B indices
    parts = zeros(Int64, length(tpt_res.sets.S))
    parts[C] = parts_spectral

    return parts
end


"""
    partition_hitting_location(tpt_res)

Partition the space minimally based on the B-hitting location of the states.

States are in the same partition if they tend to hit B with similar distributions.
"""
function partition_hitting_location(tpt_res::TPTHomogStatResult)
    # recall that rij[i, j] is the probability that state i hits B at state j
    # Of course, rij[i, X] = 0 for X not in B
    rij = tpt_res.hitting_location_distribution

    # the distance matrix between states, where "distance" is the L2 norm between their distributions
    # we exclude union(A, B)
    outside_AB = setdiff(tpt_res.sets.S, union(tpt_res.sets.A, tpt_res.sets.B))
    dist = [norm(rij[i,:] - rij[j, :]) for i in outside_AB, j in outside_AB]

    # construct partitions using kmedoids, since kmeans requires actual coordinates, not just distances
    km = kmedoids(dist, 2).assignments

    # add 0's in A/B locations
    parts = zeros(Int64, length(tpt_res.sets.S))
    parts[outside_AB] = km

    return parts
end


"""
    partition_P_plus(tpt_res)

    Partition the space minimally based on [`partition_spectral`](@ref) applied to `tpt_res.reactive_current`.
"""
function partition_P_plus(tpt_homog::TPTHomog)
    # we apply the spectral method with P_plus restricted to C and normalized accordingly
    C = tpt_homog.sets.C
    P_plus = tpt_homog.P_plus[C, C]
    P_plus = P_plus ./ sum(P_plus, dims = 2)

    tpt_homog_C = TPTHomog(P_plus, [1], [2]) # A and B are arbitrary
    parts_spectral = partition_spectral(tpt_homog_C)

    # add zeros to A/B indices
    parts = zeros(Int64, length(tpt_homog.sets.S))
    parts[C] = parts_spectral

    return parts
end


"""
    minimal_partitions(tpt_homog, tpt_res::TPTHomogStatResult)

Construct all the TPT minimal partitions and return a [`PartitionsStatResult`](@ref).
"""
function minimal_partitions(tpt_homog::TPTHomog, tpt_res::TPTHomogStatResult)
    spectral_P = partition_spectral(tpt_homog)
    spectral_f = partition_current(tpt_res)
    hitting_location = partition_hitting_location(tpt_res)
    
    return PartitionsStatResult(tpt_homog, spectral_P, spectral_f, hitting_location)
end


"""
    minimal_partitions(tpt_homog, tpt_res::TPTHomogNonStatResult)

Construct all the TPT minimal partitions and return a [`PartitionsNonStatResult`](@ref).
"""
function minimal_partitions(tpt_homog::TPTHomog, tpt_res::TPTHomogNonStatResult)
    P_plus = partition_P_plus(tpt_homog)
    
    return PartitionsNonStatResult(tpt_homog, P_plus)
end