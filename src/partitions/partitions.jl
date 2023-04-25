using Clustering:kmeans
using HDF5
using StatsBase:countmap, mode
using LinearAlgebra

struct PartitionsResult{U<:Integer}
    spectral_P::Vector{U}
    hitting_location::Vector{U}
end


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
    partition_hitting_location(tpt_res)

Partition the space minimally based on the B-hitting location of the states.

States are in the same partition if they tend to hit B with similar distributions.
"""
function partition_hitting_location(tpt_res::TPTHomogStatResult)
    # recall that rij[i, j] is the probability that state i hits B at state j
    # Of course, rij[i, X] = 0 for X not in B
    rij = tpt_res.statistics.hitting_location_distribution

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
    standardize_minimal_partition(tpt_homog, partition)

Standardize `partition` such that the partition that contains (most of) `tpt_homog.sets.A_true` is Partition 1.
"""
function standardize_minimal_partition(tpt_homog::TPTHomog, partition::Vector{<:Integer})
    # the partition that contains the most of A
    part_A = mode(partition[tpt_homog.sets.A_true])

    if part_A == 2
        partition = [part == 0 ? 0 : 3 - part for part in partition] # switches 1 to 2 and 2 to 1 while leaving 0 as 0
    end

    return partition
end


"""
    minimal_partitions(s)
"""
function minimal_partitions(tpt_homog::TPTHomog)
    # the standard spectral partitions
    spectral_P = partition_spectral(tpt_homog)


    gary_PrevP = min_part_gary_P(ulam, type = "PrevP")
    gary_frevf = min_part_gary_f(tpt, type = "frevf")
    B_hit = min_part_B(tpt)

    # Now we standardize the partitions such that the partition that contains A is Partition 1
    # This convention is also applied to Gary partition even though it's independent of A

    P_Acon =  ulam["P_open"][tpt["Ainds"], :] # chunk of P that gives P(A -> elsewhere)
    gary_PrevP = standardize_parts(P_Acon, gary_PrevP)
    gary_frevf = standardize_parts(P_Acon, gary_frevf)
    B_hit = standardize_parts(P_Acon, B_hit)

    return_dict = Dict(
        "gary_PrevP" => gary_PrevP,
        "gary_frevf" => gary_frevf,
        "B_hit" => B_hit
        )

    if h5out
        fname = "partitions" * extra_suffix * ".h5"
        rm(fname, force = true)
        fout = h5open(fname, "w")
        write_dict_to_h5(fout, "partitions", return_dict)
        close(fout)
    end
    
    return return_dict
end

