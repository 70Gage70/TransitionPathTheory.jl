using Clustering:kmeans
using HDF5
using StatsBase:countmap, mode
using LinearAlgebra

struct PartitionsResult{U<:Integer}
    spectral_P::Vector{U}
    hitting_location::Vector{U}
end



"""
    partition_spectral(tpt_homog::TPTHomog)

Partition the space minimally based on Gary Froyland's spectral method applied to the matrix `P`.
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

function partition_statistic(stat::Vector{<:Real})
    return kmeans
end


"""
Based on Gary Froyland's spectral method (not identical) with TPT probabilities.

fij: matrix = fij normalized, pi = muAB normalized
frevf: matrix = (fij normalized + fij- normalized)/2, pi = muAB normalized
"""

function min_part_gary_f(tpt; type = "frevf")
    fij = tpt["fij"] # note that these already exclude nirvana
    muAB = tpt["muABnorm"]

    all_inds = 1:length(muAB)
    good_inds = findall(x -> x != 0.0, muAB) # have to avoid A and B since there's no density there

    fij = fij[good_inds, good_inds]
    for i = 1:size(fij, 1)
        fij[i,:] = fij[i,:]/sum(fij[i,:])
    end

    muAB = muAB[good_inds]

    fij_rev = Pminus(fij, muAB)


    if type == "frevf"
        R = (fij + fij_rev)/2
    elseif type == "fij"
        R = fij
    end

    L = 2 # number of partitions
    Reigvecs = real(eigvecs(R)[:, end - L:end - 1]')
    
    # feature extraction with k means
    km = kmeans(Reigvecs, L)
    parts = zeros(Int64, length(all_inds))
    parts[good_inds] = km.assignments
 
    return parts
end


"""
Partition based on the distribution of B-hitting locations. That is, states are in the same partition
if they hit B with roughly the same distribution.
"""

function min_part_B(tpt)
    rij = tpt["rij"]
    dist = zeros(size(rij)) # distance matrix between states
    for i = 1:size(rij, 1)
        for j = 1:size(rij, 1)
            dist[i, j] = norm(rij[i,:] - rij[j, :])
        end
    end

    muAB = tpt["muABnorm"]
    good_inds = findall(x -> x != 0.0, muAB)
    dist = dist[good_inds, good_inds] # don't want to include A, B in clustering
    km = kmedoids(dist, 2) # have to use medioids since kmeans doesn't work without actual coordinates


    parts = zeros(Int64, length(muAB))
    parts[good_inds] = km.assignments

    return parts
end


"""
Standardize the partitions such that the partition that contains (most of) A is Partition 1
"""

function standardize_parts(P_Acon, parts)
    connected_to_A = Int64[]

    for i = 1:size(P_Acon, 1)
        connected_to_A = union(connected_to_A, findall(x -> x != 0.0, P_Acon[i,:])) # indices of states A hits
    end

    indcon = mode(parts[connected_to_A]) # most common part of states that A hits

    parts_new = parts

    if indcon == 2
        for i = 1:length(parts)
            if parts[i] == 1
                parts_new[i] = 2
            elseif parts[i] == 2
                parts_new[i] = 1
            end
        end
    end

    return parts_new
end


"""
    minimal_partitions(s)
"""
function minimal_partitions(tpt_homog::TPTHomog)
    # the standard spectral partitions
    spectral_P = partition_spectral(tpt_homog.P, tpt_homog.pi_stat)

    # the spectral partitions based on forward currents
    abs.(LinearAlgebra.normalize(LinearAlgebra.eigvecs(transpose(P))[:,size(P)[1]], 1))

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

