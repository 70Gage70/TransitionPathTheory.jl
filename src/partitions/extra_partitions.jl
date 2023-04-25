"""
Create the minimal partition in a data set from Ulam and TPT.
"""

using LinearAlgebra
using Clustering
using HDF5
using StatsBase:countmap, mode

####################################################################################
####################################################################################
####################################################################################


"""
Left eigenvector of P, normalized appropriately. The stationary distribution of P if it's row stochastic.
"""

function Ppi(P)
    return transpose(abs.(normalize(eigvecs(transpose(P))[:,size(P)[1]], 1)))
end
 

"""
Transition matrix of the reversed Markov chain
"""


function Pminus(P, Ppi)
    Pmin = zeros(length(Ppi), length(Ppi))
    for i = 1:length(Ppi)
        for j = 1:length(Ppi)
            if Ppi[i] > 0.0 # else: already 0.0
                Pmin[i, j] = (Ppi[j]/Ppi[i])*P[j,i]
            end
        end
    end
    
    return Pmin
end

####################################################################################

"""
Based on Gary Froyland's spectral method (not identical)

PrevP: matrix = (P + P-)/2, pi = pi
Pij: matrix = P, pi = pi
"""

function min_part_gary_P(ulam; type = "PrevP")
    P_closed = ulam["P_closed"]
    pi_closed = ulam["pi_closed"]
    P_closed_rev = Pminus(P_closed, pi_closed)


    if type == "PrevP"
        R = (P_closed + P_closed_rev)/2
    elseif type == "Pij"
        R = P_closed
    end

    # L = 2 # number of partitions
    # Reigvecs = real(eigvecs(R)[:, end - L:end - 1]')
    
    # # feature extraction with k means
    # km = kmeans(Reigvecs, L)

    Reigvecs = real(eigvecs(R)[:, end - 1:end - 1]')
    
    # feature extraction with k means
    km = kmeans(Reigvecs, 2)

    parts = km.assignments[1:end - 1] # get rid of nirvana
    
    return parts
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
Partition based on cutting the current (really J(x))
"""

function min_part_J(ulam, tpt)
    f_plus = tpt["f+"]
    centers = ulam["polys_centers"]
    fluxij = zeros(size(f_plus, 1), 2)
    
    for i = 1:size(fluxij, 1)
        coords_this = centers[i,:]
        flux = [0.0, 0.0]
        for j = 1:size(fluxij, 1)
            if j != i
                coords_next = centers[j,:]
                flux = flux + f_plus[i, j]*normalize(coords_next - coords_this)
            end
        end
    
        fluxij[i, :] = flux
    end
    
    mags = [norm(fluxij[i, :]) for i = 1:size(fluxij,1)]

    muAB = tpt["muABnorm"]
    good_inds = findall(x -> x != 0.0, muAB)
    km = kmeans(mags[good_inds]', 2) # don't want to include A, B in the clustering


    parts = zeros(Int64, length(muAB))
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
Partition based on the remaining time
"""

function min_part_trem(tpt)
    trem = tpt["t_rem"]

    muAB = tpt["muABnorm"]
    good_inds = findall(x -> x != 0.0, muAB)
    km = kmeans(trem[good_inds]', 2).assignments 

    parts = zeros(Int64, length(muAB))
    parts[good_inds] = km

    return parts
end

"""
Switches 1 to 2 and 2 to 1 in a list
"""


function swap_1_2(list)
    return map(x -> 3 - x, list)
end


"""
Standardize partition colors so that the partition that contains A is "1".

Based on which part A belongs to in Gary Part.

If A straddles multiple partitions, warn user and take largest partition to be 1.
"""


function standardize_parts(gary_parts, tpt)
    Ainds = tpt["Ainds"]
    partA = gary_parts[Ainds]
    cmp = countmap(partsA)

    if length(keys(cmp)) > 1
        display("Warning: A straddles two partitions.")
        partA = mode(partA)
    else
        partA = partA[1]
    end

    


    return 1

end


"""
All minimal partitions, returns a dictionary and outputs a .h5 file.
"""

function min_parts(ulam, tpt; extra_suffix = "")
    gary_PrevP = min_part_gary_P(ulam, type = "PrevP")
    # gary_Pij = min_part_gary_P(ulam, type = "Pij")
    gary_frevf = min_part_gary_f(tpt, type = "frevf")
    # gary_fij = min_part_gary_f(tpt, type = "fij")
    # J_cut = min_part_J(ulam, tpt)
    B_hit = min_part_B(tpt)
    # trem_cut = min_part_trem(tpt)

    # Now we standardize the partitions such that the partition that contains A is Partition 1
    # The "partition that contains A" is based on gary_PrevP


    return_dict = Dict(
        "gary_PrevP" => gary_PrevP,
        # "gary_Pij" => gary_Pij,
        "gary_frevf" => gary_frevf,
        # "gary_fij" => gary_fij,
        # "J_cut" => J_cut,
        "B_hit" => B_hit,
        # "trem_cut" => trem_cut
        )

    fname = "partitions" * extra_suffix * ".h5"
    rm(fname, force = true)
    fout = h5open(fname, "w")
    write_dict_to_h5(fout, "partitions", return_dict)
    close(fout)
    
    return return_dict
end

