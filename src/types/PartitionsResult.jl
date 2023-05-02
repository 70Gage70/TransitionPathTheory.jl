"""
    AbstractPartitionsResult    

An abstract partitions results object for homogeneous TPT.
"""
abstract type AbstractPartitionsResult end


"""
    PartitionsStatResult{U}

A container for the state partitions of stationary TPT.
"""
struct PartitionsStatResult{U<:Integer} <: AbstractPartitionsResult
    spectral_P::Vector{U}
    spectral_f::Vector{U}
    hitting_location::Vector{U}
end


"""
    PartitionsStatResult(tpt_homog, partitions...)

Construct a [`PartitionsStatResult`](@ref) standardized according to [`standardize_partitions`](@ref).
"""
function PartitionsStatResult(
    tpt_homog::TPTHomog,
    partitions::Vector{U}...) where {U <: Integer}

    # ensure all partitions are of the same length
    if length(partitions) > 1 
        for i = 1:length(partitions)-1
            l1, l2 = length(partitions[i]), length(partitions[i+1])
            @assert l1 == l2 "Partitions must be the same length (got $(l1) and $(l2))."
        end
    end

    std_parts = standardize_partitions(tpt_homog, partitions...)
    return PartitionsStatResult(std_parts...)
end


"""
PartitionsNonStatResult{U}

A container for the state partitions of nonstationary TPT.
"""
struct PartitionsNonStatResult{U<:Integer} <: AbstractPartitionsResult
    spectral_P_plus::Vector{U}
end


"""
    PartitionsNonStatResult(tpt_homog, partitions...)

Construct a [`PartitionsNonStatResult`](@ref) standardized according to [`standardize_partitions`](@ref).
"""
function PartitionsNonStatResult(
    tpt_homog::TPTHomog,
    partitions::Vector{U}...) where {U <: Integer}

    # ensure all partitions are of the same length
    if length(partitions) > 1 
        for i = 1:length(partitions)-1
            l1, l2 = length(partitions[i]), length(partitions[i+1])
            @assert l1 == l2 "Partitions must be the same length (got $(l1) and $(l2))."
        end
    end

    std_parts = standardize_partitions(tpt_homog, partitions...)
    return PartitionsNonStatResult(std_parts...)
end

"""
    standardize_partitions(tpt_homog, partitions...)

Return `partitions...` standardized to the convention that the states connected to `tpt_homog.sets.A_true.` are partition `1`.
"""
function standardize_partitions(
    tpt_homog::TPTHomog,
    partitions::Vector{U}...) where {U <: Integer}

    P = tpt_homog.P
    A_true = tpt_homog.sets.A_true
    # P[A_true, :] is the transition probabilities out of A_true, therefore
    # sum(P[A_true, :], dims = 1) is a row vector such that its i'th entry is only positive if A_true is connected via P to state i
    connected_to_A_true = findall(x -> x != 0.0, vec(sum(P[A_true, :], dims = 1))) 
    connected_to_A_true = union(connected_to_A_true, A_true)

    standardized_partitions = Vector{U}[]
    for part in partitions
        part_A = mode(filter(x -> x != 0, part[connected_to_A_true])) # 1 or 2 for minimal partition

        # if part_A == 2, the convention of `part` is backwards, so it has to be switched
        if part_A == 2
            part = [i == 0 ? 0 : 3 - i for i in part] # 3 - 1 = 2 and 3 - 2 = 1
        end

        push!(standardized_partitions, part)
    end

    return standardized_partitions
end


function Base.show(io::IO, x::PartitionsStatResult)
    print(io, "PartitionsStatResult[")
    show(io, length(x.spectral_P))
    print(" states]")
end


function Base.show(io::IO, x::PartitionsNonStatResult)
    print(io, "PartitionsNonStatResult[")
    show(io, length(x.spectral_P_plus))
    print(" states]")
end