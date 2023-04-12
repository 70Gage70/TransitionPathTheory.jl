"""
    TPTStats{T, U}

A container for all the results of TPT calculations.  
"""
struct TPTStats{T<:Real, U<:Integer}
    q_plus::Vector{T}
    q_minus::Vector{T}
    reactive_density::Vector{T}
    normalized_reactive_density::Vector{T}
    reactive_current::Matrix{T}
    forward_current::Matrix{T}
    reactive_rate_vanE::T
    transition_time_vanE::T
    remaining_time::Vector{T}
    hitting_distribution::Matrix{T}
    hitting_locations::Vector{U}
    time_cdf::Matrix{T}
    time_cdf_AB::Vector{T}
end

function Base.show(io::IO, x::TPTStats)
    print(io, "TPTStats[")
    show(io, length(x.q_plus))
    print("]")
end