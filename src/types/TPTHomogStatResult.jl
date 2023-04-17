"""
    TPTHomogStatResult{T, U}

A container for all the results of a homogeneous, stationary TPT calculations.  
"""
struct TPTHomogStatResult{T<:Real, U<:Integer}
    sets::TPTSets{U}
    pi_stationary::Vector{T}
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

function Base.show(io::IO, x::TPTHomogStatResult)
    print(io, "TPTHomogStatResult[")
    show(io, length(x.q_plus))
    print("]")
end