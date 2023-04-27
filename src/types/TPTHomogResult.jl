"""
    AbstractTPTHomogResult

An abstract statistics result object for homogeneous TPT.
"""
abstract type AbstractTPTHomogResult end

"""
    TPTHomogStatResult{T, U}

A container for all the results of a homogeneous, stationary TPT calculations.  
"""
struct TPTHomogStatResult{T<:Real, U<:Integer} <: AbstractTPTHomogResult
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
    hitting_location_distribution::Matrix{T}
    most_probable_hitting_location::Vector{U}
end

function Base.show(io::IO, x::TPTHomogStatResult)
    print(io, "TPTHomogStatResult[")
    show(io, length(x.q_plus))
    print(" states]")
end

"""
    TPTHomogNonStatResult{T, U}

A container for all the results of a homogeneous, nonstationary TPT calculations.  
"""
struct TPTHomogNonStatResult{T<:Real, U<:Integer} <: AbstractTPTHomogResult
    sets::TPTSets{U}
    horizon::U
    pi_stationary::Vector{T}
    q_plus::Vector{T}
    q_minus::Vector{T}
    # TPT statistics are time-dependent now, so each has one extra dimension relative to the stationary case.
    density::Matrix{T}
    normalized_reactive_density::Matrix{T}
    hitting_time_cdf::Matrix{T}
    hitting_time_cdf_AB::Vector{T}
    reactive_current::Array{T, 3}
    forward_current::Array{T, 3}
end

function Base.show(io::IO, x::TPTHomogNonStatResult)
    print(io, "TPTHomogNonStatResult[")
    show(io, length(x.q_plus))
    print(" states, ")
    show(io, x.horizon)
    print(" steps]")
end