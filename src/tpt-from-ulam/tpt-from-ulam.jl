using UlamMethod

"""
    remove_omega(result::Union{AbstractTPTHomogResult, AbstractPartitionsResult})

Remove the last state from indices, statistics and/or partitions of `result`. Used for applications to Ulam's method.

Note that this does not change the  actual results of any calculations, it is only for convenient reference
to TPT statistics without the omega state.
"""
function remove_omega(result::Union{AbstractTPTHomogResult, AbstractPartitionsResult})
    return remove_omega(result)
end


function remove_omega(tpt_result::TPTHomogStatResult)
    omega_tpt_result = []

    for fn in fieldnames(typeof(tpt_result))
        if fn == :sets # subtract the set [omega] from each set
            omega = [tpt_result.sets.S[end]]
            omega_sets = TPTSets([setdiff(getfield(tpt_result.sets, fn), omega) for fn in fieldnames(TPTSets)]...)
            push!(omega_tpt_result, omega_sets)
        else
            gf = getfield(tpt_result, fn)

            if fn in [:reactive_rate_vanE, :transition_time_vanE] # scalar, don't have to remove omega 
                push!(omega_tpt_result, gf)
            elseif fn in [:pi_stationary, :q_plus, :q_minus, :reactive_density, :normalized_reactive_density, :remaining_time, :most_probable_hitting_location] # vectors with omega at the end
                push!(omega_tpt_result, gf[1:end - 1])
            elseif fn in [:reactive_current, :forward_current, :hitting_location_distribution] # matrices with omega in the last row/column
                push!(omega_tpt_result, gf[1:end - 1, 1:end - 1])
            else
                @error "TPTHomogStatResult has an extra field." # shouldn't happen
            end
        end
    end
    
    omega_tpt_result = TPTHomogStatResult(omega_tpt_result...)

    return omega_tpt_result
end


function remove_omega(tpt_result::TPTHomogNonStatResult)
    omega_tpt_result = []

    for fn in fieldnames(typeof(tpt_result))
        if fn == :sets # subtract the set [omega] from each set
            omega = [tpt_result.sets.S[end]]
            omega_sets = TPTSets([setdiff(getfield(tpt_result.sets, fn), omega) for fn in fieldnames(TPTSets)]...)
            push!(omega_tpt_result, omega_sets)
        else
            gf = getfield(tpt_result, fn)

            if fn in [:horizon] # scalar, don't have to remove omega 
                push!(omega_tpt_result, gf)
            elseif fn in [:pi_stationary, :q_plus, :q_minus] # vectors with omega at the end
                push!(omega_tpt_result, gf[1:end - 1])
            elseif fn in [:hitting_time_cdf_AB] # vectors which don't contain omega
            push!(omega_tpt_result, gf)                
            elseif fn in [:density, :normalized_reactive_density, :hitting_time_cdf] # time-dependent vectors with omega in the last column
                push!(omega_tpt_result, gf[:, 1:end - 1])
            elseif fn in [:reactive_current, :forward_current] # time-dependent matrices with omega in the last row/column
                push!(omega_tpt_result, gf[:, 1:end - 1, 1:end - 1])
            else
                @error "TPTHomogNonStatResult has an extra field." # shouldn't happen
            end
        end
    end
    
    omega_tpt_result = TPTHomogNonStatResult(omega_tpt_result...)

    return omega_tpt_result
end


function remove_omega(parts_result::AbstractPartitionsResult)
    omega_parts_result = []

    for fn in fieldnames(typeof(parts_result))
        gf = getfield(parts_result, fn)
        push!(omega_parts_result, gf[1:end-1])
    end
    
    omega_parts_result = typeof(parts_result)(omega_parts_result...)

    return omega_parts_result
end


"""
    ulam_polys_to_indices(ulam::UlamResult, region_points::Matrix)

Return a list of polygon indices i such that the i'th polygon contains at least one point in `region_points`.
"""
function ulam_polys_to_indices(ulam::UlamResult, region_points::Matrix{<:Real})
    polys = ulam.polys
    inds = unique(inpoly(region_points, PolyTable(polys)).inds)
    filter!(x->x!=0, inds)

    if length(inds) == 0 
        @warn "None of the provided points is in any Ulam polygon." 
    end

    return inds
end

"""
    ulam_polys_to_indices(ulam::UlamResult, region_polygon::Matrix)

Return a list of polygon indices i such that the i'th polygon has a non-empty intersection with the polygon `region_polygon`.
"""
function ulam_polys_to_indices(ulam::UlamResult, region_polygon::UlamPolygon)
    polys = ulam.polys
    inds = [i for i in 1:length(polys) if ulam_intersects(region_polygon, polys[i])]
    filter!(x->x!=0, inds)

    if length(inds) == 0 
        @warn "The provided polygon does not intersect with any Ulam polygon." 
    end

    return inds
end

# """
#     tpt_from_ulam(s)

# Assumes UlamMethod.jl is loaded.

# Default is A is centers but B defines polygon.
# """
# function tpt_from_ulam_nirvana(
#     fin::String, 
#     domain::Matrix{<:Real}, 
#     poly_type::String, 
#     poly_number::Integer,
#     A_centers::Matrix{<:Real},
#     B_vertices::Matrix{<:Real};
#     avoid_vertices::Union{Matrix{<:Real}, Nothing} = nothing,
#     tpt_out::Union{String, Nothing} = nothing,
#     ulam_out::Union{String, Nothing} = nothing)

#     traj = UlamTrajectories(fin)
#     domain = UlamDomain(UlamPolygon(domain), poly_type = poly_type, poly_number = poly_number)
#     ulam = ulam_method(traj, domain)
    
#     A = inpoly(A_centers, PolyTable(ulam.polys)).inds
#     B = [i for i in 1:length(ulam.polys) if ulam_intersects(UlamPolygon(B_vertices), ulam.polys[i])]
    
#     if avoid_vertices !== nothing
#         avoid = [i for i in 1:length(ulam.polys) if ulam_intersects(UlamPolygon(avoid_vertices), ulam.polys[i])]
#         A = sort(union(A, avoid))
#         B = sort(union(B, avoid))
#     end

#     P = ulam.P_closed
#     omega = [size(P, 1)]
#     A = union(A, omega)
#     B = union(B, omega)
    
#     tpt_homog = TPTHomog(P, A, B)
#     tpt_res = tpt_infinite(tpt_homog)
    
#     tpt_homog, tpt_res = stats_remove_omega(tpt_homog, tpt_res)

#     if ulam_out !== nothing
#         ulam_write(ulam_out, ulam)
#     end
    
#     if tpt_out !== nothing
#         tpt_write(tpt_out, tpt_homog, tpt_res)
#     end

#     return ulam, tpt_res
# end