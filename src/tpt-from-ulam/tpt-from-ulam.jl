"""
    remove_omega(tpt_result)

Removes the last state from indices and statistics of tpt_result. Used for applications to Ulam's method.

Note that this does not change the  actual results of the TPT calculation, it is only for convenient reference
to  TPT statistics without the omega state.
"""
function remove_omega(tpt_result::TPTHomogStatResult)
    omega_tpt_result = []

    for fn in fieldnames(TPTHomogStatResult)
        if fn == :sets
            omega = [tpt_result.sets.S[end]]
            omega_sets = TPTSets([setdiff(getfield(tpt_result.sets, fn), omega) for fn in fieldnames(TPTSets)]...)
            push!(omega_tpt_result, omega_sets)
        else
            gf = getfield(tpt_result, fn)

            if typeof(gf) <: AbstractVector
                push!(omega_tpt_result, gf[1:end-1])
            elseif typeof(gf) <: AbstractMatrix
                if size(gf, 1) == size(gf, 2)
                    push!(omega_tpt_result, gf[1:end-1, 1:end-1])
                else
                    push!(omega_tpt_result, gf[1:end-1, :])
                end
            else
                push!(omega_tpt_result, gf) # quantity is a scalar, can't remove anything
            end
        end
    end
    
    omega_tpt_result = TPTHomogStatResult(omega_tpt_result...)

    return omega_tpt_result
end


"""
    tpt_from_ulam(s)

Assumes UlamMethod.jl is loaded.

Default is A is centers but B defines polygon.
"""
function tpt_from_ulam(
    fin::String, 
    domain::Matrix{<:Real}, 
    poly_type::String, 
    poly_number::Integer, 
    A_centers::Matrix{<:Real},
    B_vertices::Matrix{<:Real};
    avoid_vertices::Union{Matrix{<:Real}, Nothing} = nothing,
    tpt_out::Union{String, Nothing} = nothing,
    ulam_out::Union{String, Nothing} = nothing)

    traj = UlamTrajectories(fin)
    domain = UlamDomain(UlamPolygon(domain), poly_type = poly_type, poly_number = poly_number)
    ulam = ulam_method(traj, domain)
    
    A = inpoly(A_centers, PolyTable(ulam.polys)).inds
    B = [i for i in 1:length(ulam.polys) if ulam_intersects(UlamPolygon(B_vertices), ulam.polys[i])]
    
    if avoid_vertices !== nothing
        avoid = [i for i in 1:length(ulam.polys) if ulam_intersects(UlamPolygon(avoid_vertices), ulam.polys[i])]
        A = sort(union(A, avoid))
        B = sort(union(B, avoid))
    end

    P = ulam.P_closed
    omega = [size(P, 1)]
    A = union(A, omega)
    B = union(B, omega)
    
    tpt_homog = TPTHomog(P, A, B)
    tpt_res = tpt_infinite(tpt_homog)
    
    tpt_homog, tpt_res = stats_remove_omega(tpt_homog, tpt_res)

    if ulam_out !== nothing
        ulam_write(ulam_out, ulam)
    end
    
    if tpt_out !== nothing
        tpt_write(tpt_out, tpt_homog, tpt_res)
    end

    return ulam, tpt_res
end