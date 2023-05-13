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