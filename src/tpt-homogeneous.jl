using HDF5

include("statistics/homogeneous/remaining-time.jl")
include("statistics/homogeneous/hitting-location-distribution.jl")

"""
    tpt_stationary_statistics(tpt_homog)

Compute the stationary TPT statistics in the homogenous case.
"""
function tpt_stationary_statistics(tpt_homog::TPTHomog)
    P = tpt_homog.P
    qp = tpt_homog.q_plus
    qm = tpt_homog.q_minus
    pi_stat = tpt_homog.pi_stat
    A = tpt_homog.sets.A
    S = tpt_homog.sets.S

    # reactive density
    muAB = [qm[i]*pi_stat[i]*qp[i] for i in S]

    # reactive current
    fij = [qm[i]*pi_stat[i]*P[i, j]*qp[j] for i in S, j in S]

    # forward current
    fplusij = [max(fij[i, j] - fij[j, i], 0) for i in S, j in S]

    # vanE rate
    kAB = sum(fij[i, j] for i in A, j in S)

    # vanE time
    tAB = sum(muAB)/kAB

    # remaining time
    tAB_rem = remaining_time(tpt_homog)

    # hitting distribution
    rij, ri = hitting_location_distribution(tpt_homog)

    res = TPTHomogStatResult(
        tpt_homog.sets,
        pi_stat,
        tpt_homog.q_plus,
        tpt_homog.q_minus,
        muAB,
        muAB/sum(muAB),
        fij,
        fplusij,
        kAB,
        tAB,
        tAB_rem,
        rij,
        ri
    )

    return res
end


"""
    tpt_nonstationary_statistics(tpt_homog; horizon = 100)

Compute the nonstationary TPT statistics in the homogenous case.

### Optional Arguments
- `"horizon"`: The time step at which to cut off the calculation. Default `100`. Note that the `horizon` value does NOT enforce that trajectories leaving A hit B by that time.
"""
function tpt_nonstationary_statistics(tpt_homog::TPTHomog; horizon::Integer = 100)
    P = tpt_homog.P
    P_plus = tpt_homog.P_plus
    # qp = tpt_homog.q_plus
    # qm = tpt_homog.q_minus
    pi_stat = tpt_homog.pi_stat
    A_true = tpt_homog.sets.A_true
    B_true = tpt_homog.sets.B_true
    S = tpt_homog.sets.S
    S_plus = tpt_homog.sets.S_plus

    i0 = [i in A_true ? 1.0/length(A_true) : 0.0 for i in S] # uniform distribution supported on A_true

    density = zeros(horizon, length(S)) # increasing time is down the matrix
    density[1, :] = i0

    muABnorm = zeros(horizon, length(S))
    muABnorm[1, :] = i0

    outside_B_true = setdiff(S_plus, B_true)
    t_cdf = zeros(horizon, length(outside_B_true))
    t_cdf[1, :] = [sum(P_plus[i, j] for j in B_true) for i in outside_B_true] 

    # currents under construction
    fij = zeros(1, 1, 1)
    fplusij = zeros(1, 1, 1)

    for n = 2:horizon
        # density
        density[n, :] = transpose(P) * density[n - 1, :]  

        # normalized reactive density
        muABnorm[n, :] = transpose(P_plus) * muABnorm[n - 1, :] 

        # hitting time distribution
        t_cdf[n, :] = t_cdf[1, :] + P_plus[outside_B_true, outside_B_true] * t_cdf[n - 1, :]
        
        # currents under construction

        # reactive current
        # fij[n, :, :] = [qm[i]*pn[i]*P[i, j]*qp[j] for i in S, j in S]
        # fij[n, :, :] = zeros(length(S), length(S))

        # forward current
        # fplusij[n, :, :] = [max(fij[n, i, j] - fij[n, j, i], 0) for i in S, j in S]
        # fplusij[n, :, :] = zeros(length(S), length(S))
    end

    # create full time cdf
    # the cdf for hitting B_true at B_true should be ones
    t_cdf_full = zeros(horizon, length(S))
    t_cdf_full[:, B_true] .= 1.0
    t_cdf_full[:, outside_B_true] .= t_cdf

    # create time cdf for A to B
    t_cdf_AB = [sum(i0[i]*t_cdf_full[n, i] for i in A_true) for n = 1:horizon]

    res = TPTHomogNonStatResult(
        tpt_homog.sets,
        horizon,
        pi_stat,
        tpt_homog.q_plus,
        tpt_homog.q_minus,
        density,
        muABnorm,
        t_cdf_full,
        t_cdf_AB,
        fij,
        fplusij
    )

    return res
end


"""
    tpt_write(outfile, tpt_result; dir_name, overwrite, P_out)

Write `tpt_result` to the file `outfile`, which must be in the `.h5` format.

If `outfile` does not exist, it will be created. Results are written to the directory specified by `dir_name`.

Indices are written to the sub-directory `"indices"` and statistics to the sub-directory `"statistics"`.

### Optional Arguments
- `dir_name`: The name of the directory that `tpt_result` is written to, default `"tpt_homog"`. Directories can be nested, e.g. `"trial1/tpt"`.
- `overwrite`: If `true`, the directory `dir_name` will overwrite a directory with the same name if it exists in `outfile`. Default `false`.
"""
function tpt_write(
    outfile::String, 
    tpt_result::AbstractTPTHomogResult;
    dir_name::String = "tpt_homog", 
    overwrite::Bool = false)

    @assert outfile[end-2:end] == ".h5" "The output file must be of the form filename.h5"
    fout = h5open(outfile, "cw")

    if dir_name in keys(fout)
        if !overwrite
            @assert !(dir_name in keys(fout)) "This file already has a group with the name: $(dir_name). Pass `overwrite = true` to force a replacement."
        end

        delete_object(fout, dir_name)
    end

    tpt = create_group(fout, dir_name)
    inds = create_group(tpt, "indices")
    stats = create_group(tpt, "statistics")

    for fn in fieldnames(typeof(tpt_result))
        if fn == :sets
            for fns in fieldnames(typeof(tpt_result.sets))
                inds["$fns"] = getfield(tpt_result.sets, fns)
            end
        else
            stats["$fn"] = getfield(tpt_result, fn)
        end
    end

    close(fout)

    @info "TPT results written to $(outfile)."

    return 
end