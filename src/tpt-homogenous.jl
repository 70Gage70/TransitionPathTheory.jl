using HDF5

include("statistics/homogenous/remaining-time.jl")
include("statistics/homogenous/hitting-distribution.jl")
include("statistics/homogenous/time-cdf.jl")

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
    rij, ri = hitting_distribution(tpt_homog)

    # time cdf
    tcdf, tcdf_AB = time_cdf(tpt_homog)

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
        ri,
        tcdf,
        tcdf_AB,
    )

    return res
end

"""
    tpt_write(outfile, tpt_result)

Write the set indices and statistics from `tpt_result` to the file `outfile`, which must be in the `.h5` format.

The indices are found in `tpt_homog_stationary/indices` and the results are in `tpt_homog_stationary/statistics`.
"""
function tpt_write(outfile::String, tpt_result::TPTHomogStatResult)
    @assert outfile[end-2:end] == ".h5" "The output file must be of the form filename.h5"
    fout = h5open(outfile, "w")

    tpt = create_group(fout, "tpt_homog_stationary")
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