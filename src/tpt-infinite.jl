using HDF5

include("statistics/infinite/remaining-time.jl")
include("statistics/infinite/hitting-distribution.jl")
include("statistics/infinite/time-cdf.jl")

"""
    tpt_infinite(tpt)

Compute TPT statistics in the infinite time case.
"""
function tpt_infinite(tpt::TPTHomog)
    P = tpt.P
    qp = tpt.q_plus
    qm = tpt.q_minus
    pi_stat = tpt.pi_stat
    A = tpt.sets.A
    S = tpt.sets.S

    # reactive density
    muAB = [qm[i]*pi_stat[i]*qp[i] for i in S]

    # reactive current
    fij = [qm[i]*pi_stat[i]*P[i, j]*qp[i] for i in S, j in S]

    # forward current
    fplusij = [max(fij[i, j] - fij[j, i], 0) for i in S, j in S]

    # vanE rate
    kAB = sum(fij[i, j] for i in A, j in S)

    # vanE time
    tAB = sum(muAB)/kAB

    # remaining time
    tAB_rem = remaining_time(tpt)

    # hitting distribution
    rij, ri = hitting_distribution(tpt)

    # time cdf
    tcdf, tcdf_AB = time_cdf(tpt)

    res = TPTStats(
        tpt.q_plus,
        tpt.q_minus,
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
    tpt_write(outfile, tpt_homog, tpt_res)

Write TPT indices from `tpt_homog` and statistics from `tpt_res` to the file `outfile`, which must be in the `.h5` format.

The indices are found in `tpt_homog/inds` and the results are in `tpt_homog/stats`.
"""
function tpt_write(outfile::String, tpt_homog::TPTHomog, tpt_res::TPTStats)
    @assert outfile[end-2:end] == ".h5" "The output file must be of the form filename.h5"

    fout = h5open(outfile, "w")

    tpt = create_group(fout, "tpt_homog")

    inds = create_group(tpt, "inds")
    for fn in fieldnames(typeof(tpt_homog.sets))
        inds["$fn"] = getfield(tpt_homog.sets, fn)
    end

    stats = create_group(tpt, "stats")
    for fn in fieldnames(typeof(tpt_res))
        stats["$fn"] = getfield(tpt_res, fn)
    end

    close(fout)

    @info "TPT results written to $(outfile)."

    return 
end