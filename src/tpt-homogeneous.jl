using HDF5

include("statistics/homogeneous/remaining-time.jl")
include("statistics/homogeneous/hitting-distribution.jl")
include("statistics/homogeneous/time-cdf.jl")

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
        tcdf_AB
    )

    return res
end


"""
    tpt_nonstationary_statistics(tpt_homog)

Compute the nonstationary TPT statistics in the homogenous case.
"""
function tpt_nonstationary_statistics(tpt_homog::TPTHomog; horizon::Integer = 100)
    P = tpt_homog.P
    qp = tpt_homog.q_plus
    qm = tpt_homog.q_minus
    pi_stat = tpt_homog.pi_stat
    A_true = tpt_homog.sets.A_true
    S = tpt_homog.sets.S

    i0 = [i in A_true ? 1.0/length(A_true) : 0.0 for i in S] # uniform distribution supported on A_true

    density = Matrix{Float64}(undef, horizon, length(S)) # increasing time is down the matrix
    density[1, :] = i0

    muABnorm = Matrix{Float64}(undef, horizon, length(S))
    muABnorm[1, :] = i0
    P_plus = [P[i, j]*qp[j]/sum(P[i, k]*qp[k] for k in S) for i in S, j in S]

    # currents under construction
    fij = Array{Float64, 3}(undef, horizon, length(S), length(S))
    fplusij = Array{Float64, 3}(undef, horizon, length(S), length(S))

    for n = 2:horizon
        # density
        density[n, :] = transpose(P) * density[n - 1, :]  

        # normalized reactive density
        muABnorm[n, :] = transpose(P_plus) * muABnorm[n - 1, :] 

        # currents under construction

        # reactive current
        # fij[n, :, :] = [qm[i]*pn[i]*P[i, j]*qp[j] for i in S, j in S]
        fij[n, :, :] = zeros(length(S), length(S))

        # forward current
        # fplusij[n, :, :] = [max(fij[n, i, j] - fij[n, j, i], 0) for i in S, j in S]
        fplusij[n, :, :] = zeros(length(S), length(S))
    end

    res = TPTHomogNonStatResult(
        tpt_homog.sets,
        pi_stat,
        tpt_homog.q_plus,
        tpt_homog.q_minus,
        density,
        muABnorm,
        fij,
        fplusij
    )

    return res
end


"""
    tpt_write(outfile, tpt_result)

Write the set indices and statistics from `tpt_result` to the file `outfile`, which must be in the `.h5` format.

The indices are found in `tpt_homog/indices` and the results are in `tpt_homog/statistics`.
"""
function tpt_write(outfile::String, tpt_result::AbstractTPTHomogResult)
    @assert outfile[end-2:end] == ".h5" "The output file must be of the form filename.h5"
    fout = h5open(outfile, "w")

    tpt = create_group(fout, "tpt_homog")
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