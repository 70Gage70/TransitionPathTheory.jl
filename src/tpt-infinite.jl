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

function Pstoc(n)
    P = rand(n, n)
    for i = 1:n
        P[i,:] = P[i, :]/sum(P[i, :])
    end

    return P
end