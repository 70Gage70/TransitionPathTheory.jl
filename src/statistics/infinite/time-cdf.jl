using LinearAlgebra

"""
    time_cdf(tpt; cdf_thresh)

Compute the cdf of the reactive hitting time from any state, as well as the weighted cdf time from A.

### Optional Arguments
- `cdf_thresh`: The cdf value at which the calculation should stop, default 0.95.
"""
function time_cdf(tpt::TPTHomog; cdf_thresh = 0.95)
    @assert 0.0 < cdf_thresh < 1.0

    S = tpt.sets.S
    A_true = tpt.sets.A_true
    B_true = tpt.sets.B_true
    C_plus = tpt.sets.C_plus
    qp = tpt.q_plus
    P = tpt.P
    pi_stat = tpt.pi_stat

    # solve the linear algebra problem t = M t + b
    # denominator of M necessarily positive by the definition of C_plus
    M = [qp[j]*P[i, j]/sum(P[i, k]*qp[k] for k in S) for i in C_plus, j in C_plus]
    b = [sum(qp[j]*P[i, j]/sum(P[i, k]*qp[k] for k in S) for j in B_true) for i in C_plus]
    
    cdf_dist = zeros(length(C_plus), 1) # probability of reaching B in 0 steps is 0.0
    cdf_dist = hcat(cdf_dist, b) # probability of reaching B in 1 step is b

    # higher terms in the distribution come from hitting b with powers of M
    Mpow = I
    while minimum(cdf_dist[:, end]) < cdf_thresh
        Mpow = Mpow * M
        cdf_dist = hcat(cdf_dist, cdf_dist[:, end] + Mpow * b)
    end

    # expand cdf with zeros so results are in S
    cdf_dist_S = zeros(length(S), size(cdf_dist, 2))
    cdf_dist_S[C_plus, :] = cdf_dist

    # cdf to go from A to be is cdf_dist[i] weighted by pi[i] for i in A
    A_cdf = sum(pi_stat[i]*cdf_dist_S[i, :]/sum(pi_stat[j] for j in A_true) for i in A_true)

    return cdf_dist_S, A_cdf
end