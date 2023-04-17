using LinearAlgebra

"""
    time_cdf(tpt; cdf_thresh)

Compute the cdf of the reactive hitting time from any state, as well as the weighted cdf time from A.

### Optional Arguments
- `cdf_thresh`: The cdf value at which the calculation should stop, default 0.5.
"""
function time_cdf(tpt::TPTHomog; cdf_thresh = 0.95)
    @assert 0.0 < cdf_thresh < 1.0

    S = tpt.sets.S
    A_true = tpt.sets.A_true
    B_true = tpt.sets.B_true
    C = tpt.sets.C
    C_plus = tpt.sets.C_plus
    qp = tpt.q_plus
    P = tpt.P
    pi_stat = tpt.pi_stat

    # solve the linear algebra problem t = M t + b
    # denominator of M necessarily positive by the definition of C_plus
    M = [qp[j]*P[i, j]/sum(P[i, k]*qp[k] for k in S) for i in C_plus, j in C_plus]
    b = [sum(qp[j]*P[i, j]/sum(P[i, k]*qp[k] for k in S) for j in B_true) for i in C_plus]
    
    cdf_dist = zeros(length(S), 2) # probability of reaching B in 0 steps is 0.0
    cdf_dist[C_plus] = b # probability of reaching B in 1 step is b

    A_cdf = zeros(length(S))

    Mpow = I
    counter = 1
    while A_cdf[end] < cdf_thresh
        counter = counter + 1
        Mpow = Mpow * M
        cdf_dist = hcat(cdf_dist, cdf_dist[:, end])
        cdf_dist[C_plus, end] = cdf_dist[C_plus, end] + Mpow * b

        A_cdf = sum(pi_stat[i]*cdf_dist[i, :]/sum(pi_stat[j] for j in A_true) for i in A_true)

        if counter >= 500
            @warn "The time cdf is converging very slowly. Breaking after 500 terms."
            break
        end
    end

    # cdf to go from A to be is cdf_dist[i] weighted by pi[i] for i in A
    # A_cdf = sum(pi_stat[i]*cdf_dist[i, :]/sum(pi_stat[j] for j in A_true) for i in A_true)

    return cdf_dist, A_cdf
end