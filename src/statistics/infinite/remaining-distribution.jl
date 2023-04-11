function remaining_distribution(tpt::TPTHomog; cdf_thresh = 0.95)
    S = tpt.sets.s
    B = tpt.sets.B
    C = tpt.sets.C
    qp = tpt.q_plus
    P = tpt.P

    # define the set C_plus; the subset of states i of C such that B is reactive-reachable from i
    C_plus = [C[i] for i in S if qp[i] > 0.0]

    M = zeros(length(S), length(S))
    M[Cplus, C_plus] = [qp[j]*P[i, j]/qp[i] for i in C_plus, j in C_plus]
    
    v = zeros(length(S))
    v[C_Plus] = [sum(qp[j]*P[i, j]/qp[i] for j in B) for i in C_plus]
end