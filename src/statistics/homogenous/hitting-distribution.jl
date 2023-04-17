using LinearAlgebra

function hitting_distribution(tpt::TPTHomog) 
    S = tpt.sets.S
    B_true = tpt.sets.B_true
    C_plus = tpt.sets.C_plus
    qp = tpt.q_plus
    P = tpt.P

    # solve the linear algebra problem t = M t + b
    # denominator of M necessarily positive by the definition of C_plus
    M = I - [qp[j]*P[i, j]/sum(P[i, k]*qp[k] for k in S) for i in C_plus, j in C_plus]

    @assert det(M) != 0.0 "I - M in the hitting distribution is not invertible, check that P is well defined."

    # the probability that state i hits B at index j
    rij = zeros(length(S), length(S))
    for Bind in B_true
        b = [qp[Bind]*P[i, Bind]/sum(P[i, k]*qp[k] for k in S) for i in C_plus] # qp[Bind] = 1
        sol = M\b

        rij[C_plus, Bind] = sol
    end
    
    # trim very small (possibly negative due to floats) values
    rij[abs.(rij) .< 1e-16] .= 0.0

    # the most likely index for state i to hit
    ri = [i in C_plus ? argmax(rij[i,:]) : 0 for i in S]

    return rij, ri
end