using LinearAlgebra

"""
    remaining_time(tpt)

Compute the remaining time for each state in S.
"""
function remaining_time(tpt::TPTHomog) 
    S = tpt.sets.S
    C_plus = tpt.sets.C_plus
    qp = tpt.q_plus
    P = tpt.P

    # solve the linear algebra problem t = M t + b
    # denominator of M necessarily positive by the definition of C_plus
    M = I - [qp[j]*P[i, j]/sum(P[i, k]*qp[k] for k in S) for i in C_plus, j in C_plus]
    b = [1.0 for i in C_plus]
    
    @assert det(M) != 0.0 "I - M in the remaining time is not invertible, check that P is well defined."
    sol = M\b

    # assign the values to a vector; note that t_rem[B] = 0.0 is already handled
    t_rem = zeros(length(S))
    t_rem[C_plus] = sol
    
    # trim very small (possibly negative due to floats) values
    t_rem[abs.(t_rem) .< 1e-16] .= 0.0

    return t_rem        
end