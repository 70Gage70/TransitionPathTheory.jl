using LinearAlgebra

"""
    remaining_time(tpt)

Compute the remaining time for each state in S.
"""
function remaining_time(tpt::TPTHomog)
    S = tpt.sets.S 
    S_plus = tpt.sets.S_plus
    B_true = tpt.sets.B_true
    P_plus = tpt.P_plus

    # solve the linear algebra problem t = P_plus t + b restricted to outside_B
    outside_B = setdiff(S_plus, B_true)
    M = I - P_plus[outside_B, outside_B]
    b = [1.0 for i in outside_B]
    
    @assert det(M) != 0.0 "I - M in the remaining time is not invertible, check that P is well defined."
    sol = M\b

    # assign the values to a vector; note that t_rem[B] = 0.0 is already handled
    t_rem = zeros(length(S))
    t_rem[outside_B] = sol
    
    # trim very small (possibly negative due to floats) values
    t_rem[abs.(t_rem) .< 1e-16] .= 0.0

    return t_rem       
end