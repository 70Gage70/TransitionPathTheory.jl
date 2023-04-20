using LinearAlgebra

function hitting_location_distribution(tpt::TPTHomog)
    S = tpt.sets.S
    S_plus = tpt.sets.S_plus
    B_true = tpt.sets.B_true
    P_plus = tpt.P_plus

    # solve the linear algebra problem r = P_plus r + b restricted to outside_B
    outside_B = setdiff(S_plus, B_true)
    M = I - P_plus[outside_B, outside_B]

    @assert det(M) != 0.0 "I - M in the hitting distribution is not invertible, check that P is well defined."

    # the probability that state i hits B at index j
    rij = zeros(length(S), length(S))
    for Bind in B_true
        b = P_plus[outside_B, Bind]
        rij[outside_B, Bind] = M\b
    end

    # starting at any given B_true index, guaranteed to hit that 
    rij[B_true, B_true] = [i == j ? 1.0 : 0.0 for i = 1:length(B_true), j = 1:length(B_true)]
    
    # trim very small (possibly negative due to floats) values
    rij[abs.(rij) .< 1e-16] .= 0.0

    # the most likely index for state i to hit
    ri = [i in S_plus ? argmax(rij[i,:]) : 0 for i in S]

    return rij, ri
end