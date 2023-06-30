using LinearAlgebra

"""
    q_plus(P, A, B, C)

Compute the forward committor q_plus.

`P` is the transition matrix and `A`, `B` and `C` are subsets of `S`, the index set of `P`.
"""
function q_plus(
    P::Matrix{T}, 
    A::Vector{U}, 
    B::Vector{U}, 
    C::Vector{U}) where {T<:Real, U<:Integer}

    @assert isempty(intersect(A, B)) "A and B should not intersect in q."
    @assert isempty(intersect(C, union(A, B))) "C should not intersect with A or B."

    # solve the linear algebra problem q = P q + b on the C subspace
    M = I - P[C, C] 
    b = [sum(P[i, k] for k in B) for i in C]
    sol = M\b
    
    # assign the values to a vector; note that q[A] = 0.0 is already handled
    q = zeros(size(P, 1))
    q[B] .= 1.0
    q[C] = sol
    
    # trim very small (possibly negative due to floats) values
    q[abs.(q) .< 1e-16] .= 0.0

    return q           
end

"""
    q_minus(P_minus, A, B, C)

Compute the backward committor q_minus.

`P_minus` is the reversed transition matrix and `A`, `B` and `C` are subsets of `S`, the index set of `P_minus`.
"""
function q_minus(
    P_minus::Matrix{T}, 
    A::Vector{U}, 
    B::Vector{U}, 
    C::Vector{U}) where {T<:Real, U<:Integer}

    @assert isempty(intersect(A, B)) "A and B should not intersect in q."
    @assert isempty(intersect(C, union(A, B))) "C should not intersect with A or B."

    # The calculation is identical to q_plus, except we use P_minus instead of P 
    # and switch the roles of A and B
    M = I - P_minus[C, C] 
    b = [sum(P_minus[i, k] for k in A) for i in C]
    sol = M\b
    
    # assign the values to a vector; note that q[B] = 0.0 is already handled
    q = zeros(size(P_minus, 1))
    q[A] .= 1.0
    q[C] = sol
    
    # trim very small (possibly negative due to floats) values
    q[abs.(q) .< 1e-16] .= 0.0

    return q             
end