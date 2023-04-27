"""
    P_stoc(n)

Generate a stochastic matrix of size `n` by `n` with random entries. Used for debugging and testing.
"""
function P_stoc(n)
    P = rand(n, n)
    P = P ./ sum(P, dims = 2)

    return P
end