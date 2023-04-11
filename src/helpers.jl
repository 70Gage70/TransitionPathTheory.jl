"""
    P_stoc(outfile, tpt_homog, tpt_res)

Generate a stochastic matrix with random entries. Used for debugging and testing.
"""
function P_stoc(n)
    P = rand(n, n)
    for i = 1:n
        P[i,:] = P[i, :]/sum(P[i, :])
    end

    return P
end