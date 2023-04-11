"""
Infinite time, stationary TPT

Written by Gage Bonner September 2022
"""

###################################################################################################
###################################################################################################
###################################################################################################

using LinearAlgebra


"""
Calculates y - x, setwise. 
e.g.  mysetdiff([1,2,3,4,5], [1,2,3]) = [4,5]
e.g.  mysetdiff([1,2,3,4,5], [1,2,2,3]) = [4,5]
"""

function mysetdiff(y, x)
    xu = unique(x)
    res = Vector{eltype(y)}(undef, length(y) - length(xu))
    i = 1
    @inbounds for el in y
        el ∈ xu && continue
        res[i] = el
        i += 1
    end

    res
end

"""
Left eigenvector of P, normalized appropriately. The stationary distribution of P if it's row stochastic.
"""

function Ppi(P)
    return transpose(abs.(normalize(eigvecs(transpose(P))[:,size(P)[1]], 1)))
end
 
"""
Transition matrix of the reversed Markov chain
"""

function Pminus(P, Ppi)
    Pmin = zeros(length(Ppi), length(Ppi))
    for i = 1:length(Ppi)
        for j = 1:length(Ppi)
            if Ppi[i] > 0.0 # else: already 0.0
                Pmin[i, j] = (Ppi[j]/Ppi[i])*P[j,i]
            end
        end
    end
    
    return Pmin
end

"""
Forward committor. Essentially have linear matrix equation q = P q + b, where P is restricted to the C subspace (i.e. the complement of A U B)
and where b_i = sum P_{i, k} for k in B. Hence, construct M = I - P and so q = M\b
"""

# Subtle issue with the committor that's relevant for the case where A and B contain intersections.
#We want the intersection to have 0 forward and backward committors, therefore we check the one
# that gives 0 first in the if statements.

function qforward(ALLinds, Ainds, Binds, P)

    Cinds = mysetdiff(ALLinds, [Ainds ; Binds])

    ABint = intersect(Ainds, Binds)
    trueB = mysetdiff(Binds, ABint)

    Pc = P[Cinds, Cinds] # P restricted to C subspace
    M = I - Pc # "I" just works in Julia as an identity matrix object in LinearAlgebra
    b = [sum(P[i, k] for k in trueB) for i in Cinds]
    
    sol = M\b
    
    qplus = zeros(length(ALLinds))
    for i = 1:length(ALLinds)
        if i in Ainds
            qplus[i] = 0.0
        elseif i in Binds
            qplus[i] = 1.0
        else
            solindex = findfirst(isequal(i), Cinds)
            qplus[i] = sol[solindex]
        end
    end
    
    qplus[abs.(qplus) .< 1e-16] .= 0.0

    return qplus           
end

"""
Backward committor. Same idea as forward committor except we use Pminus instead of P
"""

function qbackward(ALLinds, Ainds, Binds, Pminus)
    Cinds = mysetdiff(ALLinds, [Ainds ; Binds]) 

    ABint = intersect(Ainds, Binds)
    trueA = mysetdiff(Ainds, ABint)

    Pc = Pminus[Cinds, Cinds] # P restricted to C subspace
    M = I - Pc # "I" just works in Julia as an identity matrix object in LinearAlgebra
    b = [sum(Pminus[i, k] for k in trueA) for i in Cinds]
    
    sol = M\b
    
    qback = zeros(length(ALLinds))
    for i = 1:length(ALLinds)
        if i in Binds
            qback[i] = 0.0
        elseif i in Ainds
            qback[i] = 1.0
        else
            solindex = findfirst(isequal(i), Cinds)
            qback[i] = sol[solindex]
        end
    end
    
    qback[abs.(qback) .< 1e-16] .= 0.0

    return qback           
end

"""
The remaining time of a reactive trajectory (second formulation). Satisfies a linear matrix equation.
"""

function t_remaining(ALLinds, Ainds, Binds, P, qplus) 
    Cplus = []

    outsideBinds = mysetdiff(ALLinds, Binds)
    for i in outsideBinds
        if sum(P[i, k]*qplus[k] for k in ALLinds) > 0.0
            push!(Cplus, i)
        end
    end

    M = zeros(length(Cplus), length(Cplus))
    for i in 1:length(Cplus)
        for j in 1:length(Cplus)
            ci = Cplus[i]
            cj = Cplus[j]
            M[i, j] = qplus[cj]*P[ci,cj]/sum(P[ci, k]*qplus[k] for k in ALLinds)
        end
    end

    M = I - M
    b = [1.0 for i = 1:length(Cplus)]
    
    sol = M\b

    t_rem = zeros(length(ALLinds)) # = 0 at B already handled
    for i = 1:length(ALLinds)
        if i in Cplus
            solindex = findfirst(isequal(i), Cplus)
            t_rem[i] = sol[solindex]
        end
    end

    tAB_rem = t_rem[Ainds[1]]

    return tAB_rem, t_rem          
end


"""
The expected first passage time to B.
Assumes that P is P_closed
"""

function t_first_passage(P, Binds) 

    ALLinds = range(1, length(P[1, :]))
    Cinds = mysetdiff(ALLinds, Binds)
    M = P[Cinds, Cinds]
       
    M = I - M
    b = [1.0 for i = 1:length(Cinds)]
    
    sol = M\b

    t_fp = zeros(length(ALLinds))
    for i = 1:length(ALLinds)
        if i in Binds
            t_fp[i] = 0.0
        else
            t_fp[i] = sol[findfirst(isequal(i), Cinds)]
        end
    end

    return t_fp           
end

function tAB_dist(ALLinds, Ainds, Binds, P, piStat, qplus, cutoff = 0.99)
    Cplus = []
    n_inds = length(ALLinds)
    
    ABint = intersect(Ainds, Binds)
    trueB = mysetdiff(Binds, ABint)
    trueA = mysetdiff(Ainds, ABint)
    
    Cinds = mysetdiff(ALLinds, [Ainds ; Binds])
    for i in Cinds
        if qplus[i] > 0.0
            push!(Cplus, i)
        end
    end
    
    M = zeros(n_inds, n_inds)
    for i in Cplus
        for j in Cplus
            M[i, j] = qplus[j]*P[i,j]/qplus[i]
        end
    end
    
    v = zeros(n_inds)
    for i in Cplus
        v[i] = sum(P[i, j]*qplus[j]/qplus[i] for j in Binds)
    end
    
    Nbar = zeros(n_inds, n_inds)
    for i in trueA
        for j in Binds
            denom = sum(P[i, k]*qplus[k] for k in Cplus)
            if denom > 0.0
                Nbar[i, j] = qplus[j]*P[i,j]/denom
            end
        end
    end
    
    N = zeros(n_inds, n_inds)
    for i in trueA
        for j in Cplus
            denom = sum(P[i, k]*qplus[k] for k in Cplus)
            if denom > 0.0
                N[i, j] = qplus[j]*P[i,j]/denom
            end
        end
    end
    
    u = zeros(n_inds)
    for i in trueA
        u[i] = piStat[i]*sum(P[i, l]*qplus[l] for l in Cplus)/sum(piStat[k]*sum(P[k, l]*qplus[l] for l in Cplus) for k in trueA)
    end
    
    n_powers = Integer(ceil(log(1 - cutoff)/log(abs(eigvals(M)[end]))))
    
    dist = zeros(n_powers)
    
    dist[1] = u' * Nbar * [1 for i = 1:n_inds]
    
    Mpow = M
    for i = 2:length(dist)
        Mpow = Mpow * M
        dist[i] = u' * N * Mpow * v
    end
    
    return dist
end

"""
The probability that state i hits B at B_j in particular (given that it is reactive)
"""

function B_hitters(ALLinds, Ainds, Binds, P, qplus)
    Cplus = []
    Cinds = mysetdiff(ALLinds, [Ainds ; Binds])
    for i in Cinds
        if qplus[i] > 0.0
            push!(Cplus, i)
        end
    end

    rij = zeros(length(ALLinds), length(ALLinds)) # jth row is B_j

    ABint = intersect(Ainds, Binds)
    trueB = mysetdiff(Binds, ABint)

    if length(trueB) == 0
        error("B is a subset of A")
    end

    for j in trueB
        M = zeros(length(Cplus), length(Cplus))
        for i in 1:length(Cplus)
            for k in 1:length(Cplus)
                ci = Cplus[i]
                ck = Cplus[k]
                M[i, k] = qplus[ck]*P[ci,ck]/qplus[ci]
            end
        end

        M = I - M
        b = [P[i, j]/qplus[i] for i in Cplus]
        
        sol = M\b

        for i = 1:length(ALLinds)
            if i in Cplus
                solindex = findfirst(isequal(i), Cplus)
                rij[i,j] = sol[solindex]
            end
        end
    end

    # ri[i] is the index of the B cell that's most likely to be hit.
    ri = zeros(Int64, length(ALLinds))
    for i in ALLinds
        am = argmax(rij[i,:])
        if am == 1
            ri[i] = 0
        else
            ri[i] = am
        end
    end

    return rij, ri

end


"""
All the statistics for infinite time tpt, returned as a dict.
    "piStat" => stationary distribution [states]
    "q+" => forward committor [states]
    "q-" => backward committor [states] 
    "muAB" => un-normalized reactive density [states] 
    "muABnorm" => normalized reactive density [states] 
    "fij" => reactive current [state i, state j] 
    "f+" => effective reactive current [state i, state j] 
    "ZAB" => probability that trajectory is reactive
    "kAB" => reactive rate [scalar] 
    "tAB" => reactive time [scalar]
    "t_rem" => remaining time [states] 
    "t_fp" => first passage time [states]
    "tABdist" => distribution of reactive hitting times of B from A [states]
"""

function tpt_infinite_stats(ALLinds, Ainds, Binds, P, piStat, P_closed = P)
    ABint = intersect(Ainds, Binds)
    trueA = mysetdiff(Ainds, ABint)
    trueB = mysetdiff(Binds, ABint)

    if length(trueA) == 0
        error("A is a subset of B.")
    elseif length(trueB) == 0
        error("B is a subset of A.")
    end

#     piStat = Ppi(P)
    Pmin = Pminus(P, piStat)
    qplus = qforward(ALLinds, Ainds, Binds, P)
    qminus = qbackward(ALLinds, Ainds, Binds, Pmin)
    
    # reactive density
    muAB = [qminus[i]*piStat[i]*qplus[i] for i = 1:length(piStat)]
    
    # normalization factor
    ZAB = sum(muAB)
    
    # normalized reactive density
    muABnorm = muAB/ZAB
    
    # reactive current
    fijAB = zeros(length(piStat), length(piStat))
    for i = 1:length(piStat)
        for j = 1:length(piStat)
            fijAB[i, j] = qminus[i]*piStat[i]*P[i, j]*qplus[j]
        end
    end
    
    # effective reactive current
    fijplus = zeros(length(piStat), length(piStat))
    for i = 1:length(piStat)
        for j = 1:length(piStat)
            fijplus[i, j] = max(fijAB[i, j] - fijAB[j, i], 0)
        end
    end
    
    # reactive rate
    kAB = sum(fijAB[i, j] for i in Ainds for j in ALLinds)
    
    # reactive time
    tAB = ZAB/kAB

    # remaining time
    tAB_rem, t_rem = t_remaining(ALLinds, Ainds, Binds, P, qplus) 

    # first passage time
    t_fp = t_first_passage(P_closed, Binds) 
    
    # tAB hitting distribution
    tAB_distribution = tAB_dist(ALLinds, Ainds, Binds, P, piStat, qplus)

    # B hitters
    rij, ri = B_hitters(ALLinds, Ainds, Binds, P, qplus) 

    tptDict = begin Dict(
            "Ainds" => Ainds,
            "Binds" => Binds,
            "piStat" => piStat, 
            "q+" => qplus, 
            "q-" => qminus, 
            "muAB" => muAB, 
            "muABnorm" => muABnorm, 
            "fij" => fijAB, 
            "f+" => fijplus, 
            "kAB" => kAB, 
            "ZAB" => ZAB,
            "tAB" => tAB,
            "tAB_rem" => tAB_rem,
            "t_rem" => t_rem,
            "t_fp" => t_fp,
            "tAB_dist" => tAB_distribution,
            "rij" => rij,
            "ri" => ri) 
    end
    
    return tptDict
end