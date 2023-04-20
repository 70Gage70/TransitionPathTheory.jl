"""
    TPTSets{U}

A container for `S`, `A`, `B` and their subsets.    
"""
struct TPTSets{U<:Integer}
    S::Vector{U}
    A::Vector{U}
    B::Vector{U}
    C::Vector{U}
    A_true::Vector{U}
    B_true::Vector{U}
    S_plus::Union{Vector{U}, Nothing}
end

"""
    TPTSets(S, A, B; S_plus)

Construct all the relevant TPT sets given `S`, `A`, `B and optionally `S_plus`.  
"""
function TPTSets(
    S::Vector{U}, 
    A::Vector{U}, 
    B::Vector{U};
    S_plus::Union{Vector{T}, Nothing} = nothing) where {T<:Real, U<:Integer}

    # Scrub any accidental repeated indices
    S = unique(S)
    A = unique(A)
    B = unique(B)    

    # ensure sets are well-defined
    @assert length(A) > 0 "A can not be empty."
    @assert length(B) > 0 "B can not be empty."
    @assert minimum(A) > 0 "Elements of A must be positive."
    @assert minimum(B) > 0 "Elements of B must be positive."
    @assert minimum(S) > 0 "Elements of S must be positive."
    @assert issubset(A, S) "A must be a subset of S."
    @assert issubset(B, S) "B must be a subset of S."
    @assert !issubset(A, B) "A can not be a subset of B."
    @assert !issubset(B, A) "B can not be a subset of A."    

    C = setdiff(S, union(A, B))
    @assert !isempty(C) "A and B can not cover the entire state space."

    # the "true" A and B are indices unique to A and B
    # note that A_true = A, B_true = B if they do not intersect
    A_true = setdiff(A, intersect(A, B))
    B_true = setdiff(B, intersect(A, B))

    if S_plus !== nothing
        @assert !isempty(intersect(A_true, S_plus)) "There are no states in A_true that can reach reactively."
        @assert issubset(B_true, S_plus) "B_true must be a subset of S_plus."

        if !issubset(A_true, S_plus)
            @info "A_true has interior states. These are ignored by TPT."
        end
    end

    return TPTSets(S, A, B, C, A_true, B_true, S_plus)
end

function Base.show(io::IO, x::TPTSets)
    print(io, "TPTSets[S(")
    show(io, length(x.S))
    print(io, "), A(")
    show(io, length(x.A))
    print(io, "), B(")
    show(io, length(x.B))
    print("), AB_int(")
    show(io, length(intersect(x.A, x.B)))
    print(")]")
end