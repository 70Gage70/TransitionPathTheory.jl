"""
    current2arrows(f, x, y; excluded = nothing)

Compute forward current vectors in a 2D physical space.

Return `(u, v)` where `u[i]` and `v[i]` are the `x` and `y` components of vector centered at \
`(x, y)` and pointing in the direction of the weighted average by `f[i,:]` of each other point.

### Arguments

- `f`: `forward_current` computed from [`stationary_statistics`](@ref).
- `x`: A `Vector` of points `[x1, x2, ...]`
- `y`: A `Vector` of points `[y1, y2, ...]`

### Optional Arguments

- `excluded`: If a `Vector{<:Integer}` is provided, rows and columns corresponding to them are excluded from `f`. When \
`excluded === nothing`, the last index of `f` is excluded automatically.
"""
function current2arrows(
    f::Matrix{<:Real}, 
    points::Vector{<:Vector{<:Real}};
    excluded::Union{Nothing, Vector{<:Integer}} = nothing)

    @argcheck size(f, 1) == size(f, 2) "f must be a square matrix"
    f_size = size(f, 1)

    if isnothing(excluded)
        idx = collect(1:f_size - 1) 
    else
        @argcheck all(map(x -> 1 <= x <= f_size, excluded)) "indices must be positive numbers in the correct range"
        idx = setdiff(1:f_size, excluded)
    end

    @argcheck length(idx) == length(points) "must have as many points as valid f indices"
 
    f_tilde = view(f, idx, idx)
    uv = [sum(f_tilde[i,j]*normalize(points[j] - points[i]) for j = setdiff(1:length(points), [i])) for i = 1:length(points)]

    return ((x -> x[1]).(uv), (x -> x[2]).(uv))
end