include("../types_loader.jl")
include("../workspaces.jl")
include("partitions.jl")

mp_stat = minimal_partitions(tpt_homog, tpt_stationary)
mp_nonstat = minimal_partitions(tpt_homog, tpt_nonstationary)