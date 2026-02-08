module BcubeCGNS
using Bcube
using HDF5
using JLD2
using SparseArrays

include("./common.jl")

# CGNS tree
include("./cgns/CGNS.jl")
using .CGNS

# HDF5 reader/writer
include("./hdf5/common.jl")
include("./hdf5/read.jl")
include("./hdf5/write.jl")

# JLD2 reader/writer
include("./jld2/common.jl")
include("./jld2/read.jl")
include("./jld2/write.jl")

# Other
include("./hdf5/eclipps_to_surf.jl")

# (Bcube) CGNS Mesh
include("./mesh/cgns-mesh.jl")

end
