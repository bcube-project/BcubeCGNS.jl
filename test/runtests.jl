using Test
using Bcube, BcubeCGNS, BcubeVTK
using SHA
using DelimitedFiles

"""
Custom way to "include" a file to print infos.
"""
function custom_include(path)
    filename = split(path, "/")[end]
    print("Running test file " * filename * "...")
    include(path)
    println("done.")
end

# Assets
const assets_path = joinpath(@__DIR__, "assets")

# This dir will be removed at the end of the tests
tempdir = mktempdir()

# Reading sha1 checksums.
#
# The sha1 of the VTU files produced by the "Read" tests depends on the Julia
# version: e.g. Julia 1.12 (GitHub CI) and Julia 1.13 (local workstation)
# produce different bytes through the WriteVTK/BcubeVTK writing chain, with
# the same inputs (this has been checked on the same machine, only switching
# the Julia version). The checksum file corresponding to the Julia version
# running the tests is therefore selected:
#   - checksums-1.12.sha1 : reference for Julia 1.12 (GitHub CI)
#   - checksums-1.13.sha1 : reference for Julia 1.13 (local workstation)
version = "$(VERSION.major).$(VERSION.minor)"
filename = "checksums-$(version)"
filepath = joinpath(@__DIR__, "$(filename).sha1")
if !isfile(filepath)
    error(
        "No sha1 checksum file for Julia $(version): expected \"$(filename).sha1\" " *
        "in the test directory. To create it, duplicate an existing " *
        "\"checksums-*.sha1\" file, run the tests, and replace its sha1s by the " *
        "actual values displayed by the failing tests.",
    )
end
f = readdlm(filepath, String)
fname2sum = Dict(r[2] => r[1] for r in eachrow(f))

@testset "BcubeGmsh.jl" begin
    custom_include("./test_read.jl")
    custom_include("./test_write.jl")
    custom_include("./test_cgns_mesh.jl")
end
