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

# This dir will be removed at the end of the tests
tempdir = mktempdir()

# Reading sha1 checksums
filename = "checksums"
f = readdlm(joinpath(@__DIR__, "$(filename).sha1"), String)
fname2sum = Dict(r[2] => r[1] for r in eachrow(f))

@testset "BcubeGmsh.jl" begin
    custom_include("./test_read.jl")
end
