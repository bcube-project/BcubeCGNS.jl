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

"""
Only valid for Lagrange
"""
function build_node_to_idof(mesh, U)
    ncomps = Bcube.get_ncomponents(U)
    node2idof = zeros(Int, nnodes(mesh), ncomps)
    fs = Bcube.get_function_space(U)
    dhl = Bcube._get_dhl(U)
    for cellInfo in Bcube.DomainIterator(CellDomain(mesh))
        shape = Bcube.shape(Bcube.celltype(cellInfo))
        icell = Bcube.cellindex(cellInfo)
        c2n = Bcube.get_nodes_index(cellInfo)
        for (ivertex, idofs_l) in enumerate(Bcube.idof_by_vertex(fs, shape))
            @assert length(idofs_l) == 1
            for icomp in 1:ncomps
                node2idof[c2n[ivertex], icomp] =
                    Bcube.get_dof(dhl, icell, icomp, idofs_l[1]) # there is only one dof per vertex with Lagrange
            end
        end
    end
    return node2idof
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
