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
    lagrange_P1_node_to_idof(mesh)

WARNING : I'm facing situations where there are mesh nodes not belonging to any element.
So the "node2idof" and "idof2node" are not strictly a permutation.
"""
function lagrange_P1_node_to_idof(mesh)
    fs = FunctionSpace(:Lagrange, 1)
    U = TrialFESpace(fs, mesh)
    dhl = Bcube._get_dhl(U)
    node2idof = zeros(Int, nnodes(mesh))
    idof2node = zeros(Int, get_ndofs(U))
    for cellInfo in Bcube.DomainIterator(CellDomain(mesh))
        shape = Bcube.shape(Bcube.celltype(cellInfo))
        icell = Bcube.cellindex(cellInfo)
        c2n = Bcube.get_nodes_index(cellInfo)
        for (ivertex_l, idofs_l) in enumerate(Bcube.idof_by_vertex(fs, shape))
            ivertex_g = c2n[ivertex_l]
            idof_g = Bcube.get_dof(dhl, icell, 1, idofs_l[1]) # there is only one dof per vertex with Lagrange
            node2idof[ivertex_g] = idof_g
            idof2node[idof_g] = ivertex_g
        end
    end
    return node2idof, idof2node
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
