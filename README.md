# BcubeCGNS.jl
Implementation of [`Bcube`](https://github.com/bcube-project/Bcube.jl) IO interface for CGNS format. Checkout the relative `Bcube` [documentation](https://bcube-project.github.io/Bcube.jl/stable/api/io/io_interface/) for more infos.

Note that all CGNS specifications are not implemented (for instance, NGON elements).

## Basic usage
```julia
using Bcube
using BcubeCGNS

mesh = rectangle_mesh(10, 20)
U = TrialFESpace(FunctionSpace(:Lagrange, 1), mesh)
u = FEFunction(U)
projection_l2!(u, PhysicalFunction(x -> sum(x)), CellDomain(mesh))

write_file("output.cgns", mesh, Dict("u" => u, "grad_u" => ∇(u)))

result = read_file("output.cgns"; varnames = "*")
mesh = result.mesh
data = result.data # Dict of variable name (String) to MeshData
```
