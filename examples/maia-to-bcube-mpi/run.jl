using Bcube, BcubeCGNS, BcubeMPI
using LinearAlgebra
using MPI
using HauntedArrays
using MPIUtils

MPI.Initialized() || MPI.Init()

function tree_to_bcube(tree_as_list, ghost_tag2part, comm_py2f)
    # Parse the input tree as list into a BcubeCGNS.CGNS.Node
    tree = BcubeCGNS.CGNS.parse_tree_as_list(tree_as_list)
    BcubeCGNS.CGNS.print_tree(tree)

    # Convert it to a Bcube.Mesh
    result = BcubeCGNS.read_tree(tree; verbose = true, varnames = "*")

    # Use additionnal infos to build the BcubeMPI.DistributedMesh
    comm = MPI.Comm(comm_py2f)
    dmesh = DistributedMesh(result.mesh, ghost_tag2part, comm)

    return (; result..., dmesh)
end

"""
    setup_simulation(bcube_io, c::AbstractVector)

Set up a simulation using mesh IO data and shared coefficient storage.

The example preserves the coefficient array in memory.
"""
function setup_simulation(bcube_io, c::AbstractVector)
    dmesh = bcube_io.dmesh

    degree = 0
    fs = FunctionSpace(:Lagrange, degree)
    U = TrialFESpace(fs, dmesh)
    V = TestFESpace(U)

    # FEFunction
    # We need to trick a bit to share in memory
    # sol = FEFunction(U, get_values(bcube_io.data["u"])) # doesn't work because not an HauntedArray
    fs = bcube_io.data[first(keys(bcube_io.data))]
    arr = get_values(fs["u"])
    ex = BcubeMPI.get_exchanger(parent(U))
    l2g = BcubeMPI.local_to_global(parent(U))
    l2p = BcubeMPI.local_to_part(parent(U))
    o2l = BcubeMPI.own_to_local(parent(U))
    harr = HauntedArray(arr, ex, l2g, l2p, o2l)
    sol = FEFunction(U, harr)

    Ω = CellDomain(dmesh)
    Γ = InteriorFaceDomain(dmesh)
    dΩ = Measure(Ω, 2 * degree + 1)

    # WARNING: only `m` is a bilinear form, the others are linear form only. The argument `_u`
    # is only here to avoid a closure on the FEFunction and to obtain one big reusable linear form `l`
    m(u, v) = ∫(u ⋅ v)dΩ # Mass matrix
    l_Ω(v, _u) = ∫((c * _u) ⋅ ∇(v))dΩ # WARNING: this is not a bilinear form, this is just to avoid to capture u

    dΓ = Measure(Γ, 2 * degree + 1)
    nΓ = get_face_normals(Γ)
    flux(u) = upwind ∘ (side⁻(u), side⁺(u), c, side⁻(nΓ))
    l_Γ(v, _u) = ∫(flux(_u) * jump(v))dΓ # WARNING: this is not a bilinear form, this is just to avoid to capture u

    # All limits do not necessarily exists on all ranks so we have to deal with the different situations
    bnd_names = Bcube.boundary_names(dmesh)

    function l_Γ_in(v, t)
        if "Xmin" in bnd_names
            Γ_in = BoundaryFaceDomain(dmesh, "Xmin")
            dΓ_in = Measure(Γ_in, 2 * degree + 1)
            nΓ_in = get_face_normals(Γ_in)
            z = zeros(spacedim(parent(dmesh)))
            bc_in = t -> z
            return ∫(side⁻(bc_in(t)) ⋅ side⁻(nΓ_in) * side⁻(v))dΓ_in
        else
            return Bcube.NullOperator()
        end
    end

    function l_Γ_out(v, _u)
        if length(bnd_names) > 0
            Γ_out = BoundaryFaceDomain(dmesh, filter(x -> x != "Xmin", bnd_names))
            dΓ_out = Measure(Γ_out, 2 * degree + 1)
            nΓ_out = get_face_normals(Γ_out)
            flux_out(u) = upwind ∘ (side⁻(u), 0.0, c, side⁻(nΓ_out))
            return ∫(flux_out(_u) * side⁻(v))dΓ_out
        else
            return Bcube.NullOperator()
        end
    end

    l(v, _u, t) = l_Ω(v, _u) - l_Γ(v, _u) - l_Γ_in(v, t) - l_Γ_out(v, _u)

    # Compute mass matrix and invert
    # We cheat here because the mass matrix is actually diagonal so we compute the inverse
    # "by hand". Otherwise, we would require PETSc
    M = assemble_bilinear(m, U, V)
    Minv = M
    Minv.nzval .= 1 ./ Minv.nzval # manual inversion...
    Minv = HauntedMatrix(Minv, harr) # SparseMatrix -> Haunted(Sparse)Matrix

    # Compute time step
    CFL = 0.5
    d = compute_min_dimcar(dmesh)
    Δt = CFL * d / norm(c)
    @one_at_a_time (@show d, Δt)

    return (; sol, Minv, l, V, Δt)
end

function upwind(ui, uj, c, nij)
    cij = c ⋅ nij
    if cij > zero(cij)
        flux = cij * ui
    else
        flux = cij * uj
    end
    flux
end

function step_forward(solver)
    Δt = solver.Δt
    Minv = solver.Minv
    l = solver.l
    V = solver.V
    sol = solver.sol
    b = assemble_linear(v -> l(v, sol, 0.0), V)
    set_dof_values!(sol, get_dof_values(sol) .+ Δt .* (Minv * b))
end

function compute_min_dimcar(dmesh)
    fs = FunctionSpace(:Lagrange, 0)
    V = TestFESpace(fs, dmesh; size = 1, isContinuous = false)

    # Define measures for cell and interior face integrations
    degquad = 2
    dΩ = Measure(CellDomain(dmesh), degquad)
    dΓ = Measure(InteriorFaceDomain(dmesh), degquad)
    dΓ_bc = Measure(BoundaryFaceDomain(dmesh), degquad)

    f1 = PhysicalFunction(x -> 1.0)
    l(v) = ∫(f1 ⋅ v)dΩ
    function l_face(v)
        ∫(side⁻(f1) ⋅ side⁻(v) + side⁺(f1) ⋅ side⁺(v))dΓ + ∫(side⁻(f1) ⋅ side⁻(v))dΓ_bc
    end

    vol = assemble_linear(l, V)
    surf = assemble_linear(l_face, V)

    return MPI.Allreduce(minimum(vol ./ surf), MPI.MIN, BcubeMPI.get_comm(dmesh))
end
