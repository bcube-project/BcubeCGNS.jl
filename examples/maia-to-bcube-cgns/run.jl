using Bcube
using BcubeCGNS

function build_bcube_mesh(tree_as_list)
    tree = BcubeCGNS.CGNS.parse_tree_as_list(tree_as_list)
    BcubeCGNS.CGNS.print_tree(tree)
    mesh = BcubeCGNS.CGNSMesh(tree)
    return mesh
end

function compute_gravity_center(mesh)
    dΩ = Measure(CellDomain(mesh), 1)
    x = sum(Bcube.compute(∫(PhysicalFunction(identity))dΩ))
    @show x
end
