using Bcube, BcubeCGNS, BcubeMPI

function build_bcube_dmesh(tree_as_list, ghost_tag2part, comm)
    # Parse the input tree as list into a BcubeCGNS.CGNS.Node
    tree = BcubeCGNS.CGNS.parse_tree_as_list(tree_as_list)
    BcubeCGNS.CGNS.print_tree(tree)

    # Convert it to a Bcube.Mesh
    result = BcubeCGNS.read_tree(tree; verbose = true)

    # Use additionnal infos to build the BcubeMPI.DistributedMesh
    dmesh = DistributedMesh(result.mesh, ghost_tag2part, comm)

    return dmesh
end