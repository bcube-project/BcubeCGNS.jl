using BcubeCGNS
using Serialization

function translate(tree_as_list)
    tree = BcubeCGNS.CGNS.parse_tree_as_list(tree_as_list)
    BcubeCGNS.CGNS.print_tree(tree)
    mesh = BcubeCGNS.CGNSMesh(tree)
end