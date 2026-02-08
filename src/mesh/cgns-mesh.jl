import Bcube: AbstractMesh
import BcubeCGNS.CGNS as CGNS
import BcubeCGNS.CGNS: CGNSBase_t, Zone_t, get_children_from_label

abstract type AbstractCGNSMesh{topoDim, spaceDim} <: AbstractMesh{topoDim, spaceDim} end

"""
`Z` is an array or tuple of zones.

Only one zone supported for now
"""
struct CGNSMesh{topoDim, spaceDim, Z} <: AbstractCGNSMesh{topoDim, spaceDim}
    zones::Z
end

Base.parent(mesh::CGNSMesh) = mesh

function CGNSMesh(cgnsTree)
    base = first(get_children_from_label(cgnsTree, CGNSBase_t))
    @show typeof(base)
    topoDim, spaceDim = CGNS.get_value(base)
    zones = get_children_from_label(base, Zone_t)
    @assert length(zones) == 1 "Only one Zone_t supported for now"
    return CGNSMesh{topoDim, spaceDim, typeof(zones)}(zones)
end