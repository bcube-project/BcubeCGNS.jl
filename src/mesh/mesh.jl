import Bcube: AbstractMesh, BcubeBackendCPUSerial
import BcubeCGNS.CGNS as CGNS
import BcubeCGNS.CGNS:
    CGNSBase_t,
    CGNSNode,
    GridCoordinates_t,
    Zone_t,
    Elements_t,
    get_children_from_label,
    get_children_from_name,
    get_child_from_label,
    get_label

abstract type AbstractCGNSMesh{topoDim, spaceDim} <: AbstractMesh{topoDim, spaceDim} end

"""
`Z` is an array or tuple of zones.

`coords_by_dim` is, for each zone, a Tuple containing the arrays (CoordinateX, CoordinateY, ...)
So coords_by_dim[z][d][i] is the `d`-th space coordinate (x, y or z) of the `i`-th node of zone `z`.

Only one zone supported for now
"""
struct CGNSMesh{topoDim, spaceDim, Z, C, B} <: AbstractCGNSMesh{topoDim, spaceDim}
    zones::Z
    coords_by_dim::C
    backend::B
end

Bcube.get_bcube_backend(mesh::CGNSMesh) = mesh.backend
Base.parent(mesh::CGNSMesh) = mesh

function Bcube.ncells(mesh::CGNSMesh{topoDim}) where {topoDim}
    function isVolElts(n)
        (n isa CGNSNode{Elements_t}) || (return false)
        return is_volumic_entity(first(CGNS.get_value(n)), topoDim)
    end
    zone = first(mesh.zones)
    elts = get_children_from_predicate(zone, isVolElts)
    function countElts(elt)
        erange = CGNS.get_value(get_child_from_name(elt, "ElementRange"))
        return erange[2] - erange[1] + 1
    end
    return mapreduce(countElts, +, elts)
end

"""
Potential performance issue here!
"""
function Bcube.get_nodes(
    mesh::CGNSMesh{topoDim, spaceDim},
    indices,
) where {topoDim, spaceDim}
    @assert length(mesh.zones) == 1 "Only one Zone_t supported for now"
    iZone = 1
    return map(indices) do inode
        t = ntuple(d -> mesh.coords_by_dim[iZone][d][inode], spaceDim)
        Bcube.Node(t)
    end
end

function CGNSMesh(cgnsTree, backend = Bcube.get_bcube_backend())
    base = first(get_children_from_label(cgnsTree, CGNSBase_t))
    topoDim, spaceDim = CGNS.get_value(base)
    zones = get_children_from_label(base, Zone_t)
    @assert length(zones) == 1 "Only one Zone_t supported for now"

    # See the struct `CGNSMesh` for more info on `coords_by_dim`
    names = ("CoordinateX", "CoordinateY", "CoordinateZ")
    coords_by_dim = map(zones) do zone
        gc = get_child_from_label(zone, GridCoordinates_t)
        return ntuple(
            d -> BcubeCGNS.CGNS.get_value(get_child_from_name(gc, names[d])),
            spaceDim,
        )
    end

    return CGNSMesh{
        topoDim,
        spaceDim,
        typeof(zones),
        typeof(coords_by_dim),
        typeof(backend),
    }(
        zones,
        coords_by_dim,
        backend,
    )
end
