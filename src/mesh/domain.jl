import BcubeCGNS.CGNS: Elements_t, get_children_from_predicate, get_child_from_name

"""
Is `ctype` necessary here since we already have it as a type parameter?
"""
struct CGNSSubDomain{CT <: Bcube.AbstractEntityType, C2N} <: Bcube.AbstractSubDomain
    ctype::CT
    c2n::C2N
    indices::Vector{Int} # local (in the subdomain) to global (in the mesh-cell-numbering, or mesh-face-numbering)
end
Bcube.get_indices(sdom::CGNSSubDomain) = sdom.indices
Bcube.get_tag(::CGNSSubDomain) = 1 # tag is not really supported for now
Bcube.get_element_type(sdom::CGNSSubDomain) = sdom.ctype

"""
`selection` are the elements, in the `eltNode`, that should populate this subdomain. Indeed the user might want only a portion of
the subdomain.
`l2g` is the local (in the subdomain) to global (in the mesh-cell-numbering, or mesh-face-numbering) correspondance
"""
function CGNSSubDomain(eltNode::CGNSNode{Elements_t}, selection, l2g)
    CT = CGNS_ENTITY_TO_BCUBE_ENTITY[first(BcubeCGNS.CGNS.get_value(eltNode))]
    ctype = CT()
    nnodes_by_elt = Bcube.nnodes(ctype)
    c2n_flat = BcubeCGNS.CGNS.get_value(get_child_from_name(eltNode, "ElementConnectivity"))
    c2n = Matrix(reshape(c2n_flat, nnodes_by_elt, :)')
    # Note : we could also decide to keep `c2n_flat` an access the nodes inside easily
    return CGNSSubDomain{CT, typeof(c2n)}(ctype, view(c2n, selection, :), l2g)
end

"""
`indices` are to be interpreted as indices in the "cell numbering", not the global CGNS element numbering.

Assomptions:
* unique Zone_t

This function seems over-complicated, it needs to be simplified/clarified
"""
function Bcube.build_subdomains_by_celltypes(
    ::Bcube.AbstractBcubeBackend,
    mesh::CGNSMesh{topoDim},
    indices,
) where {topoDim}
    # Get all volumic element nodes
    function isVolElts(n)
        (n isa CGNSNode{Elements_t}) || (return false)
        return is_volumic_entity(first(CGNS.get_value(n)), topoDim)
    end
    zone = first(mesh.zones)
    elts = get_children_from_predicate(zone, isVolElts)

    # Build the map "cell-numbering" -> "cgns-elements-numbering"
    ranges = zeros(Int, length(elts), 2)
    for (elt, r) in zip(elts, eachrow(ranges))
        r .= BcubeCGNS.CGNS.get_value(get_child_from_name(elt, "ElementRange"))
    end
    nc = mapreduce(r -> r[2] - r[1] + 1, +, eachrow(ranges))
    cell2elts = zeros(Int, nc)
    offset = 0
    for r in eachrow(ranges)
        n = r[2] - r[1] + 1
        cell2elts[(offset + 1):(offset + n)] .= collect(r[1]:r[2])
        offset += n
    end

    # For each element node, find the cells included in the input `indices`
    result = map(elts) do elt
        r = BcubeCGNS.CGNS.get_value(get_child_from_name(elt, "ElementRange"))
        n = r[2] - r[1] + 1
        l2g = Int[]
        selected = Int[]
        sizehint!(l2g, n)
        sizehint!(selected, n)
        for iglo in view(cell2elts, indices)
            if r[1] ≤ iglo ≤ r[2]
                push!(l2g, iglo)
                push!(selected, iglo - r[1] + 1)
            end
        end
        return (; selected, l2g)
    end

    # Remove elts nodes that don't contain any "indices"
    ind = findall(info -> length(info.selected) > 0, result)

    return [
        CGNSSubDomain(elt, info.selected, info.l2g) for
        (elt, info) in zip(elts[ind], result[ind])
    ]
end

function Bcube._get_index(
    domain::D,
    subdomain::CGNSSubDomain{CT},
    i::Integer,
) where {D <: Bcube.AbstractCellDomain, CT}
    # Recall that `i` designates the `i`th element of the subdomain.
    mesh = get_mesh(domain)
    c2n = subdomain.c2n[i, :]
    cnodes = get_nodes(mesh, c2n)
    icell = Bcube.get_indices(subdomain)[i]
    ctype = Bcube.get_element_type(subdomain)
    Bcube.CellInfo(icell, ctype, cnodes, c2n)
end

Bcube.CellDomain(mesh::CGNSMesh, indices) = Bcube.build_cell_domain(mesh, indices)