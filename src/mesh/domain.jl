import BcubeCGNS.CGNS: Elements_t, get_children_from_predicate, get_child_from_name
import Bcube:
    MeshConnectivity, nnodes, get_element_type, connectivities_indices, connectivities

"""
Is `ctype` necessary here since we already have it as a type parameter?
`C` is a NamedTuple of connectivities (ex "cell-to-node" or "face-to-node")
"""
struct CGNSSubDomain{T <: Bcube.AbstractEntityType, C} <: Bcube.AbstractSubDomain
    etype::T
    connectivities::C
    indices::Vector{Int} # local (in the subdomain) to global (in the mesh-cell-numbering, or mesh-face-numbering)
end
Bcube.get_indices(sdom::CGNSSubDomain) = sdom.indices
Bcube.get_tag(::CGNSSubDomain) = 1 # tag is not really supported for now
Bcube.get_element_type(sdom::CGNSSubDomain) = sdom.etype
Bcube.connectivities(sdom::CGNSSubDomain) = sdom.connectivities
Bcube.connectivities(sdom::CGNSSubDomain, c::Symbol) = Bcube.connectivities(sdom)[c]
function Bcube.connectivities_indices(sdom::CGNSSubDomain, c::Symbol)
    Bcube.indices(Bcube.connectivities(sdom, c))
end

"""
`selection` are the elements, in the `eltNode`, that should populate this subdomain. Indeed the user might want only a portion of
the subdomain.
`l2g` is the local (in the subdomain) to global (in the mesh-cell-numbering, or mesh-face-numbering) correspondance
"""
function CGNSSubDomain(eltNode::CGNSNode{Elements_t}, selection, l2g)
    # Get cell type and number of nodes per element
    CT = CGNS_ENTITY_TO_BCUBE_ENTITY[first(BcubeCGNS.CGNS.get_value(eltNode))]
    ctype = CT()
    nnodes_by_elt = nnodes(ctype)

    # Get element range
    erange = BcubeCGNS.CGNS.get_value(get_child_from_name(eltNode, "ElementRange"))
    nelts = erange[2] - erange[1] + 1

    # Get full connectivity and range
    _cell2nodes =
        BcubeCGNS.CGNS.get_value(get_child_from_name(eltNode, "ElementConnectivity"))

    # Reshape to filter more easily
    mat = reshape(_cell2nodes, nnodes_by_elt, nelts)'

    # Build Bcube connectivity from selection
    cell2nodes = Array(vec(mat[selection, :]')) # flatten
    cell2nnodes = fill(eltype(cell2nodes)(nnodes_by_elt), length(selection))
    c2n = Bcube.Connectivity(cell2nnodes, cell2nodes)

    # Build MeshConnectivities
    conn = (; c2n = MeshConnectivity(:cell, :node, c2n))

    return CGNSSubDomain{CT, typeof(conn)}(ctype, conn, l2g)
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
    subdomain::CGNSSubDomain,
    i::Integer,
) where {D <: Bcube.AbstractCellDomain}
    # Recall that `i` designates the `i`th element of the subdomain.
    # This function is very similar to Bcube._get_cellinfo, but `c2n`
    # comes from the subdomain, not the mesh and `i` # `icell`
    mesh = get_mesh(domain)
    ctype = Bcube.get_element_type(subdomain)
    n_nodes = Val(nnodes(ctype))

    c2n = connectivities_indices(subdomain, :c2n)
    _c2n = c2n[i, n_nodes]

    cnodes = get_nodes(mesh, _c2n)
    icell = Bcube.get_indices(subdomain)[i]
    Bcube.CellInfo(icell, ctype, cnodes, _c2n)
end

Bcube.CellDomain(mesh::CGNSMesh, indices) = Bcube.build_cell_domain(mesh, indices)

"""
Note : `indices` not supported for now
"""
function Bcube.InteriorFaceDomain(mesh::CGNSMesh)
    Ω = Bcube.build_cell_domain(mesh, 1:ncells(mesh))

    # Create one (face-)subdomain for each (cell-)subdomain,
    # corresponding to the internal faces of theses subdomains
    sdoms_internal_f = map(Bcube.get_subdomains(Ω)) do sdom_vol
        ctype = get_element_type(sdom_vol)
        c2n = connectivities(sdom_vol, :c2n)
        ind = Bcube.indices(c2n)
        celltypes = fill(ctype, length(ind))
        _face_types, c2f, f2c, f2n = Bcube._build_faces!(c2n, celltypes)
        connectivities = (; c2n, c2f, f2c, f2n) # do we really need c2f? I don't think so...
        ftype = first(_face_types) # they all have the same type
        indices = zeros(Int, length(_face_types)) # `indices` is not supported yet
        CGNSSubdomain{typeof(ftype), typeof(connectivities)}(
            ftype,
            connectivities,
            zeros(Int, indices),
        )
    end

    # Identify "border" faces for each face-subdomain
    all_border_faces = map(sdoms_internal_f) do subdomain
        f2c = connectivities_indices(subdomain, :f2c)
        findall(_f2c -> length(_f2c) == 1, f2c)
    end

    # Loop over all vol subdomains combination
    sdoms_join_f = CGNSSubDomain[]
    for (i, sdom_i) in enumerate(sdoms_internal_f[1:(end - 1)])
        f2n_i = connectivities_indices(sdom_i, :f2n)
        f2c_i = connectivities_indices(sdom_i, :f2c)

        # Filter
        bnd_faces = all_border_faces[i]
        bnd_f2n_i = f2n_i[bnd_faces]
        bnd_f2c_i = f2c_i[bnd_faces]

        for (j, sdom_j) in enumerate(sdoms_internal_f[i + 1, :end])
            f2n_j = connectivities_indices(sdom_j, :f2n)

            # Filter
            bnd_faces = all_border_faces[j]
            bnd_f2n_j = f2n_j[bnd_faces]
            bnd_f2c_j = f2c_j[bnd_faces]

            identified_couples = zeros(Int, length(bnd_f2n_j), 2)
            n_identifier_coupled = 0

            # Note : I use a double for loop rather than an Iterators.product
            # in order the break the inner loop as soon as a match is found
            for (bnd_fi, nodes_i) in enumerate(bnd_f2n_i)
                for (bnd_fj, nodes_j) in enumerate(bnd_f2n_j)
                    if Set(nodes_i) == Set(nodes_j)
                        n_identifier_coupled += 1
                        identified_couples[n_identifier_coupled, 1] = bnd_fi
                        identified_couples[n_identifier_coupled, 2] = bnd_fj
                        break
                    end
                end
            end

            # Create a new CGNSSubDomain for the join "i-j"
            if n_identifier_coupled > 0
                c2n_n = connectivities(sdom_i, :f2n)
                c2n_p = connectivities(sdom_j, :f2n)
                ftype = get_element_type(sdom_i)
                f2c = map(
                    couple ->
                        (first(bnd_f2c_i[couple[1]]), first(bnd_f2c_j[couple[2]])),
                    eachrow(identified_couples[1:n_identifier_coupled, :]),
                )
                f2n = map(
                    couple -> bnd_f2n_i[couple[1]],
                    eachrow(identified_couples[1:n_identifier_coupled, :]),
                )
                #TODO: need to create Connectivity and then CGNSSubDomain
            end
        end
    end

    error(
        "we will miss faces between Elements_t, we could loop over volumic subdomains, concatenate all c2n, then call `_build_faces!` and then filter by face type",
    )

    error("not implemented yet")
    # Bcube._build_faces!(c2n::MeshConnectivity, celltypes)
end