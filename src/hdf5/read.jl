function Bcube.read_file(
    ::CGNSHDF5IoHandler,
    filepath::String;
    domains = String[],
    varnames = nothing,
    topodim = 0,
    spacedim = 0,
    verbose = false,
    kwargs...,
)
    @assert length(domains) == 0 "Reading only some domains is not supported yet (but easy to implement)"

    # Open the file
    file = h5open(filepath, "r")
    root = file

    # Read (unique) CGNS base
    cgnsBase = get_cgns_base(root)

    # Reference state
    refState = read_ref_state(cgnsBase)

    # Read base dimensions (topo and space)
    dims = get_value(cgnsBase)
    topodim = topodim > 0 ? topodim : dims[1]
    zone_space_dim = dims[2]
    verbose && println("topodim = $topodim, spacedim = $(zone_space_dim)")

    # Find the list of Zone_t
    zones = get_children(cgnsBase; type = "Zone_t")
    if length(zones) == 0
        error("Could not find any Zone_t node in the file")
    elseif length(zones) > 1
        error("The file contains several Zone_t nodes, only one zone is supported for now")
    end
    zone = first(zones)

    # Read zone
    zoneCGNS = read_zone(zone, varnames, topodim, zone_space_dim, spacedim, verbose)

    # Close the file
    close(file)

    # Build Bcube Mesh
    mesh = cgns_mesh_to_bcube_mesh(zoneCGNS)

    # Build Bcube MeshCellData & MeshPointData
    if !isnothing(varnames)
        data = flow_solutions_to_bcube_data(zoneCGNS.fSols)
    else
        data = nothing
    end

    # Should we return something when pointData and/or cellData is nothing? Or remove it completely from the returned Tuple?
    return (; mesh, data, refState)
end

"""
Read a CGNS Zone node.

Return a NamedTuple with node coordinates, cell-to-type connectivity (type is a integer),
cell-to-node connectivity, boundaries (see `read_zoneBC`), and a dictionnary of flow solutions (see `read_solutions`)

-> (; coords, c2t, c2n, c2g, bcs, fSols)

* `coords` is a matrix (nnodes, nspa) of the nodes coordinates
* `c2t` is the cell-to-type connectivity
* `c2n` is the cell-to-node connectivity (flattened, this a Vector)
* `c2g` is the cell-to-globnumber connectity (i-th cell of has the c2g[i] global number in the original numbering)
* `bcs` the BC
* `fSols` a dict of FlowSolutions

Note : the `zone_space_dim` is the number of spatial dimension according to the CGNS "Zone" node; whereas
`usr_space_dim` is the number of spatial dimensions asked by the user (0 to select automatic detection).
"""
function read_zone(zone, varnames, topo_dim, zone_space_dim, usr_space_dim, verbose)
    # Preliminary check
    zoneType = get_value(get_child(zone; type = "ZoneType_t"))

    if zoneType == "Unstructured"
        read_zone_unstr(zone, varnames, topo_dim, zone_space_dim, usr_space_dim, verbose)
    elseif zoneType == "Structured"
        read_zone_struct(zone, varnames, topo_dim, zone_space_dim, usr_space_dim, verbose)
    else
        error("Unhandled ZoneType ($zoneType)")
    end
end

function read_zone_unstr(zone, varnames, topo_dim, zone_space_dim, usr_space_dim, verbose)
    # Number of elements
    nvertices, ncells, nbnd = get_value(zone)
    verbose && println("nvertices = $nvertices, ncells = $ncells")

    # Read GridCoordinates
    gridCoordinates = get_child(zone; type = "GridCoordinates_t")
    coordXNode = get_child(gridCoordinates; name = "CoordinateX")
    coords = zeros(eltype(get_value(coordXNode)), nvertices, zone_space_dim)
    suffixes = ["X", "Y", "Z"]
    for (idim, suffix) in enumerate(suffixes[1:zone_space_dim])
        node = get_child(gridCoordinates; name = "Coordinate" * suffix)
        coords[:, idim] .= get_value(node)
    end

    # Resize the `coords` array if necessary
    _space_dim =
        usr_space_dim > 0 ? usr_space_dim : compute_space_dim(topo_dim, coords; verbose)
    coords = coords[:, 1:_space_dim]

    # Read all elements
    elts = map(read_connectivity, get_children(zone; type = "Elements_t"))

    # Filter "volumic" elements to build the volumic connectivities arrays
    volumicElts = filter(elt -> is_volumic_entity(first(elt.c2t), topo_dim), elts)
    @assert length(volumicElts) > 0 "Could not find elements of topo dim = $(topo_dim) in the file"
    c2t = mapreduce(elt -> elt.c2t, vcat, volumicElts)
    c2n = mapreduce(elt -> elt.c2n, vcat, volumicElts)
    c2g = mapreduce(elt -> collect(elt.erange), vcat, volumicElts)

    # Read all BCs and then keep only the ones whose topo dim is equal to the base topo dim minus 1
    bcs = read_zoneBC(zone, elts, verbose)
    if !isnothing(bcs)
        filter!(bc -> (bc.bcdim == topo_dim - 1) || (bc.bcdim == -1), bcs)
    end

    # Read FlowSolutions
    fSols = isnothing(varnames) ? nothing : read_solutions(zone, varnames, verbose)

    return (; coords, c2t, c2n, c2g, bcs, fSols)
end

function read_zone_struct(zone, varnames, topo_dim, zone_space_dim, usr_space_dim, verbose)
    # Number of elements
    zonedims = get_value(zone)
    nvertices_struct = zonedims[:, 1]
    nvertices = prod(nvertices_struct)
    ncells_struct = zonedims[:, 2]
    ncells = prod(ncells_struct)
    verbose && println(
        "nvertices = $(nvertices_struct) ($nvertices), ncells = $(ncells_struct) ($ncells)",
    )

    # Numbering conversion struct->unstr
    # tdims = ntuple(l -> 1:nvertices_struct[l], length(nvertices_struct))
    # ijk2I = LinearIndices(tdims)
    N = length(ncells_struct)
    tdims_cell = ntuple(l -> 1:ncells_struct[l], N)
    # tdims_cell = ntuple(l -> 1:ncells_struct[N - l + 1], N)
    ijk2I = LinearIndices(Dims(nvertices_struct))
    @warn "check if vec is the right function to apply to coordinates"

    # Read GridCoordinates
    gridCoordinates = get_child(zone; type = "GridCoordinates_t")
    coordXNode = get_child(gridCoordinates; name = "CoordinateX")
    X = flatten_struct_to_unstr(get_value(coordXNode))
    coords = zeros(eltype(X), nvertices, zone_space_dim)
    suffixes = ["X", "Y", "Z"]
    for (idim, suffix) in enumerate(suffixes[1:zone_space_dim])
        node = get_child(gridCoordinates; name = "Coordinate" * suffix)
        coords[:, idim] .= flatten_struct_to_unstr(get_value(node))
    end

    # Resize the `coords` array if necessary
    _space_dim =
        usr_space_dim > 0 ? usr_space_dim : compute_space_dim(topo_dim, coords; verbose)
    coords = coords[:, 1:_space_dim]

    # Read all elements
    if topo_dim == 3
        c2n = [
            (
                ijk2I[i, j, k],
                ijk2I[i + 1, j, k],
                ijk2I[i + 1, j + 1, k],
                ijk2I[i, j + 1, k],
                ijk2I[i, j, k + 1],
                ijk2I[i + 1, j, k + 1],
                ijk2I[i + 1, j + 1, k + 1],
                ijk2I[i, j + 1, k + 1],
            ) for (i, j, k) in Iterators.product(tdims_cell...)
        ]
        c2t = fill(BCUBE_ENTITY_TO_CGNS_ENTITY[Bcube.Hexa8_t], ncells)
    elseif topo_dim == 2
        c2n = [
            (ijk2I[i, j], ijk2I[i + 1, j], ijk2I[i + 1, j + 1], ijk2I[i, j + 1]) for
            (i, j) in Iterators.product(tdims_cell...)
        ]
        c2t = fill(BCUBE_ENTITY_TO_CGNS_ENTITY[Bcube.Quad4_t], ncells)
    else
        error("Structured file with topodim = $(topo_dim) not implemented yet")
    end
    c2n = collect(Iterators.flatten(c2n))

    # Read all BCs and then keep only the ones whose topo dim is equal to the base topo dim minus 1
    bcs = read_zoneBC(zone, nothing, verbose)
    if !isnothing(bcs)
        filter!(bc -> (bc.bcdim == topo_dim - 1) || (bc.bcdim == -1), bcs)
    end

    # Read FlowSolutions
    fSols = isnothing(varnames) ? nothing : read_solutions(zone, varnames, verbose)

    return (; coords, c2t, c2n, bcs, fSols)
end

"""
Read an "Elements_t" node and returns a named Tuple of three elements:
* `erange`, the content of the `ElementRange` node
* `c2t`, the cell -> entity_type connectivity
* `c2n`, the cell -> node connectivity, flattened if `reshape = false`, as an array (nelts, nnodes_by_elt) if `reshape = true`
* `name`, only for dbg
"""
function read_connectivity(node, reshape = false)
    @assert get_cgns_type(node) == "Elements_t"

    # Build cell to (cgns) type
    code, _ = get_value(node)
    erange = read_index(get_child(node; name = "ElementRange"))
    nelts = length(erange)
    c2t = fill(code, nelts)

    # Build cell to node and reshapce
    c2n = get_value(get_child(node; name = "ElementConnectivity"))

    nnodes_by_elt = nnodes(cgns_entity_to_bcube_entity(code))
    reshape && (c2n = reshape(c2n, nelts, nnodes_by_elt))

    return (; erange, c2t, c2n, nelts, name = get_name(node))
end

"""
    read_index(node)

Read an element of type "IndexRange_t" (a "PointList" or a "ElementRange" for instance).

Returns an Array (or an iterator) of the global number of the elements.
"""
function read_index(node)
    type = get_cgns_type(node)
    if type == "IndexRange_t"
        erange = get_value(node)
        return erange[1]:erange[2]
    elseif type == "IndexArray_t"
        return vec(get_value(node))
    else
        error("Reading $type not implemented yet")
    end
end

"""
Read the "ZoneBC_t" node to build bnd connectivities.

See `read_bc` for more information of what is returned.
"""
function read_zoneBC(zone, elts, verbose)
    zoneBC = get_child(zone; type = "ZoneBC_t")

    # Premature exit if no ZoneBC is present
    isnothing(zoneBC) && (return nothing)

    # Read each BC
    bcs = map(bc -> read_bc(bc, elts, verbose), get_children(zoneBC; type = "BC_t"))
    return bcs
end

"""
    read_bc(bc, elts, verbose)

Read a CGNS BC node.

`elts` is an input corresponding to the zone connectivity "Elements" nodes. For
an structured zone, this argument is `nothing`

Return a named Tuple (bcname, bcnodes, bcdim) where bcnodes is an array of the nodes
belonging to this BC.
"""
function read_bc(bc, elts, verbose)
    # BC name
    familyName = get_child(bc; type = "FamilyName_t")
    bcname = isnothing(familyName) ? get_name(bc) : get_value(familyName)
    verbose && println("Reading BC '$bcname'")

    # BC connectivity
    bc_type = get_value(get_child(bc; type = "GridLocation_t"))
    indexRange = get_child(bc; type = "IndexRange_t")
    pointList = get_child(bc; type = "IndexArray_t")

    # BC topodim : it's not always possible to determine it, so it's negative by default
    bcdim = -1

    if bc_type in ["CellCenter", "FaceCenter"]
        if !isnothing(indexRange)
            verbose && println("GridLocation is $(bc_type) with IndexRange")

            # This is a bit complex because nothing prevents an IndexRange to span over multiples Elements_t
            erange = read_index(indexRange)

            # Allocate the array of node indices corresponding to the BC
            nelts_bc = length(erange)
            T = eltype(first(elts).c2n[1])
            bcnodes = T[]
            sizehint!(bcnodes, nelts_bc * 4) # we assume 4 nodes by elements

            nelts_found = 0

            # Loop over all the Elements_t 'nodes'
            for elt in elts
                # verbose && println("Searching for elements in Elements_t '$(elt.name)'")
                i1 = elt.erange[1]
                i2 = elt.erange[end]
                etype = cgns_entity_to_bcube_entity(first(elt.c2t))
                nnodes_by_elt = nnodes(etype)

                if i1 <= erange[1] <= i2
                    # Compute how many nodes are concerned in this Elements_t,
                    # and the offset in the connectivity
                    nelts_concerned = min(i2, erange[end]) - erange[1] + 1
                    nnodes_concerned = nelts_concerned * nnodes_by_elt
                    offset = (erange[1] - i1) * nnodes_by_elt

                    push!(bcnodes, elt.c2n[(1 + offset):(offset + nnodes_concerned)]...)

                    nelts_found += nelts_concerned
                    verbose && println("$(nelts_concerned) elts found in '$(elt.name)'")

                    bcdim = Bcube.topodim(etype)

                    # Check if we've found all the elements in this connectivity
                    (erange[end] <= i2) && break
                end

                if i1 <= erange[end] <= i2
                    # Compute how many nodes are concerned in this Elements_t,
                    # and the offset in the connectivity
                    nelts_concerned = erange[end] - max(i1, erange[1]) + 1
                    nnodes_concerned = nelts_concerned * nnodes_by_elt
                    offset = (max(i1, erange[1]) - i1) * nnodes_by_elt

                    push!(bcnodes, elt.c2n[(1 + offset):(offset + nnodes_concerned)]...)

                    nelts_found += nelts_concerned
                    verbose && println("$(nelts_concerned) elts found in '$(elt.name)'")

                    bcdim = Bcube.topodim(etype)
                end
            end

            @assert nelts_found == nelts_bc "Missing elements for BC"

            # Once we've found all the nodes, we must remove duplicates
            # Note : using sort! + unique! is much more efficient than calling "unique"
            sort!(bcnodes)
            unique!(bcnodes)

        elseif !isnothing(pointList)
            # Elements indices
            bc_elts_ind = read_index(pointList)
            sort!(bc_elts_ind)

            # Allocate the array of node indices corresponding to the BC
            nelts_bc = length(bc_elts_ind)
            T = eltype(first(elts).c2n[1])
            bcnodes = T[]
            sizehint!(bcnodes, nelts_bc * 4) # we assume 4 nodes by elements

            icurr = 1
            for elt in elts
                verbose && println("Searching for elements in Elements_t '$(elt.name)'")
                i1 = elt.erange[1]
                i2 = elt.erange[end]
                etype = cgns_entity_to_bcube_entity(first(elt.c2t))
                nnodes_by_elt = nnodes(etype)

                (bc_elts_ind[icurr] < i1) && continue
                (bc_elts_ind[icurr] > i2) && continue

                if bc_elts_ind[end] <= i2
                    iEnd = bc_elts_ind[end]
                else
                    iEnd = findfirst(i -> i > i2, view(bc_elts_ind, icurr:nelts_bc)) - 1
                end
                offset = (bc_elts_ind[icurr] - i1) * nnodes_by_elt
                push!(bcnodes, elt.c2n[(1 + offset):(nnodes_by_elt * (iEnd - i1 + 1))]...)
                icurr = iEnd + 1

                (icurr > nelts_bc) && break

                # Element-wise version (OK, but very slow)
                # while i1 <= elts_ind[icurr] <= i2
                #     offset = (elts_ind[icurr] - i1) * nnodes_by_elt
                #     push!(bcnodes, elt.c2n[(1 + offset):(offset + nnodes_by_elt)]...)
                #     (icurr == nelts_bc) && break
                #     icurr += 1
                # end
            end

            @assert icurr >= nelts_bc
        else
            error("Could not find either the PointRange nor the PointList")
        end
    elseif bc_type == "Vertex"
        # todo : I am pretty sure we can remove the "if" here
        if !isnothing(pointList)
            bcnodes = collect(read_index(pointList))
        elseif !isnothing(indexRange)
            bcnodes = collect(read_index(indexRange))
        else
            error("Could not find either the PointRange nor the PointList")
        end

        # TODO : we could try to guess `bcdim` by search the Elements_t containing
        # the points of the PointList
    else
        error("BC GridLocation '$(bc_type)' not implemented")
    end

    return (; bcname, bcnodes, bcdim)
end

"""
Read all the flow solutions in the Zone, filtering data arrays whose name is not in the `varnames` list

# TODO : check if all varnames have been found
"""
function read_solutions(zone, varnames, verbose)
    # fSols =
    #     map(fs -> read_solution(fs, varnames), get_children(zone; type = "FlowSolution_t"))

    # n_vertex_fsol = count(fs -> fs.gridLocation == "Vertex", fSols)
    # n_cell_fsol = count(fs -> fs.gridLocation == "CellCenter", fSols)

    # if verbose
    #     (n_vertex_fsol > 1) && println(
    #         "WARNING : found more than one Vertex FlowSolution, reading them all...",
    #     )
    #     (n_cell_fsol > 1) && println(
    #         "WARNING : found more than one CellCenter FlowSolution, reading them all...",
    #     )
    # end

    fSols = Dict(
        get_name(fs) => read_solution(zone, fs, varnames) for
        fs in get_children(zone; type = "FlowSolution_t")
    )

    return fSols
end

"""
Read a FlowSolution node.

Return a NamedTuple with flow solution name, grid location and array of vectors.
"""
function read_solution(zone, fs, varnames)
    # Read GridLocation : we could deal with a missing GridLocation node, by later comparing
    # the length of the DataArray to the number of cells / nodes of the zone. Let's do this
    # later.
    node = get_child(fs; type = "GridLocation_t")
    if isnothing(node)
        _nnodes, _ncells, _ = get_zone_dims(zone)
        dArray = get_child(fs; type = "DataArray_t")
        err_msg = "Could not determine GridLocation in FlowSolution '$(get_name(fs))'"
        @assert !isnothing(dArray) err_msg
        x = get_value(dArray)
        if length(x) == _nnodes
            gridLocation = "Vertex"
        elseif length(x) == _ncells
            gridLocation = "CellCenter"
        else
            error(err_msg)
        end
        @warn "Missing GridLocation in FlowSolution '$(get_name(fs))', autoset to '$gridLocation'"
    else
        gridLocation = get_value(node)
    end

    # Read variables matching asked "varnames"
    dArrays = get_children(fs; type = "DataArray_t")
    if varnames != "*"
        # filter to obtain only the desired variables names
        filter!(dArray -> get_name(dArray) in varnames, dArrays)
    end
    data = Dict(
        get_name(dArray) => flatten_struct_to_unstr(get_value(dArray)) for
        dArray in dArrays
    ) # apply `vec` to handle structured data

    # Flow solution name
    name = get_name(fs)

    return (; name, gridLocation, data)
end

"""
Convert CGNS Zone information into a Bcube `Mesh`.
"""
function cgns_mesh_to_bcube_mesh(zoneCGNS)
    nodes = [Bcube.Node(zoneCGNS.coords[i, :]) for i in 1:size(zoneCGNS.coords, 1)]
    # nodes = map(row -> Bcube.Node(row), eachrow(zoneCGNS.coords)) # problem with Node + Slice

    c2n = Int.(zoneCGNS.c2n)
    c2t = map(cgns_entity_to_bcube_entity, zoneCGNS.c2t)
    c2nnodes = map(nnodes, c2t)

    if isnothing(zoneCGNS.bcs)
        return Bcube.Mesh(nodes, c2t, Bcube.Connectivity(c2nnodes, c2n))
    else
        bc_names = Dict(i => bc.bcname for (i, bc) in enumerate(zoneCGNS.bcs))
        bc_nodes = Dict(i => Int.(bc.bcnodes) for (i, bc) in enumerate(zoneCGNS.bcs))
        return Bcube.Mesh(nodes, c2t, Bcube.Connectivity(c2nnodes, c2n); bc_names, bc_nodes)
    end
end

"""
The input `fSols` is suppose to be a dictionnary FlowSolutionName => (gridlocation, Dict(varname => array))

The output is a dictionnary FlowSolutionName => dictionnary(varname => MeshData)
"""
function flow_solutions_to_bcube_data(fSols)
    # length(fSols) == 0 && (return Dict(), Dict())

    # pointSols = filter(fs -> fs.gridLocation == "Vertex", fSols)
    # cellSols = filter(fs -> fs.gridLocation == "CellCenter", fSols)

    # pointDicts = [fs.data for fs in pointSols]
    # cellDicts = [fs.data for fs in cellSols]

    # pointDict = merge(pointDicts)
    # cellDict = merge(cellDicts)

    # pointDict = Dict(key => MeshPointData(val) for (key, val) in pointDict)
    # cellDict = Dict(key => MeshCellData(val) for (key, val) in cellDict)

    # return pointDict, cellDict

    return Dict(
        fname => Dict(
            varname =>
                t.gridLocation == "Vertex" ? MeshPointData(array) : MeshCellData(array)
            for (varname, array) in t.data
        ) for (fname, t) in fSols
    )
end

function compute_space_dim(topodim, coords, tol = 1e-15; verbose = true)
    spacedim = size(coords, 2)

    xmin, xmax = extrema(view(coords, :, 1))
    lx = xmax - xmin

    ly = lz = 0.0

    if spacedim > 1
        ymin, ymax = extrema(view(coords, :, 2))
        ly = ymax - ymin
    end

    if spacedim > 2
        zmin, zmax = extrema(view(coords, :, 3))
        lz = zmax - zmin
    end

    return Bcube._compute_space_dim(topodim, lx, ly, lz, tol, verbose)
end

"""
Return nnodes, ncells, nbnd
"""
get_zone_dims(zone) = get_value(zone)

"""
    read_ref_state(base)

Read the ReferenceState node in the base.

Return a Dict of (String) keys and associated values. If the ReferenceState node doesn't exist,
`nothing` is returned.
"""
function read_ref_state(base)
    refState = get_child(base; name = "ReferenceState", type = "ReferenceState_t")
    return isnothing(refState) ? refState : _recursive_parse(refState, 0, 2)
end

read_grid_location_child(node) = get_value(get_child(node; type = "GridLocation_t"))

function _recursive_parse(node, depth, max_depth)
    # If it's a DataArray, return the value (scalar if necessary)
    if get_cgns_type(node) == "DataArray_t"
        x = get_value(node)
        return length(x) > 1 ? x : first(x)
    end

    if (get_cgns_type(node) == "Descriptor_t") ||
       has_child(node; type = "DimensionalExponents_t")
        return get_value(node)
    end

    if depth == max_depth
        @warn "Reached max depth when parsing ReferenceState, some entries might be missing"
        return get_value(node)
    end

    Dict(get_name(child) => _recursive_parse(child, depth + 1, max_depth) for child in node)
end

flatten_struct_to_unstr(array) = vec(array)